#!/usr/bin/env python3
"""VarPep v1.6 strand-aware editing, globally filtered haplotypes: all stages in ONE command; two BAM passes for exact targeted coverage.
Exploratory research software; not a validated somatic variant caller.
"""
import argparse, collections, csv, json, re, sqlite3, sys, hashlib, multiprocessing as mp, time, os, gzip, bisect, threading
from pathlib import Path
import pysam
import subprocess, shutil, tempfile, itertools

def parse_region(value):
    """User coordinates are 1-based, inclusive: chr:start-end."""
    if not value:return None
    m=re.fullmatch(r'([^:\s]+):([\d,]+)-([\d,]+)',value.strip())
    if not m:raise ValueError('--region must be chr:start-end (1-based inclusive)')
    chrom=m.group(1);start=int(m.group(2).replace(',',''));end=int(m.group(3).replace(',',''))
    if start<1 or end<start:raise ValueError('Invalid --region coordinates')
    return chrom,start-1,end

def overlaps_region(read,region):
    if region is None:return True
    chrom,start,end=region
    return (not read.is_unmapped and read.reference_name==chrom and
            read.reference_start is not None and read.reference_end is not None and
            read.reference_start<end and read.reference_end>start)

def progress_line(stage, fraction, started, label=''):
    fraction=max(0.0,min(1.0,fraction))
    elapsed=max(0.001,time.monotonic()-started)
    eta=elapsed*(1-fraction)/fraction if fraction>0.0001 else None
    eta_txt=f'{int(eta//3600):02}:{int(eta%3600//60):02}:{int(eta%60):02}' if eta is not None else '--:--:--'
    bar='#'*int(fraction*30)+'-'*(30-int(fraction*30))
    print(f'[progress] {stage} [{bar}] {fraction*100:5.1f}%  elapsed={int(elapsed)}s ETA={eta_txt} {label}',file=sys.stderr,flush=True)


COMP=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.upper().translate(COMP)[::-1]
CODE={}
for a,rows in {'T':['FFLL','SSSS','YY**','CC*W'],'C':['LLLL','PPPP','HHQQ','RRRR'],'A':['IIIM','TTTT','NNKK','SSRR'],'G':['VVVV','AAAA','DDEE','GGGG']}.items():
    for b,aas in zip('TCAG',rows):
        for c,aa in zip('TCAG',aas): CODE[a+b+c]=aa
def translate(s): return ''.join(CODE.get(s[i:i+3],'X') for i in range(0,len(s)-2,3)).rstrip('*')
def attributes(s): return dict(re.findall(r'(\w+) "([^"]+)"',s))
def txid(s): return s.split('|',1)[0]
def normalize(s): return s.split('.',1)[0]
def read_fai(f, protein=False):
    d={}
    with open(str(f)+'.fai') as h:
        for line in h:
            k=line.split('\t',1)[0]; fields=k.split('|'); ident=fields[1] if protein and len(fields)>1 else fields[0]; d[ident]=k
    return d

def build(args):
    for f in [args.transcripts,args.proteins]:
        if not Path(f+'.fai').exists(): pysam.faidx(f)
    ntkeys=read_fai(args.transcripts); aakeys=read_fai(args.proteins,protein=True)
    con=sqlite3.connect(args.index); cur=con.cursor()
    cur.executescript('''DROP TABLE IF EXISTS transcripts; DROP TABLE IF EXISTS exons; DROP TABLE IF EXISTS exon_lookup;
      CREATE TABLE transcripts(tx TEXT PRIMARY KEY,gene TEXT,chrom TEXT,strand TEXT,ntkey TEXT,aakey TEXT,cds_start INT,cds_end INT,valid INT,reason TEXT);
      CREATE TABLE exons(tx TEXT,chrom TEXT,start INT,end INT,tx_start INT,strand TEXT);
      CREATE INDEX exon_chrom_start ON exons(chrom,start,end);
      CREATE INDEX exon_tx ON exons(tx);
      CREATE TABLE exon_lookup(id INTEGER PRIMARY KEY, chrom TEXT,tx TEXT);
      CREATE VIRTUAL TABLE exon_rtree USING rtree(id,start,end);''')
    models={}
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'): continue
            a=line.rstrip('\n').split('\t')
            if len(a)<9 or a[2]!='exon': continue
            attr=attributes(a[8]); t=attr.get('transcript_id')
            if t not in ntkeys: continue
            m=models.setdefault(t,{'gene':attr.get('gene_name',attr.get('gene_id','')),'chrom':a[0],'strand':a[6],'exons':[]})
            m['exons'].append((int(a[3])-1,int(a[4])))
    nts=pysam.FastaFile(args.transcripts); aas=pysam.FastaFile(args.proteins)
    valid=0; invalid=collections.Counter(); n=0
    for t,m in models.items():
        key=ntkeys[t]; hdr=key.split('|'); cds=re.search(r'(?:^|\|)CDS:(\d+)-(\d+)(?:\||$)',key)
        if not cds: continue
        cs,ce=int(cds[1])-1,int(cds[2]); exons=sorted(m['exons'],reverse=m['strand']=='-')
        seq=nts.fetch(key); aa_key=aakeys.get(t)
        reason='OK'
        if not aa_key: reason='no_translation'
        elif sum(e-s for s,e in exons)!=len(seq): reason='exon_length_mismatch'
        elif not (0<=cs<ce<=len(seq)): reason='invalid_CDS_bounds'
        else:
            expected=aas.fetch(aa_key).rstrip('*'); observed=translate(seq[cs:ce]);
            if observed!=expected: reason='translation_mismatch'
        ok=int(reason=='OK'); valid+=ok; invalid[reason]+=1
        cur.execute('INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?)',(t,m['gene'],m['chrom'],m['strand'],key,aa_key,cs,ce,ok,reason))
        off=0
        for s,e in exons:
            cur.execute('INSERT INTO exons VALUES (?,?,?,?,?,?)',(t,m['chrom'],s,e,off,m['strand']))
            cur.execute('INSERT INTO exon_lookup(chrom,tx) VALUES (?,?)',(m['chrom'],t)); rid=cur.lastrowid
            cur.execute('INSERT INTO exon_rtree VALUES (?,?,?)',(rid,s,e-1)); off+=e-s
        n+=1
        if n%10000==0: con.commit(); print(f'[build] {n} transcripts, {valid} valid',file=sys.stderr,flush=True)
    con.commit(); con.close(); nts.close(); aas.close()
    print(json.dumps({'indexed_transcripts':n,'validated':valid,'status_counts':dict(invalid),'index':args.index},indent=2))

def collect_events(read,reference,minbq,include_indels):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail: return []
    seq=(read.query_sequence or '').upper(); q=read.query_qualities; events=[]
    # Fetch one contiguous genomic interval, not one FASTA fetch per aligned base.
    start=read.reference_start; end=read.reference_end
    if start is None or end is None or end<=start:return []
    refseq=reference.fetch(read.reference_name,start,end).upper()
    for qp,rp in read.get_aligned_pairs(matches_only=True,with_seq=False):
        if qp is None or rp is None or not (start<=rp<end):continue
        if q is not None and q[qp]<minbq:continue
        ref=refseq[rp-start]; alt=seq[qp]
        if ref in 'ACGT' and alt in 'ACGT' and ref!=alt:events.append((rp,'SNV',ref,alt))
    if include_indels:
        qpos=0; rpos=start
        for op,n in read.cigartuples or []:
            if op in (0,7,8):qpos+=n;rpos+=n
            elif op==1:
                ins=seq[qpos:qpos+n]
                if rpos>0 and ins and set(ins)<=set('ACGT') and (q is None or min(q[qpos:qpos+n])>=minbq):
                    events.append((rpos-1,'INS','',ins))
                qpos+=n
            elif op==2:
                deleted=refseq[rpos-start:rpos-start+n]
                if deleted and set(deleted)<=set('ACGT'):events.append((rpos,'DEL',deleted,''))
                rpos+=n
            elif op==3:rpos+=n
            elif op==4:qpos+=n
    return sorted(set(events))

def exons_for(con,chrom,pos):
    return con.execute('SELECT l.tx FROM exon_rtree r JOIN exon_lookup l ON r.id=l.id WHERE l.chrom=? AND r.start<=? AND r.end>=?',(chrom,pos,pos)).fetchall()

def project(con,tx,ev):
    rows=con.execute('SELECT start,end,tx_start,strand FROM exons WHERE tx=?',(tx,)).fetchall()
    pos,typ,ref,alt=ev
    for s,e,off,strand in rows:
        if s<=pos<e:
            tpos=off+(pos-s if strand=='+' else e-1-pos)
            if typ=='SNV':return (tpos,typ,ref if strand=='+' else rc(ref),alt if strand=='+' else rc(alt))
            if typ=='DEL':
                if pos+len(ref)>e: return None # boundary-spanning indel not handled
                tpos=off+(pos-s if strand=='+' else e-pos-len(ref))
                return (tpos,typ,ref if strand=='+' else rc(ref),'')
            if typ=='INS':
                # Inserted immediately AFTER genomic base pos.
                # For negative-strand transcripts insertion precedes its projected base.
                if strand=='+':return (tpos+1,typ,'',alt)
                return (tpos,typ,'',rc(alt))
    return None

def peptides(refaa,altaa,starts,ends):
    # Only peptide windows covering a changed amino-acid interval, not all downstream peptides.
    for cls,lens in [('I',range(8,12)),('II',range(13,26))]:
        for k in lens:
            for i in range(max(0,starts-k+1),min(len(altaa)-k,ends-1)+1):
                pep=altaa[i:i+k]
                if len(pep)!=k or '*' in pep or 'X' in pep:continue
                if i+k<=len(refaa) and refaa[i:i+k]==pep:continue
                yield (cls,k,i+1,pep)

_W_PROTEIN_ONLY=False
_W_ARGS=None
_W_DB=None
_W_NT=None
_W_AA=None
_W_GENOME=None
_W_HEADER=None

def _init_worker(index,transcripts,proteins,reference,bam,min_mapq,min_baseq,include_indels,region,protein_only=False,snp_index=None,redi_index=None,exclude_known=False):
    global _W_PROTEIN_ONLY
    _W_PROTEIN_ONLY=protein_only
    global _W_ARGS,_W_DB,_W_NT,_W_AA,_W_GENOME,_W_HEADER,_W_SNP,_W_REDI,_W_EXCLUDE
    _W_ARGS=(min_mapq,min_baseq,include_indels,region)
    _W_DB=sqlite3.connect('file:'+os.path.abspath(index)+'?mode=ro',uri=True)
    _W_DB.row_factory=sqlite3.Row
    _W_NT=pysam.FastaFile(transcripts);_W_AA=pysam.FastaFile(proteins)
    _W_GENOME=pysam.FastaFile(reference)
    with pysam.AlignmentFile(bam,'rb') as handle:_W_HEADER=handle.header
    _W_EXCLUDE=exclude_known
    _W_SNP=KnownSNPAnnotator(snp_index) if exclude_known else None
    _W_REDI=EditingAnnotator(redi_index) if exclude_known else None


def _worker_protein_only(task):
    """Count transcript-specific CDS haplotypes; never aggregate noncoding events."""
    batch_id, offset, batch = task
    haplotypes = collections.Counter()
    variants = collections.Counter()
    groups = 0
    if not hasattr(_worker_protein_only, '_sites'):
        _worker_protein_only._sites = {}
        _worker_protein_only._models = {}
        _worker_protein_only._junctions = {}
    sites, models, junctions_cache = _worker_protein_only._sites, _worker_protein_only._models, _worker_protein_only._junctions
    if len(sites)>100000: sites.clear()
    if len(models)>20000: models.clear()
    if len(junctions_cache)>20000: junctions_cache.clear()
    for sam_group in batch:
        reads = []
        for sam in sam_group:
            r = pysam.AlignedSegment.fromstring(sam, _W_HEADER)
            if (not (r.is_unmapped or r.is_secondary or r.is_supplementary or
                     r.is_duplicate or r.is_qcfail) and r.mapping_quality >= _W_ARGS[0]):
                reads.append(r)
        if not reads: continue
        groups += 1
        chrom = reads[0].reference_name
        if any(r.reference_name != chrom for r in reads): continue
        observed, junctions = set(), set()
        for r in reads:
            observed.update(ev for ev in collect_events(r, _W_GENOME, _W_ARGS[1], _W_ARGS[2])
                            if _W_ARGS[3] is None or
                            (chrom == _W_ARGS[3][0] and _W_ARGS[3][1] <= ev[0] < _W_ARGS[3][2]))
            rp = r.reference_start
            for op, n in r.cigartuples or ():
                if op in (0, 2, 7, 8): rp += n
                elif op == 3: junctions.add((rp, rp+n)); rp += n
        if not observed: continue
        # Union, not intersection: noncoding events must not invalidate CDS calls.
        candidates = set()
        for ev in observed:
            key = (chrom, ev[0])
            if key not in sites:
                sites[key] = frozenset(x[0] for x in exons_for(_W_DB, chrom, ev[0]))
            candidates.update(sites[key])
        coding_union = set()
        for tx in candidates:
            if tx not in models:
                row = _W_DB.execute('SELECT cds_start,cds_end,valid FROM transcripts WHERE tx=?',(tx,)).fetchone()
                models[tx] = row
            model = models[tx]
            if model is None or not model['valid']: continue
            if tx not in junctions_cache:
                rows = _W_DB.execute('SELECT start,end FROM exons WHERE tx=? ORDER BY start',(tx,)).fetchall()
                junctions_cache[tx] = frozenset((rows[i]['end'],rows[i+1]['start']) for i in range(len(rows)-1))
            if junctions and not junctions.issubset(junctions_cache[tx]): continue
            cds_events = []
            for ev in observed:
                if tx not in sites[(chrom,ev[0])]: continue
                projected = project(_W_DB, tx, ev)
                if projected is None: continue
                tpos, typ, ref, alt = projected
                cs, ce = model['cds_start'],model['cds_end']
                if typ == 'INS':
                    inside = cs < tpos < ce
                elif typ == 'DEL':
                    inside = cs <= tpos and tpos + len(ref) <= ce
                else:
                    inside = cs <= tpos < ce
                if inside: cds_events.append(ev)
            if cds_events:
                # Enumerate singletons and pairs observed on the same fragment.
                # Do NOT store full mixed haplotypes in strict mode.
                for subset in itertools.chain(((e,) for e in sorted(set(cds_events))),
                                              itertools.combinations(sorted(set(cds_events)), 2)):
                    haplotypes[(tx,chrom,json.dumps(subset,separators=(',',':')))] += 1
        # Count every coding event once per fragment, then remove those that
        # are not in a qualifying minimal set after translation.
        for ev in set(e for tx in candidates for e in observed if tx in sites[(chrom,e[0])]):
            variants[(chrom,str(ev[0]),ev[1],ev[2],ev[3])] += 1
    return batch_id, offset, groups, variants, haplotypes


def _worker_all(task):
    """Fast first pass: count events and fragment-level observed combinations only."""
    batch_id, offset, batch = task
    variants = collections.Counter()
    haplotypes = collections.Counter()
    groups = 0
    if not hasattr(_worker_all, '_sites'):
        _worker_all._sites = {}
        _worker_all._junctions = {}
    sites = _worker_all._sites
    junctions_cache = _worker_all._junctions
    if len(sites)>100000: sites.clear()
    if len(junctions_cache)>20000: junctions_cache.clear()
    for sam_group in batch:
        reads = []
        for sam in sam_group:
            r = pysam.AlignedSegment.fromstring(sam, _W_HEADER)
            if (not (r.is_unmapped or r.is_secondary or r.is_supplementary or
                     r.is_duplicate or r.is_qcfail) and r.mapping_quality >= _W_ARGS[0]):
                reads.append(r)
        if not reads:
            continue
        groups += 1
        chrom = reads[0].reference_name
        if any(r.reference_name != chrom for r in reads):
            continue
        observed = set()
        junctions = set()
        for r in reads:
            observed.update(ev for ev in collect_events(r, _W_GENOME, _W_ARGS[1], _W_ARGS[2])
                            if _W_ARGS[3] is None or
                            (chrom == _W_ARGS[3][0] and _W_ARGS[3][1] <= ev[0] < _W_ARGS[3][2]))
            rp = r.reference_start
            for op, n in r.cigartuples or ():
                if op in (0, 2, 7, 8): rp += n
                elif op == 3:
                    junctions.add((rp, rp+n)); rp += n
        if not observed:
            continue
        events = tuple(sorted(observed))
        for ev in events:
            variants[(chrom,) + ev] += 1
        candidates = None
        for pos, typ, ref, alt in events:
            key = (chrom, pos)
            if key not in sites:
                sites[key] = frozenset(x[0] for x in exons_for(_W_DB, chrom, pos))
            candidates = sites[key] if candidates is None else candidates & sites[key]
            if not candidates: break
        if not candidates: continue
        for tx in candidates:
            if junctions:
                if tx not in junctions_cache:
                    rows = _W_DB.execute('SELECT start,end FROM exons WHERE tx=? ORDER BY start', (tx,)).fetchall()
                    junctions_cache[tx] = frozenset((rows[i]['end'], rows[i+1]['start']) for i in range(len(rows)-1))
                if not junctions.issubset(junctions_cache[tx]): continue
            # Compact JSON avoids delimiter ambiguity with indel alleles.
            haplotypes[(tx, chrom, json.dumps(events, separators=(',', ':')))] += 1
    return batch_id, offset, groups, variants, haplotypes


def _worker(task):
    return _worker_protein_only(task) if _W_PROTEIN_ONLY else _worker_all(task)

def finalize_haplotypes(db, args):
    """Translate only globally supported haplotypes; reuse original biological rules."""
    from functools import lru_cache
    idx = sqlite3.connect('file:'+os.path.abspath(args.index)+'?mode=ro', uri=True)
    idx.row_factory = sqlite3.Row
    nt = pysam.FastaFile(args.transcripts)
    aa = pysam.FastaFile(args.proteins)
    @lru_cache(maxsize=10000)
    def model(tx):
        row = idx.execute('SELECT * FROM transcripts WHERE tx=?', (tx,)).fetchone()
        if row is None or not row['valid']: return None
        exons = idx.execute('SELECT start,end,tx_start,strand FROM exons WHERE tx=?', (tx,)).fetchall()
        return row, exons, nt.fetch(row['ntkey']), aa.fetch(row['aakey']).rstrip('*')
    @lru_cache(maxsize=25000)
    def protein(tx, events_json):
        m = model(tx)
        if m is None: return None
        row, exons, seq, refaa = m
        events = json.loads(events_json)
        projected = []
        for pos, typ, ref, alt in events:
            p = None
            for e in exons:
                start, end, off, strand = e
                if not start <= pos < end: continue
                tpos = off + (pos-start if strand=='+' else end-1-pos)
                if typ=='SNV': p=(tpos,typ,ref if strand=='+' else rc(ref),alt if strand=='+' else rc(alt))
                elif typ=='DEL' and pos+len(ref)<=end:
                    p=(off+(pos-start if strand=='+' else end-pos-len(ref)),typ,ref if strand=='+' else rc(ref),'')
                elif typ=='INS': p=(tpos+1 if strand=='+' else tpos,typ,'',alt if strand=='+' else rc(alt))
                break
            if p is None: return None
            projected.append(p)
        cs,ce=row['cds_start'],row['cds_end']
        coding=sorted((p for p in projected if cs<=p[0]<(ce+1 if p[1]=='INS' else ce)),key=lambda p:p[0])
        if not coding:return None
        cds=seq[cs:ce]; edits=[]
        for pos,typ,ref,alt in coding:
            o=pos-cs
            if (typ=='DEL' and pos+len(ref)>ce) or (ref and cds[o:o+len(ref)].upper()!=ref):return None
            edits.append((o,o+len(ref),alt))
        edits.sort(reverse=True);prev=len(cds)+1
        for a,b,alt in edits:
            if b>prev:return None
            cds=cds[:a]+alt+cds[b:];prev=a
        altaa=translate(cds)
        if refaa==altaa:return None
        diffs=[i for i,(a,b) in enumerate(zip(refaa,altaa)) if a!=b]
        if not diffs and len(refaa)==len(altaa):return None
        left=min(diffs) if diffs else min(len(refaa),len(altaa))
        right=max(diffs)+1 if len(refaa)==len(altaa) else len(altaa)
        return row['gene'],row['strand'],refaa,altaa,tuple(set(peptides(refaa,altaa,left,right)))
    processed=0
    snp_filter=KnownSNPAnnotator(args.snp151_index) if args.exclude_known else None
    redi_filter=EditingAnnotator(args.rediportal_index) if args.exclude_known else None
    try:
        cursor=db.execute('SELECT tx,chrom,events,n FROM raw_haplotypes WHERE n>=?', (args.min_support,))
        hc=collections.Counter();pc=collections.Counter();contexts=[]
        for tx,chrom,events_json,n in cursor:
            result=protein(tx,events_json)
            if result is None:continue
            gene,strand,refaa,altaa,peps=result
            events=json.loads(events_json)
            label=';'.join(f'{chrom}:{pos+1}:{typ}:{ref}>{alt}' for pos,typ,ref,alt in events)
            hc[(tx,gene,chrom,strand,label)]+=n
            contexts.append((tx,label,refaa,altaa))
            selected=True
            if args.nonsynonymous_only:
                selected=all(protein(tx,json.dumps([ev],separators=(',',':'))) is not None for ev in events)
            if selected and args.exclude_known:
                selected=all(not _exclude_event_with_annotators(chrom,ev,snp_filter,redi_filter) for ev in events)
            if selected:
                for cls,k,start,pep in peps:pc[(tx,gene,cls,k,start,pep,label)]+=n
            processed+=1
            if processed%2000==0:
                _upsert(db,'haplotypes',hc);_upsert(db,'peptides',pc)
                db.executemany('INSERT OR IGNORE INTO protein_context VALUES (?,?,?,?)',contexts)
                db.commit();hc.clear();pc.clear();contexts.clear()
                print(f'[translate] {processed:,} supported haplotypes',file=sys.stderr,flush=True)
        _upsert(db,'haplotypes',hc);_upsert(db,'peptides',pc)
        db.executemany('INSERT OR IGNORE INTO protein_context VALUES (?,?,?,?)',contexts)
        db.commit()
        db.execute('CREATE TABLE IF NOT EXISTS nonsyn_events(chrom TEXT,pos INT,typ TEXT,ref TEXT,alt TEXT,PRIMARY KEY(chrom,pos,typ,ref,alt)) WITHOUT ROWID')
        # Determine protein impact for single events across all expressed transcript contexts.
        found=set()
        for tx,chrom,events_json,n in db.execute('SELECT tx,chrom,events,n FROM raw_haplotypes WHERE n>=?',(args.min_support,)):
            for ev in json.loads(events_json):
                key=(chrom,*ev)
                if key in found:continue
                if protein(tx,json.dumps([ev],separators=(',',':'))) is not None:
                    found.add(key)
        db.executemany('INSERT OR IGNORE INTO nonsyn_events VALUES (?,?,?,?,?)',found)
        db.commit()
        # Complete catalog: include all observed transcript haplotypes, even without a protein change.
        # Existing translated entries are already stored with their support counts.
        existing={(r[0],r[2],r[4]) for r in db.execute('SELECT k0,k1,k2,k3,k4 FROM haplotypes')}
        batch=collections.Counter()
        for tx,chrom,events_json,n in db.execute('SELECT tx,chrom,events,n FROM raw_haplotypes WHERE n>=?',(args.min_support,)):
            m=model(tx)
            if m is None:continue
            row=m[0]
            label=';'.join(f'{chrom}:{pos+1}:{typ}:{ref}>{alt}' for pos,typ,ref,alt in json.loads(events_json))
            if (tx,chrom,label) not in existing:
                batch[(tx,row['gene'],chrom,row['strand'],label)]+=n
            if len(batch)>=10000:_upsert(db,'haplotypes',batch);batch.clear()
        _upsert(db,'haplotypes',batch)
        db.commit()
    finally:
        idx.close();nt.close();aa.close()
        if snp_filter:snp_filter.close()
        if redi_filter:redi_filter.close()
    print(f'[translate] completed {processed:,} protein-changing haplotypes',file=sys.stderr,flush=True)

def _batches(bam_path,name_sorted,batch_size,indexed_coordinate=False,region=None):
    """Stream name-grouped fragments, yielding batch index and BAM byte offset.
    Region restriction is read-overlap based, preserving mates in a group when
    either mate overlaps the region; worker events are restricted to region.
    """
    with pysam.AlignmentFile(bam_path,'rb') as bam:
        if region and region[0] not in bam.references:
            raise ValueError('Region chromosome not present in BAM: '+region[0])
        batch=[];last=None;group=[];batch_id=0
        if indexed_coordinate:
            if not bam.has_index():raise RuntimeError('Coordinate-sorted BAM requires an index')
            reads=(r for chrom in bam.references for r in bam.fetch(chrom))
        else:reads=bam.fetch(until_eof=True)
        for r in reads:
            if name_sorted:
                if last is not None and r.query_name!=last:
                    if group and (region is None or any(overlaps_region(x,region) for x in group)):
                        batch.append([x.to_string() for x in group])
                        if len(batch)>=batch_size:
                            yield batch_id,bam.tell()>>16,batch
                            batch_id+=1;batch=[]
                    group=[]
                group.append(r);last=r.query_name
            else:
                if overlaps_region(r,region):
                    batch.append([r.to_string()])
                    if len(batch)>=batch_size:
                        yield batch_id,bam.tell()>>16,batch
                        batch_id+=1;batch=[]
        if group and (region is None or any(overlaps_region(x,region) for x in group)):
            batch.append([x.to_string() for x in group])
        if batch:yield batch_id,bam.tell()>>16,batch


def _upsert(db,table,counts):
    if not counts:return
    width={'variants':5,'haplotypes':5,'peptides':7}[table]
    cols=','.join('k'+str(i) for i in range(width))
    keys=','.join('?' for _ in range(width))
    db.executemany(f'INSERT INTO {table} ({cols},n) VALUES ({keys},?) ON CONFLICT DO UPDATE SET n=n+excluded.n',
                   (tuple(k)+(v,) for k,v in counts.items()))

def _chrom_key(chrom):
    chrom=str(chrom).strip()
    return chrom[3:] if chrom.lower().startswith('chr') else chrom

def build_rediportal_index(args):
    """Stream TABLE1_hg38_v3 into a reusable position+strand SQLite index."""
    if os.path.exists(args.out):
        raise FileExistsError(f'Index already exists: {args.out}; remove it explicitly to rebuild')
    db=sqlite3.connect(args.out)
    try:
        db.execute('PRAGMA journal_mode=DELETE')
        db.execute('CREATE TABLE sites(chrom TEXT NOT NULL,pos INTEGER NOT NULL,strand TEXT NOT NULL,PRIMARY KEY(chrom,pos,strand)) WITHOUT ROWID')
        opener=gzip.open if args.rediportal.endswith('.gz') else open
        total=0;valid=0; batch=[]
        with opener(args.rediportal,'rt',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            if not reader.fieldnames or not {'Region','Position','Strand'}.issubset(reader.fieldnames):
                raise ValueError('REDIportal must have columns Region, Position, Strand; found '+str(reader.fieldnames))
            for row in reader:
                total+=1
                chrom=_chrom_key(row['Region'])
                strand={'1':'+','0':'-','+':'+','-':'-'}.get(str(row['Strand']).strip())
                try: pos=int(row['Position'])
                except (ValueError,TypeError):continue
                if not chrom or not strand or pos<1:continue
                batch.append((chrom,pos,strand));valid+=1
                if len(batch)>=50000:
                    db.executemany('INSERT OR IGNORE INTO sites VALUES (?,?,?)',batch);db.commit();batch=[]
                    if total%500000<50000:print(f'[rediportal-index] read={total:,}',file=sys.stderr,flush=True)
            if batch:db.executemany('INSERT OR IGNORE INTO sites VALUES (?,?,?)',batch)
        db.commit()
        unique=db.execute('SELECT COUNT(*) FROM sites').fetchone()[0]
        print(json.dumps({'rows_read':total,'valid_rows':valid,'unique_sites':unique,'index':args.out},indent=2),flush=True)
    except BaseException:
        db.close()
        try:os.remove(args.out)
        except OSError:pass
        raise
    else:db.close()

class EditingAnnotator:
    def __init__(self,path):
        self.db=sqlite3.connect('file:'+os.path.abspath(path)+'?mode=ro',uri=True)
        if not self.db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='sites'").fetchone():
            raise ValueError('Not a REDIportal index; use build-rediportal-index')
        self.cache={}
    def annotate(self,chrom,pos,typ,ref,alt,transcript_strand=None):
        """A-to-I candidate classification requires an explicit transcript strand.

        A genomic SNV by itself cannot establish transcript-oriented editing.
        """
        if typ!='SNV': return 'no','not_applicable'
        ref,alt=ref.upper(),alt.upper()
        if transcript_strand not in ('+','-'):
            return 'no','transcript_strand_unresolved'
        expected=('A','G') if transcript_strand=='+' else ('T','C')
        if (ref,alt)!=expected: return 'no','other_snv'
        key=(_chrom_key(chrom),int(pos),transcript_strand)
        if key not in self.cache:
            if len(self.cache)>250000:self.cache.clear()
            self.cache[key]=bool(self.db.execute('SELECT 1 FROM sites WHERE chrom=? AND pos=? AND strand=?',key).fetchone())
        return ('yes','known_editing_candidate') if self.cache[key] else ('no','possible_novel_editing')

    def variant_annotation(self,chrom,pos,typ,ref,alt):
        """Genomic variant annotation independent of transcript assignment.

        rediportal_v3 means the same position AND A-to-I-compatible
        genomic substitution on the REDIportal strand. A separate positional
        hit is reported even when the substitution does not agree.
        Transcript-oriented editing status remains unresolved here.
        """
        if typ!='SNV':
            return ['no','not_applicable','no','']
        key=(_chrom_key(chrom),int(pos))
        strands=tuple(row[0] for row in self.db.execute(
            'SELECT strand FROM sites WHERE chrom=? AND pos=? ORDER BY strand',key))
        compatible=('+' if (ref.upper(),alt.upper())==('A','G') else
                    '-' if (ref.upper(),alt.upper())==('T','C') else None)
        return ['yes' if compatible in strands else 'no',
                'transcript_strand_unresolved',
                'yes' if strands else 'no',
                ','.join(strands)]

    def event_labels(self,events,transcript_strand):
        n_known=n_possible=n_other=n_indel=0
        for event in events.split(';'):
            try:
                chrom,pos,typ,change=event.split(':',3)
                ref,alt=change.split('>',1)
                _,status=self.annotate(chrom,int(pos),typ,ref,alt,transcript_strand)
            except (ValueError,TypeError):
                n_other+=1;continue
            if status=='known_editing_candidate':n_known+=1
            elif status=='possible_novel_editing':n_possible+=1
            elif status=='not_applicable':n_indel+=1
            else:n_other+=1
        return [n_known,n_possible,n_other,n_indel,'yes' if n_known else 'no']
    def close(self):self.db.close()

def build_snp151_index(args):
    """UCSC snp151.txt[.gz]: bin,chrom,chromStart,chromEnd,name,...

    Index genomic loci of simple SNVs (genomic single), retaining rs identifiers.
    SNP151 is a positional annotation; it does not assert an allele match.
    """
    if os.path.exists(args.out):
        raise FileExistsError('Index already exists: '+args.out)
    db=sqlite3.connect(args.out)
    try:
        db.execute('CREATE TABLE snp_sites(chrom TEXT NOT NULL,pos INTEGER NOT NULL,rsid TEXT NOT NULL,PRIMARY KEY(chrom,pos,rsid)) WITHOUT ROWID')
        db.execute('CREATE INDEX snp_sites_lookup ON snp_sites(chrom,pos)')
        opener=gzip.open if args.snp151.endswith('.gz') else open
        count=0;kept=0;batch=[]
        with opener(args.snp151,'rt') as h:
            for line in h:
                if not line.strip() or line.startswith('#'):continue
                fields=line.rstrip('\n').split('\t')
                if len(fields)<12:raise ValueError('SNP151 needs at least 12 tab-separated UCSC columns')
                count+=1
                chrom,start,end,rsid=fields[1],fields[2],fields[3],fields[4]
                # UCSC snp151 class is field 12 (index 11).
                if fields[11]!='single':continue
                try:start,end=int(start),int(end)
                except ValueError:continue
                if end!=start+1 or start<0:continue
                batch.append((_chrom_key(chrom),start+1,rsid));kept+=1
                if len(batch)>=50000:
                    db.executemany('INSERT OR IGNORE INTO snp_sites VALUES (?,?,?)',batch)
                    db.commit();batch=[]
        if batch:db.executemany('INSERT OR IGNORE INTO snp_sites VALUES (?,?,?)',batch)
        db.commit()
        print(json.dumps({'rows_read':count,'simple_snv_rows':kept,'indexed_records':db.execute('SELECT COUNT(*) FROM snp_sites').fetchone()[0],'index':args.out},indent=2))
    except BaseException:
        db.close()
        try:os.remove(args.out)
        except OSError:pass
        raise
    else:db.close()

class KnownSNPAnnotator:
    def __init__(self,path):
        self.db=sqlite3.connect('file:'+os.path.abspath(path)+'?mode=ro',uri=True)
        if not self.db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='snp_sites'").fetchone():
            raise ValueError('Not a SNP151 index; use build-snp151-index')
        self.cache={}
    def annotate(self,chrom,pos,typ):
        if typ!='SNV':return 'no',''
        key=(_chrom_key(chrom),int(pos))
        if key not in self.cache:
            if len(self.cache)>250000:self.cache.clear()
            self.cache[key]=tuple(x[0] for x in self.db.execute('SELECT rsid FROM snp_sites WHERE chrom=? AND pos=? ORDER BY rsid',key))
        rsids=self.cache[key]
        return ('yes' if rsids else 'no'),';'.join(rsids)
    def event_labels(self,events):
        found=[]
        for ev in events.split(';'):
            try:
                chrom,pos,typ,_=ev.split(':',3)
                hit,rsids=self.annotate(chrom,int(pos),typ)
                if hit=='yes':found.append(rsids)
            except (ValueError,TypeError):continue
        return [len(found),'yes' if found else 'no','|'.join(found)]
    def close(self):self.db.close()

def _exclude_event_with_annotators(chrom,ev,snp,redi):
    pos,typ,ref,alt=ev
    if typ!='SNV':return False
    if snp and snp.annotate(chrom,pos+1,typ)[0]=='yes':return True
    strand='+' if (ref,alt)==('A','G') else '-' if (ref,alt)==('T','C') else None
    return bool(redi and strand and redi.annotate(chrom,pos+1,typ,ref,alt,strand)[0]=='yes')

def _exclude_event(chrom,ev):
    pos,typ,ref,alt=ev
    if typ!='SNV':return False
    if _W_SNP and _W_SNP.annotate(chrom,pos+1,typ)[0]=='yes':return True
    # REDIportal strand must agree with A>G (+) or T>C (-).
    strand='+' if (ref,alt)==('A','G') else '-' if (ref,alt)==('T','C') else None
    return bool(_W_REDI and strand and _W_REDI.annotate(chrom,pos+1,typ,ref,alt,strand)[0]=='yes')

def protein_changes(refaa, altaa):
    """Return 1-based protein changes and positions for an entire haplotype.

    Equal-length proteins: exact substitutions. Unequal-length proteins:
    report a common-prefix/suffix-trimmed delins, avoiding false claims about
    a particular frameshift or indel mechanism from protein strings alone.
    """
    if refaa == altaa:
        return '', ''
    if len(refaa) == len(altaa):
        changes = [(f'{r}{i+1}{a}', str(i+1)) for i, (r,a) in enumerate(zip(refaa, altaa)) if r != a]
        return ';'.join(c for c,_ in changes), ';'.join(pos for _,pos in changes)
    left = 0
    while left < min(len(refaa),len(altaa)) and refaa[left] == altaa[left]:
        left += 1
    right = 0
    while (right < len(refaa)-left and right < len(altaa)-left and
           refaa[len(refaa)-right-1] == altaa[len(altaa)-right-1]):
        right += 1
    old = refaa[left:len(refaa)-right if right else len(refaa)]
    new = altaa[left:len(altaa)-right if right else len(altaa)]
    # Exact minimal protein-level replacement, not an inferred HGVS frameshift.
    start = left + 1
    end = left + len(old)
    label = f'{start}_{end}:{old or "-"}>{new or "-"}'
    return label, f'{start}-{end}' if old else str(start)


def protein_contexts(db):
    return {(tx, events): protein_changes(ref, alt)
            for tx, events, ref, alt in db.execute(
                'SELECT tx,events,refaa,altaa FROM protein_context')}


def classify_event(chrom,pos,typ,ref,alt,strand,redi,snp):
    snp_hit,snp_ids=snp.annotate(chrom,pos,typ) if snp else ('unknown','')
    if typ!='SNV':kind='indel_unclassified'
    elif snp_hit=='yes':kind='dbsnp151_position_match'
    elif not redi or strand not in ('+','-'):kind='unresolved'
    else:
        _,edit=redi.annotate(chrom,pos,typ,ref,alt,strand)
        kind={'known_editing_candidate':'known_rediportal_editing',
              'possible_novel_editing':'possible_novel_editing'}.get(edit,'other_snv_not_in_dbsnp151' if snp_hit=='no' else 'other_snv_unresolved')
    return kind,snp_ids

def event_annotations(events,strand,redi,snp):
    kinds=[];rsids=[]
    for item in events.split(';'):
        chrom,pos,typ,change=item.split(':',3)
        ref,alt=change.split('>',1)
        kind,ids=classify_event(chrom,int(pos),typ,ref,alt,strand,redi,snp)
        kinds.append(kind);rsids.append(ids)
    return ';'.join(kinds),';'.join(rsids)

def candidate_flag(chrom,pos,typ,ref,alt,redi,snp,non_syn):
    if not non_syn:return 'no'
    if typ!='SNV':return 'yes' # no SNP151/REDIportal classification for indels
    if not (redi and snp):return 'unresolved'
    if snp.annotate(chrom,pos,typ)[0]=='yes':return 'no'
    if redi.variant_annotation(chrom,pos,typ,ref,alt)[0]=='yes':return 'no'
    return 'yes'

def _snapshot(db,prefix,min_support=1,rediportal=None,snp151=None):
    strand_by_tx={tx:strand for tx,strand in db.execute('SELECT DISTINCT k0,k3 FROM haplotypes')} if rediportal is not None else {}
    changes_by_context=protein_contexts(db)
    has_nonsyn=bool(db.execute("SELECT 1 FROM sqlite_master WHERE name='nonsyn_events'").fetchone())
    for table,header in [('variants',['chrom','pos1','type','ref','alt','fragment_support']),
                         ('haplotypes',['transcript','gene','chrom','strand','events','fragment_support']),
                         ('peptides',['transcript','gene','mhc_class','length','aa_start1','peptide','events','fragment_support'])]:
        nkeys=len(header)-1
        target=prefix+'.'+table+'.tsv';tmp=target+'.tmp'
        with open(tmp,'w',newline='') as f:
            w=csv.writer(f,delimiter='\t')
            if rediportal is not None:
                header=header+(['rediportal_v3','editing_candidate','rediportal_position_match','rediportal_strands'] if table=='variants' else ['n_known_editing_sites','n_possible_novel_editing_sites','n_other_snv','n_indels','contains_known_editing'])
            if snp151 is not None:
                header=header+(['snp151_position_match','snp151_rsids'] if table=='variants' else ['n_snp151_sites','contains_snp151','snp151_rsids'])
            if table=='variants':header+=['nonsynonymous_in_any_transcript','potential_neoantigen_candidate']
            elif table=='haplotypes':header+=['event_classifications','event_rsids_ordered','contains_potential_neoantigen_event']
            if table=='peptides': header=header+['aa_change_ref_to_mut','aa_position1']
            w.writerow(header)
            for row in db.execute('SELECT * FROM '+table+' WHERE n >= ? ORDER BY '+','.join('k'+str(i) for i in range(nkeys)), (min_support,)):
                v=list(row)
                if table=='variants':
                    v[1]=int(v[1])+1
                    if rediportal is not None:v.extend(rediportal.variant_annotation(v[0],v[1],v[2],v[3],v[4]))
                elif rediportal is not None:
                    v.extend(rediportal.event_labels(v[4],v[3]) if table=='haplotypes' else rediportal.event_labels(v[6],strand_by_tx.get(v[0])))
                if snp151 is not None:
                    v.extend(snp151.annotate(v[0],v[1],v[2]) if table=='variants' else snp151.event_labels(v[4] if table=='haplotypes' else v[6]))
                if table=='variants':
                    nonsyn=bool(db.execute('SELECT 1 FROM nonsyn_events WHERE chrom=? AND pos=? AND typ=? AND ref=? AND alt=?',(v[0],v[1]-1,v[2],v[3],v[4])).fetchone()) if has_nonsyn else False
                    v.extend(['yes' if nonsyn else 'no',candidate_flag(v[0],v[1],v[2],v[3],v[4],rediportal,snp151,nonsyn)])
                elif table=='haplotypes':
                    kinds,ids=event_annotations(v[4],v[3],rediportal,snp151)
                    v.extend([kinds,ids,'yes' if has_nonsyn and any((k in ('possible_novel_editing','other_snv_not_in_dbsnp151') and bool(db.execute('SELECT 1 FROM nonsyn_events WHERE chrom=? AND pos=? AND typ=? AND ref=? AND alt=?',(chrom,int(pos)-1,typ,*change.split('>',1))).fetchone())) for k,(chrom,pos,typ,change) in zip(kinds.split(';'),(e.split(':',3) for e in (v[4] if table=='haplotypes' else events).split(';')))) else 'no'])
                if table=='peptides': v.extend(changes_by_context.get((v[0],v[6]), ('','')))
                w.writerow(v)
        os.replace(tmp,target)

def call(args):
    if args.threads<1:raise ValueError('--threads must be >=1')
    if args.name_sorted and args.indexed_coordinate:raise ValueError('--name-sorted and --indexed-coordinate are incompatible')
    if args.batch_size<1:raise ValueError('--batch-size must be >=1')
    if args.min_support<1:raise ValueError('--min-support must be >=1')
    if args.max_inflight_batches<args.threads:raise ValueError('--max-inflight-batches must be >= --threads')
    if args.commit_every<1:raise ValueError('--commit-every must be >= 1')
    prefix=args.out_prefix
    state=prefix+'.counts.sqlite'
    if os.path.exists(state):
        raise RuntimeError('Existing '+state+'; choose a fresh --out-prefix (avoids double counting)')
    db=sqlite3.connect(state)
    db.execute('PRAGMA journal_mode=WAL')
    for table,width in [('variants',5),('haplotypes',5),('peptides',7)]:
        keys=','.join('k'+str(i)+' TEXT NOT NULL' for i in range(width))
        pk=','.join('k'+str(i) for i in range(width))
        db.execute(f'CREATE TABLE {table} ({keys}, n INTEGER NOT NULL, PRIMARY KEY ({pk})) WITHOUT ROWID')
    db.execute('CREATE TABLE raw_haplotypes (tx TEXT, chrom TEXT, events TEXT, n INTEGER NOT NULL, PRIMARY KEY(tx,chrom,events)) WITHOUT ROWID')
    db.execute('CREATE TABLE protein_context (tx TEXT NOT NULL, events TEXT NOT NULL, refaa TEXT NOT NULL, altaa TEXT NOT NULL, PRIMARY KEY(tx,events)) WITHOUT ROWID')
    db.commit()
    # Initial headers only; variants are populated after protein-level filtering.
    rp=EditingAnnotator(args.rediportal_index) if args.rediportal_index else None
    snp=KnownSNPAnnotator(args.snp151_index) if args.snp151_index else None
    _snapshot(db,prefix,args.min_support,rp,snp)
    done=0;chunks=0;last=time.monotonic();started=last
    pending_v=collections.Counter();pending_h=collections.Counter()
    def flush_counts():
        if not pending_v and not pending_h:return
        _upsert(db,'variants',pending_v)
        db.executemany('INSERT INTO raw_haplotypes VALUES (?,?,?,?) ON CONFLICT DO UPDATE SET n=n+excluded.n',
                       (k+(v,) for k,v in pending_h.items()))
        db.commit();pending_v.clear();pending_h.clear()
    total_bytes=os.path.getsize(args.bam);completed_offsets={};next_completed=0;contiguous_offset=0
    ctx=mp.get_context('spawn')
    init=(args.index,args.transcripts,args.proteins,args.reference,args.bam,args.min_mapq,args.min_baseq,args.include_indels,args.region_parsed,False,args.snp151_index,args.rediportal_index,False)
    try:
        # Backpressure: the feeder may submit at most max_inflight_batches
        # tasks that have not yet been consumed and merged by the writer.
        # The generator is pulled by Pool's feeder thread; the main thread
        # releases one slot only AFTER writing the completed result to SQLite.
        slots=threading.BoundedSemaphore(args.max_inflight_batches)
        def bounded_batches():
            for batch in _batches(args.bam,args.name_sorted,args.batch_size,args.indexed_coordinate,args.region_parsed):
                slots.acquire()
                yield batch
        with ctx.Pool(processes=args.threads,initializer=_init_worker,initargs=init) as pool:
            iterator=pool.imap_unordered(_worker,bounded_batches(),chunksize=1)
            for batch_id,byte_offset,n,vc,hc in iterator:
                done+=n;chunks+=1
                completed_offsets[batch_id]=byte_offset
                while next_completed in completed_offsets:
                    contiguous_offset=completed_offsets.pop(next_completed)
                    next_completed+=1
                pending_v.update(vc);pending_h.update(hc)
                if chunks % args.commit_every == 0:
                    flush_counts()
                slots.release()
                now=time.monotonic()
                if now-last>=args.flush_seconds:
                    db.commit()
                    progress_line('1/2 calling',contiguous_offset/total_bytes,started,f'groups={done:,} batches={chunks:,}')
                    last=now
        flush_counts()
        print(json.dumps({'groups':done,'batches':chunks,'seconds':round(time.monotonic()-started,1),'mode':'fragment' if args.name_sorted else 'read','state':state},indent=2))
    except KeyboardInterrupt:
        print('[interrupt] Saved accumulated counts in '+state+'; last TSV snapshot remains available',file=sys.stderr,flush=True)
        flush_counts()
        _snapshot(db,prefix,args.min_support,rp,snp)
        raise
    finally:
        db.close()
        if rp is not None:rp.close()
        if snp is not None:snp.close()

def export(args):
    if args.min_support<1:raise ValueError('--min-support must be >=1')
    if not os.path.isfile(args.state):raise FileNotFoundError(args.state)
    db=sqlite3.connect('file:'+os.path.abspath(args.state)+'?mode=ro',uri=True)
    rp=EditingAnnotator(args.rediportal_index) if args.rediportal_index else None
    snp=KnownSNPAnnotator(args.snp151_index) if args.snp151_index else None
    try:
        _snapshot(db,args.out_prefix,args.min_support,rp,snp)
    finally:
        db.close()
        if rp is not None:rp.close()
        if snp is not None:snp.close()
    print(f'[export] min_support={args.min_support} prefix={args.out_prefix}',file=sys.stderr,flush=True)


def coverage_for_sites(bam_path, state, min_mapq, min_baseq, progress=100000, region=None):
    """Depth at variant loci: distinct QNAMEs with >=1 high-quality aligned base.
    INS: anchor genomic base; DEL: base immediately preceding deletion.
    """
    db=sqlite3.connect('file:'+os.path.abspath(state)+'?mode=ro',uri=True)
    targets=collections.defaultdict(set)
    for chrom,pos,typ,ref,alt in db.execute('SELECT k0,k1,k2,k3,k4 FROM variants'):
        pos=int(pos)
        anchor=max(0,pos-1) if typ=='DEL' else pos
        targets[chrom].add(anchor)
    db.close()
    positions={c:sorted(v) for c,v in targets.items()}
    depth=collections.Counter(); group=[]; prev=None; n=0
    stage_started=time.monotonic();bam_size=os.path.getsize(bam_path)
    def process(group):
        covered=set()
        for r in group:
            if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or
                r.is_qcfail or r.mapping_quality<min_mapq):continue
            chrom=r.reference_name; sites=positions.get(chrom)
            if not sites or r.reference_start is None or r.reference_end is None:continue
            lo=bisect.bisect_left(sites,r.reference_start)
            hi=bisect.bisect_left(sites,r.reference_end)
            if lo==hi:continue
            wanted=set(sites[lo:hi])
            qual=r.query_qualities
            for qp,rp in r.get_aligned_pairs(matches_only=True,with_seq=False):
                if rp in wanted and qp is not None and (qual is None or qual[qp]>=min_baseq):
                    covered.add((chrom,rp))
        for k in covered:depth[k]+=1
    with pysam.AlignmentFile(bam_path,'rb') as bam:
        for r in bam.fetch(until_eof=True):
            if prev is not None and r.query_name!=prev:
                process(group);group=[];n+=1
                if progress and n%progress==0:
                    progress_line('2/2 coverage',(bam.tell()>>16)/bam_size,stage_started,f'fragments={n:,} sites={sum(map(len,positions.values())):,}')
            if overlaps_region(r,region):group.append(r)
            prev=r.query_name
        if group:process(group);n+=1
    progress_line('2/2 coverage',1.0,stage_started,f'fragments={n:,}')
    return depth,n

def export_vaf(state,prefix,depth,min_support,min_vaf,rediportal_index=None,snp151_index=None):
    if not (0<=min_vaf<=1):raise ValueError('--min-vaf must be between 0 and 1')
    if min_support<1:raise ValueError('--min-support must be >= 1')
    db=sqlite3.connect('file:'+os.path.abspath(state)+'?mode=ro',uri=True)
    rp=EditingAnnotator(rediportal_index) if rediportal_index else None
    snp=KnownSNPAnnotator(snp151_index) if snp151_index else None
    event_vafs={}
    strand_by_tx={tx:strand for tx,strand in db.execute('SELECT DISTINCT k0,k3 FROM haplotypes')} if rp else {}
    changes_by_context=protein_contexts(db)
    has_nonsyn=bool(db.execute("SELECT 1 FROM sqlite_master WHERE name='nonsyn_events'").fetchone())
    try:
        for chrom,pos,typ,ref,alt,n in db.execute('SELECT k0,k1,k2,k3,k4,n FROM variants'):
            pos=int(pos);n=int(n)
            anchor=max(0,pos-1) if typ=='DEL' else pos
            den=depth.get((chrom,anchor),0)
            # Numerator may exceed anchor depth for some complex indels; mark undefined.
            vaf=n/den if den>=n and den>0 else None
            event_vafs[(chrom,pos,typ,ref,alt)]=(n,den,vaf)
        for table,headers in (
            ('variants',['chrom','pos1','type','ref','alt','fragment_support']),
            ('haplotypes',['transcript','gene','chrom','strand','events','fragment_support']),
            ('peptides',['transcript','gene','mhc_class','length','aa_start1','peptide','events','fragment_support'])):
            dest=prefix+'.'+table+'.tsv';tmp=dest+'.tmp'
            with open(tmp,'w',newline='') as f:
                w=csv.writer(f,delimiter='\t')
                extra=(['fragment_depth','allele_frequency','depth_definition'] if table=='variants' else
                       ['event_allele_frequencies','min_event_allele_frequency'])
                redi=(['rediportal_v3','editing_candidate','rediportal_position_match','rediportal_strands'] if table=='variants' else
                      ['n_known_editing_sites','n_possible_novel_editing_sites','n_other_snv','n_indels','contains_known_editing']) if rp else []
                snpcols=(['snp151_position_match','snp151_rsids'] if table=='variants' else ['n_snp151_sites','contains_snp151','snp151_rsids']) if snp else []
                proteincols=['aa_change_ref_to_mut','aa_position1'] if table=='peptides' else []
                newcols=(['nonsynonymous_in_any_transcript','potential_neoantigen_candidate'] if table=='variants' else ['event_classifications','event_rsids_ordered','contains_potential_neoantigen_event'] if table=='haplotypes' else [])
                w.writerow(headers+extra+redi+snpcols+newcols+proteincols)
                for row in db.execute('SELECT * FROM '+table+' WHERE n>=?',(min_support,)):
                    values=list(row)
                    if table=='variants':
                        chrom,pos,typ,ref,alt,n=values
                        n,den,af=event_vafs[(chrom,int(pos),typ,ref,alt)]
                        if af is None or af<min_vaf:continue
                        values[1]=int(pos)+1
                        values.extend([den,round(af,6),'aligned_base' if typ=='SNV' else 'indel_anchor_approx'])
                        if rp:values.extend(rp.variant_annotation(chrom,int(pos)+1,typ,ref,alt))
                        if snp:values.extend(snp.annotate(chrom,int(pos)+1,typ))
                    else:
                        events=values[4] if table=='haplotypes' else values[6]
                        afs=[];bad=False
                        for item in events.split(';'):
                            try:
                                chrom,pos,typ,change=item.split(':',3)
                                ref,alt=change.split('>',1)
                                _,_,af=event_vafs[(chrom,int(pos)-1,typ,ref,alt)]
                            except (ValueError,KeyError):af=None
                            if af is None or af<min_vaf:bad=True
                            afs.append(af)
                        if bad:continue
                        values.extend([';'.join(f'{x:.6f}' for x in afs),round(min(afs),6) if afs else ''])
                        if rp:values.extend(rp.event_labels(events,values[3] if table=='haplotypes' else strand_by_tx.get(values[0])))
                        if snp:values.extend(snp.event_labels(events))
                        if table=='peptides': values.extend(changes_by_context.get((values[0],values[6]), ('','')))
                    if table=='variants':
                        nonsyn=bool(db.execute('SELECT 1 FROM nonsyn_events WHERE chrom=? AND pos=? AND typ=? AND ref=? AND alt=?',(chrom,int(pos),typ,ref,alt)).fetchone()) if has_nonsyn else False
                        values.extend(['yes' if nonsyn else 'no',candidate_flag(chrom,int(pos)+1,typ,ref,alt,rp,snp,nonsyn)])
                    elif table=='haplotypes':
                        kinds,ids=event_annotations(events,values[3],rp,snp)
                        values.extend([kinds,ids,'yes' if has_nonsyn and any((k in ('possible_novel_editing','other_snv_not_in_dbsnp151') and bool(db.execute('SELECT 1 FROM nonsyn_events WHERE chrom=? AND pos=? AND typ=? AND ref=? AND alt=?',(c,int(pos)-1,t,*change.split('>',1))).fetchone())) for k,(c,pos,t,change) in zip(kinds.split(';'),(e.split(':',3) for e in events.split(';')))) else 'no'])
                    w.writerow(values)
            os.replace(tmp,dest)
            print(f'[vaf-export] {dest}',file=sys.stderr,flush=True)
    finally:
        db.close()
        if rp:rp.close()
        if snp:snp.close()


def process_peptides(state,prefix,pepsickle_bin='pepsickle',model='epitope',proteasome='I',threshold=0.5,batch_sequences=250):
    """Score full-length WT and mutant proteins with real pepsickle model.

    No fallback heuristic: missing predictor causes an explicit error. All MHC-I
    theoretical peptides remain in peptides_all.tsv; scored candidates are in
    peptides_processed.tsv. MHC-II stays in the theoretical output only.
    """
    if shutil.which(pepsickle_bin) is None and not os.path.isfile(pepsickle_bin):
        raise RuntimeError('pepsickle not found: install it in the environment or supply --pepsickle-bin')
    db=sqlite3.connect('file:'+os.path.abspath(state)+'?mode=ro',uri=True)
    try:
        contexts={(tx,ev):(ref,alt) for tx,ev,ref,alt in db.execute('SELECT tx,events,refaa,altaa FROM protein_context')}
    finally:db.close()
    source=prefix+'.peptides.tsv'
    if not os.path.exists(source):raise FileNotFoundError(source)
    with open(source,newline='') as h:
        reader=csv.DictReader(h,delimiter='\t');columns=reader.fieldnames
        rows=list(reader)
    all_dest=prefix+'.peptides_all.tsv'
    with open(all_dest+'.tmp','w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=columns,delimiter='\t');writer.writeheader();writer.writerows(rows)
    os.replace(all_dest+'.tmp',all_dest)
    # Digest only protein contexts that actually contributed retained peptide candidates.
    selected_contexts={(row['transcript'],row['events']) for row in rows if row['mhc_class']=='I'}
    contexts={k:v for k,v in contexts.items() if k in selected_contexts}
    # De-duplicate identical protein sequences before prediction.
    unique={}; sequences={}
    for key,(wt,mut) in contexts.items():
        for seq in (wt,mut):
            if not seq or any(x not in 'ACDEFGHIKLMNPQRSTVWY' for x in seq):continue
            digest=hashlib.sha256(seq.encode()).hexdigest()
            unique[digest]=seq
            sequences[seq]=digest
    if batch_sequences<1: raise ValueError('pepsickle batch size must be positive')
    scores=collections.defaultdict(dict)
    with tempfile.TemporaryDirectory(prefix='varpep_cleavage_') as tmp:
        entries=list(unique.items())
        total=(len(entries)+batch_sequences-1)//batch_sequences
        started=time.monotonic()
        for batch_i in range(total):
            chunk=entries[batch_i*batch_sequences:(batch_i+1)*batch_sequences]
            fasta=os.path.join(tmp,'proteins.fasta')
            result=os.path.join(tmp,'cleavage.tsv')
            with open(fasta,'w') as h:
                for digest,seq in chunk:h.write('>'+digest+'\n'+seq+'\n')
            cmd=[pepsickle_bin,'-f',fasta,'-o',result,'-m',model]
            if model!='epitope':cmd+=['-p',proteasome]
            print(f'[cleavage] batch {batch_i+1}/{total}; proteins={len(chunk)}; elapsed={time.monotonic()-started:.0f}s',file=sys.stderr,flush=True)
            subprocess.run(cmd,check=True)
            with open(result,newline='') as h:
                rd=csv.DictReader(h,delimiter='\t')
                expected={'position','cleav_prob','protein_id'}
                if not rd.fieldnames or not expected.issubset(set(rd.fieldnames)):
                    raise ValueError('Unexpected pepsickle output header: '+str(rd.fieldnames))
                for r in rd:
                    try:scores[r['protein_id']][int(r['position'])]=float(r['cleav_prob'])
                    except (ValueError,TypeError):continue
        print(f'[cleavage] completed {len(entries)} protein sequences in {time.monotonic()-started:.0f}s',file=sys.stderr,flush=True)
    outcols=columns+['proteasome_model','c_terminal_cleavage_prob_mut',
                     'c_terminal_cleavage_prob_wt','delta_c_terminal_cleavage_prob',
                     'processing_pass','processing_note']
    out=prefix+'.peptides_processed.tsv'
    with open(out+'.tmp','w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=outcols,delimiter='\t');w.writeheader()
        for row in rows:
            if row['mhc_class']!='I':continue
            key=(row['transcript'],row['events'])
            if key not in contexts:raise RuntimeError('Missing full protein for '+str(key))
            wt,mut=contexts[key]
            start=int(row['aa_start1']);length=int(row['length']);end=start+length-1
            if mut[start-1:end]!=row['peptide']:
                raise RuntimeError('Peptide sequence mismatch for '+str(key))
            mutscore=scores.get(sequences.get(mut,''),{}).get(end)
            # WT comparison is meaningful only when the same peptide boundary
            # maps to the same residue coordinate (no upstream indel shift).
            events=row['events'].split(';')
            indel=any(':INS:' in ev or ':DEL:' in ev for ev in events)
            wtscore=None if indel else scores.get(sequences.get(wt,''),{}).get(end)
            row.update({'proteasome_model':model+(':'+proteasome if model!='epitope' else ''),
                        'c_terminal_cleavage_prob_mut': '' if mutscore is None else f'{mutscore:.6f}',
                        'c_terminal_cleavage_prob_wt': '' if wtscore is None else f'{wtscore:.6f}',
                        'delta_c_terminal_cleavage_prob': '' if mutscore is None or wtscore is None else f'{mutscore-wtscore:.6f}',
                        'processing_pass': '' if mutscore is None else ('yes' if mutscore>=threshold else 'no'),
                        'processing_note': 'WT_boundary_not_comparable_indel' if indel else ('score_unavailable' if mutscore is None else 'C_terminal_only_not_HLA_presentation')})
            w.writerow(row)
    os.replace(out+'.tmp',out)
    print(f'[cleavage] candidates={len(rows):,}; output={out}; model={model}',file=sys.stderr,flush=True)


def main():
    p=argparse.ArgumentParser(description='VarPep v1.7: complete variant/haplotype catalogs with optional peptide-only filtering')
    sp=p.add_subparsers(dest='command',required=True)
    b=sp.add_parser('build-index')
    for x in ('gtf','transcripts','proteins','index'): b.add_argument('--'+x,required=True)
    r=sp.add_parser('build-rediportal-index')
    r.add_argument('--rediportal',required=True);r.add_argument('--out',required=True)
    sn=sp.add_parser('build-snp151-index',help='Build reusable SNP151 position/rsID index from UCSC snp151.txt[.gz]')
    sn.add_argument('--snp151',required=True);sn.add_argument('--out',required=True)
    c=sp.add_parser('call',help='One command: variant calling, periodic TSV, VAF, filtering and REDIportal annotation')
    for x in ('bam','reference','transcripts','proteins','index','out-prefix'):c.add_argument('--'+x,required=True)
    c.add_argument('--name-sorted',action='store_true',help='Required for fragment-level VAF')
    c.add_argument('--include-indels',action='store_true')
    c.add_argument('--nonsynonymous-only',action='store_true',help='Restrict peptide generation to individually protein-changing events; preserve complete variant and haplotype catalogs')
    c.add_argument('--min-mapq',type=int,default=255)
    c.add_argument('--min-baseq',type=int,default=30)
    c.add_argument('--threads',type=int,default=8)
    c.add_argument('--batch-size',type=int,default=2000)
    c.add_argument('--max-inflight-batches',type=int,default=64,help='Maximum queued/running/unmerged batches; must be >= threads')
    c.add_argument('--commit-every',type=int,default=8,help='Commit SQLite after this many completed batches')
    c.add_argument('--flush-seconds',type=float,default=60)
    c.add_argument('--min-support',type=int,default=1)
    c.add_argument('--min-vaf',type=float,default=0.0)
    c.add_argument('--rediportal-index')
    c.add_argument('--snp151-index',help='SQLite index created by build-snp151-index')
    c.add_argument('--exclude-known',action='store_true',help='Restrict peptide generation to events absent from REDIportal and SNP151; preserve complete catalogs')
    c.add_argument('--report-every',type=int,default=100000)
    c.add_argument('--pepsickle-bin',default='pepsickle',help='Executable for trained pepsickle cleavage predictor')
    c.add_argument('--proteasome-model',choices=['epitope','in-vitro','in-vitro-2'],default='epitope')
    c.add_argument('--proteasome-type',choices=['C','I'],default='I')
    c.add_argument('--cleavage-threshold',type=float,default=0.5)
    c.add_argument('--pepsickle-batch-sequences',type=int,default=250,help='Number of protein sequences per Pepsickle invocation; progress logged per batch')
    c.add_argument('--region',help='Optional genomic interval chr:start-end, 1-based inclusive (e.g. chr7:55000000-56000000)')
    args=p.parse_args()
    if args.command=='build-index':build(args);return
    if args.command=='build-rediportal-index':build_rediportal_index(args);return
    if args.command=='build-snp151-index':build_snp151_index(args);return
    if args.exclude_known and not (args.rediportal_index and args.snp151_index):
        p.error('--exclude-known requires BOTH --rediportal-index and --snp151-index')
    if not args.name_sorted: p.error('call requires --name-sorted for paired-end fragment-level VAF')
    if not 0<=args.min_vaf<=1:p.error('--min-vaf must be in [0,1]')
    if args.min_support<1:p.error('--min-support must be >=1')
    if shutil.which(args.pepsickle_bin) is None and not os.path.isfile(args.pepsickle_bin):
        p.error('Pepsickle executable not found; set --pepsickle-bin before starting calling')
    args.indexed_coordinate=False
    args.region_parsed=parse_region(args.region)
    print('[phase 1/2] Variant/haplotype/peptide calling with periodic TSV snapshots',file=sys.stderr,flush=True)
    call(args)
    state=args.out_prefix+'.counts.sqlite'
    with sqlite3.connect(state) as db:
        print('[phase 1b] Deferred translation of globally supported haplotypes',file=sys.stderr,flush=True)
        finalize_haplotypes(db,args)
    print('[phase 2/2] Fragment coverage and per-variant VAF from same BAM',file=sys.stderr,flush=True)
    depth,fragments=coverage_for_sites(args.bam,state,args.min_mapq,args.min_baseq,args.report_every,args.region_parsed)
    export_vaf(state,args.out_prefix,depth,args.min_support,args.min_vaf,args.rediportal_index,args.snp151_index)
    process_peptides(state,args.out_prefix,args.pepsickle_bin,args.proteasome_model,
                     args.proteasome_type,args.cleavage_threshold,args.pepsickle_batch_sequences)
    print(json.dumps({'status':'complete','fragments_vaf_pass':fragments,'min_support':args.min_support,'min_vaf':args.min_vaf,'prefix':args.out_prefix},indent=2),flush=True)

if __name__=='__main__':main()
