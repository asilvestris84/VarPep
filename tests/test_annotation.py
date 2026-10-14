import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import varpep

class AnnotationTests(unittest.TestCase):
    def test_rediportal_strand_and_position(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'redi.sqlite'
            with sqlite3.connect(p) as con:
                con.execute('CREATE TABLE sites(chrom TEXT,pos INTEGER,strand TEXT)')
                con.execute("INSERT INTO sites VALUES ('chr12',14882147,'-')")
            a=varpep.EditingAnnotator(str(p))
            self.assertEqual(a.variant_annotation('chr12',14882147,'SNV','T','C'),['yes','transcript_strand_unresolved','yes','-'])
            self.assertEqual(a.variant_annotation('chr12',14882147,'SNV','T','G'),['no','transcript_strand_unresolved','yes','-'])
            self.assertEqual(a.annotate('chr12',14882147,'SNV','T','C','-')[0],'yes')
            self.assertEqual(a.annotate('chr12',14882147,'SNV','T','C',None)[0],'no')
            a.close()

    def test_snp151_one_based(self):
        with tempfile.TemporaryDirectory() as d:
            t=Path(d)/'snp151.txt';p=Path(d)/'snp.sqlite'
            fields=['585','chr12','14882146','14882147','rs123','0','+','T','T','T/C','genomic','single']
            t.write_text('\t'.join(fields)+'\n')
            varpep.build_snp151_index(SimpleNamespace(snp151=str(t),out=str(p)))
            a=varpep.KnownSNPAnnotator(str(p))
            self.assertEqual(a.annotate('chr12',14882147,'SNV'),('yes','rs123'))
            self.assertEqual(a.annotate('chr12',14882146,'SNV'),('no',''))
            a.close()

if __name__=='__main__':unittest.main()
