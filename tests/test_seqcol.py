import pytest
import os
from seqcol import SeqColClient
import tempfile

DEMO_FILES = ["demo.fa", "demo2.fa", "demo3.fa", "demo4.fa", "demo5.fa"]


class TestGeneral:
    def test_no_schemas_required(self):
        """
        In contrast to the generic Henge object, SeqColClient does not
        require schemas as input, they are predefined in the constructor
        """
        assert isinstance(SeqColClient({}), SeqColClient)


class TestFastaInserting:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_insert_works(self, fasta_name, fasta_path):
        scc = SeqColClient({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        res = scc.load_fasta(f)
        assert len(res) == 2  # returns digest and list of ASDs


class TestRetrieval:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_retrieval_works(self, fasta_name, fasta_path):
        scc = SeqColClient({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = scc.load_fasta(f)
        assert scc.retrieve(d, reclimit=0) == \
               [{k: str(v) for k, v in asd.items()} for asd in asds]

    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_refget_works(self, fasta_name, fasta_path):
        scc = SeqColClient({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = scc.load_fasta(f)
        assert scc.refget(d, reclimit=0) == \
               [{k: str(v) for k, v in asd.items()} for asd in asds]


class TestBoth:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_roundtrip(self, fasta_name, fasta_path):
        """
        Verify that FASTA loaded, retrieved and decoded from RefGetHenge
        yields identical collection digest after loading
        """
        rgdb = SeqColClient({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = rgdb.load_fasta(f)
        _, path = tempfile.mkstemp(".fa")
        with open(path, "w") as fh:
            fh.write(rgdb.fasta_fmt(rgdb.refget(d)))
        d1, asds1 = rgdb.load_fasta(path)
        assert d1 == d
        os.remove(path)

