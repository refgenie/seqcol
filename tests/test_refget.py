import pytest
import os
from refget import RefGetHenge

DEMO_FILES = ["demo.fa", "demo2.fa", "demo3.fa", "demo4.fa", "demo5.fa"]


class TestGeneral:
    def test_no_schemas_required(self):
        """
        In contrast to the generic Henge object, RefGetHenge does not
        require schemas as input, they are predefined in the constructor
        """
        assert isinstance(RefGetHenge({}), RefGetHenge)


class TestFastaInserting:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_insert_works(self, fasta_name, fasta_path):
        rgdb = RefGetHenge({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        res = rgdb.load_fasta(f)
        assert len(res) == 2  # returns digest and list of ASDs


class TestRetrieval:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_insert_works(self, fasta_name, fasta_path):
        rgdb = RefGetHenge({})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = rgdb.load_fasta(f)
        assert rgdb.retrieve(d, reclimit=0) == \
               [{k: str(v) for k, v in asd.items()} for asd in asds]

