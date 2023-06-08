import json
import pytest
from seqcol import SeqColClient
from seqcol.const import *

DEMO_FILES = [
    "demo0.fa",
    "demo1.fa.gz",
    "demo2.fa",
    "demo3.fa",
    "demo4.fa",
    "demo5.fa.gz",
    "demo6.fa",
]

# Pairs of files to compare, with the "correct" compare response
COMPARE_TESTS = [
    (DEMO_FILES[1], DEMO_FILES[1], "demo_fasta/compare-1vs1.json"),
    (DEMO_FILES[0], DEMO_FILES[1], "demo_fasta/compare-0vs1.json"),
    (DEMO_FILES[5], DEMO_FILES[6], "demo_fasta/compare-5vs6.json"),
]

class TestGeneral:
    def test_no_schemas_required(self):
        """
        In contrast to the generic Henge object, SeqColClient does not
        require schemas as input, they are predefined in the constructor
        """
        assert isinstance(SeqColClient(database={}), SeqColClient)


class TestFastaInserting:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_fasta_loading_works(self, fasta_name, fasta_path):
        scc = SeqColClient(database={})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        res = scc.load_fasta(f)
        assert len(res) == 2  # returns digest and list of AnnotatedSequencesList


class TestRetrieval:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_retrieval_works(self, fasta_name, fasta_path):
        scc = SeqColClient(database={})
        f = os.path.join(fasta_path, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = scc.load_fasta(f)
        # convert integers in the dicts to strings
        # {k: str(v) if isinstance(v, int) else v for k, v in asd.items()}
        lst = [
            {k:v for k, v in asd.items()}
            for asd in asds
        ]
        assert scc.retrieve(d) == lst


class TestCompare:
    """
    Test the compare function, using demo fasta files, and pre-computed
    compare function results stored as answer files.
    """
    @pytest.mark.parametrize(["fasta1", "fasta2", "answer_file"], COMPARE_TESTS)
    def test_fasta_compare(self, fasta1, fasta2, answer_file, fasta_path):
        print(f"Fasta1: {fasta1}")
        print(f"Fasta2: {fasta2}")
        print(f"answer_file: {answer_file}")
        scc = SeqColClient(database={})
        d = scc.load_fasta_from_filepath(os.path.join(fasta_path, fasta1))
        d2 = scc.load_fasta_from_filepath(os.path.join(fasta_path, fasta2))
        with open(answer_file) as fp:
            correct_compare_response = json.load(fp)
            proposed_compare_response = scc.compare(d["SCAS"], d2["SCAS"])
            print(json.dumps(proposed_compare_response, separators=(",", ":"), ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2))
            assert proposed_compare_response == correct_compare_response
