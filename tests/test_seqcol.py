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
]

NAMES_LENGTHS_TESTS = [
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
    def test_fasta_loading_works(self, fasta_name, fa_root):
        scc = SeqColClient(database={})
        f = os.path.join(fa_root, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        res = scc.load_fasta(f)
        assert len(res) == 2  # returns digest and list of AnnotatedSequencesList


class TestRetrieval:
    @pytest.mark.parametrize("fasta_name", DEMO_FILES)
    def test_retrieval_works(self, fasta_name, fa_root):
        scc = SeqColClient(database={})
        f = os.path.join(fa_root, fasta_name)
        print("Fasta file to be loaded: {}".format(f))
        d, asds = scc.load_fasta(f)
        # convert integers in the dicts to strings
        # {k: str(v) if isinstance(v, int) else v for k, v in asd.items()}
        lst = [
            {k:v for k, v in asd.items()}
            for asd in asds
        ]
        assert scc.retrieve(d) == lst


def check_comparison(fasta1, fasta2, expected_comparison):
    print(f"Comparison: Fasta1: {fasta1} vs Fasta2: {fasta2}. Expected: {expected_comparison}")
    scc = SeqColClient(database={})
    d = scc.load_fasta_from_filepath(fasta1)
    d2 = scc.load_fasta_from_filepath(fasta2)
    with open(expected_comparison) as fp:
        correct_compare_response = json.load(fp)
        proposed_compare_response = scc.compare(d["SCAS"], d2["SCAS"])
        print(json.dumps(proposed_compare_response, separators=(",", ":"), ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2))
        assert proposed_compare_response == correct_compare_response

class TestCompare:
    """
    Test the compare function, using demo fasta files, and pre-computed
    compare function results stored as answer files.
    """
    @pytest.mark.parametrize(["fasta1", "fasta2", "answer_file"], COMPARE_TESTS)
    def test_fasta_compare(self, fasta1, fasta2, answer_file, fa_root):
        check_comparison(os.path.join(fa_root, fasta1), os.path.join(fa_root, fasta2), answer_file)


    @pytest.mark.parametrize(["fasta1", "fasta2", "answer_file"], NAMES_LENGTHS_TESTS)
    def test_names_lengths_order(self, fasta1, fasta2, answer_file, fa_root):
        """ Does the names_lengths array correctly identify order variants """
        check_comparison(os.path.join(fa_root, fasta1), os.path.join(fa_root, fasta2), answer_file)

