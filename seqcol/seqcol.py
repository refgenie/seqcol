import henge
import logging
import os
import pyfaidx
import refget

from copy import copy
from functools import reduce
from itertools import compress

from .utilities import trunc512_digest
from .const import *


_LOGGER = logging.getLogger(__name__)
henge.ITEM_TYPE = "_item_type"


class SeqColClient(refget.RefGetClient):
    """
    Extension of henge that accommodates collections of sequences.
    """

    def __init__(
        self,
        api_url_base=None,
        database={},
        schemas=None,
        henges=None,
        checksum_function=trunc512_digest,
    ):
        """
        A user interface to insert and retrieve decomposable recursive unique
        identifiers (DRUIDs).

        :param dict database: Dict-like lookup database with sequences
            and hashes
        :param dict schemas: One or more jsonschema schemas describing the
            data types stored by this Henge
        :param function(str) -> str checksum_function: Default function to
            handle the digest of the
            serialized items stored in this henge.
        """
        super(SeqColClient, self).__init__(
            api_url_base=api_url_base,
            database=database,
            schemas=schemas or INTERNAL_SCHEMAS,
            henges=henges,
            checksum_function=checksum_function,
        )
        _LOGGER.info("Initializing SeqColClient")

    def load_fasta(self, fa_file, skip_seq=False, topology_default="linear"):
        """
        Load a sequence collection into the database

        :param str fa_file: path to the FASTA file to parse and load
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths and topology
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths and topology
        :param str topology_default: the default topology assigned to
            every sequence
        """
        # TODO: any systematic way infer topology from a FASTA file?
        if topology_default not in KNOWN_TOPOS:
            raise ValueError(
                f"Invalid topology ({topology_default}). " f"Choose from: {','.join(KNOWN_TOPOS)}"
            )
        fa_object = parse_fasta(fa_file)
        aslist = []
        for k in fa_object.keys():
            seq = str(fa_object[k])
            aslist.append(
                {
                    NAME_KEY: k,
                    LEN_KEY: len(seq),
                    TOPO_KEY: topology_default,
                    SEQ_KEY: {"" if skip_seq else SEQ_KEY: seq},
                }
            )
        collection_checksum = self.insert(aslist, ASL_NAME)
        _LOGGER.debug(f"Loaded {ASL_NAME}: {aslist}")
        return collection_checksum, aslist

    def load_fasta2(self, fa_file, skip_seq=False, topology_default="linear"):
        """
        Load a sequence collection into the database

        :param str fa_file: path to the FASTA file to parse and load
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths and topology
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths and topology
        :param str topology_default: the default topology assigned to
            every sequence
        """
        # TODO: any systematic way infer topology from a FASTA file?
        _LOGGER.info("Loading fasta file...")
        fa_object = parse_fasta(fa_file)
        aslist = []
        for k in fa_object.keys():
            seq = str(fa_object[k])
            _LOGGER.info("Loading key: {k} / Length: {l}...".format(k=k, l=len(seq)))
            aslist.append(
                {
                    NAME_KEY: k,
                    LEN_KEY: len(seq),
                    TOPO_KEY: topology_default,
                    SEQ_KEY: "" if skip_seq else seq,
                }
            )
        _LOGGER.info("Inserting into database...")
        collection_checksum = self.insert(aslist, "RawSeqCol")
        _LOGGER.debug(f"Loaded {ASL_NAME}: {aslist}")
        return collection_checksum, aslist

    def compare_digests(self, digestA, digestB):
        A = self.retrieve(digestA, reclimit=1)
        B = self.retrieve(digestB, reclimit=1)
        # _LOGGER.info(A)
        # _LOGGER.info(B)
        return self.compare(A, B)

    @staticmethod
    def compare_elements(A, B):
        """
        Compare two arrays between two arrays
        """

        A_filtered = list(filter(lambda x: x in B, A))
        B_filtered = list(filter(lambda x: x in A, B))
        A_count = len(A_filtered)
        B_count = len(B_filtered)
        overlap = min(len(A_filtered), len(B_filtered))  ## counts duplicates

        if A_count + B_count < 1:
            # order match requires at least 2 matching elements
            order = None
        elif not (A_count == B_count == overlap):
            # duplicated matches means order match is undefined
            order = None
        else:
            order = A_filtered == B_filtered
        return {"a-and-b": overlap, "a-and-b-same-order": order}

    @staticmethod
    def compare(A, B):
        """
        Workhorse comparison function

        @param A Sequence collection A
        @param B Sequence collection B
        @return dict Following formal seqcol specification comparison function return value
        """
        all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
        result = {}
        return_obj = {
            "arrays": {"a-only": [], "b-only": [], "a-and-b": []},
            "elements": {
                "total": {"a": len(A["lengths"]), "b": len(B["lengths"])},
                "a-and-b": {},
                "a-and-b-same-order": {},
            },
        }

        for k in all_keys:
            _LOGGER.info(k)
            if k not in A:
                result[k] = {"flag": -1}
                return_obj["arrays"]["b-only"].append(k)
            elif k not in B:
                return_obj["arrays"]["a-only"].append(k)
            else:
                return_obj["arrays"]["a-and-b"].append(k)
                res = SeqColClient.compare_elements(A[k], B[k])
                return_obj["elements"]["a-and-b"][k] = res["a-and-b"]
                return_obj["elements"]["a-and-b-same-order"][k] = res["a-and-b-same-order"]
        return return_obj

    def retrieve(self, druid, reclimit=None, raw=False):
        try:
            return super(SeqColClient, self).retrieve(druid, reclimit, raw)
        except henge.NotFoundException as e:
            _LOGGER.debug(e)
            try:
                return self.refget(druid)
            except Exception as e:
                _LOGGER.debug(e)
                raise e
                return henge.NotFoundException(
                    "{} not found in database, or in refget.".format(druid)
                )

    def load_fasta_from_refgenie(self, rgc, refgenie_key):
        """
        @param rgc RefGenConf object
        @param refgenie_key key of genome to load
        @param scc SeqColClient object to load into
        """
        filepath = rgc.seek(refgenie_key, "fasta")
        return self.load_fasta_from_filepath(filepath)

    def load_fasta_from_filepath(self, filepath):
        """
        @param filepath Path to fasta file
        @param sc
        """
        fa_object = parse_fasta(filepath)
        SCAS = fasta_obj_to_seqcol(fa_object, digest_function=self.checksum_function)
        digest = self.insert(SCAS, "SeqColArraySet", reclimit=1)
        return {
            "fa_file": filepath,
            "fa_object": fa_object,
            "SCAS": SCAS,
            "digest": digest,
        }

    def load_multiple_fastas(self, fasta_dict):
        """
        Wrapper for load_fasta_from_filepath

        @param fasta_list
        """
        results = {}
        for name, path in fasta_dict.items():
            print(f"Processing fasta '{name}'' at path '{path}'...")
            results[name] = self.load_fasta_from_filepath(path)
        return results


# Static functions below (these don't require a database)


def explain_flag(flag):
    """Explains a compare flag"""
    print(f"Flag: {flag}\nBinary: {bin(flag)}\n")
    for e in range(0, 13):
        if flag & 2**e:
            print(FLAGS[2**e])


def parse_fasta(fa_file) -> pyfaidx.Fasta:
    """
    Read in a gzipped or not gzipped FASTA file
    """
    try:
        return pyfaidx.Fasta(fa_file)
    except pyfaidx.UnsupportedCompressionFormat:
        # pyfaidx can handle bgzip but not gzip; so we just hack it here and
        # gunzip the file into a temporary one and read it in not to interfere
        # with the original one.
        from gzip import open as gzopen
        from shutil import copyfileobj
        from tempfile import NamedTemporaryFile

        with gzopen(fa_file, "rt") as f_in, NamedTemporaryFile(mode="w+t", suffix=".fa") as f_out:
            f_out.writelines(f_in.read())
            f_out.seek(0)
            return pyfaidx.Fasta(f_out.name)


from typing import Callable


def fasta_to_seqcol(fa_file_path: str) -> dict:
    """Given a fasta, return a canonical seqcol object"""
    fa_obj = parse_fasta(fa_file_path)
    return fasta_obj_to_seqcol(fa_obj)


def fasta_obj_to_seqcol(
    fa_object: pyfaidx.Fasta,
    verbose: bool = True,
    digest_function: Callable[[str], str] = henge.md5,
) -> dict:
    """
    Given a fasta object, return a CSC (Canonical Sequence Collection object)
    """
    # CSC = SeqColArraySet
    # Or maybe should be "Level 1 SC"

    CSC = {"lengths": [], "names": [], "sequences": [], "sorted_name_length_pairs": []}
    seqs = fa_object.keys()
    nseqs = len(seqs)
    print(f"Found {nseqs} chromosomes")
    i = 1
    for k in fa_object.keys():
        if verbose:
            print(f"Processing ({i} of {nseqs}) {k}...")
        seq = str(fa_object[k])
        seq_length = len(seq)
        seq_name = fa_object[k].name
        seq_digest = digest_function(seq.upper())
        snlp = {"length": seq_length, "name": seq_name}  # sorted_name_length_pairs
        snlp_digest = digest_function(henge.canonical_str(snlp))
        CSC["lengths"].append(seq_length)
        CSC["names"].append(seq_name)
        CSC["sorted_name_length_pairs"].append(snlp_digest)
        CSC["sequences"].append(seq_digest)
        i += 1
    CSC["sorted_name_length_pairs"].sort()
    return CSC


def build_sorted_name_length_pairs(obj: dict, digest_function):
    """Builds the sorted_name_length_pairs attribute, which corresponds to the coordinate system"""
    sorted_name_length_pairs = []
    for i in range(len(obj["names"])):
        sorted_name_length_pairs.append({"length": obj["lengths"][i], "name": obj["names"][i]})
    nl_digests = []
    for i in range(len(sorted_name_length_pairs)):
        nl_digests.append(digest_function(henge.canonical_str(sorted_name_length_pairs[i])))

    nl_digests.sort()
    return nl_digests
