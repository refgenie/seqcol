import henge
import logging
import logmuse
import os
import pyfaidx
import refget

from copy import copy
from functools import reduce
from itertools import compress

from .hash_functions import trunc512_digest
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
                f"Invalid topology ({topology_default}). "
                f"Choose from: {','.join(KNOWN_TOPOS)}"
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
        return self.compat_all(A, B)

    @staticmethod
    def compat(A, B):
        """
        New compatibility function for array-based data model.
        """

        lenA = len(A)
        lenB = len(B)
        dupeA = lenA - len(dict.fromkeys(A))
        dupeB = lenB - len(dict.fromkeys(B))
        ainb = [x in B for x in A]
        bina = [x in A for x in B]
        sum_ainb = sum(ainb)
        if sum_ainb > 1:
            order = list(compress(B, bina)) == list(compress(A, ainb))
        else:
            order = False

        result = {
            "a-in-b": sum_ainb,
            "b-in-a":  sum(bina),
            "a-total": lenA,
            "b-total": lenB,
            "a-duplicated": dupeA,
            "b-duplicated": dupeB,
            "order-match": order
        }
        return result

    @staticmethod
    def compatoverlap(A, B):
        A_filtered = list(filter(lambda x: x in B, A))
        B_filtered = list(filter(lambda x: x in A, B))
        A_count = len(A_filtered)
        B_count = len(B_filtered)
        overlap = min(len(A_filtered), len(B_filtered)) ## counts duplicates

        if A_count + B_count < 1:
            # order match requires at least 2 matching elements
            order = None
        elif not (A_count == B_count == overlap):
            # duplicated matches means order match is undefined
            order = None
        else:
            order = (A_filtered == B_filtered)
        return { "overlap": overlap, "order-match": order }

    @staticmethod
    def compat_all(A, B):
        all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
        result = {}
        # flipped_format = {
        #     "overlap": {},
        #     "a-total": len(A["lengths"]),
        #     "b-total": len(B["lengths"]),
        #     "order-match": {},
        #     "only-in-a": [],
        #     "only-in-b": [],
        # }
        new_format = {
            "arrays": {
                "a-only": [],
                "b-only": [],
                "a-and-b": []
            },
            "elements": {
                "total": {
                    "a": len(A["lengths"]),
                    "b": len(B["lengths"])
                },
                "overlap": {},
                "order-match": {},
            }
        }

        for k in all_keys:
            _LOGGER.info(k)
            if k not in A:
                result[k] = {"flag": -1}
                new_format["arrays"]["b-only"].append(k)
            elif k not in B:
                new_format["arrays"]["a-only"].append(k)
            else:
                new_format["arrays"]["a-and-b"].append(k)
                res = SeqColClient.compatoverlap(A[k], B[k])
                new_format["elements"]["overlap"][k] = res["overlap"]
                new_format["elements"]["order-match"][k] = res["order-match"]
        return new_format

    @staticmethod
    def compat_all_old(A, B):
        all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
        result = {}
        flipped_format = {
            "a-in-b": {},
            "b-in-a": {},
            "a-total": {},
            "b-total": {},
            "a-duplicated": {},
            "b-duplicated": {},
            "order-match": [],
            "only-in-a": [],
            "only-in-b": [],
        }
        for k in all_keys:
            _LOGGER.info(k)
            if k not in A:
                result[k] = {"flag": -1}
                flipped_format["only-in-b"].append(k)
            elif k not in B:
                flipped_format["only-in-a"].append(k)
            else:
                v = SeqColClient.compat(A[k], B[k])
                result[k] = v
                if "a-in-b" in v:
                    flipped_format["a-in-b"][k] = v['a-in-b']
                if "b-in-a":
                    flipped_format["b-in-a"][k] = v['b-in-a']
                if "a-total" in v:
                    flipped_format["a-total"][k] = v['a-total']
                if "b-total" in v:
                    flipped_format["b-total"][k] = v['b-total']
                if "a-duplicated" in v:
                    flipped_format["a-duplicated"][k] = v['a-duplicated']
                if "b-duplicated" in v:
                    flipped_format["b-duplicated"][k] = v['b-duplicated']
                if "order-match" in v:
                    flipped_format["order-match"].append(k)

        # result = {
        #     "any-elements-shared": any(ainb),
        #     "all-a-in-b": all(ainb),
        #     "all-b-in-a": all(bina),
        #     "order-match": order,
        #     "flag": flag
        # }

        return flipped_format

    @staticmethod
    def compare_asds(asdA, asdB, explain=False):
        """
        Compare Annotated Sequence Digests (ASDs) -- digested sequences and `data

        :param str asdA: ASD for first sequence collection to compare.
        :param str asdB: ASD for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """

        def _xp(prop, lst):
            """Extract property from a list of dicts"""
            return list(map(lambda x: x[prop], lst))

        def _index(x, lst):
            """Find an index of a sequence element in a list of dicts"""
            try:
                return _xp(SEQ_KEY, lst).index(x)
            except:
                return None

        def _get_common_content(lstA, lstB):
            """
            Find the intersection between two list of dicts with sequences
            """
            return list(
                filter(None.__ne__, [_index(x, lstB) for x in _xp(SEQ_KEY, lstA)])
            )

        # Not ideal, but we expect these to return lists, but if the item was
        # singular only a dict is returned
        if not isinstance(asdA, list):
            asdA = [asdA]
        if not isinstance(asdB, list):
            asdB = [asdB]

        ainb = [x in _xp(SEQ_KEY, asdB) for x in _xp(SEQ_KEY, asdA)]
        bina = [x in _xp(SEQ_KEY, asdA) for x in _xp(SEQ_KEY, asdB)]

        return_flag = 0  # initialize
        if sum(ainb) > 1:
            ordA = _get_common_content(asdA, asdB)
            if ordA == sorted(ordA):
                return_flag += CONTENT_A_ORDER
        if sum(bina) > 1:
            ordB = _get_common_content(asdB, asdA)
            if ordB == sorted(ordB):
                return_flag += CONTENT_B_ORDER

        ainb_len = [x in _xp(LEN_KEY, asdB) for x in _xp(LEN_KEY, asdA)]
        bina_len = [x in _xp(LEN_KEY, asdA) for x in _xp(LEN_KEY, asdB)]

        ainb_name = [x in _xp(NAME_KEY, asdB) for x in _xp(NAME_KEY, asdA)]
        bina_name = [x in _xp(NAME_KEY, asdA) for x in _xp(NAME_KEY, asdB)]

        ainb_topo = [x in _xp(TOPO_KEY, asdB) for x in _xp(TOPO_KEY, asdA)]
        bina_topo = [x in _xp(TOPO_KEY, asdA) for x in _xp(TOPO_KEY, asdB)]

        if all(ainb):
            return_flag += CONTENT_ALL_A_IN_B
        if all(bina):
            return_flag += CONTENT_ALL_B_IN_A

        if all(ainb_name):
            return_flag += NAMES_ALL_A_IN_B
        if all(bina_name):
            return_flag += NAMES_ALL_B_IN_A

        if all(ainb_topo):
            return_flag += TOPO_ALL_A_IN_B
        if all(bina_topo):
            return_flag += TOPO_ALL_B_IN_A

        if all(ainb_len):
            return_flag += LENGTHS_ALL_A_IN_B
        if all(bina_len):
            return_flag += LENGTHS_ALL_B_IN_A

        if explain:
            explain_flag(return_flag)
        return return_flag

    def compare(self, digestA, digestB, explain=False):
        """
        Given two collection checksums in the database, provide some information
        about how they are related.

        :param str digestA: Digest for first sequence collection to compare.
        :param str digestB: Digest for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """
        typeA = self.database[digestA + henge.ITEM_TYPE]
        typeB = self.database[digestB + henge.ITEM_TYPE]

        if typeA != typeB:
            _LOGGER.error(
                f"Can't compare objects of different types: " f"{typeA} vs {typeB}"
            )

        asdA = self.retrieve(digestA, reclimit=1)
        asdB = self.retrieve(digestB, reclimit=1)
        return self.compare_asds(asdA, asdB, explain=explain)

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
        SCAS = fasta_to_scas(fa_object)
        digest = self.insert(SCAS, "SeqColArraySet", reclimit=1)
        return {
          "fa_file": filepath,
          "fa_object": fa_object,
          "SCAS": SCAS,
          "digest": digest
        }


# Static functions below (these don't require a database)


def explain_flag(flag):
    """Explains a compare flag"""
    print(f"Flag: {flag}\nBinary: {bin(flag)}\n")
    for e in range(0, 13):
        if flag & 2 ** e:
            print(FLAGS[2 ** e])


def parse_fasta(fa_file):
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

        with gzopen(fa_file, "rt") as f_in, NamedTemporaryFile(
            mode="w+t", suffix=".fa"
        ) as f_out:
            f_out.writelines(f_in.read())
            f_out.seek(0)
            return pyfaidx.Fasta(f_out.name)

# static 
def fasta_to_scas(fa_object, verbose=True):
    """
    Given a fasta object, return a SCAS (Sequence Colletion Array Set)
    """
    # SCAS = SeqColArraySet
    # Or maybe should be "Level 1 SC"
    SCAS = {"lengths": [] , "names": [], "sequences": []}
    seqs = fa_object.keys()
    nseqs = len(seqs)
    print(f"Found {nseqs} chromosomes")
    i=1
    for k in fa_object.keys():
        if verbose:
            print(f"Processing ({i} of {nseqs}) {k}...")
        seq = str(fa_object[k])
        digest = henge.md5(seq.upper())
        SCAS["lengths"].append(str(len(seq)))
        SCAS["names"].append(fa_object[k].name)
        SCAS["sequences"].append(digest)
        i += 1
    return SCAS
