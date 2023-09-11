import henge
import logging
import refget

from copy import copy
from functools import reduce
from itertools import compress

from .const import *
from .seqcol import *
from .utilities import sha512t24u_digest


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
        checksum_function=sha512t24u_digest,
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
        return compare(A, B)

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
