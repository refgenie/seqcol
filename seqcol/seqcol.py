import os
import pyfaidx
import logging
from .hash_functions import trunc512_digest

import logmuse
import henge

from .const import *

from copy import copy


_LOGGER = logging.getLogger(__name__)
henge.ITEM_TYPE = "_item_type"


class SeqColClient(henge.Henge):
    """
    Extension of henge that accommodates collections of sequences.
    """

    def __init__(self, database, schemas=None, henges=None,
                 checksum_function=trunc512_digest):
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
            database=database, schemas=schemas or INTERNAL_SCHEMAS,
            henges=henges, checksum_function=checksum_function
        )

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
            raise ValueError(f"Invalid topology ({topology_default}). "
                             f"Choose from: {','.join(KNOWN_TOPOS)}")
        fa_object = parse_fasta(fa_file)
        aslist = []
        for k in fa_object.keys():
            seq = str(fa_object[k])
            aslist.append(
                {NAME_KEY: k, LEN_KEY: len(seq), TOPO_KEY: topology_default,
                 SEQ_KEY: {"" if skip_seq else SEQ_KEY: seq}}
            )
        collection_checksum = self.insert(aslist, ASL_NAME)
        _LOGGER.debug(f"Loaded {ASL_NAME}: {aslist}")
        return collection_checksum, aslist

    @staticmethod
    def compare_asds(asdA, asdB, explain=False):
        """
        Compare Annotated Sequence Digests (ASDs) -- digested sequences and metadata

        :param str asdA: ASD for first sequence collection to compare.
        :param str asdB: ASD for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """

        def _xp(prop, lst):
            """ Extract property from a list of dicts """
            return list(map(lambda x: x[prop], lst))

        def _index(x, lst):
            """ Find an index of a sequence element in a list of dicts """
            try:
                return _xp(SEQ_KEY, lst).index(x)
            except:
                return None

        def _get_common_content(lstA, lstB):
            """
            Find the intersection between two list of dicts with sequences
            """
            return list(filter(None.__ne__,
                               [_index(x, lstB) for x in _xp(SEQ_KEY, lstA)]))

        # Not ideal, but we expect these to return lists, but if the item was
        # singular only a dict is returned
        if not isinstance(asdA, list):
            asdA = [asdA]
        if not isinstance(asdB, list):
            asdB = [asdB]

        ainb = [x in _xp(SEQ_KEY, asdB) for x in
                _xp(SEQ_KEY, asdA)]
        bina = [x in _xp(SEQ_KEY, asdA) for x in
                _xp(SEQ_KEY, asdB)]

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
            _LOGGER.error(f"Can't compare objects of different types: "
                          f"{typeA} vs {typeB}")

        asdA = self.retrieve(digestA, reclimit=1)
        asdB = self.retrieve(digestB, reclimit=1)
        return self.compare_asds(asdA, asdB, explain=explain)


# Static functions below (these don't require a database)

def explain_flag(flag):
    """ Explains a compare flag """
    print(f"Flag: {flag}\nBinary: {bin(flag)}\n")
    for e in range(0, 13):
        if flag & 2**e:
            print(FLAGS[2**e])


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
        with gzopen(fa_file, 'rt') as f_in, \
                NamedTemporaryFile(mode='w+t', suffix=".fa") as f_out:
            f_out.writelines(f_in.read())
            f_out.seek(0)
            return pyfaidx.Fasta(f_out.name)
