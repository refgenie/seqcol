from copy import copy
import henge
import logmuse
import os

import pyfaidx
import logging
_LOGGER = logging.getLogger(__name__)
henge.ITEM_TYPE = "_item_type"

from .const import *

SCHEMA_FILEPATH = os.path.join(
        os.path.dirname(__file__),
        "schemas")


class RefGetHenge(henge.Henge):
    """
    Extension of henge that accommodates refget sequences.
    """

    def __init__(self, database, schemas=None, henges=None, checksum_function=henge.md5):
        """
        A user interface to insert and retrieve decomposable recursive unique
        identifiers (DRUIDs).

        :param dict database: Dict-like lookup database with sequences and hashes.
        :param dict schemas: One or more jsonschema schemas describing the
            data types stored by this Henge
        :param function(str) -> str checksum_function: Default function to handle the digest of the
            serialized items stored in this henge.
        """
        def _load_schema(name):
            return henge.load_yaml(os.path.join(SCHEMA_FILEPATH, name))

        # These are the item types that this henge can understand.
        if not schemas:
            schemas = {
                "sequence": _load_schema("sequence.yaml"),
                "ASD": _load_schema("annotated_sequence_digest.yaml"),
                "ASDList": _load_schema("ASDList.yaml"),
                "ACDList": _load_schema("ACDList.yaml"),
                "ACD": _load_schema("annotated_collection_digest.yaml")
            }
        super(RefGetHenge, self).__init__(database, schemas, henges=henges,
                                          checksum_function=checksum_function)

    def refget(self, digest, reclimit=None, postprocess=None):
        item_type = self.database[digest + henge.ITEM_TYPE]
        full_data = self.retrieve(digest, reclimit=reclimit)
        if not postprocess:
            return full_data

        if postprocess == "simplify":
            if item_type == "sequence":
                return full_data['sequence']
            elif item_type == "asd":
                asdlist = {}
                for x in full_data:
                    asdlist[x['name']] = x['sequence_digest']['sequence']
                return asdlist
        elif postprocess == "fasta":
            if item_type == "sequence":
                raise Exception("can't postprocess a sequence into fasta")
            elif item_type == "asd":
                asdlist = {}
                for x in full_data:
                    asdlist[x['name']] = x['sequence_digest']['sequence']                
                return self.fasta_fmt(asdlist)
        else:
            raise NotImplementedError(
                "This postprocessing mode is not implemented")

    def fasta_fmt(self, content):
        """
        Given a content dict return by refget for a sequence collection,
        convert it to a string that can be printed as a fasta file.
        """
        return "\n".join(
            ["\n".join([">" + x["name"], x["sequence_digest"]["sequence"]])
             for x in content])

    def load_seq(self, seq):
        checksum = self.insert({'sequence': seq}, "sequence")
        _LOGGER.debug("Loaded {}".format(checksum))
        return checksum

    def load_fasta(self, fa_file, lengths_only=False):
        """
        Calculates checksums and loads each sequence in a fasta file into the
        database, and loads a level 2 collection checksum representing the
        entire collection into the database.
        """
        fa_object = parse_fasta(fa_file)
        asdlist = []
        for k in fa_object.keys():
            seq = str(fa_object[k])
            if lengths_only:
                seq_digest = ""
            else:
                seq_digest = self.load_seq(seq)
            asdlist.append({'name': k,
                          'length': len(seq), 
                          'topology': 'linear',
                          'sequence_digest': seq_digest})

        _LOGGER.debug(asdlist)
        collection_checksum = self.insert(asdlist, 'ASDList')
        return collection_checksum, asdlist

    def load_seqset(self, seqset):
        """
        Convert a 'seqset', which is a dict with names as sequence names and
        values as sequences, into the 'asdlist' required for henge insert.
        """
        seqset_new = copy(seqset)
        for k, v in seqset.items():
            if isinstance(v, str):
                seq = v
                v = {'sequence': seq}
            if 'length' not in v.keys():
                if 'sequence' not in v.keys():
                    _LOGGER.warning(
                        "Each sequence must have either length or a sequence.")
                else:
                    v['length'] = len(v['sequence'])
            if 'sequence' in v.keys():
                v['sequence_digest'] = self.load_seq(seq)
                del v['sequence']
            if 'name' not in v.keys():
                v['name'] = k
            if 'toplogy' not in v.keys():
                v['toplogy'] = 'linear'

            seqset_new[k] = v

        collection_checksum = self.insert(list(seqset_new.values()), 'ASDList')
        return collection_checksum, seqset_new

    def compare_asds(self, asdA, asdB, explain=False):
        """
        Compare Annotated Sequence Digests (ASDs) -- digested sequences and metadata

        :param str asdA: ASD for first sequence collection to compare.
        :param str asdB: ASD for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """
        # Not ideal, but we expect these to return lists, but if the item was
        # singular only a dict is returned
        if not isinstance(asdA, list):
            asdA = [asdA]
        if not isinstance(asdB, list):
            asdB = [asdB]

        def xp(prop, lst):
            """ Extract property """
            return list(map(lambda x: x[prop], lst))

        ainb = [x in xp('sequence_digest', asdB) for x in
                xp('sequence_digest', asdA)]
        bina = [x in xp('sequence_digest', asdA) for x in
                xp('sequence_digest', asdB)]

        def index(x, lst):
            try:
                return xp('sequence_digest', lst).index(x)
            except:
                return None

        return_flag = 0  # initialize
        if sum(ainb) > 1:
            ordA = list(filter(None.__ne__, [index(x, asdB) for x in
                                             xp('sequence_digest', asdA)]))
            if ordA == sorted(ordA):
                return_flag += CONTENT_A_ORDER
        if sum(bina) > 1:
            ordB = list(filter(None.__ne__, [index(x, asdA) for x in
                                             xp('sequence_digest', asdB)]))
            if ordB == sorted(ordB):
                return_flag += CONTENT_B_ORDER

        ainb_len = [x in xp('length', asdB) for x in xp('length', asdA)]
        bina_len = [x in xp('length', asdA) for x in xp('length', asdB)]

        ainb_name = [x in xp('name', asdB) for x in xp('name', asdA)]
        bina_name = [x in xp('name', asdA) for x in xp('name', asdB)]

        if all(ainb):
            return_flag += CONTENT_ALL_A_IN_B

        if all(bina):
            return_flag += CONTENT_ALL_B_IN_A

        if all(ainb_name):
            return_flag += NAMES_ALL_A_IN_B
        if all(bina_name):
            return_flag += NAMES_ALL_B_IN_A

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
            _LOGGER.error("Can't compare objects of different types: {} vs {}".
                          format(typeA, typeB))

        asdA = self.refget(digestA, reclimit=1)
        asdB = self.refget(digestB, reclimit=1)
        return self.compare_asds(asdA, asdB, explain=explain)


# Static functions below (these don't require a database)

def explain_flag(flag):
    """ Explains a compare flag """
    print("Flag: {}\nBinary: {}\n".format(flag, bin(flag)))
    for e in range(0, 13):
        if flag & 2**e:
            print(FLAGS[2**e])


def parse_fasta(fa_file):
    _LOGGER.debug("Hashing {}".format(fa_file))
    try:
        fa_object = pyfaidx.Fasta(fa_file)
    except pyfaidx.UnsupportedCompressionFormat:
        # pyfaidx can handle bgzip but not gzip; so we just hack it here and
        # unzip the file for checksumming, then rezip it for the rest of the
        # asset build.
        # TODO: streamline this to avoid repeated compress/decompress
        os.system("gunzip {}".format(fa_file))
        fa_file_unzipped = fa_file.replace(".gz", "")
        fa_object = pyfaidx.Fasta(fa_file_unzipped)
        os.system("gzip {}".format(fa_file_unzipped))
    return fa_object