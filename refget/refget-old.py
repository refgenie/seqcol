""" Computing configuration representation """

import redis
import logging
import pyfaidx
import pymongo
pymongo.Connection = lambda host, port, **kwargs: pymongo.MongoClient(host=host, port=port)
import mongodict


from collections import OrderedDict
from . import trunc512_digest
from . import __version__

_LOGGER = logging.getLogger(__name__)

DELIM_LVL1 = u"\u00B7"
DELIM_LVL2 = u"\u2016"
DELIM_LVL3 = u"\u2980"

class RedisDict(redis.Redis):
    """
    Dict-like interface to a redis back-end
    """
    def __delitem__(self, key):
        self.set(key, None)

    def __getitem__(self, key):
        result = self.get(key)
        if result:
            return result.decode()
        else:
            raise KeyError(key)
    
    def __setitem__(self, key, value):
        self.set(key, value)

class MongoDict(mongodict.MongoDict):
    pass


class RefDB(object):

    def __init__(self, database, checksum_function=trunc512_digest):
        """
        :@param database: Dict-like lookup database with sequences and hashes.
        """
        self.database = database
        self.checksum_function = checksum_function


    def refget(self, checksum, lookup_table=None, reclimit=None):
        """
        Recursive refget lookup implementation

        :param str checksum: 
        :param int reclimit: Limit the number of recursion layers in the refget lookup
        :return str: Sequence corresponding to given checksum
        """
        if not lookup_table:
            lookup_table = self.database

        try:
            result = lookup_table[checksum]
        except KeyError:
            return "Not found"

        if (DELIM_LVL2 not in result):
            # Base case
            return result
        else:
            # Recursive case
            if isinstance(reclimit, int):
                reclimit = reclimit - 1

            content = OrderedDict()
            for unit in result.split(DELIM_LVL3):
                name, laseq = unit.split(DELIM_LVL2)
                try:
                    length, seq = laseq.split(DELIM_LVL1)
                except ValueError:
                    length = None
                    seq = laseq
                if (isinstance(reclimit, int) and reclimit == 0):
                    content[name] = {
                        'length': length,
                        'seq': seq,
                        }
                else:
                    content[name] = {
                        'length': length,
                        'seq': self.refget(seq, lookup_table, reclimit)
                        }
            return content

    def fasta_fmt(self, content):
        """
        Given a content dict return by refget for a sequence collection,
        convert it to a string that can be printed as a fasta file.
        """
        return "\n".join(["\n".join(
            [">" + name, seq]) for name, seq in content.items()])


    def load_seq(self, seq, checksum_function=None):
        if not checksum_function:
            checksum_function = self.checksum_function
        checksum = checksum_function(seq)
        self.database[checksum] = seq
        _LOGGER.info("Loaded {}".format(checksum))

        return checksum


    def load_fasta(self, fa_file, checksum_function=None):
        """
        Calculates checksums and loads each sequence in a fasta file into the
        database, and loads a level 2 collection checksum representing the
        entire collection into the database.
        """
        if not checksum_function:
            checksum_function = self.checksum_function        
        fa_object = parse_fasta(fa_file)
        content_checksums = {}
        for k in fa_object.keys():
            seq = str(fa_object[k])
            content_checksums[k] = {'length': len(seq), 'seq': self.load_seq(seq)}
        # Produce a length-annotated seq
        #collection_string = ";".join([":".join(i) for i in content_checksums.items()])
        collection_string = DELIM_LVL3.join([("{}" + DELIM_LVL2 + "{}" + 
                                DELIM_LVL1 + "{}").format(name, value["length"], value["seq"]) 
                                for name, value in content_checksums.items()])
        _LOGGER.info("collection_string: {}".format(collection_string))
        collection_checksum = self.load_seq(collection_string)
       
        return collection_checksum, content_checksums


    def compare(self, checksumA, checksumB):
        """
        Given two collection checksums in the database, provide some information
        about how they are related.

        """
        contents1 = self.refget(checksumA, reclimit=1)
        contents2 = self.refget(checksumB, reclimit=1)

        ainb = [x in contents2.values() for x in contents1.values()]
        bina = [x in contents1.values() for x in contents2.values()]

        ainb_length = [x['length'] for x in contents1.values()]
        # [[name, val['length']] for name, val in content1.items()]


        if all(ainb):
            if all(bina):
                names_check = [x in contents2.keys() for x in contents1.keys()]
                if all(names_check):
                    res = "Sequence-level identical, order mismatch"
                else:
                    res = "Sequence-level identical, names mismatch"
            else:
                res = "A is a sequence-level subset of B"
        elif any(ainb):
            if all(bina):
                res = "B is a sequence-level subset of A"
            else:
                res = "A and B share some sequences"
        else:
            res = "No sequences shared"



        return res


# Static functions below (these don't require a database)

def parse_fasta(fa_file):
    _LOGGER.info("Hashing {}".format(fa_file))
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

def fasta_checksum(fa_file):
    """
    Just calculate checksum of fasta file without loading it.
    """
    fa_object = parse_fasta(fa_file)
    content_checksums = {}
    for k in fa_object.keys():
        content_checksums[k] = self.checksum_function(str(fa_object[k]))
    collection_string = ";".join([":".join(i) for i in content_checksums.items()])
    collection_checksum = self.load_seq(collection_string)
    return collection_checksum, content_checksums