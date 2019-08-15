""" Computing configuration representation """

import redis
import logging
import pyfaidx

from collections import OrderedDict
from . import trunc512_digest
from . import __version__

_LOGGER = logging.getLogger(__name__)

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

        if (':' not in result):
            return result
        else:
            if isinstance(reclimit, int):
                reclimit = reclimit - 1

            content = OrderedDict()
            for unit in result.split(";"):
                name, seq = unit.split(':')
                if (isinstance(reclimit, int) and reclimit == 0):
                    content[name] = seq
                else:
                    content[name] = self.refget(seq, lookup_table, reclimit)
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


    def fasta_checksum(self, fa_file):
        """
        Just calculate checksum of fasta file without loading it.
        """
        fa_object = pyfaidx.Fasta(fa_file)
        content_checksums = {}
        for k in fa_object.keys():
            content_checksums[k] = self.checksum_function(str(fa_object[k]))
        collection_string = ";".join([":".join(i) for i in content_checksums.items()])
        collection_checksum = self.load_seq(collection_string)
        return collection_checksum, content_checksums


    def load_fasta(self, fa_file, checksum_function=None):
        """
        Calculates checksums and loads each sequence in a fasta file into the
        database, and loads a level 2 collection checksum representing the
        entire collection into the database.
        """
        if not checksum_function:
            checksum_function = self.checksum_function        
        fa_object = pyfaidx.Fasta(fa_file)
        content_checksums = {}
        for k in fa_object.keys():
            content_checksums[k] = self.load_seq(str(fa_object[k]))
        collection_string = ";".join([":".join(i) for i in content_checksums.items()])
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
