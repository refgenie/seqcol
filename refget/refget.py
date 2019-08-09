""" Computing configuration representation """

import redis
import logging
import pyfaidx

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

        if (':' not in result) or (isinstance(reclimit, int) and reclimit == 1):
            return result
        else:
            if isinstance(reclimit, int):
                reclimit = reclimit - 1
            return "\n".join(["\n".join(
                [">" + seq.split(':')[0], self.refget(
                    seq.split(':')[1], lookup_table)]) for seq in result.split(";")])



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
        fa_object = pyfaidx.Fasta(fa_file)
        content = {}
        for k in fa_object.keys():
            content[k] = self.load_seq(str(fa_object[k]))
            # content[k.encode()] = checksum_function(str(fa_object[k]).encode()).hexdigest().encode()
        collection_string = ";".join([":".join(i) for i in content.items()])
        collection_checksum = self.load_seq(collection_string)
       
        return collection_checksum, content


