import henge
import logmuse
import mongodict
import yaml


henge.ITEM_TYPE = "_item_type"

class RefGetHenge(henge.Henge):
    """
    Extension of henge that accommodates refget sequences.
    """
    def refget(self, digest, reclimit=None, return_fmt=None):
        item_type = self.database[digest + henge.ITEM_TYPE]

        if item_type == "sequence":
            return self.retrieve(digest)['sequence']
        elif item_type == "asd":
            asdlist = [{x['name']: x['sequence_digest']['sequence']} for x in self.retrieve(digest)]
            if return_fmt == "fasta":
                return self.fasta_fmt(asdlist)
            else:
                return asdlist

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

        checksum = self.insert({'sequence': seq}, "sequence")
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
        asdlist = {}
        for k in fa_object.keys():
            seq = str(fa_object[k])
            asdlist[k] = {'name': k,
                          'length': len(seq), 
                          'toplogy': 'linear',
                          'sequence_digest': self.load_seq(seq)}

        collection_checksum = self.insert(asdlist, 'asd')
        return collection_checksum, asdlist

    def load_seqset(self, seqset):
        """
        Convert a 'seqset', which is a dict with names as sequence names and
        values as sequences, into the 'asdlist' required for henge insert.
        """
        seqset_new = copy(seqset)
        for k, v in seqset:
            if 'length' not in v.keys():
                if 'sequence' not in v.keys():
                    _LOGGER.error("Each sequence must have either length or a sequence.")
                else:
                    v['length'] = len(v['sequence'])
            if 'sequence' in i.keys():
                v['sequence_digest'] = self.load_seq(seq)
            if 'name' not in v.keys():
                v['name'] = k

            seqset_new[k] = v

        collection_checksum = self.insert(seqset_new, 'asd')
        return collection_checksum, seqset_new


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



logmuse.init_logger("henge", "DEBUG", devmode=True)
backend = henge.MongoDict(host='localhost', port=27017, database='my_dict',
                        collection='store')
# schemas = {"sequence": yaml.safe_load(seq_schema), "asd": yaml.safe_load(asd),
            # "acd": yaml.safe_load(acd)}

def load_yaml(filename):
    with open(filename) as f:
        return yaml.safe_load(f)

schemas = {"sequence": load_yaml("sequence.yaml"), "asd": load_yaml("annotated_sequence_digest.yaml"),
            "acd": load_yaml("annotated_collection_digest.yaml")}

h = RefGetHenge(backend, schemas=schemas)

item_seq1 = {'sequence': "TCGA"}
item_seq2 = {'sequence': "TCGATCGATCGATCGA"}
item_seq3 = {'sequence': "GGAA"}
item_seq4 = {'sequence': "CGGCCCGGCGC"}

druids1 = h.insert([item_seq1], "sequence")
druids2 = h.insert([item_seq2], "sequence")
druids3 = h.insert([item_seq3], "sequence")
druids4 = h.insert([item_seq4], "sequence")


asd1 = {'sequence_digest': druids1,
        'name': "chr1",
        'length': 10, 
        'topology': "linear"}
asd2 = {'sequence_digest': druids2,
        'name': "chr2",
        'length': 20, 
        'topology': "linear"}
asd3 = {'sequence_digest': druids3,
        'name': "chr3",
        'length': 30, 
        'topology': "circular"}
asd4 = {'sequence_digest': druids4,
        'name': "chr4::mod",
        'length': 40, 
        'topology': "linear"}                

druidasd1 = h.insert([asd1, asd2], "asd")
druidasd2 = h.insert([asd3, asd4], "asd")

acd1 = {'collection_digest': druidasd1,
        'name': "fasta1"} 
acd2 = {'collection_digest': druidasd2,
        'name': "fasta2"} 

druidacd = h.insert([acd1, acd2], "acd")    

h.retrieve(druidacd, reclimit=1)
h.retrieve(druidacd, reclimit=2)
h.retrieve(druidacd, reclimit=0)
druidacd

h.retrieve(druids1)

h.refget(druids1)

h.retrieve(druidasd1)
h.refget(druidasd1)

h.show()

