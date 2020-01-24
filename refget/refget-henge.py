from copy import copy
import henge
import logmuse
import mongodict
import yaml

import pyfaidx
import logging
_LOGGER = logging.getLogger(__name__)
henge.ITEM_TYPE = "_item_type"

class RefGetHenge(henge.Henge):
    """
    Extension of henge that accommodates refget sequences.
    """
    def refget(self, digest, reclimit=None, postprocess=None):
        item_type = self.database[digest + henge.ITEM_TYPE]

        full_data = self.retrieve(digest, reclimit=reclimit)
        if not postprocess:
            return full_data
    
            full_data = self.retrieve(digest)
        if postprocess == "simplify":
            if item_type == "sequence":
                return full_data['sequence']
            elif item_type == "asd":
                asdlist = {}
                for x in full_data:
                    asdlist[x['name']] = x['sequence_digest']['sequence']
                return asdlist

        if postprocess == "fasta":
            if item_type == "sequence":
                raise Exception("can't postprocess a sequence into fasta")
            elif item_type == "asd":
                asdlist = {}
                for x in full_data:
                    asdlist[x['name']] = x['sequence_digest']['sequence']                
                return self.fasta_fmt(asdlist)

        _LOGGER.error("Not implemented.")

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

        checksum = self.insert([{'sequence': seq}], "sequence")
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
        asdlist = []
        for k in fa_object.keys():
            seq = str(fa_object[k])
            asdlist.append({'name': k,
                          'length': len(seq), 
                          'topology': 'linear',
                          'sequence_digest': self.load_seq(seq)})

        print(asdlist)
        _LOGGER.info(asdlist)
        collection_checksum = self.insert(asdlist, 'asd')
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
                    _LOGGER.error("Each sequence must have either length or a sequence.")
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

        collection_checksum = self.insert(list(seqset_new.values()), 'asd')
        return collection_checksum, seqset_new


CONTENT_ALL_A_IN_B = 2**0
CONTENT_ALL_B_IN_A = 2**1
LENGTHS_ALL_A_IN_B = 2**2
LENGTHS_ALL_B_IN_A = 2**3
NAMES_ALL_A_IN_B = 2**4
NAMES_ALL_B_IN_A = 2**5
TOPO_ALL_A_IN_B = 2**6
TOPO_ALL_B_IN_A = 2**7
CONTENT_ANY_SHARED = 2**8
LENGTHS_ANY_SHARED = 2**9
NAMES_ANY_SHARED = 2**10
CONTENT_A_ORDER = 2**11
CONTENT_B_ORDER = 2**12

FLAGS = {}
FLAGS[CONTENT_ALL_A_IN_B] = "CONTENT_ALL_A_IN_B"
FLAGS[CONTENT_ALL_B_IN_A] = "CONTENT_ALL_B_IN_A"
FLAGS[LENGTHS_ALL_A_IN_B] = "LENGTHS_ALL_A_IN_B"
FLAGS[LENGTHS_ALL_B_IN_A] = "LENGTHS_ALL_B_IN_A"
FLAGS[NAMES_ALL_A_IN_B] = "NAMES_ALL_A_IN_B"
FLAGS[NAMES_ALL_B_IN_A] = "NAMES_ALL_B_IN_A"
FLAGS[TOPO_ALL_A_IN_B] = "TOPO_ALL_A_IN_B"
FLAGS[TOPO_ALL_B_IN_A] = "TOPO_ALL_B_IN_A"
FLAGS[CONTENT_ANY_SHARED] = "CONTENT_ANY_SHARED"
FLAGS[LENGTHS_ANY_SHARED] = "LENGTHS_ANY_SHARED"
FLAGS[NAMES_ANY_SHARED] = "NAMES_ANY_SHARED"
FLAGS[CONTENT_A_ORDER] = "CONTENT_A_ORDER"
FLAGS[CONTENT_B_ORDER] = "CONTENT_B_ORDER"

def explain(flag):
    print(bin(flag))
    for e in range(0,13):
        if flag & 2**e:
            print(FLAGS[2**e])



def compare(rgdb, digestA, digestB):
    """
    Given two collection checksums in the database, provide some information
    about how they are related.
    """
    typeA = rgdb.database[digestA + henge.ITEM_TYPE]
    typeB = rgdb.database[digestB + henge.ITEM_TYPE]

    if typeA != typeB:
        _LOGGER.error("Can't compare objects of different types: {} vs {}".format(typeA, typeB))

    asdA = rgdb.refget(digestA, reclimit=1)
    asdB = rgdb.refget(digestB, reclimit=1)

    def xp(prop, lst):
        """ Extract property """
        return list(map(lambda x: x[prop], lst))

    ainb = [x in xp('sequence_digest', asdB) for x in xp('sequence_digest', asdA)]
    bina = [x in xp('sequence_digest', asdA) for x in xp('sequence_digest', asdB)]
    
    def index(x, lst):
        try:
            return xp('sequence_digest', lst).index(x)
        except:
            return None

    return_flag = 0  # initialize
    if sum(ainb) > 1:
        ordA = list(filter(None.__ne__, [index(x, asdB) for x in xp('sequence_digest', asdA)]))
        if (ordA == sorted(ordA)):
            return_flag += CONTENT_A_ORDER
    if sum(bina) > 1:        
        ordB = list(filter(None.__ne__, [index(x, asdA) for x in xp('sequence_digest', asdB)]))
        if (ordB == sorted(ordB)):
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

    _LOGGER.info(explain(return_flag))
    return return_flag

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
rgdb = h

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


h.load_seq("TCGATTTT")

fa_file = "../demo_fasta/demo2.fa"
checksum2, content2 = h.load_fasta(fa_file)

h.refget(checksum2)
print(h.refget(checksum2, postprocess="fasta"))
h.refget(checksum2, postprocess="simplify")
h.database[checksum2]

h.refget(h.refget(checksum2, reclimit=1)[0]['sequence_digest'])

h.retrieve(druids1)
h.refget(druids1)
h.refget(druids1, postprocess="simplify")



fa_file = "../demo_fasta/demo.fa"
checksum, content = h.load_fasta(fa_file)

h.refget(checksum)
h.refget(checksum2)

compare(h, checksum, checksum2)

ss = [{'name': "chr1",
        'length': 10, 
        'topology': "linear"},
        {'name': "chr2",
        'length': 20, 
        'topology': "linear"}]

ss




druid_ss = h.insert(ss, "asd")

h.retrieve(druid_ss)
h.refget(druid_ss)

h.refget(checksum)
h.refget(checksum2)

f = compare(h, checksum, druid_ss)
explain(f)
f = compare(h, checksum, checksum2)
explain(f)
f = compare(h, checksum2, druid_ss)
explain(f)

digestA = checksum
digestB = druid_ss
digestB = checksum2
bin(40)
bin(38)



lres = {'chr1': {'length': 10},
      "chr2": {'length': 20}}



h.load_seqset(lres)


h.insert(lres, 'asd')
