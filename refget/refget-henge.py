import henge
import logmuse
import mongodict
import yaml

class RefGetHenge(henge.Henge):
    """
    Extension of henge that accommodates refget sequences.
    """
    def refget(self, digest):
        if self.database[digest + "_item_type"] == "sequence":
            return self.retrieve(digest)['sequence']
        elif self.database[digest + "_item_type"] == "asd":
            return[{x['name']: x['sequence_digest']['sequence']} for x in self.retrieve(digest)]



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

