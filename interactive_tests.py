import seqcol
from seqcol import SeqColClient


scc = SeqColClient({}, schemas=["seqcol/schemas/SeqColArraySet.yaml"])
scc

fa_file = "demo_fasta/demo.fa.gz"
fa_object = seqcol.parse_fasta(fa_file)

skip_seq = False
aslist = []
names = []
lengths = []
sequences = []
for k in fa_object.keys():
    seq = str(fa_object[k])
    names.append(k)
    lengths.append(str(len(seq)))
    sequences.append(seq)

array_set = {"names": names, "lengths": lengths, "sequences": sequences}

collection_checksum = scc.insert(array_set, "SeqColArraySet")
collection_checksum

scc.retrieve("d229d5c16b3a1b3788f01aa439f01e682ba84bc9935ad08a")

scc.retrieve("d229d5c16b3a1b3788f01aa439f01e682ba84bc9935ad08a", reclimit=1)


from henge import md5

names = []
lengths = []
seq_digests = []
for k in fa_object.keys():
    seq = str(fa_object[k])
    names.append(k)
    lengths.append(str(len(seq)))
    seq_digests.append(scc.checksum_function(seq))


array_set = {"names": names, "lengths": lengths, "sequences": seq_digests}
array_set
collection_checksum = scc.insert(array_set, "SeqColArraySet", reclimit=1)

scc.database = {}
scc.retrieve("d229d5c16b3a1b3788f01aa439f01e682ba84bc9935ad08a", reclimit=1)
scc.retrieve("d229d5c16b3a1b3788f01aa439f01e682ba84bc9935ad08a", reclimit=2)
scc.database[collection_checksum]
scc.checksum_function


scc.database["d229d5c16b3a1b3788f01aa439f01e682ba84bc9935ad08a"]

scc.database["a99fe2e92875099c84f73b20ef8e7dd2f2d12f063383bed0"]
scc.database["ca82b053295b6f49923d0b2cedb83de49c6be59688c3dfd9"]


import os

os.getcwd()
