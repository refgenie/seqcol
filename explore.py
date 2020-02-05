import yaml
import jsonschema

s = yaml.safe_load(open("ASDList.yaml"))
asd = yaml.safe_load(open("annotated_sequence_digest.yaml"))


ss = [{'name': "chr1",
        'length': 10, 
        'topology': "linear"},
        {'name': "chr2",
        'length': 20, 
        'topology': "linear"}]

jsonschema.validate(ss, s)

jsonschema.validate({'name': "chr1",
        'length': 10, 
        'topology': "linear"}, asd)

