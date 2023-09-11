import binascii
import hashlib
import json
import os

from jsonschema import Draft7Validator
from typing import Optional
from yacman import load_yaml

from .const import SeqCol
from .exceptions import *


# Refget digests from published seqcol v1.0 protocol
# Retrieved July 2019
# http://samtools.github.io/hts-specs/refget.html
def trunc512_digest(seq, offset=24) -> str:
    digest = hashlib.sha512(seq.encode()).digest()
    hex_digest = binascii.hexlify(digest[:offset])
    return hex_digest.decode()

def sha512t24u_digest(seq: str, offset: int = 24) -> str:
    """ GA4GH digest function """
    digest = hashlib.sha512(seq.encode()).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:offset])
    return tdigest_b64us.decode("ascii")


def canonical_str(item: dict) -> str:
    """Convert a dict into a canonical string representation"""
    return json.dumps(
        item, separators=(",", ":"), ensure_ascii=False, allow_nan=False, sort_keys=True
    )


def print_csc(csc: dict) -> str:
    """Convenience function to pretty-print a canonical sequence collection"""
    return print(json.dumps(csc, indent=2))


# Simple true/false validation
def validate_seqcol_bool(seqcol_obj: SeqCol, schema=None) -> bool:
    schema_path = os.path.join(os.path.dirname(__file__), "schemas", "seqcol.yaml")
    schema = load_yaml(schema_path)
    validator = Draft7Validator(schema)
    return validator.is_valid(seqcol_obj)


# Get errors if invalid (use this one)
# Get the errors with exception.errors
def validate_seqcol(seqcol_obj: SeqCol, schema=None) -> Optional[dict]:
    schema_path = os.path.join(os.path.dirname(__file__), "schemas", "seqcol.yaml")
    schema = load_yaml(schema_path)
    validator = Draft7Validator(schema)
    if not validator.is_valid(seqcol_obj):
        errors = sorted(validator.iter_errors(seqcol_obj), key=lambda e: e.path)
        raise InvalidSeqColError("Validation failed", errors)
    return True
