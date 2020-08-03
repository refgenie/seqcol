# Refget digests from published seqcol v1.0 protocol
# Retrieved July 2019
# http://samtools.github.io/hts-specs/refget.html

import hashlib
import binascii


def trunc512_digest(seq, offset=24):
    digest = hashlib.sha512(seq.encode()).digest()
    hex_digest = binascii.hexlify(digest[:offset])
    return hex_digest.decode()
