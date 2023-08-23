import binascii
import hashlib
import json


# Refget digests from published seqcol v1.0 protocol
# Retrieved July 2019
# http://samtools.github.io/hts-specs/refget.html
def trunc512_digest(seq, offset=24) -> str:
    """GA4GH digest algorithm"""
    digest = hashlib.sha512(seq.encode()).digest()
    hex_digest = binascii.hexlify(digest[:offset])
    return hex_digest.decode()


def canonical_str(item: dict) -> str:
    """Convert a dict into a canonical string representation"""
    return json.dumps(
        item, separators=(",", ":"), ensure_ascii=False, allow_nan=False, sort_keys=True
    )


def print_csc(csc: dict) -> str:
    """Convenience function to pretty-print a canonical sequence collection"""
    return print(json.dumps(csc, indent=2))


def seqcol_digest(seqcol_obj: dict, schema: dict = None) -> str:
    """
    Given a canonical sequence collection, compute its digest.

    :param dict seqcol_obj: Dictionary representation of a canonical sequence collection object
    :param dict schema: Schema defining the inherent attributes to digest
    :return str: The sequence collection digest
    """

    seqcol_obj2 = {}
    for attribute in seqcol_obj:
        seqcol_obj2[attribute] = canonical_str(seqcol_obj[attribute])

    # Step 1a: Remove any non-inherent attributes,
    # so that only the inherent attributes contribute to the digest.
    if schema:
        seqcol_obj2_filtered = {}
        for k in schema["inherent"]:
            seqcol_obj2_filtered[k] = seqcol_obj2[k]
    else:
        seqcol_obj2_filtered = seqcol_obj2

    # Step 2: Apply RFC-8785 to canonicalize the value
    # associated with each attribute individually.
    seqcol_obj2 = {}
    for attribute in seqcol_obj:
        seqcol_obj2[attribute] = canonical_str(seqcol_obj[attribute])
    # seqcol_obj2  # visualize the result

    # Step 3: Digest each canonicalized attribute value
    # using the GA4GH digest algorithm.

    seqcol_obj3 = {}
    for attribute in seqcol_obj2:
        seqcol_obj3[attribute] = trunc512_digest(seqcol_obj2[attribute])
    # print(json.dumps(seqcol_obj3, indent=2))  # visualize the result

    # Step 4: Apply RFC-8785 again to canonicalize the JSON
    # of new seqcol object representation.

    seqcol_obj4 = canonical_str(seqcol_obj3)
    # seqcol_obj4  # visualize the result

    # Step 5: Digest the final canonical representation again.
    seqcol_digest = trunc512_digest(seqcol_obj4)
    return seqcol_digest
