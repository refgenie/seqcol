import henge
import logging
import pyfaidx

from typing import Callable

from .utilities import *
from .const import *


_LOGGER = logging.getLogger(__name__)
henge.ITEM_TYPE = "_item_type"


def explain_flag(flag):
    """Explains a compare flag"""
    print(f"Flag: {flag}\nBinary: {bin(flag)}\n")
    for e in range(0, 13):
        if flag & 2**e:
            print(FLAGS[2**e])

def fasta_to_digest(fa_file_path: str) -> str:
    """Given a fasta, return a digest"""
    seqcol_obj = fasta_to_seqcol(fa_file_path)
    return seqcol_digest(seqcol_obj)


def parse_fasta(fa_file) -> pyfaidx.Fasta:
    """
    Read in a gzipped or not gzipped FASTA file
    """
    try:
        return pyfaidx.Fasta(fa_file)
    except pyfaidx.UnsupportedCompressionFormat:
        # pyfaidx can handle bgzip but not gzip; so we just hack it here and
        # gunzip the file into a temporary one and read it in not to interfere
        # with the original one.
        from gzip import open as gzopen
        from tempfile import NamedTemporaryFile

        with gzopen(fa_file, "rt") as f_in, NamedTemporaryFile(mode="w+t", suffix=".fa") as f_out:
            f_out.writelines(f_in.read())
            f_out.seek(0)
            return pyfaidx.Fasta(f_out.name)


def fasta_to_seqcol(fa_file_path: str) -> dict:
    """Given a fasta, return a canonical seqcol object"""
    fa_obj = parse_fasta(fa_file_path)
    return fasta_obj_to_seqcol(fa_obj)


def fasta_obj_to_seqcol(
    fa_object: pyfaidx.Fasta,
    verbose: bool = True,
    digest_function: Callable[[str], str] = sha512t24u_digest,
) -> dict:
    """
    Given a fasta object, return a CSC (Canonical Sequence Collection object)
    """
    # CSC = SeqColArraySet
    # Or maybe should be "Level 1 SC"

    CSC = {"lengths": [], "names": [], "sequences": [], "sorted_name_length_pairs": []}
    seqs = fa_object.keys()
    nseqs = len(seqs)
    print(f"Found {nseqs} chromosomes")
    i = 1
    for k in fa_object.keys():
        if verbose:
            print(f"Processing ({i} of {nseqs}) {k}...")
        seq = str(fa_object[k])
        seq_length = len(seq)
        seq_name = fa_object[k].name
        seq_digest = digest_function(seq.upper())
        snlp = {"length": seq_length, "name": seq_name}  # sorted_name_length_pairs
        snlp_digest = digest_function(canonical_str(snlp))
        CSC["lengths"].append(seq_length)
        CSC["names"].append(seq_name)
        CSC["sorted_name_length_pairs"].append(snlp_digest)
        CSC["sequences"].append(seq_digest)
        i += 1
    CSC["sorted_name_length_pairs"].sort()
    return CSC


def build_sorted_name_length_pairs(obj: dict, digest_function):
    """Builds the sorted_name_length_pairs attribute, which corresponds to the coordinate system"""
    sorted_name_length_pairs = []
    for i in range(len(obj["names"])):
        sorted_name_length_pairs.append({"length": obj["lengths"][i], "name": obj["names"][i]})
    nl_digests = []
    for i in range(len(sorted_name_length_pairs)):
        nl_digests.append(digest_function(canonical_str(sorted_name_length_pairs[i])))

    nl_digests.sort()
    return nl_digests


def compare_seqcols(A: SeqCol, B: SeqCol):
    """
    Workhorse comparison function

    @param A Sequence collection A
    @param B Sequence collection B
    @return dict Following formal seqcol specification comparison function return value
    """
    validate_seqcol(A)  # First ensure these are the right structure
    validate_seqcol(B)

    all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
    result = {}
    return_obj = {
        "arrays": {"a-only": [], "b-only": [], "a-and-b": []},
        "elements": {
            "total": {"a": len(A["lengths"]), "b": len(B["lengths"])},
            "a-and-b": {},
            "a-and-b-same-order": {},
        },
    }

    for k in all_keys:
        _LOGGER.info(k)
        if k not in A:
            result[k] = {"flag": -1}
            return_obj["arrays"]["b-only"].append(k)
        elif k not in B:
            return_obj["arrays"]["a-only"].append(k)
        else:
            return_obj["arrays"]["a-and-b"].append(k)
            res = _compare_elements(A[k], B[k])
            return_obj["elements"]["a-and-b"][k] = res["a-and-b"]
            return_obj["elements"]["a-and-b-same-order"][k] = res["a-and-b-same-order"]
    return return_obj


def _compare_elements(A: list, B: list):
    """
    Compare elements between two arrays. Helper function for individual elements used by workhorse compare_seqcols function
    """

    A_filtered = list(filter(lambda x: x in B, A))
    B_filtered = list(filter(lambda x: x in A, B))
    A_count = len(A_filtered)
    B_count = len(B_filtered)
    overlap = min(len(A_filtered), len(B_filtered))  # counts duplicates

    if A_count + B_count < 1:
        # order match requires at least 2 matching elements
        order = None
    elif not (A_count == B_count == overlap):
        # duplicated matches means order match is undefined
        order = None
    else:
        order = A_filtered == B_filtered
    return {"a-and-b": overlap, "a-and-b-same-order": order}


def seqcol_digest(seqcol_obj: SeqCol, schema: dict = None) -> str:
    """
    Given a canonical sequence collection, compute its digest.

    :param dict seqcol_obj: Dictionary representation of a canonical sequence collection object
    :param dict schema: Schema defining the inherent attributes to digest
    :return str: The sequence collection digest
    """

    validate_seqcol(seqcol_obj)
    # Step 1a: Remove any non-inherent attributes,
    # so that only the inherent attributes contribute to the digest.
    seqcol_obj2 = {}
    if schema:
        for k in schema["inherent"]:
            # Step 2: Apply RFC-8785 to canonicalize the value
            # associated with each attribute individually.
            seqcol_obj2[k] = canonical_str(seqcol_obj[k])
    else:  # no schema provided, so assume all attributes are inherent
        for k in seqcol_obj:
            seqcol_obj2[k] = canonical_str(seqcol_obj[k])
    # Step 3: Digest each canonicalized attribute value
    # using the GA4GH digest algorithm.

    seqcol_obj3 = {}
    for attribute in seqcol_obj2:
        seqcol_obj3[attribute] = sha512t24u_digest(seqcol_obj2[attribute])
    # print(json.dumps(seqcol_obj3, indent=2))  # visualize the result

    # Step 4: Apply RFC-8785 again to canonicalize the JSON
    # of new seqcol object representation.

    seqcol_obj4 = canonical_str(seqcol_obj3)

    # Step 5: Digest the final canonical representation again.
    seqcol_digest = sha512t24u_digest(seqcol_obj4)
    return seqcol_digest
