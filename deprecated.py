from itertools import compress
from functools import reduce

import seqcol

sc1 = {"names": ["chr1", "chr2"], "sequences": ["ay89fw", "we4f9x"]}

sc2 = {"names": ["1", "2", "3"], "sequences": ["ay89fw", "we4f9x", "3n20xk2"]}

sc3 = {"names": ["2", "3", "1"], "sequences": ["we4f9x", "3n20xk2", "ay89fw"]}

sc4 = {"names": ["chr1", "3", "2"], "sequences": ["zyt2fw", "9snm23k", "fsd2x3"]}

sc5 = {
    "names": ["chr1", "3", "2"],
    "sequences": ["zyt2fw", "9snm23k", "fsd2x3"],
    "topologies": ["circular", "linear", "linear"],
}









   @staticmethod
    def compat_all_old(A, B):
        all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
        result = {}
        flipped_format = {
            "a-in-b": {},
            "b-in-a": {},
            "a-total": {},
            "b-total": {},
            "a-duplicated": {},
            "b-duplicated": {},
            "order-match": [],
            "only-in-a": [],
            "only-in-b": [],
        }
        for k in all_keys:
            _LOGGER.info(k)
            if k not in A:
                result[k] = {"flag": -1}
                flipped_format["only-in-b"].append(k)
            elif k not in B:
                flipped_format["only-in-a"].append(k)
            else:
                v = SeqColClient.compat(A[k], B[k])
                result[k] = v
                if "a-in-b" in v:
                    flipped_format["a-in-b"][k] = v['a-in-b']
                if "b-in-a":
                    flipped_format["b-in-a"][k] = v['b-in-a']
                if "a-total" in v:
                    flipped_format["a-total"][k] = v['a-total']
                if "b-total" in v:
                    flipped_format["b-total"][k] = v['b-total']
                if "a-duplicated" in v:
                    flipped_format["a-duplicated"][k] = v['a-duplicated']
                if "b-duplicated" in v:
                    flipped_format["b-duplicated"][k] = v['b-duplicated']
                if "order-match" in v:
                    flipped_format["order-match"].append(k)

        # result = {
        #     "any-elements-shared": any(ainb),
        #     "all-a-in-b": all(ainb),
        #     "all-b-in-a": all(bina),
        #     "order-match": order,
        #     "flag": flag
        # }

        return flipped_format


    def compare_digests_old(self, digestA, digestB, explain=False):
        """
        Given two collection checksums in the database, provide some information
        about how they are related.

        :param str digestA: Digest for first sequence collection to compare.
        :param str digestB: Digest for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """
        typeA = self.database[digestA + henge.ITEM_TYPE]
        typeB = self.database[digestB + henge.ITEM_TYPE]

        if typeA != typeB:
            _LOGGER.error(
                f"Can't compare objects of different types: " f"{typeA} vs {typeB}"
            )

        asdA = self.retrieve(digestA, reclimit=1)
        asdB = self.retrieve(digestB, reclimit=1)
        return self.compare_asds(asdA, asdB, explain=explain)


   @staticmethod
    def compare_asds(asdA, asdB, explain=False):
        """
        Compare Annotated Sequence Digests (ASDs) -- digested sequences and `data

        :param str asdA: ASD for first sequence collection to compare.
        :param str asdB: ASD for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """

        def _xp(prop, lst):
            """Extract property from a list of dicts"""
            return list(map(lambda x: x[prop], lst))

        def _index(x, lst):
            """Find an index of a sequence element in a list of dicts"""
            try:
                return _xp(SEQ_KEY, lst).index(x)
            except:
                return None

        def _get_common_content(lstA, lstB):
            """
            Find the intersection between two list of dicts with sequences
            """
            return list(
                filter(None.__ne__, [_index(x, lstB) for x in _xp(SEQ_KEY, lstA)])
            )

        # Not ideal, but we expect these to return lists, but if the item was
        # singular only a dict is returned
        if not isinstance(asdA, list):
            asdA = [asdA]
        if not isinstance(asdB, list):
            asdB = [asdB]

        ainb = [x in _xp(SEQ_KEY, asdB) for x in _xp(SEQ_KEY, asdA)]
        bina = [x in _xp(SEQ_KEY, asdA) for x in _xp(SEQ_KEY, asdB)]

        return_flag = 0  # initialize
        if sum(ainb) > 1:
            ordA = _get_common_content(asdA, asdB)
            if ordA == sorted(ordA):
                return_flag += CONTENT_A_ORDER
        if sum(bina) > 1:
            ordB = _get_common_content(asdB, asdA)
            if ordB == sorted(ordB):
                return_flag += CONTENT_B_ORDER

        ainb_len = [x in _xp(LEN_KEY, asdB) for x in _xp(LEN_KEY, asdA)]
        bina_len = [x in _xp(LEN_KEY, asdA) for x in _xp(LEN_KEY, asdB)]

        ainb_name = [x in _xp(NAME_KEY, asdB) for x in _xp(NAME_KEY, asdA)]
        bina_name = [x in _xp(NAME_KEY, asdA) for x in _xp(NAME_KEY, asdB)]

        ainb_topo = [x in _xp(TOPO_KEY, asdB) for x in _xp(TOPO_KEY, asdA)]
        bina_topo = [x in _xp(TOPO_KEY, asdA) for x in _xp(TOPO_KEY, asdB)]

        if all(ainb):
            return_flag += CONTENT_ALL_A_IN_B
        if all(bina):
            return_flag += CONTENT_ALL_B_IN_A

        if all(ainb_name):
            return_flag += NAMES_ALL_A_IN_B
        if all(bina_name):
            return_flag += NAMES_ALL_B_IN_A

        if all(ainb_topo):
            return_flag += TOPO_ALL_A_IN_B
        if all(bina_topo):
            return_flag += TOPO_ALL_B_IN_A

        if all(ainb_len):
            return_flag += LENGTHS_ALL_A_IN_B
        if all(bina_len):
            return_flag += LENGTHS_ALL_B_IN_A

        if explain:
            explain_flag(return_flag)
        return return_flag


    @staticmethod
    def compat(A, B):
        """
        New compatibility function for array-based data model.
        """

        lenA = len(A)
        lenB = len(B)
        dupeA = lenA - len(dict.fromkeys(A))
        dupeB = lenB - len(dict.fromkeys(B))
        ainb = [x in B for x in A]
        bina = [x in A for x in B]
        sum_ainb = sum(ainb)
        if sum_ainb > 1:
            order = list(compress(B, bina)) == list(compress(A, ainb))
        else:
            order = False

        result = {
            "a-in-b": sum_ainb,
            "b-in-a":  sum(bina),
            "a-total": lenA,
            "b-total": lenB,
            "a-duplicated": dupeA,
            "b-duplicated": dupeB,
            "order-match": order
        }
        return result


def compat(A, B):
    ainb = [x in B for x in A]
    bina = [x in A for x in B]
    if any(ainb):
        order = list(compress(B, bina)) == list(compress(A, ainb))
    else:
        order = False

    any(ainb)

    flag = 0
    flag += 2 if all(ainb) else 0
    flag += 4 if all(bina) else 0
    flag += 8 if order else 0
    flag += 1 if any(ainb) else 0

    return flag


# New compat function that adds true/false
def compat(A, B):
    ainb = [x in B for x in A]
    bina = [x in A for x in B]
    if any(ainb):
        order = list(compress(B, bina)) == list(compress(A, ainb))
    else:
        order = False

    any(ainb)

    flag = 0
    flag += 2 if all(ainb) else 0
    flag += 4 if all(bina) else 0
    flag += 8 if order else 0
    flag += 1 if any(ainb) else 0
    result = {
        "any-elements-shared": any(ainb),
        "a-subset-of-b": all(ainb),
        "b-subset-of-a": all(bina),
        "order-match": order,
        "flag": flag,
    }
    return result


# For each array:
# - any-elements-shared (1) (0001)
# - all-a-in-b        (2) (0010)
# - all-b-in-a        (4) (0100)
# - order-match       (8) (1000)

# no match: 0000 = 0
# one or more shared elements: 0001 = 1
# all-a-in-b (a is a subset of b, in different order) = 0011 = 3
# all-b-in-a (b is a subset of a, in different order) = 0101 = 5
# same content, different order: 0111 = 7
# a is a subset of b, in same order: 1011 = 11
# b is a subset of a, in same order: 1101 = 13
# identity: 1111 = 15


compat(sc1["sequences"], sc2["sequences"])
compat(sc1["names"], sc2["names"])

compat(sc3["sequences"], sc2["sequences"])
compat(sc3["names"], sc2["names"])

compat(sc1["sequences"], sc3["sequences"])


def compat_all(A, B):
    all_keys = list(A.keys()) + list(set(B.keys()) - set(list(A.keys())))
    result = {}
    for k in all_keys:
        if k not in A or k not in B:
            result[k] = {"flag": -1}
        else:
            result[k] = compat(A[k], B[k])
    # result["all"] = reduce(lambda x,y: x['flag']&y['flag'], list(result.values()))
    return result


compat_all(sc1, sc2)
compat_all(sc3, sc2)
compat_all(sc1, sc3)
compat_all(sc1, sc5)
