"""
make_block.py - write a caDNAno (square-lattice) JSON for blocks of parallel helices.

The mechanical test structures need no scaffold routing: every helix segment carries its
own "scaffold" strand, and all helix-to-helix connections are staple crossovers at the
canonical square-lattice positions (pairs 7|8, 15|16, 23|24, 31|32 in each 32-bp period,
cycling through the four neighbours). Staples are then nicked mid-segment into pieces of
about 32-40 nt, as an autobreak would do.

    helices = {(row, col): (start, end)}                    # one segment per helix, end exclusive
    helices = {(row, col): [(s1, e1), (s2, e2)]}            # or several disjoint segments
    write_cadnano(helices, "blade.json")
"""
import json

PERIOD = 32
# crossover pair (p | p+1) with p = PAIR_LOW[k], to neighbour direction k
PAIR_LOW = {0: 31, 1: 23, 2: 15, 3: 7}


def neighbours(r, c):
    """caDNAno square-lattice neighbour order p0..p3 (parity dependent)."""
    if (r + c) % 2 == 0:
        return [(r, c + 1), (r + 1, c), (r, c - 1), (r - 1, c)]
    return [(r, c - 1), (r - 1, c), (r, c + 1), (r + 1, c)]


def assign_numbers(coords):
    num, even, odd = {}, 0, 1
    for rc in sorted(coords):
        if (rc[0] + rc[1]) % 2 == 0:
            num[rc] = even
            even += 2
        else:
            num[rc] = odd
            odd += 2
    return num


def normalise(helices):
    out = {}
    for rc, segs in helices.items():
        if isinstance(segs, tuple) and len(segs) == 2 and isinstance(segs[0], int):
            segs = [segs]
        out[rc] = sorted(tuple(s) for s in segs)
    return out


def build(helices, margin=3, min_len=32, max_len=56, skip_xover=None, verbose=True):
    """
    margin: minimum number of base pairs between a crossover and a segment end
    skip_xover: optional callable (rcA, rcB, p) -> True to omit that crossover
    """
    helices = normalise(helices)
    num = assign_numbers(helices.keys())
    array_len = PERIOD * ((max(e for segs in helices.values() for s, e in segs) + PERIOD - 1) // PERIOD)

    scaf_next, stap_next = {}, {}
    for rc, segs in helices.items():
        even = (rc[0] + rc[1]) % 2 == 0
        for s, e in segs:
            for i in range(s, e):
                if even:      # scaffold 5'->3' towards increasing index; staple decreasing
                    if i + 1 < e:
                        scaf_next[(rc, i)] = (rc, i + 1)
                    if i - 1 >= s:
                        stap_next[(rc, i)] = (rc, i - 1)
                else:
                    if i - 1 >= s:
                        scaf_next[(rc, i)] = (rc, i - 1)
                    if i + 1 < e:
                        stap_next[(rc, i)] = (rc, i + 1)

    xovers = []
    for rc, segs in helices.items():
        if (rc[0] + rc[1]) % 2 != 0:
            continue                                  # handle each pair once, from the even helix
        for k, nb in enumerate(neighbours(*rc)):
            if nb not in helices:
                continue
            for s, e in segs:
                for s2, e2 in helices[nb]:
                    lo, hi = max(s, s2), min(e, e2)
                    p = PAIR_LOW[k]
                    while p + 1 < hi:
                        if p - margin + 1 >= lo and p + 1 + margin <= hi:
                            if skip_xover is None or not skip_xover(rc, nb, p):
                                # even helix E (staple runs down), odd helix O (staple runs up)
                                stap_next[(rc, p + 1)] = (nb, p + 1)
                                stap_next[(nb, p)] = (rc, p)
                                xovers.append((rc, nb, p))
                        p += PERIOD
    if verbose:
        print(f"helices: {len(helices)}   staple crossover pairs: {len(xovers)}")

    def is_cut_site(a, b):
        """legal nick: same helix, between indices m | m+1 with m % 8 == 3 (mid-segment)."""
        if a[0] != b[0]:
            return False
        m = min(a[1], b[1])
        return m % 8 == 3

    def strands_of(nxt):
        prev = {v: k for k, v in nxt.items()}
        nodes = set(nxt) | set(prev)
        for rc, segs in helices.items():
            for s, e in segs:
                nodes |= {(rc, i) for i in range(s, e)}
        seen, out = set(), []
        for n in sorted(nodes):
            if n in prev or n in seen:
                continue
            path = [n]
            seen.add(n)
            while path[-1] in nxt:
                path.append(nxt[path[-1]])
                seen.add(path[-1])
            out.append((path, False))
        for n in sorted(nodes):                       # whatever is left is circular
            if n in seen:
                continue
            path = [n]
            seen.add(n)
            while nxt[path[-1]] != n:
                path.append(nxt[path[-1]])
                seen.add(path[-1])
            out.append((path, True))
        return out

    # open circular staples at a legal site, then cut long staples
    for path, circ in strands_of(stap_next):
        if not circ:
            continue
        L = len(path)
        for j in range(L):
            a, b = path[j], path[(j + 1) % L]
            if is_cut_site(a, b):
                del stap_next[a]
                break
        else:
            raise RuntimeError("circular staple with no legal nick site")

    for path, circ in strands_of(stap_next):
        assert not circ
        since = 0
        for j in range(len(path) - 1):
            since += 1
            remaining = len(path) - (j + 1)
            if since >= min_len and remaining >= 16 and is_cut_site(path[j], path[j + 1]):
                del stap_next[path[j]]
                since = 0

    staples = [p for p, c in strands_of(stap_next)]
    lens = sorted(len(p) for p in staples)
    if verbose:
        print(f"staples: {len(staples)}   length min/median/max: {lens[0]}/{lens[len(lens)//2]}/{lens[-1]}")
        too_long = [l for l in lens if l > max_len]
        if too_long:
            print(f"  note: {len(too_long)} staples longer than {max_len} nt: {too_long}")
    return helices, num, array_len, scaf_next, stap_next, staples, xovers


def write_cadnano(helices, path, **kw):
    helices, num, n, scaf_next, stap_next, staples, xovers = build(helices, **kw)
    scaf_prev = {v: k for k, v in scaf_next.items()}
    stap_prev = {v: k for k, v in stap_next.items()}

    def cell(rc, i, nxt, prv):
        if not any(s <= i < e for s, e in helices[rc]):
            return [-1, -1, -1, -1]
        p5 = prv.get((rc, i))
        p3 = nxt.get((rc, i))
        return [num[p5[0]] if p5 else -1, p5[1] if p5 else -1,
                num[p3[0]] if p3 else -1, p3[1] if p3 else -1]

    vstrands = []
    for rc in sorted(helices, key=lambda x: num[x]):
        vstrands.append({
            "row": rc[0], "col": rc[1], "num": num[rc],
            "scaf": [cell(rc, i, scaf_next, scaf_prev) for i in range(n)],
            "stap": [cell(rc, i, stap_next, stap_prev) for i in range(n)],
            "loop": [0] * n, "skip": [0] * n,
            "scafLoop": [], "stapLoop": [], "stap_colors": [],
        })
    with open(path, "w") as f:
        json.dump({"name": path, "vstrands": vstrands}, f)
    return num, staples, xovers


if __name__ == "__main__":
    import sys
    rows, cols, length = (int(x) for x in sys.argv[1:4])
    out = sys.argv[4]
    helices = {(r, c): (0, length) for r in range(rows) for c in range(cols)}
    write_cadnano(helices, out)
