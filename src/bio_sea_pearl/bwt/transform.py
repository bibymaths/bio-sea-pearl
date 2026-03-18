#!/usr/bin/env python3
import argparse
import sys
import pickle

def suffix_array(s: str) -> list[int]:
    n = len(s)
    rank = [ord(c) for c in s]
    sa = list(range(n))
    tmp = [0]*n
    k = 1
    while True:
        sa.sort(key=lambda i: (rank[i], rank[i+k] if i+k<n else -1))
        tmp[sa[0]] = 0
        for i in range(1,n):
            prev, curr = sa[i-1], sa[i]
            prev_key = (rank[prev], rank[prev+k] if prev+k<n else -1)
            curr_key = (rank[curr], rank[curr+k] if curr+k<n else -1)
            tmp[curr] = tmp[prev] + (curr_key != prev_key)
        rank, tmp = tmp, rank
        if rank[sa[-1]] == n-1:
            break
        k <<= 1
    return sa

def bwt(s: str, sentinel: str='$') -> str:
    s = s + sentinel
    sa = suffix_array(s)
    return ''.join(s[i-1] if i>0 else s[-1] for i in sa)

class FMIndex:
    def __init__(self, text: str, sentinel: str='$'):
        # add the sentinel to your text
        self._text = text + sentinel

        # build & save the suffix array on the full text
        self.sa = suffix_array(self._text)

        # build the BWT directly from that SA
        self.bwt = ''.join(
            self._text[i-1] if i>0 else self._text[-1]
            for i in self.sa
        )
        # build C and Occ exactly as before
        self.sentinel = sentinel
        self._build_c()
        self._build_occ()


    def _build_c(self):
        # C[c] = total number of chars in text < c
        counts = {}
        for ch in self.bwt:
            counts[ch] = counts.get(ch, 0) + 1
        self.C = {}
        total = 0
        for ch in sorted(counts):
            self.C[ch] = total
            total += counts[ch]

    def _build_occ(self):
        # Occ[ch][i] = number of occurrences of ch in BWT[0:i]
        self.Occ = {}
        # initialize zero-row
        for ch in self.C:
            self.Occ[ch] = [0]
        for i, ch in enumerate(self.bwt, start=1):
            for c in self.Occ:
                self.Occ[c].append(self.Occ[c][-1] + (1 if c==ch else 0))

    def backward_search(self, pattern: str) -> (int,int):
        """Return the half-open interval [l,r) in the suffix array
           where suffixes start with pattern."""
        l = 0
        r = len(self.bwt)
        for ch in reversed(pattern):
            if ch not in self.C:
                return 0, 0
            l = self.C[ch] + self.Occ[ch][l]
            r = self.C[ch] + self.Occ[ch][r]
            if l >= r:
                return 0, 0
        return l, r

def parse_fasta(fp):
    header = None
    seq_lines = []
    for line in fp:
        line = line.rstrip()
        if not line: continue
        if line.startswith('>'):
            if header:
                yield header, ''.join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        yield header, ''.join(seq_lines)

def main():
    p = argparse.ArgumentParser(
        description="Build BWT + FM-index for each sequence in a FASTA file."
    )
    p.add_argument("fasta", help="Input FASTA file (‘-’ for stdin)")
    p.add_argument("-s","--sentinel",metavar="CHAR",
                   help="Sentinel (not in sequences)", default="$")
    args = p.parse_args()

    inf = sys.stdin if args.fasta=="-" else open(args.fasta)
    for header, seq in parse_fasta(inf):
        idx = FMIndex(seq, sentinel=args.sentinel)
        # save the FMIndex object
        fname = f"{header}.fmidx"
        with open(fname, "wb") as f:
            pickle.dump(idx, f)
        print(f"[{header}] BWT length={len(idx.bwt)}  C={idx.C}  (saved → {fname})",
              file=sys.stderr)

    if inf is not sys.stdin:
        inf.close()

if __name__=="__main__":
    main()
