from transform import FMIndex
import pickle

# Load the pickled FMIndex instance
with open("bwt.fmidx", "rb") as fh:
    idx: FMIndex = pickle.load(fh)

# Now idx is your FMIndex—call its methods:
pattern = "ACGTAG"
l, r = idx.backward_search(pattern)
positions = idx.sa[l:r]

print(f"Found {r-l} hit{'s' if r-l != 1 else ''} for {pattern!r}:")
for pos in positions:
    print(pos)