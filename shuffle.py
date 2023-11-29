import numpy as np
import struct
import os

def shuffle(fn):
    d = np.fromfile("data/" + fn + "_200M_uint64", dtype=np.uint64)[1:]
    np.random.shuffle(d)
    with open("data/" + fn + "_200M_uint64_shuffle", "wb") as f:
        f.write(struct.pack("Q", len(d)))
        print(len(d))
        d.tofile(f)

shuffle("books")
shuffle("fb")
shuffle("osm_cellids")
shuffle("wiki_ts")

