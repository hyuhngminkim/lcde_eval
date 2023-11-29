import numpy as np
import struct
import os

def downsample(fn):
    d = np.fromfile("data/" + fn + "_200M_uint64", dtype=np.uint64)[1:]
    offset = 10
    nd = d[::offset]
    with open("data/" + fn + "_20M_uint64", "wb") as f:
        f.write(struct.pack("Q", len(nd)))
        print(len(nd))
        nd.tofile(f)

downsample("books")
downsample("fb")
downsample("osm_cellids")
downsample("wiki_ts")

