#!/usr/bin/env python3

import msgpack
import bz2
import numpy as np
from argparse import ArgumentParser
from pprint import pprint

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('path')
    args = parser.parse_args()

    path = args.path
    if path.endswith(".bz2"):
        fp = bz2.open(path, "rb")
    else:
        fp = open(path, "rb")

    unpacker = msgpack.Unpacker(fp, raw=False)
    obj = unpacker.unpack()
    fp.close()

    scores = []
    num_queries = len(obj)

    for dists in obj:
        scores += [dist[1] for dist in dists]

    scores.sort()
    print('Top-1:   ', scores[num_queries - 1])
    print('Top-10:  ', scores[num_queries * 10 - 1])
    print('Top-100: ', scores[num_queries * 100 - 1])
    print('Top-1000:', scores[num_queries * 1000 - 1])
