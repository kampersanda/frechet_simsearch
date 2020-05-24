#!/usr/bin/env python3

import os
import sys
import msgpack
import json
import bz2
from argparse import ArgumentParser


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
    for obj in unpacker:
        json.dump(obj, sys.stdout, indent=2)
        sys.stdout.write('\n')
    fp.close()
