#!/usr/bin/env python3

import os
import sys
import msgpack
import pandas as pd
import numpy as np
import random
from pprint import pprint
from argparse import ArgumentParser


def pack(path, metadata, objs):
    with open(path, 'wb') as f:
        packer = msgpack.Packer(encoding='utf-8')
        f.write(packer.pack(metadata))
        for obj in objs:
            f.write(packer.pack(obj))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('query_size', type=int, default=1000)
    args = parser.parse_args()

    path = args.path
    unpacker = msgpack.Unpacker(open(path, 'rb'), encoding='utf-8')
    metadata = unpacker.unpack()

    objs = [obj for obj in unpacker]
    random.shuffle(objs)

    base, ext = os.path.splitext(path)
    pack(f'{base}_query{ext}', metadata, objs[:args.query_size])
    pack(f'{base}_base{ext}', metadata, objs[args.query_size:])
