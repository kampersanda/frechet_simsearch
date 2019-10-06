#!/usr/bin/env python3

import os
import sys
import msgpack
from glob import glob
from argparse import ArgumentParser


def read_giscup(input_dir):
    curves = []
    filenames = glob(input_dir + '/*.dat')
    for i, filename in enumerate(filenames):
        points = []
        for row in open(filename).readlines()[1:]:
            items = row.rstrip().split(' ')
            points.append([float(items[0]), float(items[1])])
        curves.append({
            'points': points,
            'id': i
        })
    return curves


def write_curves(curves, path):
    with open(path, 'wb') as fp:
        msgpack.pack({'dimensions': '2'}, fp)
        for curve in curves:
            msgpack.pack(curve, fp)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_dir')
    parser.add_argument('output_path')
    args = parser.parse_args()

    write_curves(read_giscup(args.input_dir), args.output_path)
