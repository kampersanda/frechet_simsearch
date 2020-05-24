#!/usr/bin/env python3

import msgpack
import numpy as np
from argparse import ArgumentParser


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('curve_path')
    args = parser.parse_args()

    curve_path = args.curve_path

    unpacker = msgpack.Unpacker(open(curve_path, 'rb'), raw=False)
    metadata = unpacker.unpack()

    curves = [curve for curve in unpacker]
    lengths = np.array([len(curve['points']) for curve in curves])

    num_curves = len(curves)
    num_points = np.sum(lengths)
    min_length = np.min(lengths)
    max_length = np.max(lengths)
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)

    print('num_curves:', num_curves)
    print('num_points:', num_points)
    print('min_length:', min_length)
    print('max_length:', max_length)
    print('mean_length:', mean_length)
    print('median_length:', median_length)
