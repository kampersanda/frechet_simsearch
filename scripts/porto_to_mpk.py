#!/usr/bin/env python3

import os
import sys
import re
import msgpack
import pandas as pd
from argparse import ArgumentParser


def deg_to_ms(deg):
    d = int(deg)
    md = abs(deg - d) * 60
    m = int(md)
    sd = (md - m) * 60
    return 1000 * (d * 3600 + m * 60 + sd)


def parse_data(path):
    curves = []
    df = pd.read_csv(path, usecols=['POLYLINE'])
    repatter = re.compile(r'([\-0-9\.]+),\s*([\-0-9\.]+)')
    for i, (_, row) in enumerate(df.iterrows()):
        polyline = row['POLYLINE']
        results = repatter.findall(polyline)
        if not results:
            print(f'skip empty row of id {i}')
            continue
        points = []
        for result in results:
            lng = float(result[0])  # keido
            lat = float(result[1])  # ido
            points.append((deg_to_ms(lng), deg_to_ms(lat)))
        curves.append({'id': i, 'points': points})
    return curves


def write_curves(curves, path):
    with open(path, 'wb') as fp:
        msgpack.pack({'dimensions': '2'}, fp)
        for curve in curves:
            msgpack.pack(curve, fp)
    n_curves = len(curves)
    print(f'wrote {n_curves} curves')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_path')
    parser.add_argument('output_path')
    args = parser.parse_args()

    curves = parse_data(args.input_path)
    write_curves(curves, args.output_path)
