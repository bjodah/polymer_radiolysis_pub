#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json


def main(result_file):
    with open(result_file) as fh:
        lines = fh.read().splitlines()
        keys = lines[0].strip('#').split()[1:]
        values = map(float, lines[-1].split()[1:])
    print(json.dumps({k: v for k, v in zip(keys, values)}), end='')

if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
