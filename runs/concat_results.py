#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def main(*paths):
    first = True
    t_accu = 0.0
    for path in paths:
        with open(path) as fh:
            header = fh.readline() + fh.readline()
            if first:
                print(header, end='')
            for line in fh:
                spl = line.split()
                spl0 = spl[0]
                nspl0 = len(spl0)
                fspl0 = float(spl0)
                print(('%23.17e' % (fspl0 + t_accu)) + line[nspl0:], end='')
            t_accu += fspl0
        first = False


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
