#!/usr/bin/env python3

import numpy as np
import os

def main(r1="results1.txt", r2="results2.txt", outfile='delta.txt'):
    if os.path.exists(outfile):
        raise ValueError("files exists: %s" % outfile)
    with open(r1) as ifh:
        h1 = ifh.readline()
        info = ifh.readline()[1:].replace(' ', '').replace(',,', ',')
        info = eval('dict(%s)' % info)
    with open(r2) as ifh:
        h2 = ifh.readline()
        if h1 != h2:
            raise ValueError("headers differs")
    d1 = np.genfromtxt(r1)
    d2 = np.genfromtxt(r2)
    if not np.all(d1[:, 0] == d2[:, 0]):
        raise ValueError("Not identical output times")
    do = d1 - d2
    do[:, 0] = d1[:, 0]
    header = h1.lstrip('#') + '%s' % ', '.join('%s=%s' % kv for kv in info.items())
    np.savetxt(outfile, do, header=header)

if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
