#!/usr/bin/env python3
"""state_follow.py -- follow IQ-TREE's .state while it is written and keep only the alignment columns.

With -fconst, IQ-TREE writes for every internal node the alignment columns followed by every constant site (millions of rows per node). This script runs beside the --ancestral run, reads the growing file and writes the comment / header lines and every row with Site <= the number of alignment columns to --out, so the huge raw file can be deleted afterwards.
Standard library only.
"""
import argparse
import os
import sys
import time


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--path', required=True, help='.state written by IQ-TREE')
    ap.add_argument('--pid', type=int, required=True, help='PID of the IQ-TREE --ancestral process')
    ap.add_argument('--columns', required=True, help='L1 columns table (header + one row per alignment column)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    with open(a.columns) as fh:
        ncols = sum(1 for _ in fh) - 1
    print(f'[{time.strftime("%F %T")}] >>> .state: keeping the rows of the {ncols} alignment columns of every internal node', flush=True)

    while not os.path.exists(a.path):
        if not alive(a.pid):
            time.sleep(5)
            if os.path.exists(a.path):
                break
            sys.exit(f'[ERROR] IQ-TREE exited without writing {a.path}')
        time.sleep(5)

    buf = b''
    with open(a.path, 'rb') as fin, open(a.out, 'wb') as fout:
        writer_done = False
        while True:
            chunk = fin.read(16 << 20)
            if not chunk:
                if writer_done:
                    break
                if not alive(a.pid):
                    writer_done = True     # one more pass drains anything written just before exit
                else:
                    time.sleep(2)
                continue
            buf += chunk
            lines = buf.split(b'\n')
            buf = lines.pop()
            for line in lines:
                if not line or line[0] == 35 or line.startswith(b'Node\t'):   # 35 = '#'
                    fout.write(line + b'\n')
                    continue
                f = line.split(b'\t', 2)
                if int(f[1]) <= ncols:
                    fout.write(line + b'\n')
        if buf.strip():
            sys.exit(f'[ERROR] {a.path} ends with an incomplete line')
    print(f'[{time.strftime("%F %T")}] >>> .state read; alignment-column rows written to {a.out}', flush=True)


if __name__ == '__main__':
    main()
