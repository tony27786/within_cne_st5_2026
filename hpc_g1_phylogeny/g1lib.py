"""Shared helpers for the G1/G2 scripts (standard library; numpy only where noted).

BED intervals are 0-based half-open and VCF POS is 1-based; sequences are handled upper-case (the reference is soft-masked).
"""
import bisect
import csv
import gzip
import os
import sys

BASES = b'ACGT'


def die(msg, code=1):
    print(f'[ERROR] {msg}', file=sys.stderr, flush=True)
    sys.exit(code)


def log(msg):
    print(msg, flush=True)


def open_text(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode.replace('t', '') if 'b' in mode else mode)


def read_tsv(path):
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, cols):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_samples(path):
    """Sample table written by write_sample_table (sample, group, role, callset), file order."""
    rows = read_tsv(path)
    for k in ('sample', 'group', 'role', 'callset'):
        if rows and k not in rows[0]:
            die(f'{path}: missing column {k}')
    return rows


def read_genome(path):
    """contig<TAB>length (genome.txt or a .fai); returns [(contig, length)] in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2 and f[0]:
                out.append((f[0], int(f[1])))
    return out


def read_names(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


# BED intervals

def merge(ivs):
    ivs = sorted(ivs)
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


class Intervals:
    """Merged interval sets per contig with O(log n) membership."""

    def __init__(self, by_contig=None):
        self.starts, self.ends = {}, {}
        for c, ivs in (by_contig or {}).items():
            m = merge(ivs)
            self.starts[c] = [s for s, _ in m]
            self.ends[c] = [e for _, e in m]

    @classmethod
    def from_bed(cls, path):
        d = {}
        if path and path != 'none':
            with open_text(path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith(('#', 'track', 'browser')):
                        continue
                    f = line.split('\t')
                    d.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        return cls(d)

    def contains(self, contig, pos0):
        st = self.starts.get(contig)
        if not st:
            return False
        i = bisect.bisect_right(st, pos0) - 1
        return i >= 0 and pos0 < self.ends[contig][i]

    def total(self):
        return sum(e - s for c in self.starts for s, e in zip(self.starts[c], self.ends[c]))

    def intervals(self, contig):
        return list(zip(self.starts.get(contig, []), self.ends.get(contig, [])))

    def contigs(self):
        return list(self.starts)


def write_bed(path, rows):
    """rows: iterable of (contig, start, end[, extra...]); written as given."""
    with open(path, 'w') as fh:
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')


# FASTA with .fai

class Faidx:
    """Whole-contig access to a FASTA indexed by samtools faidx; returns upper-case bytes without newlines."""

    def __init__(self, path):
        self.path = path
        self.index = {}
        fai = path + '.fai'
        if not os.path.exists(fai):
            die(f'missing index {fai} (samtools faidx)')
        with open(fai) as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split('\t')[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))

    def length(self, name):
        return self.index[name][0]

    def fetch(self, name):
        length, offset, lb, lw = self.index[name]
        nlines = (length + lb - 1) // lb if lb else 0
        nbytes = length + nlines * (lw - lb)
        with open(self.path, 'rb') as fh:
            fh.seek(offset)
            raw = fh.read(nbytes)
        seq = raw.replace(b'\n', b'').replace(b'\r', b'')[:length].upper()
        if len(seq) != length:
            die(f'{self.path}: contig {name} has {len(seq)} bases, index says {length}')
        return seq


def fasta_lengths(path):
    """{contig: length} from a .fai."""
    out = {}
    with open(path + '.fai') as fh:
        for line in fh:
            f = line.split('\t')
            out[f[0]] = int(f[1])
    return out


def read_fasta_alignment(path):
    """Small/medium alignment FASTA -> (names in order, {name: bytes upper})."""
    names, seqs, cur, buf = [], {}, None, []
    with open_text(path, 'rb') as fh:
        for line in fh:
            line = line.rstrip(b'\r\n')
            if line.startswith(b'>'):
                if cur is not None:
                    seqs[cur] = b''.join(buf).upper()
                cur = line[1:].split()[0].decode()
                names.append(cur)
                buf = []
            else:
                buf.append(line)
    if cur is not None:
        seqs[cur] = b''.join(buf).upper()
    return names, seqs


def write_fasta(path, names, seqs, width=0):
    with open(path, 'wb') as fh:
        for n in names:
            s = seqs[n]
            fh.write(b'>' + n.encode() + b'\n')
            if width:
                for i in range(0, len(s), width):
                    fh.write(s[i:i + width] + b'\n')
            else:
                fh.write(s + b'\n')


# VCF

def vcf_positions(path):
    """[(contig, pos1)] of every record, in file order."""
    out = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 2)
            out.append((f[0], int(f[1])))
    return out


def gt_allele(gt):
    """Haploid GT string -> allele index (int) or None when missing; diploid-looking calls are rejected."""
    if gt in ('.', './.', '.|.', ''):
        return None
    if '/' in gt or '|' in gt:
        die(f'non-haploid genotype {gt!r}')
    return int(gt)


def pair_key(a, b):
    return (a, b) if a <= b else (b, a)
