#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def find_chloroplast_irs(seq):
    # Finds the longest pair of inverted repeats
    _ir = [['inverted'], ('inverted',)]
    rep_regs = [f for f in seq.features
                if f.type == 'repeat_region' and f.qualifiers.get('rpt_type') in _ir]

    if not rep_regs:
        return

    assert len(rep_regs) == 2, rep_regs
    if any(not f.location for f in rep_regs):
        # Some BioPython versions can have problem with detecting circular genomes
        return

    # Repair features with location of type <~start..ira.start>
    half_length = len(seq) // 2
    for f in rep_regs:
        if len(f) >= half_length:
            loc = f.location
            f.location = CompoundLocation([FeatureLocation(loc.end, len(seq), strand=1),
                                           FeatureLocation(0, loc.start + 1, strand=1)])

    # Check order of regions
    ira, irb = rep_regs
    ira = int(ira.location.parts[0].start), int(ira.location.parts[-1].end)
    irb = int(irb.location.parts[0].start), int(irb.location.parts[-1].end)
    diff_1 = (irb[0] - ira[1]) % len(seq)
    diff_2 = (ira[0] - irb[1]) % len(seq)
    return (ira, irb) if diff_1 < diff_2 else (irb, ira)


def identify(filename):
    return find_chloroplast_irs(SeqIO.read(filename, 'genbank'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Extracts IRs location from GeSeq annotation.

Sequence has to be annotated with GeSeq web application:
https://chlorobox.mpimp-golm.mpg.de/geseq.html""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='Sequence filename, in GenBank format')
    params = parser.parse_args()
    if irs := identify(params.filename):
        print(*irs[0], *irs[1])
