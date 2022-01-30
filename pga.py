#!/usr/bin/env python3

import sys
import os.path
import shutil
import tempfile
import subprocess
from Bio import SeqIO

"""
Note:
Script sometimes crashes on working with temp files.
Something like it removes temp files before calling blast.
Workaround: repeat a call it more times of that kind of error.
"""


def remove_directory(_dir, create):
    # assert not _dir.endswith("/") or _dir.endswith("\\"), _dir
    _dir = os.path.normpath(_dir)

    if os.path.isdir(_dir):
        if sys.platform == "win32":
            temp_path = _dir + "_"

            if os.path.exists(temp_path):
                remove_directory(temp_path, False)

            try:
                os.renames(_dir, temp_path)
            except OSError as exception:
                if exception.errno != errno.ENOENT:
                    raise
            else:
                shutil.rmtree(temp_path)
        else:
            shutil.rmtree(_dir)

    if create:
        os.makedirs(_dir)


def _loc(ir, seq_len):
    ir = ir.location.parts
    assert len(ir) <= 2, ir
    if ir[0].strand > 0:
        return _loc_check(int(ir[0].start), int(ir[-1].end), seq_len)
    return _loc_check(int(ir[-1].start), int(ir[0].end), seq_len)


def _loc_check(start, end, seq_length):
    # Some (probably only) complement regions can be wrongly oriented
    if (end - start) % seq_length < seq_length // 2:
        return start, end


def identify(seq_filename, leave_tmp_file=False):
    # Tmp directories
    tmp_d = tempfile.gettempdir()
    reference = os.path.join(tmp_d, 'pga', 'reference')
    target = os.path.join(tmp_d, 'pga', 'target')
    out = os.path.join(tmp_d, 'pga', 'out')
    #
    for d in (reference, target, out):
        remove_directory(d, True)

    # Remove reference features to make annotation faster
    seq_type = 'fasta' if any(seq_filename.endswith(e) for e in ('.fa', '.fas', '.fasta', '.fs')) else 'genbank'
    seq_rec = SeqIO.read(seq_filename, seq_type)
    seq_rec.annotations['molecule_type'] = 'DNA'
    seq_rec.features = []  # seq_rec.features[:8]
    SeqIO.write([seq_rec], os.path.join(reference, os.path.basename(seq_filename)), 'genbank')
    SeqIO.convert(seq_filename, 'genbank', os.path.join(target, 'x.fa'), 'fasta')

    res_irs = None

    # Run script
    # Note: script can randomly crash. It is run more times and outputs are checked!
    for _ in range(5):
        with tempfile.TemporaryFile() as stderr:
            subprocess.run(['PGA.pl', '-r', reference, '-t', target, '-o', out],
                           check=True,
                           stdout=subprocess.DEVNULL, stderr=stderr)
            stderr.seek(0)
            err = stderr.read().decode('utf-8')

        # Read IRs
        seq = SeqIO.read(os.path.join(out, 'x.gb'), 'genbank')
        if rep_regs := [f for f in seq.features
                        if f.type == 'repeat_region' and f.qualifiers.get('rpt_type', ('inverted',))[0] == 'inverted']:
            if len(rep_regs) == 2:
                ira, irb = rep_regs
                seq_len = len(seq.seq)
                if (ira := _loc(ira, seq_len)) and (irb := _loc(irb, seq_len)):
                    res_irs = ira, irb
                break

        if 'readline() on closed filehandle' not in err:
            break

    if not leave_tmp_file:
        for d in (reference, target, out):
            remove_directory(d, False)

    return res_irs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Runs PGA annotaiton tool and extracts IRs location from the result.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='Sequence filename, in fasta or GenBank format')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')

    params = parser.parse_args()
    if irs := identify(params.filename, leave_tmp_file=params.leave_tmp_file):
        print(*irs[0], *irs[1])
