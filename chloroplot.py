#!/usr/bin/env python3

import os.path
import subprocess
import tempfile

"""
Python wrapper method around chloroplot.R wrapper :-)
"""


def convert_to_fasta(seq_filename):
    # Note: seq_filename can be in fasta or genbank format
    if any(seq_filename.endswith(e) for e in ('.fa', '.fas', '.fasta', '.fs')):
        return seq_filename
    # Convert file into fasta format
    from Bio import SeqIO
    fa_filename = os.path.join(tempfile.gettempdir(), f'tmp_chloroplot.fa')
    SeqIO.convert(seq_filename, 'genbank', fa_filename, 'fasta')
    return fa_filename


def identify(seq_filename, print_chloroplot_output=False, leave_tmp_file=False):
    fa_filename = convert_to_fasta(seq_filename)

    _dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(_dir, 'chloroplot_resources', 'chloroplot.R')
    try:
        result = subprocess.run(['Rscript', r_script, fa_filename],
                                check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return

    output = result.stdout.decode('utf-8')
    if print_chloroplot_output:
        print(output)
    if fa_filename != seq_filename:
        if leave_tmp_file:
            print(f'Note: tmp file {fa_filename} is not removed!')
        else:
            os.remove(fa_filename)

    return parse_output(output.split('\n'))


def parse_output(lines):
    # Format of output is:
    # $ir_table
    #    chr  start    end name       text center
    # 1 chr1      0  81884  LSC LSC: 81934  40917
    # 2 chr1  81884 107546  IRA IRA: 25662  94715
    # 3 chr1 107546 124801  SSC SSC: 17255 116174
    # 4 chr1 124801 150460  IRB IRB: 25659 137630
    # 5 chr1 150460 150510  LSC            150485
    #
    # $indel_table
    #    mismatch_type position string    col
    # 1         insert    82161      C  green
    # 2         insert    82162      C  green
    # 3         insert    82163      C  green
    # 4         delete   150187      D yellow
    # 5        replace   104414      G    red
    # 6        replace   104420      G    red
    # 7        replace   104421      C    red
    # ...

    iras, irbs = [], []
    seq_length = 0
    for line in lines:
        fields = line.split()
        if len(fields) < 5:
            continue
        if fields[4] in ('IRA', 'IRB', 'LSC', 'SSC'):
            seq_length = int(fields[3])
            if fields[4] == 'IRA':
                iras.append((int(fields[2]), int(fields[3])))
            if fields[4] == 'IRB':
                irbs.append((int(fields[2]), int(fields[3])))
    if iras and irbs:
        ira, irb = _loc(iras, seq_length), _loc(irbs, seq_length)
        if ira and irb:
            return (ira, irb) if ((irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length) else (irb, ira)


def _loc(ir_parts, seq_length):
    assert 1 <= len(ir_parts) <= 2, ir_parts
    if len(ir_parts) == 1:
        return _loc_check(*ir_parts[0], seq_length)
    # Of type:
    #    chr  start    end name       text center
    # 1 chr1      0      3  IRB IRB: 24645 139014
    # 2 chr1      3  83717  LSC LSC: 83714  41860
    # 3 chr1  83717 108362  IRA IRA: 24645  96040
    # 4 chr1 108362 126691  SSC SSC: 18329 117526
    # 5 chr1 126691 151333  IRB            139012

    if ir_parts[0][1] < 0:
        # ??? NC_028349
        #   chr  start    end name       text center
        # 1 chr1      0   -242  IRB IRB: 25926 146143
        # 2 chr1   -242  88174  LSC LSC: 88416  43966
        # 3 chr1  88174 114236  IRA IRA: 26062 101205
        # 4 chr1 114236 133179  SSC SSC: 18943 123708
        # 5 chr1 133179 159347  IRB            146263
        return _loc_check(ir_parts[1][0], (ir_parts[1][1] + ir_parts[0][1]), seq_length)

    if ir_parts[-1][0] > ir_parts[-1][1]:
        # ??? NC_057051
        #    chr  start    end name        text center
        # 1 chr1      0   8826  IRB    IRB: 106   8773
        # 2 chr1   8826  14076  SSC   SSC: 5250  11451
        # 3 chr1  14076  14182  IRA    IRA: 106  14129
        # 4 chr1  14182 173226  LSC LSC: 159044  93704
        # 5 chr1 173226 164506  IRB             168866
        return tuple(ir_parts[0])

    return _loc_check(ir_parts[1][0], ir_parts[0][1], seq_length)


def _loc_check(start, end, seq_length):
    if (end - start) % seq_length < seq_length // 2:
        return start, end


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Run Chlorolpot IR detection on sequence data and extracts IRs location from the result.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='Sequence filename, in fasta or GenBank format')
    parser.add_argument('-C', '--print-chloroplot-output', action='store_true', help='Print Chloroplot R output.')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')
    params = parser.parse_args()
    if irs := identify(params.filename,
                       print_chloroplot_output=params.print_chloroplot_output,
                       leave_tmp_file=params.leave_tmp_file):
        print(*irs[0], *irs[1])
