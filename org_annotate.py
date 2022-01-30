#!/usr/bin/env python3

import sys
import os.path
import shutil
import tempfile
import subprocess
from Bio import SeqIO

"""
Implementation (copy) of ORG.Annotate IR detection method.

ORG.Annotate repository: https://git.metabarcoding.org/org-asm/org-annotate
IR detection is implemented by scripts in detectors/normalize/lib directory.

Main script is lookforIR.lib.sh.
Detection is done in 3 steps.
 1. blast sequence to database of LSC and SSC sequences to find possible IR positions.
    Database is located in directory data/ir.
    Database contains 45 LSC and 72 SSC. Max accession number is NC_026103
 2. repseek on sequence.
    repseek is located in directory src/repseek(/repseek-2014.09)
 3. Choose IR location between blasted data and repseek results.
    Done with script selectIR.py.

Note: (optimization) implementation swaps steps 1 and 2, since process is terminated if repseek don't find any IR.
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


def convert_to_fasta(seq_filename, tmp_d):
    # Note: seq_filename can be in fasta or genbank format
    if any(seq_filename.endswith(e) for e in ('.fa', '.fas', '.fasta', '.fs')):
        return seq_filename
    # Convert file into fasta format
    from Bio import SeqIO
    fa_filename = os.path.join(tmp_d, 'query.fa')
    SeqIO.convert(seq_filename, 'genbank', fa_filename, 'fasta')
    return fa_filename


def identify(seq_filename, print_repseek_output=False, print_blast_output=False, leave_tmp_file=False):
    # Locate repseek executable
    if not (repseek_exe := os.environ.get('REPSEEK_EXE')):
        # Chech is repseek executable and on the PATH
        if not (repseek_exe := shutil.which('repseek')):
            raise OSError('Executable repseek not found! Set REPSEEK_EXE environment command.')

    # Filenames
    tmp_d = os.path.join(tempfile.gettempdir(), 'org_annotate')
    remove_directory(tmp_d, True)
    repseek_results = os.path.join(tmp_d, 'repseek_results')
    blastn_results = os.path.join(tmp_d, 'blastn_results')
    oa_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ORGAnnotate_resources')

    # Create fasta file
    fa_filename = convert_to_fasta(seq_filename, tmp_d)

    # 1. repseek
    # repseek -c -p 0.001 -i ${QUERY} 2>> /dev/null > ${REPEATS}
    # nrepeat="$(wc -l ${REPEATS} | awk '{print $1}')"
    # nrepeat != 0
    try:
        result = subprocess.run([repseek_exe, '-c', '-p', '0.001', '-1', fa_filename],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        raise
        # return
    output = result.stdout.decode('utf-8')
    if print_repseek_output:
        print(output)
    if not output:
        return
    with open(repseek_results, 'w') as _out:
        _out.write(output)

    # 2. blastn
    # (comment)
    # Blast columns:
    #   query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    # We keep blast matches if :
    #   The match is longer than 1000
    #   The identity is higher than 80%
    #
    # The match file has the following format:
    #   LSC/SSC  begin end  same_strand=1/diff_strand=0
    #
    # (code)
    # blastn -db ${SCDB} \
    #        -query ${QUERY} \
    #        -outfmt 6 \
    #        -max_target_seqs 10000 | \
    #   awk '($4 > 100) && ($3>80) { \
    #              SAME=(($7 < $8) && ($9 < $10)) || (($7 > $8) && ($9 > $10)); \
    #              if ($7 < $8) \
    #                 {print substr($2,1,3),$7,$8,SAME}  \
    #              else \
    #                 {print substr($2,1,3),$8,$7,SAME}}' | \
    #   sort -nk 2 > ${MATCHES}
    try:
        result = subprocess.run(['blastn', '-db', os.path.join(oa_dir, 'SC_RefDB'),
                                 '-query', fa_filename,
                                 '-outfmt', '6',
                                 '-max_target_seqs', '10000'],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        raise
        # return
    output = result.stdout.decode('utf-8')
    if print_blast_output:
        print(output)
    res = []
    for line in output.split('\n'):
        fs = line.split()
        if fs and int(fs[3]) > 100 and float(fs[2]) > 80:
            f7, f8, f9, f10 = map(int, fs[6:10])
            same = int(((f7 < f8) and (f9 < f10)) or ((f7 > f8) and (f9 > f10)))
            res.append((fs[1][:3], f7, f8, same) if f7 < f8 else (fs[1][:3], f8, f7, same))
    with open(blastn_results, 'w') as _out:
        _out.write('\n'.join(' '.join(map(str, fs)) for fs in sorted(res, key=lambda x: x[1])))

    # 3.
    # SELECTIR = detectors/normalize/lib/selectIR.py
    # IR=( $(${SELECTIR} ${MATCHES} ${REPEATS}) )
    try:
        result = subprocess.run([os.path.join(oa_dir, 'selectIR.py'), blastn_results, repseek_results],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return
    output = result.stdout.decode('utf-8')
    res_irs = None
    if output:
        fs = output.split('\n')[0].split()
        ira_s, ira_l, irb_s, irb_l = map(int, fs[4:8])
        ira_s -= 1
        irb_s -= 1
        res_irs = (ira_s, ira_s + ira_l), (irb_s, irb_s + irb_l)

    if not leave_tmp_file:
        remove_directory(tmp_d, False)

    return res_irs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Run ORG.Annotate IR detection on sequence data and extracts IRs location from the result.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='Sequence filename, in fasta or GenBank format')
    parser.add_argument('-R', '--print-repseek-output', action='store_true', help='Print repseek output.')
    parser.add_argument('-B', '--print-blast-output', action='store_true', help='Print Blastn output.')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')
    params = parser.parse_args()
    if irs := identify(params.filename,
                       print_repseek_output=params.print_repseek_output,
                       print_blast_output=params.print_blast_output,
                       leave_tmp_file=params.leave_tmp_file):
        print(*irs[0], *irs[1])
