#!/usr/bin/env python3

import sys
import os.path
import shutil
import tempfile
import subprocess

# Notes:
#  * It is not possible to run more instances of this script since all of them
#    would use same tmp directory: join(tempfile.gettempdir(), 'plann')!
#  * Script uses reference genome, in genbank format, to annotate genes.
#    That genome has to have at least one annotated gene, since script crashes if blast result does not exist.


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


def convert_to_fasta(seq_filename, fa_filename):
    # Note: seq_filename can be in fasta or genbank format
    if any(seq_filename.endswith(e) for e in ('.fa', '.fas', '.fasta', '.fs')):
        # ToDo: possible speedup, remove features except one, or create dummy one as in convert_to_genbank()
        shutil.copy(seq_filename, fa_filename)
    else:
        # Convert file into fasta format
        from Bio import SeqIO
        SeqIO.convert(seq_filename, 'genbank', fa_filename, 'fasta')


def convert_to_genbank(seq_filename, gb_filename):
    # Note: seq_filename can be in fasta or genbank format
    if any(seq_filename.endswith(e) for e in ('.gb', '.genbank', '.gbs')):
        shutil.copy(seq_filename, gb_filename)
    else:
        # Convert file into fasta format
        from Bio import SeqIO, SeqFeature
        from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
        seq_rec = SeqIO.read(seq_filename, 'fasta')
        seq_rec.annotations['molecule_type'] = 'DNA'

        # Add some dummy gene to prevent crash because missing blast results
        loc = FeatureLocation(ExactPosition(36), ExactPosition(137))
        seq_rec.features.append(SeqFeature(loc, type='gene', qualifiers=dict(gene='rps19')))

        SeqIO.write([seq_rec], gb_filename, 'genbank')


def identify(seq_filename, leave_tmp_file=False):
    # Find script location
    if not (plann_script := os.environ.get('PLANN_SCRIPT')):
        # Chech is plann.pl executable and on the PATH
        if p_pl := shutil.which('plann.pl'):
            plann_script = os.path.realpath(p_pl)

    if not plann_script:
        raise OSError('Plann script not found! Set PLANN_SCRIPT environment command.')

    # Tmp directories
    tmp_plann = os.path.join(tempfile.gettempdir(), 'plann')
    remove_directory(tmp_plann, True)

    fa_filename = os.path.join(tmp_plann, 'x.fa')
    gb_filename = os.path.join(tmp_plann, 'x.gb')
    convert_to_fasta(seq_filename, fa_filename)
    convert_to_genbank(seq_filename, gb_filename)

    res_irs = None
    try:
        cmd = ['perl', plann_script, '-reference', gb_filename, '-fasta', fa_filename, '-out', 'plann_result']
        result = subprocess.run(cmd, cwd=tmp_plann, check=True,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        tbl_filename = os.path.join(tmp_plann, 'plann_result.tbl')
        if os.path.isfile(tbl_filename):
            with open(tbl_filename, 'r') as _in:
                in_repeat = False
                irs = []
                for line in _in.readlines():
                    if line[0].startswith('>'):
                        continue
                    if not line[0].isspace():  # line.startswith('\t'):
                        # Feature start. 2-3 fields, first two are location
                        fs = line.split()
                        in_repeat = (len(fs) == 3) and (fs[2] == 'repeat_region')
                        if in_repeat:
                            irs.append((int(fs[0]), int(fs[1])))
                    elif in_repeat:
                        fs = line.split()
                        if len(fs) == 2 and fs[0] == 'rpt_type' and fs[1] != 'inverted':
                            irs.pop()
                            in_repeat = False  # To make sure only one repeat is removed
                #
                if len(irs) == 2:
                    ira, irb = irs
                    # Our indexing
                    ira = (ira[0] - 1, ira[1])
                    irb = (irb[0] - 1, irb[1])
                    # Find sequence length
                    with open(os.path.join(tmp_plann, 'plann_result.fsa'), 'r') as _fsa:
                        seq_length = len(list(_fsa.readlines())[1]) - 1
                    # Check IR order
                    res_irs = (ira, irb) if (irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length else (irb, ira)

    except subprocess.CalledProcessError:
        print(seq_filename)
        raise

    if not leave_tmp_file:
        remove_directory(tmp_plann, False)

    return res_irs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Runs plann annotaiton tool and extracts IRs location from the result.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='Sequence filename, in fasta or GenBank format')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')
    params = parser.parse_args()
    if irs := identify(params.filename, leave_tmp_file=params.leave_tmp_file):
        print(*irs[0], *irs[1])
