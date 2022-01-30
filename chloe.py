#!/usr/bin/env python3

import os
import urllib3
import json
import time
import tempfile
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def get_url_json(url):
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    if data := response.data.decode('utf-8'):
        return json.loads(data)


def download_url(url, filename):
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    with open(filename, 'wb') as f:
        f.write(response.data)


def identify(accession, leave_tmp_file=False, num_retries=3, num_progress=8):
    accession = os.path.basename(accession).split('.', 1)[0]  # To be sure, and to work *.gb
    num_tries = 0
    while num_tries < num_retries:
        if data := get_url_json(f'https://chloe.plantenergy.edu.au/annotate-ncbi?ncid={accession}&force_circular=true'):
            status = data['status'].upper()
            num_p = 0
            while status in ('PENDING', 'PROGRESS') and num_p < num_progress:
                time.sleep(3)
                if not (data := get_url_json('https://chloe.plantenergy.edu.au' + data['task_url'])):
                    continue
                status = data['status'].upper()
                num_p += 1

            if status != 'SUCCESS':
                num_tries += 1
                continue

            gb_filename = os.path.join(tempfile.gettempdir(), f'{accession}.gb')
            download_url('https://chloe.plantenergy.edu.au' + data['download_gbff'], gb_filename)
            irs = find_chloroplast_irs(SeqIO.read(gb_filename, 'genbank'))
            if not leave_tmp_file:
                os.remove(gb_filename)

            return irs


def find_chloroplast_irs(seq):
    rep_regs = [f for f in seq.features if f.type == 'repeat_region']

    if not rep_regs:
        return

    if len(rep_regs) == 1:
        # Chloe annotates IRs as one composite feature
        loc = rep_regs[0].location
        if isinstance(loc, CompoundLocation):
            strands = [[], []]  # +1, -1 strand
            for p in loc.parts:
                strands[int(p.strand < 0)].append(p)
            assert all(len(s) >= 1 for s in strands), strands
            ira = strands[0][0] if len(strands[0]) == 1 else CompoundLocation(strands[0])
            # Note: take a care about order of parts in strand -1
            #       join(83185..108842,complement(126069..151707),complement(1..19))
            irb = strands[1][0] if len(strands[1]) == 1 else CompoundLocation(strands[1][::-1])
            #
            ira = int(ira.parts[0].start), int(ira.parts[-1].end)
            irb = int(irb.parts[0].start), int(irb.parts[-1].end)
            diff_1 = (irb[0] - ira[1]) % len(seq)
            diff_2 = (ira[0] - irb[1]) % len(seq)
            return (ira, irb) if diff_1 < diff_2 else (irb, ira)
        return

    assert len(rep_regs) == 2, rep_regs

    # Repair features with location of type <~start..ira.start>
    half_length = len(seq) // 2
    for f in rep_regs:
        if len(f) >= half_length:
            loc = f.location
            f.location = CompoundLocation([FeatureLocation(loc.end, len(seq), strand=1),
                                           FeatureLocation(0, loc.start + 1, strand=1)])
        elif isinstance(f.location, CompoundLocation) and \
                len(f.location.parts) == 2 and \
                all(p.strand == -1 for p in f.location.parts) and \
                f.location.parts[1].start == 0:  # Note: BioPyhton changes order parts if location is complement!
            # Chloe can have annotation of format:
            #   repeat_region complement(join(1..24302,168220..173587))
            # It should be: complement(join(168220..173587,1..24302))
            f.location.parts = f.location.parts[::-1]

    ira, irb = rep_regs
    ira = int(ira.location.parts[0].start), int(ira.location.parts[-1].end)
    irb = int(irb.location.parts[0].start), int(irb.location.parts[-1].end)
    diff_1 = (irb[0] - ira[1]) % len(seq)
    diff_2 = (ira[0] - irb[1]) % len(seq)
    return (ira, irb) if diff_1 < diff_2 else (irb, ira)


def ir_loc(ir_parts):
    assert 1 <= len(ir_parts) <= 2, ir_parts
    p1 = ir_parts[0]
    if len(ir_parts) == 1:
        return [int(p1.start), int(p1.end)]
    if p1.start == 0:
        return [int(ir_parts[1].start), int(p1.end)]
    assert ir_parts[1].start == 0, ir_parts
    return [int(p1.start), int(ir_parts[1].end)]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
Fetch annotation from Chloe web application and extracts IRs location from the result.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('accession', help='Accession number')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')
    parser.add_argument('-R', '--num-retries', type=int, default=3, help='Maximal number of retries')
    parser.add_argument('-P', '--num-progress', type=int, default=8, help='Maximal number of progress checks')
    params = parser.parse_args()
    if irs := identify(params.accession,
                       leave_tmp_file=params.leave_tmp_file,
                       num_retries=params.num_retries,
                       num_progress=params.num_progress):
        print(*irs[0], *irs[1])
