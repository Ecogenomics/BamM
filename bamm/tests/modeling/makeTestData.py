#!/usr/bin/env python
###############################################################################
#
# makeTestReads.py - make reads for use in unit testing
#
###############################################################################
#                                                                             #
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014,2015"
__credits__ = ["Michael Imelfort"]
__license__ = "LGPLv3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"

###############################################################################

# system imports
import argparse
import sys
import os

import numpy as np
np.seterr(all='raise')

import json

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# these numbers are to do with the scaling used in the model diagram.
# messy but too deep to remove now
mult = 25
readLen = 3 * mult

types = ['PE', 'MP', 'UP']
cov_types = ['counts', 'cmean', 'pmean', 'opmean', 'tpmean', 'pmedian']
ext_flags = {'mix_bams': [True, False],
             'mix_groups': [True, False],
             'mix_reads': [True, False]}

# len/mult, seq, pileup, starts
contigs = {'A':[17, '', {'PE':[], 'MP':[], 'UP':[]}, {'PE':0, 'MP':0, 'UP':0}],
           'B':[15, '', {'PE':[], 'MP':[], 'UP':[]}, {'PE':0, 'MP':0, 'UP':0}],
           'C':[17, '', {'PE':[], 'MP':[], 'UP':[]}, {'PE':0, 'MP':0, 'UP':0}],
           'E':[16, '', {'PE':[], 'MP':[], 'UP':[]}, {'PE':0, 'MP':0, 'UP':0}],
           'Z':[17, '', {'PE':[], 'MP':[], 'UP':[]}, {'PE':0, 'MP':0, 'UP':0}]}
           
           
# the reads we make
pe= []
mp = []
up = []

# work out what the coverages output should look like
coverages = {}
for ct in cov_types:
    coverages[ct] = {'A':{'PE':0., 'MP':0., 'UP':0.},
                     'B':{'PE':0., 'MP':0., 'UP':0.},
                     'C':{'PE':0., 'MP':0., 'UP':0.},
                     'E':{'PE':0., 'MP':0., 'UP':0.}}

# work out what the links output should look like
links = []

# bin information
cid_2_grp = {'Z':None}
cid_2_grp_comb = {'Z':None}

# predict read titles
reads = {}
read_tracker = {}

# important filenames
contigs_filename = "contigs.fa"
results_filename = "predicted_outputs.json"
out_fnames = {'PE': ["pe.1.fa", "pe.2.fa"], # coupled
              'MP': "mp.fa",                # interleaved
              'UP': "up.fa"}                # singles

bam_fnames = {'PE': "contigs.pe.1",         # the .1 remains
              'MP': "contigs.mp",
              'UP': "contigs.up"}
              
bad_contigs_filename = "contigs_bad.fa"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def revcom(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)[::-1]

def cutMers(seq, k):
    mers = {}
    maxx = len(seq) - k + 1
    for i in range(maxx):
        mer = seq[i:i+k]
        rmer = revcom(mer)
        mers[mer] = 1
        mers[rmer] = 1
    return (mers.keys(), maxx*2)

def getAndCheckSeq(seq, s_len, c_len, mers):
    tmp_mers = []
    tmp_seq = ""
    maxx = s_len - c_len
    done = False
    start = 0
    while not done:
        start = np.random.randint(maxx)
        tmp_seq = seq[start: start+c_len]
        (tmp_mers, made_mers) = cutMers(tmp_seq, readLen)
        if len(tmp_mers) == made_mers:
            # no duplicates within the contig
            found_mer = False
            for mer in tmp_mers:
                if mer in mers:
                    found_mer = True
                    break
            if not found_mer:
                done = True
    return tmp_seq

def squishLink(link):
    return ",".join([str(i) for i in link])

def doWork( args ):
    """ Main wrapper"""
    # load groups
    for group_name in args.groups:
        try:
            with open(group_name, "r") as fh:
                for line in fh:
                    cid = line.rstrip()
                    # make sure to remove the path first
                    cid_2_grp[cid] = os.path.split(group_name)[1]
                    cid_2_grp_comb[cid] = "allGroups"
        except:
            print "Error opening file:", group_name, sys.exc_info()[0]
            raise
        

    # parse fasta file
    seq = []
    import mimetypes
    all_open = open
    try:
        # handle gzipped files
        mime = mimetypes.guess_type(args.fasta)
        if mime[1] == 'gzip':
            import gzip
            all_open = gzip.open
    except:
        print "Error when guessing refs file mimetype"
        raise
    try:
        with all_open(args.fasta, "r") as fh:
            for line in fh:
                if line[0] == '>':
                    continue
                seq.append(line.rstrip())
    except:
        print "Error opening file:", args.fasta, sys.exc_info()[0]
        raise
    seq = "".join(seq)
    s_len = len(seq)

    # now make the contigs
    mers = {}
    for cid in contigs.keys():
        c_len = contigs[cid][0]*mult
        contigs[cid][1] = getAndCheckSeq(seq, s_len, c_len, mers)
        for type in types:
            contigs[cid][2][type] = np.array([0.]*c_len)

    # parse readkey file
    try:
        with open(args.readkey, "r") as fh:
            for line in fh:
                if line[0] == "#" or len(line) == 0:
                    continue
                [rid, con_1, pos1, con_2, pos2, type] = line.rstrip().split("\t")
                pos1 = int(pos1)*mult
                pos2 = int(pos2)*mult
                reads[rid] = [con_1, pos1, con_2, pos2, type]
                c1_seq = contigs[con_1][1][pos1:pos1+readLen]
                c2_seq = contigs[con_2][1][pos2:pos2+readLen]

                if type == 'PE':
                    c2_seq = revcom(c2_seq)
                    pe.append((rid, c1_seq, c2_seq))
                    rev = [0, 1]
                    header_f_string1 = "%s_PE"
                    header_f_string2 = "%s_PE"
                    fname = bam_fnames['PE']
                elif type == 'MP':
                    c1_seq = revcom(c1_seq)
                    mp.append((rid, c1_seq, c2_seq))
                    rev = [1, 0]
                    header_f_string1 = "%s_MP"
                    header_f_string2 = "%s_MP"
                    fname = bam_fnames['MP']
                else:
                    up.append((rid, c1_seq))
                    header_f_string1 = "%s_UP"


                # update starts and pileups here
                contigs[con_1][2][type][pos1:pos1+readLen] += 1
                contigs[con_1][3][type] += 1
                # one read for unpaired
                if type != 'UP':

                    contigs[con_2][2][type][pos2:pos2+readLen] += 1
                    contigs[con_2][3][type] += 1

                    # now links
                    if con_1 != con_2 and con_1 != 'Z' and con_2 != 'Z':
                        if con_1 < con_2:
                            link = [con_1,
                                    con_2,
                                    contigs[con_1][0]*mult,
                                    pos1,
                                    rev[0],
                                    contigs[con_2][0]*mult,
                                    pos2,
                                    rev[1],
                                    fname]
                        else:
                            link = [con_2,
                                    con_1,
                                    contigs[con_2][0]*mult,
                                    pos2,
                                    rev[1],
                                    contigs[con_1][0]*mult,
                                    pos1,
                                    rev[0],
                                    fname]
                        links.append(link)
    except:
        print "Error opening file:", args.readkey, sys.exc_info()[0]
        raise

#----------------------------
# Predict outputs

    # counts
    cov_type = 'counts'
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        for type in types:
            coverages[cov_type][cid][type] = contigs[cid][3][type]

    # cmean
    cov_type = 'cmean'
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        for type in types:
            coverages[cov_type][cid][type] = float(contigs[cid][3][type]) / \
                                             float(contigs[cid][0]*mult)
    # pmean
    cov_type = 'pmean'
    pmean_cov = []
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        for type in types:
            coverages[cov_type][cid][type] = np.mean(contigs[cid][2][type])

    # opmean
    cov_type = 'opmean'
    opmean_cov = []
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        c_len = contigs[cid][0]*mult
        for type in types:
            m = np.mean(contigs[cid][2][type])
            s = np.std(contigs[cid][2][type])
            u = m + s
            l = m - s
            L = 0
            S = 0
            for c in range(c_len):
                val = contigs[cid][2][type][c]
                if val >= l and val <= u:
                    L += 1
                    S += val
            if L > 0:
                coverages[cov_type][cid][type] = float(S)/float(L)
            else:
                coverages[cov_type][cid][type] = 0.

    # tpmean
    cov_type = 'tpmean'
    tpmean_cov = []
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        c_len = contigs[cid][0]*mult
        for type in types:
            from_end = int(c_len/10)
            L = 0
            S = 0
            contigs[cid][2][type].sort()
            for c in range(from_end, c_len-from_end):
                L += 1
                S += contigs[cid][2][type][c]
            if L > 0:
                coverages[cov_type][cid][type] = float(S)/float(L)
            else:
                coverages[cov_type][cid][type] = 0.

    # pmedian
    cov_type = 'pmedian'
    pmedian_cov = []
    for cid in contigs.keys():
        if cid == 'Z':
            continue
        for type in types:
            coverages[cov_type][cid][type] = np.median(contigs[cid][2][type])

    # read Ids
    # the desired output is:
    # filename    read_id    pair info
    # so we need to determine filenames and extensions
    for i_opt in ['', '--interleave']:
        read_tracker[i_opt] = {}
        for b_opt in ext_flags['mix_bams']:
            read_tracker[i_opt][b_opt] = {}
            for r_opt in ext_flags['mix_reads']:
                read_tracker[i_opt][b_opt][r_opt] = {}
                for g_opt in ext_flags['mix_groups']:
                    read_tracker[i_opt][b_opt][r_opt][g_opt] = {}
                    if b_opt:
                        b_suffix = {'PE':"allMapped",
                                    'MP':"allMapped",
                                    'UP':"allMapped"}
                    else:
                        b_suffix = bam_fnames

                    if r_opt:
                        r_suffix = {1:".allReads",
                                    2:".allReads",
                                    0:".allReads"}
                    else:
                        if i_opt == '':
                            r_suffix = {1:".1",
                                        2:".2",
                                        0:".unpairedReads"}
                        else:
                            r_suffix = {1:".pairedReads",
                                        2:".pairedReads",
                                        0:".unpairedReads"}


                    tmp_cid_2_grp = cid_2_grp
                    if g_opt:
                        tmp_cid_2_grp = cid_2_grp_comb

                    for read in reads:
                        fn_2 = ''
                        fn_1 = ''
                        pair_str_1 = ''
                        pair_str_2 = ''
                        hdr_1 = ''
                        hdr_2 = ''
                        con_1 = reads[read][0]
                        con_2 = reads[read][2]
                        type = reads[read][4]
                        if type == 'UP':
                            # unpaired reads are always so
                            if con_1 != 'Z':
                                fn_1 = "%s.%s%s.list.gz" % (b_suffix[type],
                                                          tmp_cid_2_grp[con_1],
                                                          r_suffix[0])
                                pair_str_1 = "UR_NM_NG"
                                hdr_1 = "g_%s;p_%s;b_%s;c_%s;r_%s_UP" %(cid_2_grp[con_1],
                                                                        pair_str_1,
                                                                        bam_fnames[type],
                                                                        con_1,
                                                                        read)
                                try:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_1].append(hdr_1)
                                except KeyError:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_1] = [hdr_1]
                        else:
                            # paired reads may be unpaired in a group sense
                            grp_1 = tmp_cid_2_grp[con_1]
                            grp_2 = tmp_cid_2_grp[con_2]

                            if grp_1 != grp_2:
                                rs1 = 0
                                rs2 = 0
                                pg = "_UG"
                            else:
                                rs1 = 1
                                rs2 = 2
                                pg = "_PG"

                            if con_1 == 'Z' or con_2 == 'Z':
                                pair_str_1 = "PR_UM_NG"
                                pair_str_2 = "PR_UM_NG"
                            else:
                                pair_str_1 = "PR_PM"+pg
                                pair_str_2 = "PR_PM"+pg

                            if con_1 != 'Z':
                                hdr_1 = "g_%s;p_%s;b_%s;c_%s;r_%s_%s" %(cid_2_grp[con_1],
                                                                          pair_str_1,
                                                                          bam_fnames[type],
                                                                          con_1,
                                                                          read,
                                                                          type)
                                fn_1 = "%s.%s%s.list.gz" % (b_suffix[type],
                                                           grp_1,
                                                           r_suffix[rs1])
                                try:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_1].append(hdr_1)
                                except KeyError:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_1] = [hdr_1]

                            if con_2 != 'Z':
                                hdr_2 = "g_%s;p_%s;b_%s;c_%s;r_%s_%s" %(cid_2_grp[con_2],
                                                                          pair_str_2,
                                                                          bam_fnames[type],
                                                                          con_2,
                                                                          read,
                                                                          type)
                                fn_2 = "%s.%s%s.list.gz" % (b_suffix[type],
                                                           grp_2,
                                                           r_suffix[rs2])

                                try:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_2].append(hdr_2)
                                except KeyError:
                                    read_tracker[i_opt][b_opt][r_opt][g_opt][fn_2] = [hdr_2]
#----------------------------
# Write the files

    with open(os.path.join(args.outdir, results_filename), "w") as r_fh:
        r_fh.write(json.dumps({"links": links,
                               "coverages": coverages,
                               "extracts":read_tracker
                               })
                   )

    with open(os.path.join(args.outdir, contigs_filename), "w") as con_fh:
        for cid in contigs.keys():
            if cid != 'Z':
                con_fh.write(">%s\n%s\n" % (cid, contigs[cid][1]))

    with open(os.path.join(args.outdir, out_fnames['PE'][0]), "w") as pe1_fh:
        pe2_fh = open(os.path.join(args.outdir, out_fnames['PE'][1]), "w")
        for (rid, r1, r2) in pe:
            pe1_fh.write(">%s_PE\n%s\n" % (rid, r1))
            pe2_fh.write(">%s_PE\n%s\n" % (rid, r2))
        pe2_fh.close()

    with open(os.path.join(args.outdir, out_fnames['MP']), "w") as mp_fh:
        for (rid, r1, r2) in mp:
            mp_fh.write(">%s_MP\n%s\n>%s_MP\n%s\n" % (rid, r1, rid, r2))

    with open(os.path.join(args.outdir, out_fnames['UP']), "w") as up_fh:
        for (rid, r1) in up:
            up_fh.write(">%s_UP\n%s\n" % (rid, r1))
            
    if args.bad:        
        with open(os.path.join(args.outdir, bad_contigs_filename), "w") as bad_con_fh:
            for cid in contigs.keys():
                if cid != 'Z':
                    bad_con_fh.write(">%s\n%s\n" % ('AA', contigs[cid][1]))

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--readkey', help="describes how reads should be made")
    parser.add_argument('-f', '--fasta', help="fasta file to cut contigs from")
    parser.add_argument('-g', '--groups', nargs='+', help="groups files")
    parser.add_argument('--bad', action="store_true", help="simulate invalid reference input")
    parser.add_argument('-o',
                        '--outdir',
                        default='.',
                        help="where to write output files to")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
