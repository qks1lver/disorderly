#!/usr/bin/env python3

# Imports
import os
import argparse
import numpy as np
from datetime import datetime
from random import sample
from multiprocessing import Pool, cpu_count
from functools import partial

# Constants
_amino_acids_ = 'ARNDCQEGHILKMFPSTWYVX'
_aa_set_ = set(_amino_acids_)
_blank_ = [0 for _ in range(len(_amino_acids_) + 1)]
_verbose_ = False
_p_test_ = '../test/'
_p_data_ = '../data/'
_cpu_ = cpu_count()

# Functions
def composition(sequence):

    # Compute properties
    sequence = sequence.upper()
    query_length = len(sequence)
    mass = 1 / query_length

    query_unique = list(set(sequence) - _aa_set_)
    if query_unique:
        print('There are special residues in sequence: %s (ignored during comparison)' % ','.join(query_unique))
        sequence = ''.join([i for i in sequence if i not in query_unique])

    # Compute composition
    comp = _blank_.copy()
    comp[0] = query_length
    for e in sequence:
        comp[1 + _amino_acids_.index(e)] += mass

    return comp

def build_db(p_fasta, p_db=''):

    if not p_db:
        p_db = p_fasta + '.disorderdb'

    if _verbose_:
        print('Building database from FASTA: %s ...' % p_fasta)

    n_seqs = 0
    db_comps = []
    with open(p_fasta, 'r') as fi, open(p_db, 'w+') as fo:
        seq = ''
        header = ''
        for l in fi:
            if l.startswith('>'):
                if seq:
                    comp = composition(seq)
                    _ = fo.write('%s\t%d\t%s\n' % (header, comp[0], '\t'.join(['%.6f' % i for i in comp[1:]])))
                    db_comps.append([header] + comp)
                    n_seqs += 1

                seq = ''
                header = l[1:].strip().replace('\t', ' ')

            else:
                seq += l.strip()

        if seq:
            comp = composition(seq)
            _ = fo.write('%s\t%d\t%s\n' % (header, comp[0], '\t'.join(['%.6f' % i for i in comp[1:]])))
            db_comps.append([header] + comp)
            n_seqs += 1

    if _verbose_:
        print('Generated database for %d sequences at %s' % (n_seqs, p_db))

    return p_db, db_comps

def read_db(p_db):

    db_comps = []
    with open(p_db, 'r') as f:
        for l in f:
            tmp = l.split('\t')
            comp = [tmp[0]]
            comp.append(int(tmp[1]))
            comp += [float(i) for i in tmp[2:]]
            db_comps.append(comp)

    return db_comps

def search(p_query, db_comps, p_out=''):

    if not p_out:
        stamp = 'search-%s-%s.csv' % (datetime.now().strftime('%Y%m%d%H%M%S'), ''.join(sample('ABCDEF', 4)))
        p_file, _  = os.path.splitext(p_query)
        p_out = p_file + '_' + stamp

    with open(p_out, 'w+') as f:
        _ = f.write('Queries,Hits,Distances\n')

    query_seqs = read_fasta(p_query)

    if _verbose_:
        print('Searching on %d CPUs ...' % _cpu_)

    for q_header, q_seq in query_seqs.items():

        query_comp = composition(q_seq)

        with Pool(processes=_cpu_) as pool:
            res = pool.map(partial(_compare_, q_comp=query_comp), db_comps)

        candidates = {h:d for h,d in res if h}

        candidates_sorted = sorted(candidates.items(), key=lambda x: x[1])

        with open(p_out, 'a') as f:

            for h,d in candidates_sorted:
                _ = f.write('%s,%s,%.4f\n' % (q_header, h, d))

    if _verbose_:
        print('Search complete, results saved as %s' % p_out)

    return p_out

def _compare_(t_db_comp, q_comp):

    header = ''
    dist = 1000.

    if q_comp[0] == t_db_comp[1]:
        header = t_db_comp[0]

        # Compute Euclidean distance
        dist = np.linalg.norm(np.array(q_comp[1:]) - np.array(t_db_comp[2:]))

    return header, dist

def read_fasta(p_fasta):

    sequences = {}

    n_seqs = 0
    with open(p_fasta, 'r') as f:
        seq = ''
        header = ''
        for l in f:
            if l.startswith('>'):
                if seq:
                    sequences[header] = seq
                    n_seqs += 1

                seq = ''
                header = l[1:].strip().replace('\t',' ')

            else:
                seq += l.strip()

        if seq:
            sequences[header] = seq
            n_seqs += 1

    if _verbose_:
        print('Found %d sequences in %s' % (n_seqs, p_fasta))

    return sequences

# Run
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Search similar proteins by composition.')
    parser.add_argument('-i', dest='p_query', help='Query FASTA')
    parser.add_argument('-fb', dest='p_db_fasta', default='', help='FASTA of database')
    parser.add_argument('-o', dest='p_out', default='', help='Output file name/path')
    parser.add_argument('-db', dest='p_db', default='', help='Database name/path')
    parser.add_argument('-v', dest='verbose', action='store_true', help='To verbose')

    args = parser.parse_args()

    _verbose_ = args.verbose

    if args.p_db_fasta:
        _, db_comps = build_db(args.p_db_fasta)
    else:
        db_comps = read_db(args.p_db)

    search(args.p_query, db_comps, args.p_out)
