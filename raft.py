#!bin/python
# -*- coding: utf-8 -*-
"""
Based on the algorithm described in:

Fischer, K. F. & Marqusee, S. A rapid test for identification of autonomous
folding units in proteins. Journal of molecular biology 302, 701–12 (2000).
"""

from __future__ import division
import argparse
import os
import sys
from Bio.PDB import PDBParser
from numpy.linalg import norm
import numpy as np

def _is_contacting(res_a, res_b, distance = 4):
    return any(norm(atom_a.get_coord() - atom_b.get_coord()) < distance
               for atom_a in res_a.get_list()
               for atom_b in res_b.get_list())

def _get_num_contacts(set_a, set_b):
    return sum(1 for i, j in zip(set_a, set_b)
               if i != j and _is_contacting(i, j))

def _get_residues(structure):
    return [residue
            for model in structure.get_list()
            for chain in model.get_list()
            for residue in chain.get_list()]

def _raft_score(matrix, start, end):
    intra_contacts = matrix[start:end, start:end].sum()
    inter_contacts = matrix[start:end, 0:start].sum() + \
      matrix[end:, start:end].sum()

    val = (intra_contacts - inter_contacts) / (end - start)

    return val

def _contact_matrix(structure):
    residues = _get_residues(structure)
    n_res = len(residues)
    contacts = np.zeros((n_res, n_res), dtype = np.int)

    for index_a in xrange(n_res):
        for index_b in xrange(index_a + 1, n_res):
            if _is_contacting(residues[index_a], residues[index_b]):
                contacts[index_a, index_b] = 1

    return contacts

def _raft_scores(structure):
    print "Building contact matrix"
    matrix = _contact_matrix(structure)
    print "Done"

    scores = np.zeros(matrix.shape)

    for index_a in range(matrix.shape[0]):
        for index_b in range(index_a + 1, matrix.shape[1]):
            scores[index_a, index_b] = _raft_score(matrix, index_a, index_b)

    return scores

def _create_parser():
    parser = argparse.ArgumentParser(
        description = "Process structures with RAFT")
    parser.add_argument("file", nargs = "+",
                        help = "Files containing PDB structure")

    return parser

def main(args):
    parser = _create_parser()
    args = parser.parse_args(args)
    for f in args.file:
        print("{}:".format(f))

        with open(f) as pdb_input:
            p = PDBParser(PERMISSIVE = 1)
            structure_id = ""
            structure = p.get_structure(structure_id, pdb_input)
            print "  struct", structure
            scores = _raft_scores(structure)
            print "  scores", scores

if __name__ == "__main__":
    main(sys.argv[1:])
