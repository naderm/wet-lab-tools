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

def _raft_score(residues, segment):
    seg_residues = [i for i in residues if i.id[1] in segment]
    other_residues = [i for i in residues if i.id[1] not in segment]

    intra_contacts = _get_num_contacts(seg_residues, seg_residues)
    inter_contacts = _get_num_contacts(seg_residues, other_residues)

    val = (intra_contacts - inter_contacts) / len(segment)

    return val

def _raft_scores(structure):
    residues = _get_residues(structure)
    start = min(residues).id[1]
    end = max(residues).id[1]

    return [[_raft_score(residues, range(i, j))
             for j in range(i + 2, end + 1)]
             for i in range(start, end - 2)]

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
            print structure
            print _raft_scores(structure)

if __name__ == "__main__":
    main(sys.argv[1:])
