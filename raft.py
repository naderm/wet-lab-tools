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

def _is_contacting(res_a, res_b, distance = 4):
    return any(abs(atom_a.get_coord() - atom_b.get_coord()) < distance
               for atom_a in res_a.get_list()
               for atom_b in res_b.get_list())

def _get_num_contacts(set_a, set_b):
    return len(1 for i, j in zip(set_a, set_b)
               if i != j and _is_contacting(i, j))

def _get_residues(structure):
    return [res
            for model in structure.get_list()
            for chain in model.get_list()
            for residue in chain.get_list()]

def raft(struct, segment):
    residues = _get_residues(struct)

    seg_residues = [i for i in residues if i.id in segment]
    other_residues = [i for i in residues if i.id not in segment]

    intra_contacts = _get_num_contacts(seg_residues, seg_residues)
    inter_contacts = _get_num_contacts(seg_residues, other_residues)

    val = (intra_contacts - inter_contacts) / len(segment)

    return val

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
            struct = p.get_structure(structure_id, pdb_input)
            print struct

if __name__ == "__main__":
    main(sys.argv[1:])
