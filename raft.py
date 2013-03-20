#!bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import argparse
import os
import sys

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
