#!/usr/bin/env python2


import os
import sys
import argparse
from glob import iglob

MATCH=3
MISMATCH=-3
GAP=-2


class sw_cell:
    def __init__ (self, score):
        self.score = score
        self.list  = []



def smith_waterman(seq1, seq2, match=MATCH,mismatch=MISMATCH,gap=GAP):
    pass
            


            




def params():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Imposta i parametri
    parser.add_argument("seq1",nargs="?", help="First sequence to align", type=str)
    parser.add_argument("seq2",nargs="?", help="Second sequence to align", type=str)
    parser.add_argument("-i", "--input",required=False,help="input sequences separed with tabs", type=str)
    parser.add_argument("-m", "--match", help="score of a match", type=float, default=MATCH)
    parser.add_argument("-s", "--mismatch", help="score of a mismatch", type=float, default=MISMATCH)
    parser.add_argument("-g","--gap", help="score of a gap", type=float, default=GAP)
    parser.add_argument("--minlength", help="minimum length accepted", type=int, default=1)
    parser.add_argument("--minscore", help="minimum score accepted", type=float, default=1)
    parser.add_argument("--numresult", help="number of alignaments to return", type=int, default=100)

    return parser.parse_args()

def check_params(args):

    args.seqs = []

    if args.input is not None:
        if os.path.exists(args.input):
            if ps.path.isfile(args.input):
                
                size = len(args.seqs)
                args.seqs += [tuple(line.strip().split()) for line in open(args.input) if len(line.strip().split()) == 2]
                
                if len(args.seqs) <= size:
                    print("the file could be empty or wrongly formatted: \nseq1\tseq2\nseq3\tseq4\nseq5\tseq6")
            else:
                print ("input is not a file")
        else:
            print("input file does not exists")

    if (args.seq1 is None) or (args.seq2 is None):        
        print("either seq1 or seq2 is not provided")
        if not args.seqs:
            print("you did not give me anything...")
            sys.exit(1)
    else:
        args.seqs.append((args.seq1, args.seq2))

    if args.match <= 0:
        print("match must be > 0")
        sys.exit(1)
    
    if args.mismatch > 0:
        print("mismatch must be <= 0")
        sys.exit(1)

    if args.gap > 0:
        print("gap must be <= 0")
        sys.exit(1)

    if args.minscore is not None:
        if args.minscore <= 0:
            print("negative minscore is not ok")
            sys.exit(1)

    return args

    
if __name__ == "__main__":

    # Carica i parametri
    args = params()

    # Controlla i parametri
    args = check_params (args)

    for seqs in args.seqs:
        smith_waterman (seqs[0], seqs[1], args.match, args.mismatch, args.gap)

    
    