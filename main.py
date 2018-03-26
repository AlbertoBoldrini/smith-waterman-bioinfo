#!/usr/bin/env python


import os
import sys
import argparse
import smith_waterman
from glob import iglob


__author__  = 'Alberto Boldrini'
__email__   = 'alberto.boldrini@studenti.unitn.it'
__version__ = '0.01'



MATCH_SCORE = 3
MISMATCH_SCORE = -3
GAP_OPEN_SCORE = -2
GAP_EXTEND_SCORE = -2




# Parses the command line arguments of this program
def parse_args ():

    # Creates an object to parse the command line parameters
    parser = argparse.ArgumentParser(description="This program apply pure smith-waterman algorithm to a list of sequences.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Explains to the object what parameters it should expect

    # Positional parameters
    parser.add_argument("sequences", nargs="*", help="Sequences to align", type=str)

    # Optional parameters
    parser.add_argument("-i", "--input",      help="input file", type=str)
    parser.add_argument("-m", "--match",      help="match score", type=float, default=MATCH_SCORE)
    parser.add_argument("-t", "--mismatch",   help="mismatch score", type=float, default=MISMATCH_SCORE)
    parser.add_argument("-g", "--gap_open",   help="gap open score", type=float, default=GAP_OPEN_SCORE)
    parser.add_argument("-e", "--gap_extend", help="gap extend score", type=float, default=GAP_EXTEND_SCORE)
    parser.add_argument("-l", "--minlength",  help="min length accepted", type=int, default=1)
    parser.add_argument("-s", "--minscore",   help="min score accepted", type=float, default=0)
    parser.add_argument("-n", "--numresult",  help="max alignments returned", type=int, default=10)

    # Parses and returns the object containing the params
    return parser.parse_args()

# Checks the correctness of the received parameters and loads files
def check_args (args):

    # Inserts in the object args a list of pair sequences to align.
    # They come from the command line and from the input file 
    args.seqs = []

    # Checks if the input param has been specified
    if args.input is not None:

        # Checks if exists something with that name
        if os.path.exists(args.input):

            # Checks if it is a file or something else
            if os.path.isfile(args.input):
                
                # Reads the file and splits the sequences for each line
                sequences = [line.strip().split() for line in open(args.input)]

                for pair in sequences:

                    if len(pair) != 2:
                        sys.exit('The file \'' + args.input + '\' contains one or more lines not correctly formatted.\n'
                                 'The file must be formatted in this way:\n\n'
                                 'seq1\tseq2\nseq3\tseq4\nseq5\tseq6')

                # Inserts the loaded sequences in the args structure
                args.seqs += sequences
            
            else:
                sys.exit('The path \'' + args.input + '\' is not a file.')
        
        else:
            sys.exit('The file \'' + args.input + '\' does not exists.')


    # Checks that the number of positional parameter is even
    if len(args.sequences) % 2 != 0:
        sys.exit('An odd number of sequences has been passed in the command line.')

    # Pairs the sequences two by two
    for i in range(len(args.sequences) >> 1):
        args.seqs.append ([args.sequences[i], args.sequences[i+1]])


    # Checks correctness of the values of parameters
    if args.match <= 0:
        sys.exit('The match score must be > 0.')
    
    if args.mismatch > 0:
        sys.exit('The mismatch score must be <= 0.')

    if args.gap_open > 0:
        sys.exit('The gap_open score must be <= 0.')

    if args.gap_extend > 0:
        sys.exit('The gap_extend score must be <= 0.')

    return args

    
if __name__ == "__main__":

    # Parses command line parameters
    args = parse_args()

    # Checks parameters and loads the input file
    args = check_args (args)

    # Aligns each pair of sequences
    for pair in args.seqs:

        # Print the sequences
        print ("Aligns of the sequences: ")
        print (pair[0])
        print (pair[1],"\n\n")

        # Starts the smith-waterman algorithm
        result = smith_waterman.smith_waterman (pair[0], pair[1], args.match, args.mismatch,
                                                args.gap_open, args.gap_extend)


        # Next prints the alignments resulted
        print ("Results:")

        # Counts the number of alignments already written to the stdout
        number_alignment_written = 0

        # Print a list of string for each score
        for pos in result.best_scores:

            # Fetches the score for this alignment
            score = result.matrix[pos.x][pos.y].score

            # Checks the minscore option
            if score < args.minscore:
                break
            
            # Show the score for the next list of alignments
            print ("Score:", score)

            # Format the alignment in a list of equivalent strings
            sw_f_list = smith_waterman.format_sw_alignment (result, pos)

            for sw_f in sw_f_list:

                # Checks the minlength option
                if len(sw_f.sequence1) < args.minlength:
                    continue
                
                # Prints the alignment
                print(sw_f.sequence1)
                print(sw_f.sequence2)
                print("   ")

                # Counts the new alignment written
                number_alignment_written += 1

                # Stops because the max number of alignments hab been reached
                if number_alignment_written >= args.numresult:
                    break

            # Stops because the max number of alignments hab been reached
            if number_alignment_written >= args.numresult:
                break
    