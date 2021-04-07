import argparse
import sys


# Main script
parser = argparse.ArgumentParser(description="Use this tool to execute the desired algorithm for sequence alignment "
                                             "on the sequences informed on the input file")
parser.add_argument("algorithm", help="Valid options: global (using Needleman-Wunsch) and local (using Smith-Waterman)")
parser.add_argument("input_file", help="Path to the file containing the scores for gap, match and mismatch, each on a "
                                       "separate line in this order, then the sequences as triples with identifier, "
                                       "description and sequence, each on a different line and in this order")

args = parser.parse_args()

if args.algorithm == "global":
    raise Exception("Method not implemented")
elif args.algorithm == "local":
    raise Exception("Method not implemented")
else:
    print("Invalid algorithm!", file=sys.stderr)
