import argparse
import sys


class Scores:
    def __init__(self, gap, match, mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch


class SequenceDefinition:
    def __init__(self, identifier, description, sequence):
        self.identifier = identifier
        self.description = description
        self.sequence = sequence


def read_input_file(input_file_path):
    read_sequence_list = []

    with open(input_file_path) as file:
        file_contents = file.readlines()

    read_scores = Scores(int(file_contents[0]), int(file_contents[1]), int(file_contents[2]))

    for i in range(3, len(file_contents), 3):
        read_sequence_list.append(SequenceDefinition(file_contents[i], file_contents[i+1], file_contents[i+2]))

    return read_scores, read_sequence_list


# Main script
parser = argparse.ArgumentParser(description="Use this tool to execute the desired algorithm for sequence alignment "
                                             "on the sequences informed on the input file")
parser.add_argument("algorithm", help="Valid options: global (using Needleman-Wunsch) and local (using Smith-Waterman)")
parser.add_argument("input_file", help="Path to the file containing the scores for gap, match and mismatch, each on a "
                                       "separate line in this order, then the sequences as triples with identifier, "
                                       "description and sequence, each on a different line and in this order")

args = parser.parse_args()

if args.algorithm == "global":
    scores, sequence_list = read_input_file(args.input_file)
    raise Exception("Method not fully implemented")
elif args.algorithm == "local":
    raise Exception("Method not implemented")
else:
    print("Invalid algorithm!", file=sys.stderr)
