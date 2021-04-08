import argparse
import sys


# Type definitions
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


class NeedlemanWunschAligner:
    def __init__(self, nw_scores, nw_sequence_list):
        self.scores = nw_scores
        self.sequence_definition_list = nw_sequence_list

    def run_alignment(self):
        # Takes the sequences in pairs and runs the algorithm
        for i in range(len(self.sequence_definition_list)):
            for j in range(i + 1, len(self.sequence_definition_list)):
                self.__build_similarity_matrix(self.sequence_definition_list[i].sequence,
                                               self.sequence_definition_list[j].sequence)

        raise Exception("Method not fully implemented")

    def __build_similarity_matrix(self, sequence1, sequence2):
        num_rows = len(sequence2) + 1
        num_columns = len(sequence1) + 1
        self.similarity_matrix = [[0 for _ in range(num_rows)] for _ in range(num_columns)]

        for i in range(1, num_rows):
            self.similarity_matrix[i][0] = self.scores.gap * i
        for j in range(1, num_columns):
            self.similarity_matrix[0][j] = self.scores.gap * j

        for i in range(1, num_rows):
            for j in range(1, num_columns):
                s = self.scores.match if sequence1[i - 1] == sequence2[j - 1] else self.scores.mismatch
                match = self.similarity_matrix[i - 1][j - 1] + s
                delete = self.similarity_matrix[i - 1][j] + self.scores.gap
                insert = self.similarity_matrix[i][j - 1] + self.scores.gap
                self.similarity_matrix[i][j] = max(match, delete, insert)

        # Print matrix
        for row in self.similarity_matrix:
            print(' '.join(map("{:3d}".format, row)))
        raise Exception("Method not fully implemented")


# Helper functions
def read_input_file(input_file_path):
    read_sequence_list = []

    with open(input_file_path) as file:
        file_contents = file.read().splitlines()

    read_scores = Scores(int(file_contents[0]), int(file_contents[1]), int(file_contents[2]))

    for i in range(3, len(file_contents), 3):
        read_sequence_list.append(SequenceDefinition(file_contents[i], file_contents[i + 1], file_contents[i + 2]))

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
    aligner = NeedlemanWunschAligner(scores, sequence_list)
    aligner.run_alignment()
elif args.algorithm == "local":
    raise Exception("Method not implemented")
else:
    print("Invalid algorithm!", file=sys.stderr)
