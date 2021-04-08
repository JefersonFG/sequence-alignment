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


class Direction(enumerate):
    none = 0
    diagonal = 1
    up = 2
    left = 3


class NeedlemanWunschAligner:
    def __init__(self, nw_scores, nw_sequence_list):
        self.scores = nw_scores
        self.sequence_definition_list = nw_sequence_list

    def run_alignment(self):
        # Takes the sequences in pairs and runs the algorithm
        for i in range(len(self.sequence_definition_list)):
            for j in range(i + 1, len(self.sequence_definition_list)):
                sequence_definition1 = self.sequence_definition_list[i]
                sequence_definition2 = self.sequence_definition_list[j]

                self.__build_similarity_matrix(sequence_definition1.sequence, sequence_definition2.sequence)
                final_score, alignment1, alignment2 = self.__compute_optimal_sequence(sequence_definition1.sequence,
                                                                                      sequence_definition2.sequence)

                # Print results
                print(f"Alignment for sequences {sequence_definition1.identifier} and {sequence_definition2.identifier}"
                      f":\nAlignment 1: {alignment1}\nAlignment 2: {alignment2}\nScore: {final_score}")

    def __build_similarity_matrix(self, sequence1, sequence2):
        num_rows = len(sequence1) + 1
        num_columns = len(sequence2) + 1
        self.similarity_matrix = [[(0, Direction.none) for _ in range(num_columns)] for _ in range(num_rows)]

        for i in range(1, num_rows):
            self.similarity_matrix[i][0] = (self.scores.gap * i, Direction.none)
        for j in range(1, num_columns):
            self.similarity_matrix[0][j] = (self.scores.gap * j, Direction.none)

        for i in range(1, num_rows):
            for j in range(1, num_columns):
                s = self.scores.match if sequence1[i - 1] == sequence2[j - 1] else self.scores.mismatch
                match = self.similarity_matrix[i - 1][j - 1][0] + s
                delete = self.similarity_matrix[i - 1][j][0] + self.scores.gap
                insert = self.similarity_matrix[i][j - 1][0] + self.scores.gap

                score = max(match, delete, insert)

                if score == match:
                    direction = Direction.diagonal
                elif score == delete:
                    direction = Direction.left
                elif score == insert:
                    direction = Direction.up
                else:
                    raise Exception("Failed to determine direction")

                self.similarity_matrix[i][j] = (score, direction)

    def __compute_optimal_sequence(self, sequence1, sequence2):
        num_rows = len(sequence1) + 1
        num_columns = len(sequence2) + 1
        i = num_rows - 1
        j = num_columns - 1

        alignment1 = ""
        alignment2 = ""

        final_score, _ = self.similarity_matrix[i][j]

        while i > 0 or j > 0:
            _, current_direction = self.similarity_matrix[i][j]

            if current_direction == Direction.diagonal:
                alignment1 = sequence1[i - 1] + alignment1
                alignment2 = sequence2[j - 1] + alignment2
                i -= 1
                j -= 1
            elif current_direction == Direction.left:
                alignment1 = sequence1[i - 1] + alignment1
                alignment2 = '-' + alignment2
                i -= 1
            else:
                alignment1 = '-' + alignment1
                alignment2 = sequence2[j - 1] + alignment2
                j -= 1

        return final_score, alignment1, alignment2


class SmithWatermanAligner:
    def __init__(self, nw_scores, nw_sequence_list):
        self.scores = nw_scores
        self.sequence_definition_list = nw_sequence_list

    def run_alignment(self):
        # Takes the sequences in pairs and runs the algorithm
        for i in range(len(self.sequence_definition_list)):
            for j in range(i + 1, len(self.sequence_definition_list)):
                sequence_definition1 = self.sequence_definition_list[i]
                sequence_definition2 = self.sequence_definition_list[j]

                self.__build_score_matrix(sequence_definition1.sequence, sequence_definition2.sequence)
                final_alignment_list = self.__compute_optimal_sequence(sequence_definition1.sequence,
                                                                       sequence_definition2.sequence)

                # Print results
                for score, alignment1, alignment2 in final_alignment_list:
                    print(f"Alignment for sequences {sequence_definition1.identifier} and "
                          f"{sequence_definition2.identifier}:\nAlignment 1: {alignment1}\n"
                          f"Alignment 2: {alignment2}\nScore: {score}")

    def __build_score_matrix(self, sequence1, sequence2):
        num_rows = len(sequence1) + 1
        num_columns = len(sequence2) + 1
        self.score_matrix = [[(0, Direction.none) for _ in range(num_columns)] for _ in range(num_rows)]
        self.max_score_list = []

        for i in range(1, num_rows):
            self.score_matrix[i][0] = (0, Direction.none)
        for j in range(1, num_columns):
            self.score_matrix[0][j] = (0, Direction.none)

        for i in range(1, num_rows):
            for j in range(1, num_columns):
                s = self.scores.match if sequence1[i - 1] == sequence2[j - 1] else self.scores.mismatch
                alignment = self.score_matrix[i - 1][j - 1][0] + s
                end_gap_1 = self.score_matrix[i - 1][j][0] + self.scores.gap
                end_gap_2 = self.score_matrix[i][j - 1][0] + self.scores.gap
                no_similarity = 0

                score = max(alignment, end_gap_1, end_gap_2, no_similarity)

                if not self.max_score_list:
                    self.max_score_list.append((score, i, j))
                elif score == self.max_score_list[0][0]:
                    self.max_score_list.append((score, i, j))
                elif score > self.max_score_list[0][0]:
                    del self.max_score_list[:]
                    self.max_score_list.append((score, i, j))

                if score == alignment:
                    direction = Direction.diagonal
                elif score == end_gap_1:
                    direction = Direction.left
                elif score == end_gap_2:
                    direction = Direction.up
                elif score == no_similarity:
                    direction = Direction.none
                else:
                    raise Exception("Failed to determine direction")

                self.score_matrix[i][j] = (score, direction)

    def __compute_optimal_sequence(self, sequence1, sequence2):
        alignment_list = []
        final_score = self.max_score_list[0][0]

        for _, i, j in self.max_score_list:
            alignment1 = ""
            alignment2 = ""
            current_score, current_direction = self.score_matrix[i][j]

            while current_score != 0:
                if current_direction == Direction.diagonal:
                    alignment1 = sequence1[i - 1] + alignment1
                    alignment2 = sequence2[j - 1] + alignment2
                    i -= 1
                    j -= 1
                elif current_direction == Direction.left:
                    alignment1 = sequence1[i - 1] + alignment1
                    alignment2 = '-' + alignment2
                    i -= 1
                else:
                    alignment1 = '-' + alignment1
                    alignment2 = sequence2[j - 1] + alignment2
                    j -= 1

                current_score, current_direction = self.score_matrix[i][j]

            alignment_list.append((final_score, alignment1, alignment2))

        return alignment_list


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
    scores, sequence_list = read_input_file(args.input_file)
    aligner = SmithWatermanAligner(scores, sequence_list)
    aligner.run_alignment()
else:
    print("Invalid algorithm!", file=sys.stderr)
