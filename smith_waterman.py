#!/usr/bin/env python2


import sys

# Position in the Smith-Waterman matrix
class sw_pos:
    def __init__ (self, x, y):
        self.x = x
        self.y = y

# A cell of the Smith-Waterman matrix
class sw_cell:
    def __init__ (self, score):
        self.score = score
        self.paths = []

# The result of the Smith-Waterman alignment
class sw_result:
    def __init__ (self, seq1, seq2, matrix, best_scores):
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrix = matrix
        self.best_scores = best_scores

# Computes the alignment returning a sw_result object
def smith_waterman (seq1, seq2, match_score, mismatch_score, gap_open_score, gap_extend_score):

    # The matrix has one more column ans one more row
    w = len(seq1) + 1 
    h = len(seq2) + 1

    # Creates the matrix with sw_cell objects 
    m = [ [sw_cell(0) for y in range(h)] for x in range(w)] 

    # Creates a list that contains the coordinate of best scores
    best_scores = [] 

    # For each cell compute the score
    for x in range(1,w):
        for y in range(1,h):

            # Looks for possible "horizontal" gaps
            for xg in range(x):

                # Computes the score
                score = m[xg][y].score + gap_open_score + gap_extend_score * (x - xg - 1)

                # Is score better then previous?
                if score > m[x][y].score:
                    m[x][y].score = score           # Sets the new maximum score
                    m[x][y].paths = [sw_pos(xg, y)] # Sets the path found

                # Is score equivalent?
                elif score == m[x][y].score:
                    m[x][y].paths.append (sw_pos(xg, y)) # Adds the equivalent path

            # Looks for possible "vertical" gaps
            for yg in range(y):
                
                # Computes the score
                score = m[x][yg].score + gap_open_score + gap_extend_score * (y - yg - 1)

                # Is score better then previous?
                if score > m[x][y].score:
                    m[x][y].score = score           # Sets the new maximum score
                    m[x][y].paths = [sw_pos(x, yg)] # Sets the path found

                # Is score equivalent?
                elif score == m[x][y].score:
                    m[x][y].paths.append (sw_pos(x, yg)) # Adds the equivalent path
            

            # Diagonal movement (the -1 in the seq is for the extra row and column)

            # Match case
            if seq1[x-1] == seq2[y-1]:

                # Computes the score
                score = m[x-1][y-1].score + match_score

                # Is score better or equal then previous?
                if score >= m[x][y].score:
                    m[x][y].score = score              # Sets the new maximum score
                    m[x][y].paths = [sw_pos(x-1, y-1)] # Sets the path found
                    best_scores.append(sw_pos(x,y))    # It ends with a match, is a best score!

            # Mismatch case
            else:

                # Computes the score
                score = m[x-1][y-1].score + mismatch_score

                # Is score better then previous?
                if score > m[x][y].score:
                    m[x][y].score = score           # Sets the new maximum score
                    m[x][y].paths = [sw_pos(x-1, y-1)] # Sets the path found

                # Is score equivalent?
                elif score == m[x][y].score:
                    m[x][y].paths.append (sw_pos(x-1, y-1)) # Adds the equivalent path            

    # Sorts the best_scores list using the score as key
    best_scores.sort(key=lambda p: m[p.x][p.y].score, reverse=True)

    # Returns the matrix and the best scores
    return sw_result(seq1, seq2, m, best_scores)

# Returns a list of formatted equivalent alignments that start from pos in the matrix
def format_sw_alignment (result, pos): 

    # Defines an object for the stack
    class stack_item:
        def __init__ (self, pos, seq1, seq2):
            self.pos  = pos
            self.seq1 = seq1
            self.seq2 = seq2 

    # Creates a stack for all paths. Inserts the initial path
    stack = [stack_item (pos, "", "")]

    # The output of this function
    output = []

    while len(stack) > 0:

        # Processes the last element of the stack
        item = stack.pop ()

        # For each equivalent branch, it formats the alignment
        for dest in result.matrix[item.pos.x][item.pos.y].paths:

            # Copies the stack item
            processed = stack_item (sw_pos (item.pos.x, item.pos.y), item.seq1, item.seq2)

            # It moves from processed.pos to dest adding character in the seq variables
            while processed.pos.x > dest.x or processed.pos.y > dest.y:

                # Adds the character from seq1
                if processed.pos.x > dest.x:
                    processed.seq1 = result.seq1[processed.pos.x-1] + processed.seq1
                    processed.pos.x -= 1

                else:
                    processed.seq1 = "-" + processed.seq1

                # Adds the character from seq2
                if processed.pos.y > dest.y:
                    processed.seq2 = result.seq2[processed.pos.y-1] + processed.seq2
                    processed.pos.y -= 1
                
                else:
                    processed.seq2 = "-" + processed.seq2

            # If the score is positive the processed item is reinserted in the stack
            if result.matrix[processed.pos.x][processed.pos.y].score > 0:
                stack.append(processed) 

            # Otherwise is put in the output
            else:
                output.append(processed)

    return output

def print_matrix (result):

    # The matrix has one more column ans one more row
    w = len(result.seq1) + 1 
    h = len(result.seq2) + 1

    sys.stdout.write ("\t")

    # Prints the seq2
    for x in range(1,w):
        sys.stdout.write (result.seq1[x-1]+"\t")
    
    # Prints the matrix
    for y in range(1,h):

        # Writes the letter for this line
        sys.stdout.write ("\n" + result.seq2[y-1] + "\t")

        # Write numbers in this line
        for x in range(1,w):
            sys.stdout.write (str(result.matrix[x][y].score) + "\t")

    # Closes the line
    sys.stdout.write ("\n")