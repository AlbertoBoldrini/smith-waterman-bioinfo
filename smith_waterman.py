#!/usr/bin/env python2

# Position in the Smith-Waterman matrx
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
    def __init__ (self, sequence1, sequence2, matrix, best_scores):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.matrix = matrix
        self.best_scores = best_scores

# A formatted alignment with '-' for gaps
class sw_formatted:
    def __init__ (self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2

    def __add__ (self, other):
        return sw_formatted (self.sequence1 + other.sequence1, self.sequence2 + other.sequence2)


# Computes the alignment returning a sw_result object
def smith_waterman (sequence1, sequence2, match_score, mismatch_score, gap_open_score, gap_extend_score):

    # The matrix has one more column ans one more row
    w = len(sequence1) + 1 
    h = len(sequence2) + 1

    # Creates the matrix with sw_cell objects 
    m = [ [sw_cell(0) for y in range(h)] for x in range(w)] 

    # Creates a list that contains the coordinate of best scores
    best_scores = [] 

    # For each cell compute the score
    for x in range(1,w):
        for y in range(1,h):

            # Diagonal movement (the -1 in the sequence is for the extra row and column)
            if sequence1[x-1] == sequence2[y-1]:
                score = m[x-1][y-1].score + match_score
            else:
                score = m[x-1][y-1].score + mismatch_score

            # Is score better then zero?
            if score > 0:
                m[x][y].score = score              # Sets the new maximum score
                m[x][y].paths = [sw_pos(x-1, y-1)] # Sets the path found

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
            
            # If the score is positive inserts it in the best scores list
            if m[x][y].score > 0:
                best_scores.append(sw_pos(x,y))

    # Sorts the best_scores list using the score as key
    best_scores.sort(key=lambda p: m[p.x][p.y].score, reverse=True)

    # Returns the matrix and the best scores
    return sw_result(sequence1, sequence2, m, best_scores)

# Returns a list of formatted equivalent alignments that start from start_pos in the matrix
def format_sw_alignment (result, start_pos):

    # Fetches the current cell
    cell = result.matrix[start_pos.x][start_pos.y]

    # If the cell has a negative score
    # the null alignent is returned
    if cell.score <= 0:
        return [sw_formatted("", "")]   

    # The list of equivalent alignents that end in p1
    output = []

    # For each equivalent branch, it formats the alignment
    for p2 in cell.paths:

        # The contribution of this subpath in the alignemnt
        step = sw_formatted("", "")

        # Copies start position
        p1 = sw_pos (start_pos.x, start_pos.y)

        # It moves from p1 to p2 adding character in the step variable
        while p1.x > p2.x or p1.y > p2.y:

            # Adds the character from sequence1
            if p1.x > p2.x:
                step.sequence1 = result.sequence1[p1.x-1] + step.sequence1
                p1.x -= 1

            else:
                step.sequence1 = "-" + step.sequence1

            # Adds the character from sequence2
            if p1.y > p2.y:
                step.sequence2 = result.sequence2[p1.y-1] + step.sequence2
                p1.y -= 1
            
            else:
                step.sequence2 = "-" + step.sequence2

        # For each equivalent path that starts from p2 is appended the 
        # part of alignment of this subpath
        output += [sw_f + step for sw_f in format_sw_alignment (result, p2)]

    # Returns the list of equivalent alignents that end in p1
    return output


        


# Zona test

''' result = smith_waterman ("CRCCTGGGGAGTRCRG", "CAACGAGCGCAACCCT", 3, -3, -2, -2)



for pos in result.best_scores:

    print ("Score:", result.matrix[pos.x][pos.y].score)

    sw_f_list = format_sw_alignment (result, pos)

    for sw_f in sw_f_list:
        print(sw_f.sequence1)
        print(sw_f.sequence2)
        print("   ")

    break
 '''