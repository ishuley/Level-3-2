from fractions import *
from functools import reduce
import copy
import math

"""
PatrickJMT on YouTube: https://www.youtube.com/watch?v=uvYTGEZQTEs&list=PLANMHOrJaFxPMQCMYcYqwOCYlreFswAKP was instrumental in helping me to wrap my head around the application of  absorbing markov chains

"""


# converts values to percentages, expressed as fractions
# must call before get_absorbing_form() or it will misinterpret '1' for this program's purposes
def get_transition_form(m):
    sol = m[:]
    counter = 0
    for state in sol:
        for val in state:
            counter += val
        for i, val in enumerate(state):
            if val == 0:
                continue
            state[i] = Fraction(val, counter)
        counter = 0
    return sol


# change terminal stages of a matrix into absorbing,
# since not moving is the same as moving to the same place infinitely
def get_absorbing_form(m):
    sol = m[:]
    for i, state in enumerate(sol):
        for item in state:
            if item != 0:
                break
        else:
            sol[i][i] = Fraction(1)
    return sol


# create a key to tell sorted() that we want this list of lists to be sorted by the number of non-zero elements in each list
def nonzero_count(state):
    nonzero_counter = 0
    for item in state:
        if item != 0:
            nonzero_counter += 1
    return nonzero_counter


# credit: https://www.programiz.com/python-programming/examples/multiply-matrix
def multiply_matrices(X, Y):
    return [[sum(a*b for a, b in zip(X_row, Y_col)) for Y_col in zip(*Y)] for X_row in X]


# influenced by: https://www.geeksforgeeks.org/program-subtraction-matices/
def subtract_matrices(m1, m2):
    sol = []
    for seg_index in range(len(m1)):
        temp = []
        for item_index in range(len(m1[seg_index])):
            temp.append(m1[seg_index][item_index] - m2[seg_index][item_index])
        sol.append(temp)
    return sol


# F = inverse of (I-Q)
# solution matrix is (F*R)
def get_sol_matrix(r, q):
    i = get_identity_matrix(len(q))
    i_minus_q = subtract_matrices(i, q)
    f = get_inverse(i_minus_q)
    fr = multiply_matrices(f, r)
    return fr


def get_identity_matrix(n):
    sol = []
    for i in range(n):
        temp = []
        for j in range(n):
            if i == j:
                temp.append(Fraction(1))
            else:
                temp.append(Fraction(0))
        sol.append(temp)
    return sol



# get_inverse() is not my work, credit: https://elonen.iki.fi/code/misc-notes/python-gaussj/
def get_inverse(mat):
    m = mat[:]
    """
    return the inv of the matrix M
    """
    # clone the matrix and append the identity matrix
    # [int(i==j) for j in range_M] is nothing but the i(th row of the identity matrix
    m2 = [row[:]+[int(i == j) for j in range(len(m))]
          for i, row in enumerate(m)]
    # extract the appended matrix (kind of m2[m:,...]
    return [row[len(m[0]):] for row in m2] if gauss_jordan(m2) else None


# gauss_jordan() is not my work, credit: https://elonen.iki.fi/code/misc-notes/python-gaussj/
# this algorithm is awe inspiring
def gauss_jordan(m, eps=1.0/(10**10)):
    """Puts given matrix (2D array) into the Reduced Row Echelon Form.
    Returns True if successful, False if 'm' is singular.
    NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
    Written by Jarno Elonen in April 2005, released into Public Domain"""
    (h, w) = (len(m), len(m[0]))
    for y in range(0, h):
        maxrow = y
        for y2 in range(y+1, h):    # Find max pivot
            if abs(m[y2][y]) > abs(m[maxrow][y]):
                maxrow = y2
        (m[y], m[maxrow]) = (m[maxrow], m[y])
        if abs(m[y][y]) <= eps:     # Singular?
            return False
        for y2 in range(y+1, h):    # Eliminate column y
            c = m[y2][y] / m[y][y]
            for x in range(y, w):
                m[y2][x] -= m[y][x] * c
    for y in range(h-1, 0-1, -1):  # Backsubstitute
        c = m[y][y]
        for y2 in range(0, y):
            for x in range(w-1, y-1, -1):
                m[y2][x] -= m[y][x] * m[y2][y] / c
        m[y][y] /= c
        for x in range(h, w):       # Normalize row y
            m[y][x] /= c
    return True


def get_num_absorbing(std_m):
    sol = 0
    for seg in std_m:
        if nonzero_count(seg) == 1:
            sol += 1
    return sol


# put the matrix into standard form (move the stages around)
def get_standard_form(m):
    # return the list sorted by the key
    sorted_form = sorted(m, key=nonzero_count)
    lowest_absorbing_index = None  # will correspond to number of columns to move
    for seg in sorted_form:
        if nonzero_count(seg) == 1:
            for i, val in enumerate(seg):
                if val == Fraction(1):
                    if lowest_absorbing_index == None or lowest_absorbing_index > i:
                        lowest_absorbing_index = i
    for seg in sorted_form:
        for _ in range(lowest_absorbing_index):
            seg.append(seg.pop(0))
    return sorted_form


# get subsets R and Q from matrix in standard form, matrix sent must be in standard form
def get_split_standard_form(std_m):

    r, q = [], []
    num_absorbing = get_num_absorbing(std_m)
    for state in std_m[:]:
        if nonzero_count(state) > 1:
            r.append(state[:num_absorbing])
            q.append(state[num_absorbing:])
    return r, q


# credit: TakingItCasual retrieved from
# https://stackoverflow.com/questions/37237954/calculate-the-lcm-of-a-list-of-given-numbers-in-python
def get_lcm(denoms):
    return reduce(lambda a,b: a*b // math.gcd(a,b), denoms)


def solution(m):
    if len(m) == 1:
        return [1, 1]
    m = get_transition_form(m)
    m = get_absorbing_form(m)
    # convert everything to a fraction if it isn't already (should only be zeros that aren't)
    # this is primarily for readability while testing
    for i, seg in enumerate(m):
        for j, val in enumerate(seg):
            if val == 0:
                m[i][j] = Fraction(0)
    m = get_standard_form(m)
    r, q = get_split_standard_form(m)

    m = get_sol_matrix(r, q)

    m = m[0]

    # format and return result
    denoms = [val.denominator for val in m]

    lcm = get_lcm(denoms)

    sol = []
    for frac in m:
        numerator = 0
        if frac.numerator != 0:
            numerator = (frac.numerator*lcm)//frac.denominator
        else:
            numerator = 0
        sol.append(numerator)
    sol.append(lcm)

    return sol

#############################T TEST CASES ############################################################
###################################################################################################
###################################################################################################
##############################################

print(solution([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]]))


print(solution([
    [0, 1, 0, 0, 0, 1], 
    [4, 0, 0, 3, 2, 0], 
    [0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0]]))


print(solution([[1, 2, 3, 0, 0, 0],
        [4, 5, 6, 0, 0, 0],
        [7, 8, 9, 1, 0, 0],
        [0, 0, 0, 0, 1, 2],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]]))


print(solution([[0]]))


print(solution([[0, 0, 12, 0, 15, 0, 0, 0, 1, 8],
        [0, 0, 60, 0, 0, 7, 13, 0, 0, 0],
        [0, 15, 0, 8, 7, 0, 0, 1, 9, 0],
        [23, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [37, 35, 0, 0, 0, 0, 3, 21, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), "WRONG")


print(solution([[0, 7, 0, 17, 0, 1, 0, 5, 0, 2],
        [0, 0, 29, 0, 28, 0, 3, 0, 16, 0],
        [0, 3, 0, 0, 0, 1, 0, 0, 0, 0],
        [48, 0, 3, 0, 0, 0, 17, 0, 0, 0],
        [0, 6, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), "WRONG!")


print(solution([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]))

print(solution([[1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), "WRONG!")

print(solution([[0, 86, 61, 189, 0, 18, 12, 33, 66, 39],
        [0, 0, 2, 0, 0, 1, 0, 0, 0, 0],
        [15, 187, 0, 0, 18, 23, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), "WRONG!")

print(solution([[0, 0, 0, 0, 3, 5, 0, 0, 0, 2],
        [0, 0, 4, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 4, 4, 0, 0, 0, 1, 1],
        [13, 0, 0, 0, 0, 0, 2, 0, 0, 0],
        [0, 1, 8, 7, 0, 0, 0, 1, 3, 0],
        [1, 7, 0, 0, 0, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), "WRONG!")
