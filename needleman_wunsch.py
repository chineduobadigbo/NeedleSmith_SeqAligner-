import numpy as np
from ReadScoringMatrix import read_scoring_matrix


def needleman_wunsch_linear(
    sequ_one: str, sequ_two: str, scoring_matrix_file: str, gap_penalty: int
):
    """
    Computes global alignments using the Needleman-Wunsch algorithm with a linear gap penalty.

    Args:
        sequ_one (str): The first sequence for alignment.
        sequ_two (str): The second sequence for alignment.
        scoring_matrix_file (str): Path to the scoring matrix file.
        gap_penalty (int): Gap penalty to apply for insertions or deletions.

    Returns:
        tuple: A tuple containing:
            - score (float): The optimal alignment score.
            - alignments (list of tuples): A list of all optimal alignments,
              each represented as a tuple of two strings (aligned sequence 1, aligned sequence 2).
    """
    x = sequ_one
    y = sequ_two
    n = len(x)
    m = len(y)
    d = gap_penalty
    s = read_scoring_matrix(scoring_matrix_file)
    S = np.zeros([n + 1, m + 1])
    B = np.empty([n + 1, m + 1], dtype=object)
    S[0, 0] = 0
    for j in range(1, m + 1):
        S[0, j] = -j * d
        B[0, j] = [2]  # initialize first column with arrows pointing up
    for i in range(1, n + 1):
        S[i, 0] = -i * d
        B[i, 0] = [1]  # initialize first row with arrows pointing left
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a = np.array(
                (
                    S[i - 1, j - 1] + s[x[i - 1]][y[j - 1]],
                    S[i - 1, j] - d,
                    S[i, j - 1] - d,
                )
            )
            S[i, j] = max(a)
            direction = list(np.where(a == np.max(a))[0])
            B[i, j] = direction
    score = S[n, m]

    alignments = backtrack(B, x, y, n, m)
    return score, alignments


def backtrack(B, x, y, i, j):
    """
    Performs backtracking to reconstruct optimal alignments from the scoring matrix
    and backtracking directions in the Needleman-Wunsch algorithm.

    Args:
        B (np.array): Backtracking matrix indicating directions (0 for diagonal,
                      1 for up, 2 for left).
        x (str): The first sequence.
        y (str): The second sequence.
        i (int): The current row index in the backtracking matrix.
        j (int): The current column index in the backtracking matrix.

    Returns:
        list of tuples: A list of alignments, where each alignment is a tuple of two strings
                        (aligned sequence 1, aligned sequence 2).
    """
    if (i, j) == (0, 0):
        # terminate recursion
        return [("", "")]

    alignments = []
    for direction in B[i, j]:
        if direction == 0:  # Match/mismatch
            sub_alignments = backtrack(B, x, y, i - 1, j - 1)
            for align in sub_alignments:
                alignments.append((align[0] + x[i - 1], align[1] + y[j - 1]))
        elif direction == 1:  # Gap in y
            sub_alignments = backtrack(B, x, y, i - 1, j)
            for align in sub_alignments:
                alignments.append((align[0] + x[i - 1], align[1] + "-"))
        elif direction == 2:  # Gap in x
            sub_alignments = backtrack(B, x, y, i, j - 1)
            for align in sub_alignments:
                alignments.append((align[0] + "-", align[1] + y[j - 1]))
    return alignments
