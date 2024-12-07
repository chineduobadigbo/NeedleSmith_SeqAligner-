from ReadScoringMatrix import read_scoring_matrix
import numpy as np


def smith_waterman_affine(
    sequ_one: str,
    sequ_two: str,
    scoring_matrix_file: str,
    gap_penalty: int,
    extend_penalty: int,
):
    """
    Computes local alignments using the Smith-Waterman algorithm with affine gap penalties.

    Args:
        sequ_one (str): The first sequence for alignment.
        sequ_two (str): The second sequence for alignment.
        scoring_matrix_file (str): Path to the scoring matrix file.
        gap_penalty (int): Gap opening penalty.
        extend_penalty (int): Gap extension penalty.

    Returns:
        tuple: A tuple containing:
            - score (float): The highest alignment score.
            - alignments (list of tuples): A list of all optimal alignments,
              each represented as a tuple of two strings (aligned sequence 1, aligned sequence 2).
    """
    x = sequ_one
    y = sequ_two
    n = len(x)
    m = len(y)
    d = gap_penalty
    e = extend_penalty
    s = read_scoring_matrix(scoring_matrix_file)
    Gd = np.zeros([n + 1, m + 1])
    Gx = np.zeros([n + 1, m + 1])
    Gy = np.zeros([n + 1, m + 1])
    Bx = np.empty([n + 1, m + 1], dtype=object)
    Bx.fill([])
    By = np.empty([n + 1, m + 1], dtype=object)
    By.fill([])
    Bd = np.empty([n + 1, m + 1], dtype=object)
    Bd.fill([])
    # initialize first row and column
    Gd[0, 0] = 0

    Gx[0, 0] = -np.inf
    Gy[0, 0] = -np.inf
    for j in range(1, m + 1):
        Gx[0, j] = -d - (j - 1) * e
        Gy[0, j] = -np.inf
        Gd[0, j] = -np.inf

    for i in range(1, n + 1):
        Gy[i, 0] = -d - (i - 1) * e
        Gx[i, 0] = -np.inf
        Gd[i, 0] = -np.inf

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # get similarity value for given amino acid pair
            # using scoring matrix
            matrix_val = s[x[i - 1]][y[j - 1]]

            # start with Gd
            # check which matrrix has highest value in upper left diagonal
            # and add scoring matrix value
            Gd[i, j] = (
                max(
                    mVal(Gd, (i - 1, j - 1)),
                    mVal(Gy, (i - 1, j - 1)),
                    mVal(Gx, (i - 1, j - 1)),
                )
                + matrix_val
            )
            # Now store origin matrix in Gd matrix
            # Initialize origin with empty list because multiple max values can exist
            Bd[i, j] = []
            if Gd[i, j] == mVal(Gd, (i - 1, j - 1)) + matrix_val:
                Bd[i, j].append("d")
            if Gd[i, j] == mVal(Gy, (i - 1, j - 1)) + matrix_val:
                Bd[i, j].append("y")
            if Gd[i, j] == mVal(Gx, (i - 1, j - 1)) + matrix_val:
                Bd[i, j].append("x")

            # Gx value is Max val of left value of Gx and Gd
            # If Gd: open gap
            # If Gx: extend gap
            Gx[i, j] = max(mVal(Gd, (i, j - 1)) - d, mVal(Gx, (i, j - 1)) - e)
            # Initialize arrow with empty list
            # check which matrix had max val and store in Bx
            Bx[i, j] = []
            if Gx[i, j] == mVal(Gd, (i, j - 1)) - d:
                Bx[i, j].append("d")
            if Gx[i, j] == mVal(Gx, (i, j - 1)) - e:
                Bx[i, j].append("x")

            # Gy value is Max val of up value of Gy and Gd
            # If Gd: open gap
            # If Gy: extend gap
            Gy[i, j] = max(mVal(Gd, (i - 1, j)) - d, mVal(Gy, (i - 1, j)) - e)
            # Initialize arrow with empty list
            # check which matrix had max val and store in By
            By[i, j] = []
            if Gy[i, j] == mVal(Gd, (i - 1, j)) - d:
                By[i, j].append("d")
            if Gy[i, j] == mVal(Gy, (i - 1, j)) - e:
                By[i, j].append("y")

    score = np.max([Gd.max(), Gx.max(), Gy.max()])
    # stores which matrix/matrices contained the maximum value
    mats = []
    # stores the respective index of the max val for a given matrix
    indeces = []
    if score == Gd.max():
        mats.append("d")
        # get indices as tuple
        indeces.append(np.unravel_index(Gd.argmax(), Gd.shape))
    if score == Gx.max():
        mats.append("x")
        indeces.append(np.unravel_index(Gx.argmax(), Gx.shape))
    if score == Gy.max():
        mats.append("y")
        indeces.append(np.unravel_index(Gy.argmax(), Gy.shape))
    # create dictionary with G and B matrices
    B = {"d": Bd, "x": Bx, "y": By}
    G = {"d": Gd, "x": Gx, "y": Gy}
    alignments = []
    # loop over all matrices which had the max value
    # and backtrack starting from respective index
    for mat, index in zip(mats, indeces):
        alignments += backtrack(B, G, x, y, mat, *index)
    return score, alignments


def mVal(matrix: np.array, index: tuple):
    """
    Safely retrieves the value at a given index in the scoring matrix,
    ensuring it returns 0 for negative indices.

    Args:
        matrix (np.array): The scoring matrix to retrieve values from.
        index (tuple): A tuple specifying the row and column index.

    Returns:
        float: The value at the specified index or 0 if the index is out of bounds.
    """
    return max(0, matrix[index])


def backtrack(B: dict, G: dict, x: str, y: str, mat: str, i: int, j: int):
    """
    Performs backtracking to reconstruct optimal alignments from the scoring matrices.

    Args:
        B (dict): A dictionary containing backtracking matrices for each direction ('d', 'x', 'y').
        G (dict): A dictionary containing scoring matrices for each direction ('d', 'x', 'y').
        x (str): The first sequence.
        y (str): The second sequence.
        mat (str): The matrix ('d', 'x', 'y') to start backtracking from.
        i (int): The row index in the matrix.
        j (int): The column index in the matrix.

    Returns:
        list of tuples: A list of alignments, where each alignment is a tuple of two strings
                        (aligned sequence 1, aligned sequence 2).
    """
    # get direction to next entry
    # and alignment character based on current matrix
    if mat == "d":
        next_index = (i - 1, j - 1)  # matrix d always goes diagonas
        align_char_x = x[i - 1]
        align_char_y = y[j - 1]
    elif mat == "x":
        next_index = (i, j - 1)  # matrix x always left
        align_char_x = "-"  # gap in x
        align_char_y = y[j - 1]
    else:
        next_index = (i - 1, j)  # matrix y always up
        align_char_x = x[i - 1]
        align_char_y = "-"  # gap in y

    if G[B[mat][i, j][0]][next_index] < 0:
        # if next entry is 0* terminate here
        return [(align_char_x, align_char_y)]

    alignments = []
    for next_mat in B[mat][i, j]:
        # loop over all matrices night might have had the max val
        sub_alignments = backtrack(B, G, x, y, next_mat, *next_index)
        # return the alignments recursively
        for align in sub_alignments:
            # append alignments to current chars
            alignments.append((align[0] + align_char_x, align[1] + align_char_y))
    if G[mat][i, j] == 0:
        # if 0
        # add new alignments that ends here
        alignments.append(("", ""))
    return alignments
