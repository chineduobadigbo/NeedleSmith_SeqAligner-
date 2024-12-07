import re


def read_scoring_matrix(filename):
    """
    Reads a scoring matrix from a file in BLAST format and parses it into a nested dictionary.

    The BLAST format includes amino acid letters as column and row headers with scores
    represented as integers. The matrix is returned as a dictionary where each key
    corresponds to a row and maps to another dictionary of column scores.

    Args:
        filename (str): Path to the scoring matrix file in BLAST format.

    Returns:
        dict: A nested dictionary representing the scoring matrix. The outer dictionary's
              keys are the row names, and the inner dictionary's keys are the column names.

    Raises:
        IOError: If the file cannot be opened.
        KeyError: If there are duplicate row or column captions, or more columns than expected.
        ValueError: If non-numerical symbols are found in the matrix.
        IndexError: If the number of scores in a row does not match the number of columns.

    Notes:
        - The scoring matrix file must be in BLAST format.
        - The file must not contain duplicate amino acid headers (columns or rows).
    """
    try:
        with open(filename, "r", encoding="UTF-8") as input_file:
            column_captions = " "  # Amino acid letters + *
            while re.match(r"(^\s*$|^\s*#)", column_captions):
                # Ignore comment lines and empty lines in the file.
                # First line with column captions breaks the loop.
                column_captions = input_file.readline()

            # Save the column captions as keys
            colnames = re.findall(r"[A-IK-NP-Z\*]", column_captions)
            del column_captions
            if len(colnames) > 24:  # Sanity check, no more than 23 AAs + *.
                raise KeyError("More columns in the BLAST file than reasonable")
            if len(colnames) != len(set(colnames)):  # Distinct column captions.
                raise KeyError("Duplicate column captions.")

            matrix = dict()  # Contains the return value.

            for line in input_file:
                rowname = line[0]
                if rowname in matrix:  # Every rowname should be new.
                    raise KeyError("Duplicate row captions.")
                else:
                    matrix[rowname] = dict()
                if re.match("[a-zA-Z]", line[1:]):
                    raise ValueError("Non-numerical symbol in the matrix.")
                vals = re.findall(r"-?\d+", line[1:])
                if len(vals) != len(colnames):
                    # Every column caption should have its value in every row
                    raise IndexError("Incorrect number of values")
                for i, v in enumerate(vals):
                    matrix[rowname][colnames[i]] = int(v)

            return matrix
    except IOError as e:
        print(filename + ": " + e.strerror)
        raise
