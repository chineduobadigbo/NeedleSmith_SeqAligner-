def is_valid_amino_acid_sequence(sequence):
    """
    Checks if the given sequence contains only valid amino acid characters.

    Args:
        sequence (str): The amino acid sequence to validate.

    Returns:
        bool: True if the sequence is valid, False otherwise.
    """
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    return all(char in valid_amino_acids for char in sequence.upper())
