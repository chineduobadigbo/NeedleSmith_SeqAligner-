import argparse
from smith_waterman import smith_waterman_affine
from needleman_wunsch import needleman_wunsch_linear


def main():
    parser = argparse.ArgumentParser(
        "Calculates optimal alignments of amino acid pairs and their scores."
    )
    parser.add_argument("sequ_one", type=str, help="First sequence for comparison.")
    parser.add_argument(
        "sequ_two",
        type=str,
        help="Second sequence for comparison.",
    )
    parser.add_argument(
        "scoring_matrix_file",
        type=str,
        help="Filename fo scoring matrix to use.",
    )
    parser.add_argument(
        "gap_penalty",
        type=int,
        help="Gap Penalty to use.",
    )
    parser.add_argument(
        "--mode",
        choices=["local", "global"],
        required=True,
        help="Alignment mode: 'local' or 'global'.",
    )

    parser.add_argument(
        "--e",
        type=int,
        help="Gap extension penalty to use (required for local mode).",
    )
    args = parser.parse_args()
    # Validate extend_penalty for local mode
    if args.mode == "local" and args.e is None:
        parser.error("The --e argument is required when using local mode.")
    elif args.mode == "global" and args.e is not None:
        parser.error("The --e argument should not be provided in global mode.")

    if args.mode == "global":
        score, alignments = needleman_wunsch_linear(
            args.sequ_one, args.sequ_two, args.scoring_matrix_file, args.gap_penalty
        )
    else:
        score, alignments = smith_waterman_affine(
            args.sequ_one,
            args.sequ_two,
            args.scoring_matrix_file,
            args.gap_penalty,
            args.e,
        )
    print(f"Score: {int(score)}")
    print(f"Found {len(alignments)} optimal {args.mode} alignments:")
    for alignment in alignments:
        print(f"{alignment[0]}")
        print(f"{alignment[1]}")
        print()


if __name__ == "__main__":
    main()
