import argparse
from typing import Optional, Sequence

from cholla_chem.main import resolve_compounds_to_smiles
from cholla_chem.utils.file_utils import read_names_from_file, write_results


def build_parser() -> argparse.ArgumentParser:
    """
    Build a parser for the CLI.

    Returns:
        argparse.ArgumentParser: A parser object
    """
    p = argparse.ArgumentParser(prog="cholla-chem")
    p.add_argument("names", nargs="*", type=str, help="Chemical names to resolve")
    p.add_argument("--input", "-i", type=str, help="Text file with one name per line")
    p.add_argument(
        "--input-format",
        type=str,
        default=None,
        choices=["txt", "csv", "tsv"],
        help="Optional override for input format (otherwise inferred from file extension)",
    )
    p.add_argument(
        "--input-column",
        type=str,
        default="name",
        help="Column name to read from CSV/TSV inputs (default: name)",
    )
    p.add_argument(
        "--output", "-o", type=str, help="Write results to a file (default: stdout)"
    )
    p.add_argument(
        "--output-format",
        type=str,
        default=None,
        choices=["json", "csv", "tsv", "smi", "txt"],
        help="Optional override for output format (otherwise inferred from output file extension)",
    )
    p.add_argument(
        "--smiles-selection-mode",
        default="weighted",
        type=str,
        help="Smiles selection mode",
    )
    p.add_argument(
        "--detailed-name-dict",
        default=False,
        type=bool,
        help="Whether to return a detailed name dictionary",
    )
    p.add_argument("--batch-size", default=500, type=int, help="Batch size")
    p.add_argument(
        "--normalize-unicode",
        default=True,
        type=bool,
        help="Whether to normalize unicode",
    )
    p.add_argument(
        "--split-names-to-solve",
        default=True,
        type=bool,
        help="Whether to split names to solve",
    )
    p.add_argument(
        "--resolve-peptide-shorthand",
        default=True,
        type=bool,
        help="Whether to resolve peptide shorthand",
    )
    p.add_argument(
        "--attempt-name-correction",
        default=True,
        type=bool,
        help="Whether to attempt name correction",
    )
    p.add_argument(
        "--internet-connection-available",
        default=True,
        type=bool,
        help="Whether an internet connection is available",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Main entry point for the CLI.

    Resolve compound names to SMILES strings.

    Args:
        argv: Optional sequence of command-line arguments

    Returns:
        int: Exit code (0 for success, non-zero for failure)
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    names = list(args.names)
    if args.input:
        names.extend(
            read_names_from_file(
                args.input,
                input_format=args.input_format,
                input_column=args.input_column,
            )
        )

    if not names:
        parser.error("Provide names as arguments or via --input")

    if args.output_format != "json" and args.detailed_name_dict:
        parser.error("--detailed-name-dict can only be used with JSON output format")

    results = resolve_compounds_to_smiles(
        names,
        smiles_selection_mode=args.smiles_selection_mode,
        detailed_name_dict=args.detailed_name_dict,
        batch_size=args.batch_size,
        normalize_unicode=args.normalize_unicode,
        split_names_to_solve=args.split_names_to_solve,
        resolve_peptide_shorthand=args.resolve_peptide_shorthand,
        attempt_name_correction=args.attempt_name_correction,
        internet_connection_available=args.internet_connection_available,
    )

    write_results(
        results,
        output_path=args.output,
        output_format=args.output_format,
    )

    return 0
