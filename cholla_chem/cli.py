import argparse
from typing import Optional, Sequence

from cholla_chem.main import resolve_compounds_to_smiles


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cholla-chem")
    p.add_argument("names", nargs="*", type=str, help="Chemical names to resolve")
    p.add_argument("--input", "-i", type=str, help="Text file with one name per line")
    p.add_argument(
        "--output", "-o", type=str, help="Write results to a file (default: stdout)"
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
    parser = build_parser()
    args = parser.parse_args(argv)

    names = list(args.names)
    if args.input:
        with open(args.input, "r", encoding="utf-8") as f:
            names.extend([line.strip() for line in f if line.strip()])

    if not names:
        parser.error("Provide names as arguments or via --input")

    # You’ll adapt this call to match your function signature
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

    out_text = "\n".join(f"{k}\t{v}" for k, v in results.items())

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(out_text + "\n")
    else:
        print(out_text)

    return 0
