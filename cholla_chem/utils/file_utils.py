import csv
import json
from pathlib import Path
from typing import Dict, List


def _infer_input_format(path: str) -> str:
    """
    Infer the input format of a file from its extension.

    Supported formats are CSV, TSV, and txt.

    Args:
        path: Path to the file

    Returns:
        str: One of "csv", "tsv", or "txt"
    """
    ext = Path(path).suffix.lower()
    if ext == ".csv":
        return "csv"
    if ext == ".tsv":
        return "tsv"
    if ext == ".txt":
        return "txt"
    return "txt"


def _infer_output_format(path: str) -> str:
    """
    Infer the output format of a file from its extension.

    Supported formats are JSON, CSV, and TSV.

    Args:
        path: Path to the file

    Returns:
        str: One of "json", "csv", "tsv"

    Raises:
        ValueError: If the extension is not recognized
    """
    ext = Path(path).suffix.lower()
    if ext == ".json":
        return "json"
    if ext == ".csv":
        return "csv"
    if ext == ".tsv":
        return "tsv"
    if ext == ".smi":
        return "smi"
    if ext == ".txt":
        return "txt"
    raise ValueError(
        f"Unsupported output extension: {ext} (use .json, .csv, .tsv, .smi, .txt)"
    )


def read_names_from_file(
    path: str,
    *,
    input_format: str | None = None,
    input_column: str = "name",
    encoding: str = "utf-8",
) -> List[str]:
    """
    Reads a list of chemical names from a file.

    Supports the following input formats:
    - txt: a plain text file with one name per line
    - csv: a comma-separated value file with a header row
        containing a column like "name"
    - tsv: a tab-separated value file with a header row
        containing a column like "name"

    Args:
        path: Path to the file
        input_format: Optional format of the file (txt, csv, tsv)
        input_column: Name of the column containing the names (default: "name")
        encoding: Encoding of the file (default: "utf-8")

    Returns:
        List[str]: A list of chemical names read from the file
    """
    fmt = (input_format or _infer_input_format(path)).lower()

    if fmt == "txt":
        with open(path, "r", encoding=encoding, newline="") as f:
            return [line.strip() for line in f if line.strip()]

    out: List[str] = []

    if fmt in ("csv", "tsv"):
        delimiter = "," if fmt == "csv" else "\t"
        with open(path, "r", encoding=encoding, newline="") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            if not reader.fieldnames:
                raise ValueError(
                    f"{path} has no header row; expected a column like '{input_column}'"
                )
            if input_column not in reader.fieldnames:
                raise ValueError(
                    f"{path} is missing column '{input_column}'. Columns: {reader.fieldnames}"
                )
            for row in reader:
                val = (row.get(input_column) or "").strip()
                if val:
                    out.append(val)
            return out

    raise ValueError(f"Unsupported input format: {fmt} (use txt, csv, tsv)")


def write_results(
    results: Dict,
    *,
    output_path: str | None,
    output_format: str | None = None,
    encoding: str = "utf-8",
) -> None:
    """
    Write results of chemical name resolution to a file.

    Args:
        results: A dictionary mapping compound names to their SMILES representations.
        output_path: The path to the output file. If not provided, print results to stdout.
        output_format: The format of the output file. Supported formats are json, csv, tsv.
        encoding: The encoding of the output file. Default is utf-8.

    Returns:
        None
    """
    if not output_path:
        out_text = "\n".join(f"{k}\t{v}" for k, v in results.items())
        print(out_text)
        return

    fmt = (output_format or _infer_output_format(output_path)).lower()

    if fmt == "json":
        with open(output_path, "w", encoding=encoding) as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
            f.write("\n")
        return

    if fmt == "smi":
        with open(output_path, "w", encoding=encoding) as f:
            for name, smiles in results.items():
                f.write(f"{name}\t{smiles}\n")
        return

    if fmt == "txt":
        with open(output_path, "w", encoding=encoding) as f:
            for name, smiles in results.items():
                f.write(f"{name}\t{smiles}\n")
        return

    if fmt in ("csv", "tsv"):
        delimiter = "," if fmt == "csv" else "\t"
        with open(output_path, "w", encoding=encoding, newline="") as f:
            w = csv.DictWriter(f, fieldnames=["name", "smiles"], delimiter=delimiter)
            w.writeheader()
            for name, smiles in results.items():
                w.writerow({"name": name, "smiles": smiles})
        return

    raise ValueError(f"Unsupported output format: {fmt}")
