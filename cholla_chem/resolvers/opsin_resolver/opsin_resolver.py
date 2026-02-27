import os
import subprocess
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from importlib import resources
from typing import Dict, Iterator, List, Sequence, Tuple, Union

from cholla_chem.utils.logging_config import logger

try:
    # python < 3.8
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


@dataclass(frozen=True)
class OpsinResult:
    outputs: List[str]
    errors: List[str]
    returncode: int


@contextmanager
def opsin_jar_path() -> Iterator[str]:
    jar_name = "opsin-cli-2.8.0-jar-with-dependencies.jar"

    try:
        # py >= 3.9
        jar_resource = resources.files("cholla_chem.datafiles").joinpath(jar_name)
        with resources.as_file(jar_resource) as jar_path:
            yield str(jar_path)
    except Exception:
        # py < 3.9 fallback
        with resources.path("cholla_chem.datafiles", jar_name) as jar_path:
            yield str(jar_path)


def run_opsin(
    chemical_name: Union[str, Sequence[str]],
    output_format: Literal[
        "SMILES",
        "ExtendedSMILES",
        "InChI",
        "StdInChI",
        "StdInChIKey",
    ] = "SMILES",
    allow_acid: bool = False,
    allow_radicals: bool = False,
    allow_bad_stereo: bool = False,
    wildcard_radicals: bool = False,
    failure_analysis: bool = False,
) -> OpsinResult:
    with opsin_jar_path() as jar_fpath:
        arg_list = ["java", "-jar", jar_fpath]

        format_args = {
            "SMILES": "-osmi",
            "ExtendedSMILES": "-oextendedsmiles",
            "InChI": "-oinchi",
            "StdInChI": "-ostdinchi",
            "StdInChIKey": "-ostdinchikey",
        }
        try:
            arg_list.append(format_args[output_format])
        except KeyError:
            raise RuntimeError("Output format {:s} is invalid.".format(output_format))

        if isinstance(chemical_name, str):
            names: List[str] = [chemical_name]
        else:
            names = list(chemical_name)

        sanitized_names = [s.replace("\n", "") for s in names]

        tmp_fpath = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".txt", delete=False, encoding="utf-8"
            ) as f:
                tmp_fpath = f.name
                f.write("\n".join(sanitized_names) + "\n")

            if allow_acid:
                arg_list.append("-a")
            if allow_radicals:
                arg_list.append("-r")
            if allow_bad_stereo:
                arg_list.append("-s")
            if wildcard_radicals:
                arg_list.append("-w")
            if failure_analysis:
                arg_list.append("-detailedFailureAnalysis")

            arg_list.append(tmp_fpath)

            result = subprocess.run(
                arg_list,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
            )

            if result.returncode != 0 and len(result.stdout) == 0:
                msg = "\n".join([s for s in result.stderr if s]).strip()
                if not msg:
                    msg = "OPSIN failed with non-zero return code but stderr was empty."
                n = len(sanitized_names)
                return OpsinResult(
                    outputs=[""] * n,
                    errors=[msg] * n,
                    returncode=result.returncode,
                )

            stdout_lines = (result.stdout or "").replace("\r", "").split("\n")
            stderr_lines = (result.stderr or "").replace("\r", "").split("\n")

            if stdout_lines and stdout_lines[-1] == "":
                stdout_lines = stdout_lines[:-1]
            if stderr_lines and stderr_lines[-1] == "":
                stderr_lines = stderr_lines[:-1]

            outputs = stdout_lines

            errors: List[str] = []
            err_i = 0
            for out in outputs:
                if out != "":
                    errors.append("")
                else:
                    errors.append(
                        stderr_lines[err_i] if err_i < len(stderr_lines) else ""
                    )
                    err_i += 1

            if len(outputs) != len(sanitized_names):
                logger.warning(
                    f"OPSIN output length mismatch: outputs ({len(outputs)}), inputs ({len(sanitized_names)})"
                )
                outputs = (outputs + [""] * len(sanitized_names))[
                    : len(sanitized_names)
                ]
                errors = (errors + [""] * len(sanitized_names))[: len(sanitized_names)]

            return OpsinResult(
                outputs=outputs, errors=errors, returncode=result.returncode
            )

        except Exception as e:
            logger.exception(f"Unexpected error occurred: {e}")
            n = len(names)
            if len(stdout_lines) == 0 and result.returncode != 0:
                msg = "\n".join([s for s in stderr_lines if s]).strip()
                if not msg:
                    msg = "OPSIN failed with non-zero return code but stderr was empty."
                return OpsinResult(
                    outputs=[""] * n, errors=[msg] * n, returncode=result.returncode
                )
            else:
                return OpsinResult(
                    outputs=[""] * n,
                    errors=[f"Unexpected error occurred: {e}"] * n,
                    returncode=1,
                )

        finally:
            if tmp_fpath is not None:
                try:
                    os.remove(tmp_fpath)
                except FileNotFoundError:
                    pass


def name_to_smiles_opsin(
    compound_name_list: List[str],
    allow_acid: bool = False,
    allow_radicals: bool = True,
    allow_bad_stereo: bool = False,
    wildcard_radicals: bool = False,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Convert a list of chemical names to their corresponding SMILES representations using OPSIN.

    Args:
    compound_name_list (List[str]): A list of IUPAC or common chemical names to be converted.
    allow_acid (bool): If True, allow interpretation of acids.
    allow_radicals (bool): If True, enable radical interpretation.
    allow_bad_stereo (bool): If True, allow OPSIN to ignore uninterpretable stereochem.
    wildcard_radicals (bool): If True, output radicals as wildcards.

    Returns:
    Tuple[Dict[str, str], Dict[str, str]]:
        - First dict: Mapping from successfully converted chemical names to their SMILES strings.
        - Second dict: Mapping from chemical names that failed conversion to their error messages.

    Note:
        This function uses the `py2opsin` library to interface with OPSIN.
        Newline characters in input names are stripped to avoid CLI parsing issues.
    """
    opsin_name_dict: Dict[str, str] = {}
    failure_message_dict: Dict[str, str] = {}

    # Strip newlines to prevent CLI parsing issues
    sanitized_names = [
        compound_name.replace("\n", "") for compound_name in compound_name_list
    ]

    result = run_opsin(
        chemical_name=sanitized_names,
        output_format="SMILES",
        failure_analysis=True,
        allow_acid=allow_acid,
        allow_radicals=allow_radicals,
        allow_bad_stereo=allow_bad_stereo,
        wildcard_radicals=wildcard_radicals,
    )

    smiles_strings = result.outputs
    failure_messages = result.errors

    if len(smiles_strings) != len(compound_name_list) or len(failure_messages) != len(
        compound_name_list
    ):
        logger.warning(
            f"Mismatching lengths: "
            f"smiles_strings ({len(smiles_strings)}), "
            f"compound_name_list ({len(compound_name_list)}), "
            f"failure_messages ({len(failure_messages)})"
        )
        return {}, {}

    for compound_name, smiles, msg in zip(
        compound_name_list, smiles_strings, failure_messages
    ):
        if smiles:
            opsin_name_dict[compound_name] = smiles
        if msg:
            failure_message_dict[compound_name] = msg

    return opsin_name_dict, failure_message_dict
