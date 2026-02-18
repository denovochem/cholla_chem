## Basic CLI usage:
cholla-chem resolves one or more chemical names to SMILES from the command line. You can pass names directly as positional arguments, or provide an input text file with one name per line. Names are split on whitespace, so names with spaces should be quoted.

Resolve a chemical name:

```cholla-chem "aspirin"```

Resolve multiple chemical names:

```cholla-chem "aspirin" "2-acetyloxybenzoic acid" "NaH"```

Read names from a file (one per line):

```cholla-chem --input names.txt```

Write tab-separated output (name<TAB>smiles) to a file instead of stdout:

```cholla-chem --input names.txt --output results.tsv```

<br>

## Options:
- ```--smiles-selection-mode```: SMILES selection mode when resolvers disagree. Default: weighted. (str)

- ```--detailed-name-dict```: Whether to return a more detailed structure instead of just the selected SMILES. Default: False. (bool)

- ```--batch-size```: Batch size for resolution. Default: 500. (int)

- ```--normalize-unicode```: Whether to normalize unicode in inputs. Default: True. (bool)

- ```--split-names-to-solve```: Whether to split names on common delimiters to resolve mixtures. Default: True. (bool)

- ```--resolve-peptide-shorthand```: Whether to expand/resolve peptide shorthand. Default: True. (bool)

- ```--attempt-name-correction```: Whether to attempt name correction (OCR/typos/pagination errors, etc.). Default: True. (bool)

- ```--internet-connection-available```: Whether to allow internet-backed resolvers (where applicable). Default: True. (bool)