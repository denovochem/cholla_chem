# cholla_chem
[![version](https://img.shields.io/github/v/release/denovochem/cholla_chem)](https://github.com/denovochem/cholla_chem/releases)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://gitHub.com/denovochem/cholla_chem/graphs/commit-activity)
[![License](https://img.shields.io/pypi/l/PubChemPy)](https://github.com/denovochem/cholla_chem/blob/main/LICENSE)
[![Run Tests](https://github.com/denovochem/cholla_chem/actions/workflows/tests.yml/badge.svg)](https://github.com/denovochem/cholla_chem/actions/workflows/tests.yml)
[![Build Docs](https://github.com/denovochem/cholla_chem/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/denovochem/cholla_chem/actions/workflows/pages/pages-build-deployment)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/denovochem/cholla_chem/blob/main/examples/example_notebook.ipynb)

This library is used for performant, comprehensive, and customizable name-to-SMILES conversions. 

This library can use the following existing name-to-SMILES resolvers:
- [py2opsin](https://github.com/denovochem/py2opsin)
- [PubChemPy](https://github.com/denovochem/PubChemPy)
- [CIRpy](https://github.com/mcs07/CIRpy)
- [ChemSpiPy](https://github.com/mcs07/ChemSpiPy)


This library also implements the following new resolvers:
- Manually curated dataset of common names not correctly resolved by other resolvers (e.g. 'NaH')
- Structural formula resolver (e.g. 'CH3CH2CH2COOH')
- Inorganic shorthand resolver (e.g. '[Cp*RhCl2]2')


The following string editing/manipulation strategies may be applied to compounds to assist with name-to-SMILES resolution:
- string sanitization for special characters and mojibake.
- Name correction for OCR errors, typos, pagination errors, etc. 
- Splitting compounds on common delimiters (useful for mixtures of compounds, e.g. 'BH3•THF')
- Peptide shorthand expansion (e.g. 'cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)' -> 'cyclo(l-aspartyl-l-arginyl-l-valyl-l-tyrosyl-l-isoleucyl-l-histidyl-l-prolyl-l-phenylalanyl)')


When resolvers disagree on the SMILES for a given compound, a variety of SMILES selection methods can be employed to determine the "best" SMILES for a given compound name. See the documentation for more details.

## Installation

Install cholla_chem with pip directly from this repo:

```shell
pip install git+https://github.com/denovochem/cholla_chem.git
```

## Basic usage
Resolve chemical names to SMILES by passing a string or a list of strings:
```pycon
from cholla_chem import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(['aspirin'])

"{'aspirin': 'CC(=O)Oc1ccccc1C(=O)O'}"
```

See detailed information including which resolver returned which SMILES with detailed_name_dict=True:
```pycon
from cholla_chem import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(
    ['2-acetyloxybenzoic acid'], 
    detailed_name_dict=True
)

"{'2-acetyloxybenzoic acid': {
    'SMILES': 'CC(=O)Oc1ccccc1C(=O)O',
    'SMILES_source': ['pubchem_default', 'opsin_default'],
    'SMILES_dict': {
        'CC(=O)Oc1ccccc1C(=O)O': ['pubchem_default', 'opsin_default']
    },
    'additional_info': {}
}}"
```

## Advanced usage
Many aspects of the name-to-SMILES resolution process can be customized, including the resolvers that are used, the configuration of those resolvers, and the strategy used to pick the best SMILES.

In this example, we resolve chemical names with OPSIN, PubChem, and CIRPy, and use a custom consensus weighting approach to pick the best SMILES:
```pycon
from cholla_chem import (
    OpsinNameResolver,
    PubChemNameResolver,
    CIRpyNameResolver,
    resolve_compounds_to_smiles,
)

opsin_resolver = OpsinNameResolver(
    resolver_name='opsin', 
    resolver_weight=4
)
pubchem_resolver =  PubChemNameResolver(
    resolver_name='pubchem', 
    resolver_weight=3
)
cirpy_resolver = CIRpyNameResolver(
    resolver_name='cirpy', 
    resolver_weight=2
)

resolved_smiles = resolve_compounds_to_smiles(
    ['2-acetyloxybenzoic acid'],
    [opsin_resolver, pubchem_resolver, cirpy_resolver],
    smiles_selection_mode='weighted',
    detailed_name_dict=True
)

"{'2-acetyloxybenzoic acid': {
    'SMILES': 'CC(=O)Oc1ccccc1C(=O)O',
    'SMILES_source': ['opsin', 'pubchem', 'cirpy'],
    'SMILES_dict': {
        'CC(=O)Oc1ccccc1C(=O)O': ['opsin', 'pubchem', 'cirpy']
    },
    'additional_info': {}
}}"
```

## Command line interface

The package can be used as a command line tool. The command line interface can solve single chemical names or read from a file.

```bash
cholla-chem "aspirin"
```

```bash
cholla-chem --input names.txt --output results.tsv
```

```bash
cholla-chem --help
```

See documentation for more details. 

## Documentation
Full documentation is available [here](https://denovochem.github.io/cholla_chem/)

## Contributing

- Feature ideas and bug reports are welcome on the Issue Tracker.
- Fork the [source code](https://github.com/denovochem/cholla_chem) on GitHub, make changes and file a pull request.

## License

cholla_chem is licensed under the [MIT license](https://github.com/denovochem/cholla_chem/blob/main/LICENSE).