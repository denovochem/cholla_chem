# placeholder_name
[![PyPI Version](https://img.shields.io/pypi/v/PubChemPy?logo=python&logoColor=%23ffffff)](https://pypi.python.org/pypi/PubChemPy)
[![License](https://img.shields.io/pypi/l/PubChemPy)](https://github.com/denovochem/name_to_smiles/blob/main/LICENSE)
[![Tests](https://img.shields.io/github/actions/workflow/status/mcs07/pubchempy/test.yml?logo=github&logoColor=%23ffffff&label=tests)](https://github.com/mcs07/PubChemPy/actions/workflows/test.yml)
[![Docs](https://img.shields.io/readthedocs/pubchempy?logo=readthedocs&logoColor=%23ffffff)](https://denovochem.github.io/name_to_smiles/)

This library is used for performant, comprehensive, and customizable name-to-SMILES conversions. 

This library uses the following existing name-to-SMILES resolvers:
- OPSIN
- PubChem
- CIRPy

This library implements the following new resolvers and resolution strategies:
- Manual database of common names not resolved by other resolvers
- Structural formula resolver (e.g. CH3CH2CH2COOH)
- Peptide shorthand resolver (e.g. Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)
- Mixtures of compounds (e.g. H₂O•THF)

## Installation

Install placeholder_name with pip directly from this repo:

```shell
pip install git+https://github.com/denovochem/name_to_smiles.git
```

## Basic usage
Resolve chemical names to SMILES by passing a string or a list of strings:
```pycon
from placeholder_name import resolve_compounds_to_smiles
resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'])
''
```
See detailed information including which resolver returned which SMILES with detailed_name_dict=True:
```pycon
from placeholder_name import resolve_compounds_to_smiles
resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'], detailed_name_dict=True)
''
```

## Advanced usage
Many aspects of the name-to-SMILES resolution process can be customized, including the resolvers that are used, the configuration of those resolvers, and the strategy used to pick the best SMILES.

In this example, we resolve chemical names with OPSIN, PubChem, and CIRPy, and use a custom consensus weighting approach to pick the best SMILES:
```pycon
from placeholder_name import resolve_compounds_to_smiles
from placeholder_name import OpsinNameResolver, PubChemNameResolver, CIRPyNameResolver

opsin_resolver = OpsinNameResolver('opsin', resolver_weight=4)
pubchem_resolver =  PubChemNameResolver('pubchem', resolver_weight=3)
cirpy_resolver = CIRPyNameResolver('cirpy', resolver_weight=2)

resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'], [opsin_resolver, pubchem_resolver, cirpy_resolver], smiles_selection_mode='weighted', detailed_name_dict=True)
''
```

## Documentation
Full documentation is availible at 

## Contributing

- Feature ideas and bug reports are welcome on the Issue Tracker.
- Fork the [source code](https://github.com/denovochem/name_to_smiles) on GitHub, make changes and file a pull request.

## License

PubChemPy is licensed under the [MIT license](https://github.com/denovochem/name_to_smiles/blob/main/LICENSE).
