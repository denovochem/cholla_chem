## Resolvers

placeholder_name uses a variety of resolvers to convert natural language chemical names to SMILES. Resolvers can be initialized and passed to the function resolve_compounds_to_smiles to customize how compounds are resolved to SMILES:

```
from placeholder_name import OpsinNameResolver, PubChemNameResolver, CIRPyNameResolver

opsin_resolver = OpsinNameResolver('opsin', resolver_weight=4)
pubchem_resolver =  PubChemNameResolver('pubchem', resolver_weight=3)
cirpy_resolver = CIRPyNameResolver('cirpy', resolver_weight=2)

resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'], [opsin_resolver, pubchem_resolver, cirpy_resolver])
```

## OpsinNameResolver
This resolver uses a fork of the [py2opsin](https://github.com/csnbritt/py2opsin) library that returns the error message from OPSIN if a name cannot be resolved. This resolver can be configured with the following arguments:

Arguments:
- allow_acid (bool, optional): Allow interpretation of acids. Defaults to False.
- allow_radicals (bool, optional): Enable radical interpretation. Defaults to False.
- allow_bad_stereo (bool, optional): Allow OPSIN to ignore uninterpreatable stereochem. Defaults to False.
- wildcard_radicals (bool, optional): Output radicals as wildcards. Defaults to False.

Default weight for 'weighted' SMILES selection: 3

```
from placeholder_name import OpsinNameResolver
opsin_resolver = OpsinNameResolver(
                                    'opsin',
                                    allow_acid=False,
                                    allow_radicals: True,
                                    allow_bad_stereo: False,
                                    wildcard_radicals: False
                                )

resolved_smiles = resolve_compounds_to_smiles(['2-acetyloxybenzoic acid'], [opsin_resolver])
```

## PubChemNameResolver
This resolver uses a fork of the [PubChemPy](https://github.com/csnbritt/PubChemPy) library which implements batching with the Power User Gateway XML schema to significantly speedup SMILES resolutions.

Default weight for 'weighted' SMILES selection: 2

```
from placeholder_name import PubChemNameResolver
pubchem_resolver = PubChemNameResolver('pubchem')

resolved_smiles = resolve_compounds_to_smiles(['acetone'], [pubchem_resolver])
```

## CIRPyNameResolver
This resolver uses the python library [CIRpy](https://github.com/mcs07/CIRpy), a Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.

Default weight for 'weighted' SMILES selection: 1

```
from placeholder_name import CIRPyNameResolver
cirpy_resolver = CIRPyNameResolver('cirpy')

resolved_smiles = resolve_compounds_to_smiles(['acetone'], [cirpy_resolver])
```

## ManualNameResolver 
This resolver uses a dataset of manually curated names and their corresponding SMILES, especially focused on common names that are incorrectly resolved by other resolvers (e.g. 'NAH').

Default weight for 'weighted' SMILES selection: 10

```
from placeholder_name import ManualNameResolver
manual_resolver = ManualNameResolver('manual')

resolved_smiles = resolve_compounds_to_smiles(['NAH'], [manual_resolver])
```

ManualNameResolver can be initialized with a custom dictionary mapping chemical names to SMILES:
```
from placeholder_name import ManualNameResolver

custom_name_dict = {'Foobar': 'c1ccccc1'}
manual_resolver = ManualNameResolver('manual', provided_name_dict=custom_name_dict)

resolved_smiles = resolve_compounds_to_smiles(['Foobar'], [manual_resolver])
```

## PeptideNameResolver
This resolver converts shorthand peptide names (e.g. 'cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)') to an IUPAC-like name, then attempts to resolve to SMILES with py2opsin.

Default weight for 'weighted' SMILES selection: 3

```
from placeholder_name import PeptideNameResolver
peptide_shorthand_resolver = PeptideNameResolver('peptide')

resolved_smiles = resolve_compounds_to_smiles(['cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)'], [peptide_shorthand_resolver])
```

## StructuralFormulaNameResolver
This resolver converts simple structural chemical formulas (e.g. 'CH3CH2CH2COOH') to SMILES.

Default weight for 'weighted' SMILES selection: 2

```
from placeholder_name import StructuralFormulaNameResolver
structural_formula_resolver = StructuralFormulaNameResolver('structural_formula')

resolved_smiles = resolve_compounds_to_smiles(['CH3CH2CH2COOH'], [structural_formula_resolver])
```