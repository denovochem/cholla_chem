# Name manipulation

This package contains utilities that automatically attempt to transform or correct input names if the initial name-to-SMILES conversion fails. At the moment there are three major capabilities:

- Splitting on delimiters (e.g. 'BH3•THF' → 'BH3' and 'THF')
- Peptide shorthand expansion (e.g. `Boc-Ala-Gly-OMe` → `tert-butoxycarbonyl-l-alanyl-glycine methyl ester`)
- Chemical name correction (primarily aimed at OCR/typo artifacts, with validation via OPSIN)

Each of these name manipulations is enabled by default, but can be disabled by setting the corresponding argument to False.

## Splitting on delimiters

Split chemical names on common delimiters and then solve each component separately.

```pycon
from cholla_chem import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(compounds_list=['BH3•THF'], split_names_to_solve=True)

"{'BH3•THF': 'B.C1CCOC1'}"
```


## Peptide shorthand expansion

Applied when a chemical name looks like peptide shorthand, typically containing amino-acid abbreviations delimited by hyphens, such as:

- `Ala-Gly-Ser`
- `Boc-Asp-Lys(Boc)-OMe`
- `cyclo(Ala-Gly-Ser)` or `cyclo[Ala-Gly-Ser]`

```pycon
from cholla_chem import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(compounds_list=['cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)'], resolve_peptide_shorthand=True)

"{'cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)': 'BN1[C@@H](CC(=O)O)C(=O)N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC2=CC=C(C=C2)O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CC2=CNC=N2)C(=O)N2[C@@H](CCC2)C(=O)N[C@@H](CC2=CC=CC=C2)C1=O'}"
```


## Chemical name correction

Applied when a chemical name looks like it might have OCR/typo artifacts, with validation via OPSIN.

- Examples of artifacts:
  - `3,4-dihydro- 3-hydroxy-4-oxo-l,2,3-benzotriazine` → `3,4-dihydro- 3-hydroxy-4-oxo-1,2,3-benzotriazine`
  - `ethyl perfiuorobutyrate` → `ethyl perfluorobutyrate`
  - `l-mercapto-Z-thiapropane` → `1-mercapto-2-thiapropane`

```pycon
from cholla_chem import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(compounds_list=['l-mercapto-Z-thiapropane'], attempt_name_correction=True)

"{'l-mercapto-Z-thiapropane': 'SCSC'}"
```
