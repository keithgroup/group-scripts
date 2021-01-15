# qcjson-creator

This Python script creates Quantum Chemistry JSON, QCJSON, files for computational chemistry output files using the [QCSchema](https://github.com/MolSSI/QCSchema).

Reproducibility is a cornerstone of science, and computational chemistry is not different.
In decreasing level of transparency and reproducibility in computational chemistry, researchers should be able to provide source codes, raw output files, input files, or some sort of machine readable database (XML and JSON). 
Since commercial software is commonly used, this incentivizes providing output files and a machine readable database.

This is an easy-to-use script for creating JSON databases of computational chemistry output files.

## Requirements

* cclib
* numpy

## Use


### Options

* `--overwrite`: Will overwrite JSON file if it already exists.
* `--save_dir`: Path to directory to store JSON files; for example, "`--save_dir json_dir`".
* `--name`: What to name the cumulative JSON file; for example, '`--name example-data`'.
Default is 'data'.
* `--remove`: Does not include file paths that have these strings in their path; for example, '`--remove bad bp86`' will not include any files with 'bad' or 'bp86' in the pathway.

### Example

The following command produces a `data.json` file from output files in the  `example-calculations` directory.

```
$ qcjson-creator.py . -p
QCJSON creator v0.1.1
Written by Alex M. Maldonado (@aalexmmaldonado)
Energies and distances are Hartrees and Angstroms

Looking for output files in ./
Found 6 output files

Making QCJSON for 4h2o.abc0.iter2-orca.opt-mp2.def2tzvp
Making QCJSON for 12h2o.su.etal-orca.sp-mp2.def2tzvp.frozencore
Making QCJSON for neb-bare.0-orca.sp-bp86.def2tzvp.smd
Making QCJSON for 5h2o.abc0.iter1.mp2.md.300k.iter1.mol0,1,2,3,4-orca.engrad-mp2.def2tzvp
Making QCJSON for neb-bare.0-orca.sp-wb97x.def2tzvp
Making QCJSON for neb-bare.0-orca.sp-bp86.def2svp

0 file(s) encountered errors and not written
```

## QCJSON Format

Custom modifications are made to [QCSchema](https://github.com/MolSSI/QCSchema) to meet the immediate needs of our research.

- [QCSchema](https://github.com/MolSSI/QCSchema) do not support multiple configurations in a single JSON
file: one file per structure and calculation. This is inherently incompatible
with geometry optimizations, trajectories, or any other iterative procedure.
At of the time of writing (2021-01-03), there has been no consensus of how to manage these files.
Our immediate solution is to list each [QCSchema](https://github.com/MolSSI/QCSchema) inside a list (i.e., [ { }, { }, { }, ... ] ). Furthermore, only iterations that have
the same topology are supported (no checks are provided).
- Other information with keywords not defined in the [QCSchema](https://github.com/MolSSI/QCSchema).

### Keys

Here are definitions of all [QCSchema](https://github.com/MolSSI/QCSchema) and custom keys used in this script.
Not all QCJSONs will have every key all the time; some are method dependent like the frozen core approximation.
All custom keys are marked with a *, and the Python type is specified.

* ``"schema_name"``: str, specifies the type of QCJSON. Always equal to "qc_schema_output".
* ``"schema_version"``: int, specifies the [QCSchema](https://github.com/MolSSI/QCSchema) version. Current version is 1.
* ``"qcjson_creator_version"``*: str, specifies the qcjson-creator script version.
* ``"provenance"``: dict, a brief description of the program, version, and routine used to generate the output. Can include more detailed information such as computation time, processor information, and host location.
    * ``"creator"``: str, name of the QC package.
    * ``"version"``: str, version of the QC package.
* ``"molecule"``: dict, a full description of the overall molecule(s) (e.g., its geometry, fragments, and charges).
    * ``"geometry"``: list [[x_float, y_float, z_float], [x_float, y_float, z_float], ...], vector of XYZ coordinates of the atoms of length equal to the number of atoms.
    * ``"symbols"``: list [str, str, ...], atom symbols in title case (e.g., "H" and "Na").
* ``"molecular_charge"``: int, the overall charge of the molecule.
* ``"molecular_multiplicity"``: int, The overall multiplicity of the molecule.
* ``"name"``: ``str``, Desired name of the molecule or system (currently is just the file name).
* ``"atomic_numbers"``: list [int, int, ...], atomic numbers, nuclear charge for atoms.
  Ghostedness should be indicated through ‘real’ field, not zeros here.
* ``"driver"``: str, What are you looking to calculate: energy, gradient, Hessian, or property.
  Note, we implement other driver options: optimization.
* ``"return_result"``: str, The “primary” return of a given computation.
  For energy, gradient, and Hessian quantities these are either single numbers or a array representing the derivative quantity.
* ``"model"``: dict, The overall mathematical model we are using for our calculation.
  Another way to think about this is the largest superset that still obtains roughly the same result.
  In QM, this is the Hamiltonian (HF, DFT, …) combined with the overall basis of the calculation.
  An example in QM would be HF/STO-3G or B3LYP/6-311G**.
  Custom basis sets can be handled with custom keywords.
    * ``"method"``: str, Primary method or level of theory.
    * ``"basis"``: str, Overall basis set.
* ``"keywords"``: dict, Various tunable parameters for the calculation.
  These vary widely, depending on the basis and model chemistry.
  These represent the keywords of individual programs currently.
    * ``"scf_convergence_tolerance"``*: str, Program-specific SCF convergence criteria.
      For example, "tight".
    * ``"dispersion"``*: str, If empirical dispersion corrections are used (with DFT), this specifies the method.
    * ``"scf_grid_level"``*: int, Specifies program-specific integration grid level for the scf cycle.
    * ``"final_grid_level"``*: int, Specifies program-specific integration grid level for the final energy calculation.
    * ``"implicit_solvent"``*: str, The implicit solvent model if used.
      For example, "SMD" or "CPCM".
    * ``"solvent_name"``*: str, Name of the solvent.
      For example, "water".
    * ``"rij_approximation"``*: bool, If the resolution of identity approximation for the Coulomb integrals (RI-J) is used.
    * ``"cosx_approximation"``*: bool, If the chain-of-spheres integration approximation to the exchange term (COSX) is used.
* ``"properties"``: dict, A list of valid quantum chemistry properties tracked by the schema.
    * ``"calcinfo_nbasis"``: int, The number of basis functions for the computation.

## Examples

In this directory, there are the following example ORCA 4.2.0 calculation output files with corresponding QCJSONs:

- `4h2o.abc0.iter2-orca.opt-mp2.def2tzvp.out`: An optimization using MP2/def2-TZVP with the default frozen core approximation.
- `neb-bare.0-orca.sp-bp86.def2svp.out`: A single-point energy calculation using BP86-D3BJ/def2-TZVP.
- `neb-bare.0-orca.sp-bp86.def2tzvp.smd.out`: A single-point energy calculation using BP86-D3BJ/def2-TZVP with a water implicit solvent model (SMD).
- `neb-bare.0-orca.sp-wb97x.def2tzvp.out`: A single-point energy calculation using $\omega$B97X-D3BJ/def2-TZVP.
- `5h2o.abc0.iter1.mp2.md.300k.iter1.mol0,1,2,3,4-orca.engrad-mp2.def2tzvp.out`: An ORCA job running multiple energy+gradient calculations on different configurations of the same system.

## Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2021-01-06
### Fixed
- cclib version requirement.
- Recursive option would save in current directory and not in the same directory
  of the output file.
- get_json would incorrectly catch KeyboardInterrupt exception.

## [0.1.0] - 2021-01-05
### Added
- Custom parser for ORCA information such as integration grid, scf energy
  contributions (e.g., one-electron and two-electron energies), MP2 correlation
  energies, RI approximations.
- Debug option to raise errors instead of skipping over files.
- Alpha and beta electron HOMO and LUMO information.
- 'return_energy' property regardless of driver.

### Changed
- Nest iterations into a list instead of having int labels.
- Standardized getting SCF, MP, and CC energies from cclib.
- Requires outfile path to initialize json classes.
- Write each JSON file directly after parsing instead of all at the end. That
  way if the script crashes the proceeding JSON files are already written.

## [0.0.1] - 2021-01-03
### Added
- Initial release!