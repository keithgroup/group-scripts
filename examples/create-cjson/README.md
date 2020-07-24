# create-cjson

This Python script creates a chemical JSON, cjson, data base for computational chemistry output files.

Reproducibility is a cornerstone of science, and computational chemistry is not different.
In decreasing level of transparency and reproducibility in computational chemistry, researchers should be able to provide source codes, raw output files, input files, or some sort of machine readable database (XML and JSON). 
Since commercial software is commonly used, this incentivises providing output files and a machine readable database.

This is an easy-to-use script for creating a JSON database of computational chemistry output files.
We employ a [chemical JSON](https://github.com/OpenChemistry/chemicaljson) format with additional information.

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
* `--parse_name`: Naming convention for specifically formatted naming and organizational scheme.
Calculations should be organized in nested directories by 'system', 'structure', 'program', and 'functional'. Then the file names should have the format of "system-structure-program.calctype-functional.basis.options".
For example, a single point energy calculation of a four water cluster with B3LYP/def2-TZVP and SMD continuum solvent model in ORCA would be stored as 'h2o/4/orca/b3lyp/h2o-4-orca.sp-b3lyp.def2tzvp.smd.out'.

### Example

The following command produces a `data.json` file from output files in the  `example-calculations` directory.

```
$ create-cjson example-calculations/ --save_dir example-cjson/ --overwrite
Writing example-p-orca.sp-wb97x.def2tzvp chemicalJSON file ...
Writing example-r-orca.sp-wb97x.def2tzvp chemicalJSON file ...
Writing example-ts-orca.sp-wb97x.def2tzvp chemicalJSON file ...
```

## Chemical JSON Format

### Organizational keys

* `dir_name`: directory names containing the output files.
Note that all subdirectories are nested keys, so if a file is stored in dir1/dir2/outfile_name.out then the chemical JSON file is stored under `['dir1']['dir2']`.
    * `outfile_name`: is the the file name of the output file.
    Note that this can be parsed using the `--parse_name` option if using the specific format of `system-structure-program.calctype-functional.basis.options`.
    The chemical JSON file for the example above would be under `['dir1']['dir2']['outfile_name']`.

### Keys

* `atoms`: Contains information describing the identity and position of all atoms in the system.
    * `coords`: Cartesian coordinates in Angstrom units.
    * `elements`:
* `chemical json`:
* `convergence`:
    * `scf`:
        * `targets`:
        * `values`:
* `input parameters`:
    * `basis`:
    * `dispersion`:
    * `functional`:
    * `options`:
    * `package`:
    * `task`:
    * `theory`:
    * `version`:
* `name`:
* `properties`:
    * `charge`:
    * `energy`:
    * `multiplicity`:
    * `number of atoms`:
    * `total dipole moment`:

### Example

```
{
    "dir_name": {
        "outfile_name": {
            "atoms": {
                "coords": {
                    "3d": [[1.208995, -0.056447, 2.319079],
                           [0.097037, -0.742342, 2.374182],
                           [2.038902, -0.132529, 1.43526],
                           [-0.542349, -0.041626, 5.719692],
                           [-0.575633, -1.229482, 6.069369],
                           [0.640706, 0.363733, 5.678614],
                           [-1.192228, 0.66647, 6.518487],
                           [1.341532, 0.636241, 3.19142],
                           [-1.20216, 0.053744, 4.329507],
                           [-0.465006, -0.449687, 3.231394],
                           [-1.366461, 0.987555, 4.156071]]
                },
                "elements": {
                    "atom count": 11,
                    "heavy atom count": 5,
                    "number": [6, 8, 8, 5, 1, 1, 1, 1, 8, 1, 1]
                }
            },
            "chemical json": 0,
            "convergence": {
                "scf": {
                    "targets": {
                        "energy change": 1e-08,
                        "max density change": 1e-07,
                        "rms density change": 5e-09
                    },
                    "values": {
                        "energy change": [
                            0.0, -0.027866510673, -0.016345801165,
                            -0.010261106153, -0.006958896803, -0.0161478416,
                            -3.98515e-05, -1.11236e-05, -1.4657e-06, -1.3176e-06,
                            -1.02e-08, -1.88e-08, -3.3043e-10
                        ],
                        "max density change": [
                            0.02126995, 0.01698518,
                            0.01289865, 0.01096771, 0.02845554, 0.000969,
                            0.000923, 0.001173, 0.000386, 4.4e-05, 1.9e-05,
                            9e-06, 4.6722e-06
                        ],
                        "rms density change": [
                            0.00062129, 0.00054238, 0.00039006, 0.00027431,
                            0.00068531, 0.002937, 0.001966, 0.001431, 0.000487,
                            7.4e-05, 3.2e-05, 7e-06, 1.3269e-07
                        ]
                    }
                }
            },
            "input parameters": {
                "basis": "def2-tzvp",
                "dispersion": "d3bj",
                "functional": "wb97x-d3bj",
                "options": ["grid4", "tightscf"],
                "package": "orca",
                "task": "energy",
                "theory": "dft",
                "version": "4.2.0"
            },
            "name": "example-p-orca.sp-wb97x.def2tzvp",
            "properties": {
                "charge": -1,
                "energy": {
                    "alpha": {
                        "gap": 11.7065011168203,
                        "homo": -5.0561202447554505
                    },
                    "beta": {
                        "gap": 11.7065011168203,
                        "homo": -5.0561202447554505
                    },
                    "dispersion": -0.2696803907577486,
                    "energy units": "eV",
                    "scf": -7959.819841441988,
                    "total": -7960.089521832745
                },
                "multiplicity": 1,
                "number of atoms": 11,
                "total dipole moment": 4.656141491506972
            }
        }
    }
}
```

