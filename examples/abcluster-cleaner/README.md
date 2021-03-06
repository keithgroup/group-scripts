# abcluster-cleaner

We often use ABCluster to generate low-energy clusters containing solutes and/or
solvent molecules. It leaves multiple files for the same structure, and when
water is included it still lists the lone pairs.

This script will remove all duplicate structures files (.cluster and .gjf;
leaving the .xyz files) and remove all lone pair entries.

## Requirements

* NumPy
  
## cli help

```text
usage: abcluster-cleaner.py [-h] [--decimals [decimals]] [-r] [--no_remove]
                            [--exclude EXCLUDE [EXCLUDE ...]] [--include INCLUDE [INCLUDE ...]]
                            [lm_dir]

Partitions structure(s) for energy+gradient calculations.

positional arguments:
  lm_dir                Path to local minima directory from ABCluster.

optional arguments:
  -h, --help            show this help message and exit
  --decimals [decimals]
                        Number of decimal points to write in xyz files. Defaults to 8.
  -r, --recursive       Recursively clean XYZ files.
  --no_remove           Do not remove .gjf or .cluster files.
  --exclude EXCLUDE [EXCLUDE ...]
                        Ignore paths that contain at least one of these words.
  --include INCLUDE [INCLUDE ...]
                        Only include files that contain all of these words.
```

## Examples

If inside the `original` directory, you could clean the XYZ files by running the following command.

```text
abcluster-cleaner.py
```

The results are shown in the `cleaned` directory.