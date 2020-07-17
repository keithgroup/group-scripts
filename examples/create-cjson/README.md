# create-cjson

This Python script creates a chemical JSON, cjson, data base for computational chemistry output files.

Reproducibility is a cornerstone of science, and computational chemistry is not different.
In decreasing level of transparency and reproducibility in computational chemistry, researchers should be able to provide source codes, raw output files, input files, or some sort of machine readable database (XML and JSON). 
Since commercial software is commonly used, this incentivises providing output files and a machine readable database.

This is an easy-to-use script for creating a JSON database of computational chemistry output files.
We employ a [chemical JSON](https://github.com/OpenChemistry/chemicaljson) format with additional information.

## Requirements


## Chemical JSON Format


## Use



## Example

```
python create-cjson.py example-calculations/ --save_dir example-cjson/
Writing out-example-calculation-p chemicalJSON file ...
Writing out-example-calculation-r chemicalJSON file ...
Writing out-example-calculation-ts chemicalJSON file ...
```

