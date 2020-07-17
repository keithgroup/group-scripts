# submit-calcs

This Python script finds all calculation submission files and executes them.
Designed for use with the [slurm workload manager](https://slurm.schedmd.com/overview.html) and uses the command `sbatch <file_name>` to submit jobs.

## Use

The only required argument is the path to a folder containing the desired calculations to submit.
```
python submit-calcs.py path/to/calculation/dir/
```
or
```
submit-calcs path/to/calculation/dir/
```

Submissions scripts are found by searching for files that contain the string specified by the `--submit_string` argument in all nested directories; the default is `slurm`.

## Example

```
$ submit-calcs.py example-calculations/ --submit_string slurm
Submitting job in example-calculation-r folder
Submitting job in example-calculation-ts folder
Submitting job in example-calculation-p folder
```