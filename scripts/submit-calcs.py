#!/usr/bin/env python3


# MIT License
# 
# Copyright (c) 2020, Alex M. Maldonado
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import os
import subprocess
import argparse


def get_files(path, expression, recursive=True):
    """Returns paths to all files in a given directory that matches a provided
    expression in the file name. Commonly used to find all files of a certain
    type, e.g. output or xyz files.
    
    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in 'path'.
    recursive :obj:`bool`, optional
        Recursively find all files in all subdirectories.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    if path[-1] != '/':
        path += '/'
    if recursive:
        all_files = []
        for (dirpath, _, filenames) in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != '/':
                    dirpath += '/'
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files


def main():
    
    parser = argparse.ArgumentParser(
        description='Submit calculations using the slurm command "sbatch"'
    )
    parser.add_argument(
        'dirs', metavar='dirs', nargs='+', default='',
        help='directories containing calculations to be submitted'
    )
    parser.add_argument(
        '--submit_string', nargs=1, default='slurm',
        help='search string for submission scripts, defaults to slurm'
    )
    args = parser.parse_args()

    # Submits calculations in every direction
    for directory in args.dirs:
        os.chdir(directory)
        all_jobs_list = get_files(directory, args.submit_string[0])
        for job_path in all_jobs_list:

            slurm_file = job_path.split('/')[-1]
            job_directory = '/'.join(job_path.split('/')[:-1])

            print('Submitting job in ' + job_path.split('/')[-2] + ' folder')
            os.chdir(job_directory)
            bash_command = 'sbatch ' + slurm_file
            process = subprocess.Popen(
                bash_command.split(), stdout=subprocess.PIPE
            )
            _, error = process.communicate()

            if error is not None:
                print('Job not submitted successfully. Error:')
                print(str(error))


if __name__ == "__main__":
    main()