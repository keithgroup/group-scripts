#!/usr/bin/env python3

# MIT License
# 
# Copyright (c) 2021, Alex M. Maldonado
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

__version__ = '0.0.1'

"""
# mbGDML Data Set Creator

Creates a mbGDML data set (https://keithgroup.github.io/mbGDML/tutorials/datasets.html)
from QCJSONs using a modified QCSchema (https://github.com/MolSSI/QCSchema) format
(https://github.com/keithgroup/group-scripts/blob/master/scripts/qcjson-creator.py).

## Requirements
- numpy
- mbgdml
"""

import os
import argparse
import json
import numpy as np
from mbgdml.data import mbGDMLDataset

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
                files.append(os.path.abspath(f))
    return files

def check_jsons(jsons):
    """Several checks to confirm all JSONs will be compatible.

    Parameters
    ----------
    jsons : :obj:`list` [:obj:`dict`]
        All loaded JSONs to be combined into data set.
    """
    
    # Data to check.
    driver = 'gradient'
    atom_nums = jsons[0]['atomic_numbers']

    # Loop through all jsons and assert checks.
    for i in jsons[1:]:
        assert i['driver'] == driver
        assert i['atomic_numbers'] == atom_nums

def get_json_data(json_dict):
    """

    Parameters
    ----------
    json_dict : :obj:`dict`
        QCJSON dictionary.

    Returns
    -------
    :obj:`list` [:obj:`int`]
        Atomic numbers of all atoms
    :obj:`list` [:obj:`float`]
        Cartesian coordinates.
    :obj:`float`
        Total energy of the structure
    :obj:`list` [:obj:`float`]
        Gradients of the structure.
    """
    atom_nums = json_dict['atomic_numbers']
    coords = json_dict['molecule']['geometry']
    energy = json_dict['properties']['return_energy']
    gradient = json_dict['return_result']
    return atom_nums, coords, energy, gradient

def dataset_from_jsons(json_files, dataset_name):
    """

    Parameters
    ----------
    json_files : :obj:`list`
        List of loaded JSON files. Each item could be a list of dictionaries or
        a dictionary.
    """
    theory = ''
    atom_nums = []
    coords = []
    energies = []
    gradients = []
    for json_dict in json_files:
        # Iterative jobs are stored as a list of QCJSON files.
        if isinstance(json_dict, list):
            for i in json_dict:
                json_data = get_json_data(i)
                if len(atom_nums) == 0:
                    atom_nums = json_data[0]
                else:
                    assert atom_nums == json_data[0]
                coords.append(json_data[1])
                energies.append(json_data[2])
                gradients.append(json_data[3])
                method = i['model']['method']
                basis_set = i['model']['basis']
                dict_theory = f'{method}/{basis_set}'
                if theory == '':
                    theory = dict_theory
                else:
                    assert theory.lower() == dict_theory.lower()
        elif isinstance(json_dict, dict):
            json_data = get_json_data(json_dict)
            if len(atom_nums) == 0:
                atom_nums = json_data[0]
            else:
                assert atom_nums == json_data[0]
            coords.append(json_data[1])
            energies.append(json_data[2])
            gradients.append(json_data[3])
        else:
            raise TypeError

    z = np.array(atom_nums)
    R = np.array(coords)
    E = np.array(energies)
    F = np.negative(np.array(gradients))

    gdml_dataset = mbGDMLDataset()
    gdml_dataset.name = dataset_name
    gdml_dataset.z = z
    gdml_dataset.R = R
    gdml_dataset.r_unit = 'Angstrom'
    gdml_dataset.E = E
    gdml_dataset.e_unit = 'hartree'
    gdml_dataset.convertE('kcal/mol')
    gdml_dataset.F = F
    gdml_dataset.convertF('hartree', 'Angstrom', 'kcal/mol', 'Angstrom')
    gdml_dataset.theory = theory

    return gdml_dataset


def main():

    parser = argparse.ArgumentParser(
        description='Creates NumPy npz data set using the mbGDML format from '
                    'QCJSONs from the qcjson-creator.py script. '
                    'Should only be used for systems containing the same atoms.'
    )
    parser.add_argument(
        'search_dir', metavar='search_dir', type=str, nargs='?', default='.',
        help='Path to start searching for QCJSONs. Defaults to current direcotry.'
    )
    parser.add_argument(
        '--name', metavar='name', type=str, nargs='?', default='dataset',
        help='File name for the npz data set. Defaults to dataset.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save npz data set. Defaults to current directory.'
    )
    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Recursively search for QCJSONs'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite npz data set.'
    )
    parser.add_argument(
        '--remove', metavar='remove_files', nargs='+', default='',
        help='Ignore paths that contain at least one of these words.'
    )

    args = parser.parse_args()

    print('Making data sets from QCJSONs')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)\n')

    search_string = '.json'

    # Ensures paths end in a '/'.
    search_dir = args.search_dir
    if search_dir[-1] != '/':
        search_dir += '/'
    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    
    # Checks to see if data set already exists.
    if os.path.isfile(f'{save_dir}{args.name}.npz') and not args.overwrite:
        print(f'{save_dir}{args.name}.npz already exists and overwrite is False.\n')
        raise FileExistsError

    # Finds all JSON files.
    all_json_paths = get_files(search_dir, search_string, recursive=args.recursive)
    print(f'Found {len(all_json_paths)} QCJSONs')
    
    # Removes files that match any of the words in the remove list.
    if len(args.remove) > 0:
        start_number = len(all_json_paths)
        for trigger in args.remove:
            all_json_paths = [i for i in all_json_paths if trigger not in i]
        end_number = len(all_json_paths)
        print(f'Removed {start_number-end_number} file(s) for a total of {end_number} QCJSONs')

    # Loads all JSON files.
    all_jsons = []
    for path in all_json_paths:
        with open(path, 'r') as f:
            all_jsons.append(json.load(f))

    
    # Data checks before proceeding.
    print('Checking compatibility of all QCJSONs ... ', end="")
    check_jsons(all_jsons)
    print('passed')
    
    print(f'Creating the data set ... ', end="")
    gdml_dataset = dataset_from_jsons(all_jsons, args.name)
    gdml_dataset.save(
        gdml_dataset.name, gdml_dataset.dataset, save_dir
    )
    print('done\n')

    print(f'Your data set is: {save_dir}{args.name}.npz')

if __name__ == "__main__":
    main()