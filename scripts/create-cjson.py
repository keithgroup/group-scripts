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
from functools import reduce
import argparse
import cclib
from cclib.io.cjsonwriter import NumpyAwareJSONEncoder
import numpy as np
import json

### Hard coded information ###

# Level of theories categorized as "dft" or "ab initio"
theories = {
    'dft': ['pbe', 'bp86', 'b3lyp', 'wb97x-d3bj', 'wb97xd'],
    'ab initio': ['rhf', 'hf', 'uhf', 'mp2', 'dlpno-ccsd(t)', 'ccsd(t)']
}

# ORCA keywords D3 and D3BJ in version 4.2.0 and later both mean D3BJ.
orca_dispersion = ['d3', 'd3bj']

### Utility functions ###

def standard_dir(path):
    """Ensures directory path ends in '/'
    
    Parameters
    ----------
    path : str
        Path to directory with or without '/' at the end.
    
    Returns
    -------
    str
        Path to directory with '/' at the end.
    """

    if path[-1] != '/':
        path += '/'
    
    return path


def get_files(path, expression):
    """Returns paths to all files in a given directory that matches a provided
    expression in the file name. Commonly used to find all files of a certain
    type, e.g. output or xyz files.
    
    Parameters
    ----------
    path : str
        Specifies the directory to search.
    expression : str
        Expression to be tested against all file names in 'path'.
    
    Returns
    -------
    list
        All absolute paths to files matching the provided expression.
    """
    if path[-1] != '/':
        path += '/'
    
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
    for file in all_files:
        if expression in file:
            files.append(file)

    return files


def read_cjson(cjson_path):
    """Read JSON file.
    
    Parameters
    ----------
    cjson_path : str
        Path to json file.
    
    Returns
    -------
    dict
        Contents of JSON file.
    """

    with open(cjson_path, 'r') as reader:
        cjson_dict = json.loads(reader.readline())
    
    return cjson_dict


### Writing individual json files ###

def parse_dispersion(out_file_path, package):
    """Parses and returns dispersion in eV from output file.

    Parameters
    ----------
    out_file_path: str
        Path to output file.
    package: str
        Computational chemistry package of the output file. Only 'orca' and
        'gaussian' are implemented.

    Returns
    -------
    float
        Parsed dispersion energy in eV.
    """

    if package.lower() == 'orca':
        with open(out_file_path, 'r') as out_file:
            for line in out_file:
                if 'Dispersion correction' in line and 'DFT' not in line:
                    line_split = line.split()
                    dispersion_energy = float(line_split[-1])
                    break
    elif package.lower() == 'gaussian':
        return 0
    
    dispersion_energy = cclib.parser.utils.convertor(
        dispersion_energy, 'hartree', 'eV'
    )
    
    return dispersion_energy


def add_dispersion(cjson_dict, dispersion_energy):
    """Adds dispersion information to cjson dict.

    Cclib does not parse dispersion energy consistently across all
    computational chemistry programs. This is planned for a future release
    (post 1.6.3), but this is a quick fix to add dispersion energies
    consistently.

    Parameters
    ----------
    cjson_dict : dict
        Dictionary of chemicalJSON format.
    dispersion_energy : float
        Dispersion energy of structure in eV.
    
    Returns
    -------
    dict
        ChemicalJSON dict with added dispersion information.
    """

    scf_energy = cjson_dict['properties']['energy']['scf']
    total_energy = scf_energy + dispersion_energy

    cjson_dict['properties']['energy']['dispersion'] = dispersion_energy
    cjson_dict['properties']['energy']['total'] = total_energy

    return cjson_dict


def add_filename(cjson_dict, name):
    """Add filename to chemicalJSON dict.

    Parameters
    ----------
    cjson_dict : dict
        Dictionary of chemicalJSON format.
    name : str
        File name to add to cjson_dict.

    Returns
    -------
    dict
        ChemicalJSON dict with added dispersion information.
    """

    cjson_dict['name'] = name

    return cjson_dict


def modify_cjson(cjson_dict, ccdata, filename):
    """Modifies cjson dict produced by cclib into a custom format.

    Cclib provides a framework to produce chemicalJSON files. This function
    renames, organizes, and modifies some information.

    Parameters
    ----------
    cjson_dict : dict
        ChemicalJSON dict produced by cclib.
    ccdata : cclib.ccdata
        Cclib ccdata object used to generate cjson_dict.
    filename : str
        Name of the computational chemistry file to add to cjson_dict.
    
    Returns
    -------
    dict
        Modified chemicalJSON dict.
    """

    # Removes unwanted data
    cjson_dict['atoms'].pop('orbitals', None)
    cjson_dict['atoms'].pop('mass', None)
    cjson_dict['atoms'].pop('core electrons', None)
    cjson_dict['properties'].pop('orbitals', None)
    cjson_dict['properties'].pop('partial charges', None)
    cjson_dict['optimization']['scf'].pop('scf energies', None)

    # Modifies data
    cjson_dict['atoms']['coords']['3d'] = np.reshape(
        np.array(cjson_dict['atoms']['coords']['3d']),
        (cjson_dict['atoms']['elements']['atom count'], 3)
    )
    if ccdata.metadata['package'] == 'ORCA':
        cjson_dict['optimization']['scf']['targets'] = {
            'energy change': cjson_dict['optimization']['scf']['targets'][0][0],
            'max density change': cjson_dict['optimization']['scf']['targets'][0][1],
            'rms density change': cjson_dict['optimization']['scf']['targets'][0][2]
        }
        cjson_dict['optimization']['scf']['values'] = {
            'energy change': cjson_dict['optimization']['scf']['values'][0][:,0],
            'max density change': cjson_dict['optimization']['scf']['values'][0][:,1],
            'rms density change': cjson_dict['optimization']['scf']['values'][0][:,2]
        }
        cjson_dict['convergence'] = cjson_dict.pop('optimization')
    cjson_dict['properties']['energy']['scf'] = cjson_dict['properties']['energy'].pop('total')
    
    # Adds data
    cjson_dict['properties']['energy']['energy units'] = 'eV'
    cjson_dict['input parameters'] = {}
    cjson_dict['input parameters']['package'] = ccdata.metadata['package'].lower()
    cjson_dict['input parameters']['version'] = ccdata.metadata['package_version'].lower()

    # Adds calculation information to cjson parsed from ccdata.
    if ccdata.metadata['package'] == 'ORCA':
        avail_keywords = ccdata.metadata['keywords']
        uncat_keywords = []
        for keyword in avail_keywords:
            keyword = keyword.lower()
            if keyword in theories['dft']:
                cjson_dict['input parameters']['theory'] = 'dft'
                cjson_dict['input parameters']['functional'] = keyword
                if 'd3bj' in keyword:
                    cjson_dict['input parameters']['dispersion'] = 'd3bj'
            elif keyword in theories['ab initio']:
                cjson_dict['input parameters']['theory'] = keyword
            elif 'def2' in keyword:
                cjson_dict['input parameters']['basis'] = keyword
            elif 'sp' == keyword:
                cjson_dict['input parameters']['task'] = 'energy'
            elif 'opt' == keyword:
                cjson_dict['input parameters']['task'] = 'optimization'
            elif 'freq' == keyword:
                cjson_dict['input parameters']['task'] = 'frequency'
            elif 'cpcm(water)' in keyword:
                if ccdata.metadata['package'] == 'ORCA':
                    if 'smd' in ccdata.metadata['input_file_contents'].lower():
                        cjson_dict['input parameters']['solvent model'] = 'SMD'
                        cjson_dict['input parameters']['solvent'] = 'water'
                    else:
                        cjson_dict['input parameters']['solvent model'] = 'CPCM'
                        cjson_dict['input parameters']['solvent'] = 'water'
            elif keyword in orca_dispersion:
                cjson_dict['input parameters']['dispersion'] = 'd3bj'
            else:
                uncat_keywords.append(keyword)
        cjson_dict['input parameters']['options'] = uncat_keywords
    elif ccdata.metadata['package'] == 'Gaussian':
        cjson_dict['input parameters']['theory'] = ccdata.metadata['methods'][0].lower()
        cjson_dict['input parameters']['functional'] = ccdata.metadata['functional'].lower()
        cjson_dict['input parameters']['basis'] = ccdata.metadata['basis_set'].lower()
        cjson_dict['input parameters']['task'] = 'energy'
        cjson_dict['input parameters']['dispersion'] = 'd3'  # Check

        # Solvent information
        if 'smd' in filename.lower():
            cjson_dict['input parameters']['solvent model'] = 'SMD'
            cjson_dict['input parameters']['solvent'] = 'water'
        elif 'cpcm' in filename.lower():
            cjson_dict['input parameters']['solvent model'] = 'CPCM'
            cjson_dict['input parameters']['solvent'] = 'water'
        elif 'ipcm' in filename.lower():
            cjson_dict['input parameters']['solvent model'] = 'IPCM'
            cjson_dict['input parameters']['solvent'] = 'water'

    # Checks that essential information was added correctly.
    if 'theory' not in cjson_dict['input parameters']:
        raise KeyError('No theory was identified.')
    if 'basis' not in cjson_dict['input parameters']:
        raise KeyError('No basis set was identified.')
    if 'task' not in cjson_dict['input parameters']:
        raise KeyError('No calculation type was identified.')

    return cjson_dict


def create_cjson_dict(out_path, file_name):
    """Creates single chemical JSON format dictionary.

    Parameters
    ----------
    out_path : str
        Path to computational chemistry output file parsable by cclib.
    file_name : str
        Name of the computational chemistry output file.
    
    Returns
    -------
    dict
        Chemical JSON dictionary.
    """

    ccdata = cclib.io.ccread(out_path)
    cclib_cjson = cclib.io.CJSONWriter(ccdata)
    cclib_cjson_dict = cclib_cjson.as_dict()
    cjson_dict = modify_cjson(cclib_cjson_dict, ccdata, file_name)
    cjson_dict = add_filename(cjson_dict, file_name)

    # Adds dispersion information if necessary
    if 'dispersion' in cjson_dict['input parameters']:
        disp_energy = parse_dispersion(
            out_path, cjson_dict['input parameters']['package']
        )
        if disp_energy != 0:  # ORCA
            cjson_dict = add_dispersion(cjson_dict, disp_energy)
        else:  # Gaussian
            cjson_dict['properties']['energy']['total'] = cjson_dict['properties']['energy'].pop('scf')
    
    # Renames ab initio 'scf' energy to 'total' energy (no dispersion)
    if cjson_dict['input parameters']['theory'] in theories['ab initio']:
        cjson_dict['properties']['energy']['total'] = cjson_dict['properties']['energy'].pop('scf')

    return cjson_dict


def write_all_cjson(output_dir, save_dir, remove, overwrite):
    """ Driver function to write all chemical JSON files.

    Parameters
    ----------
    output_dir : str
        Path to directory containing computational chemistry output files
        parsable by cclib.
    save_dir : str
        Path to directory where chemical JSON files will be saved.
    remove : list
        Files with paths containing strings in this list will not be included.
    overwrite : bool
        Chemical json files will be overwritten.
    """

    outfile_paths = get_files(output_dir, 'out')
    for string in remove:
        outfile_paths = [i for i in outfile_paths if string not in i]

    # Creates chemicalJSON files
    for out_path in outfile_paths:
        # Simple logging.
        file_name = '.'.join(out_path.split('/')[-1].split('.')[:-1])
        print(f'Writing {file_name} chemicalJSON file ...')

        # Creates the cjson dict
        cjson_dict = create_cjson_dict(out_path, file_name)

        # Write cjson
        org_dirs = standard_dir('/'.join(out_path.split('/')[:-1]).replace(output_dir, ''))

        json_dir = f'{save_dir}{org_dirs}/'
        json_path = f'{json_dir}{file_name}.json'
        
        os.makedirs(json_dir, exist_ok=True)
        
        if not os.path.isfile(json_path) or overwrite:
            with open(json_path, 'w') as cjson:
                json.dump(
                    cjson_dict, cjson, cls=NumpyAwareJSONEncoder,
                    sort_keys=True
                )
        else:
            print('File already exists.')


### Writing cumulative json file ###

def combined_cjson_name(cjson_dict, parse_name):
    """Custom naming function for cumulative chemical JSON file keys.

    Chemical JSON files typically named by
    "system-structure-program.calctype-functional.basis.options".
    Usually, "basis.options" (e.g., def2svp.cpcm) is used as the
    key for the individual chemical JSON dictionary.

    Parameters
    ----------
    cjson_dict : dict
        cjson dictionary that must have 'name' key.
    
    Returns
    -------
    str
        Key for individual chemical JSON dictionary.
    """
    
    json_name = cjson_dict['name']

    if parse_name:
        json_name = '.'.join(json_name.split('-')[-1].split('.')[1:])

    return json_name


def combine_cjson(save_dir, name, parse_name):
    """Writes cumulative chemical JSON file containing all files.

    Parameters
    ----------
    save_dir : str
        Path to directory where chemical JSON files will be saved.
    name : str
        File name for cumulative chemical JSON file.
    """
    all_cjson_dict = {}

    start = save_dir.rfind(os.sep) + 1
    for path, _, files in os.walk(save_dir):

        files = [i for i in files if '.json' in i]
        folders = path[start:].split(os.sep)
        try:
            cjson_files = {}
            for cjson_file in files:
                cjson_path = path + '/' + cjson_file
                cjson_dict = read_cjson(cjson_path)
                cjson_name = combined_cjson_name(cjson_dict, parse_name)
                cjson_files[cjson_name] = cjson_dict
            parent = reduce(dict.get, folders[:-1], all_cjson_dict)
            parent[folders[-1]] = cjson_files
        except KeyError:
            # Will not include JSON files without 'name' key.
            # Useful for using --overwrite when cumulative file is present.
            pass
    
    with open(f'{save_dir}{name}.json', 'w') as json_file:
        json.dump(all_cjson_dict, json_file)


### Main ###

def main():
    
    parser = argparse.ArgumentParser(
        description='Create chemical json files from output files. Files must '
                    'have "out" somewhere in the file name.'
    )
    parser.add_argument(
        'output_dir', metavar='output_dir', type=str, nargs=1,
        help='path to directory with computational chemistry output files'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs=1, default='.',
        help='path to save cjson file'
    )
    parser.add_argument(
        '--remove', metavar='remove_files', nargs='+', default='',
        help='ignore file paths matching the words in this string'
    )
    parser.add_argument(
        '--name', metavar='name', nargs=1, default='data',
        help='name of final combined file'
    )
    parser.add_argument(
        '--overwrite', action='store_true', help='overwrite json files'
    )
    parser.add_argument(
        '--parse_name', action='store_true',
        help='Rename JSON file from parsed information from output file name'
        ' in the format '
        'system-structure-program.calctype-functional.basis.options.out'
        'Output files should be organized in directories by functional '
        'with basis subdirectories.'
    )
    args = parser.parse_args()
    
    output_dir = standard_dir(args.output_dir[0])
    save_dir = standard_dir(args.save_dir[0])
    if not os.path.isdir(output_dir):
        raise ValueError(f'{output_dir} is not a valid directory.')
    if not os.path.isdir(save_dir):
        raise ValueError(f'{save_dir} is not a valid directory.')

    write_all_cjson(output_dir, save_dir, args.remove, args.overwrite)
    combine_cjson(save_dir, args.name, args.parse_name)


if __name__ == "__main__":
    main()
