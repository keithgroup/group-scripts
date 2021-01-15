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
# ABCluster Cleaner

We often use ABCluster to generate low-energy clusters containing solutes and/or
solvent molecules. It leaves multiple files for the same structure, and when
water is included it still lists the lone pairs.

This script will remove all duplicate structures files (.cluster and .gjf;
leaving the .xyz files) and remove all lone pair entries.

## Requirements
* NumPy
"""

import os
import argparse
import numpy as np
from periodictable import elements

### Global Information ###

file_extensions = ['.xyz', '.cluster', '.gjf']

element_to_z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
    'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31,'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
    'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
    'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
    'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93,
    'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106,
    'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Uuq': 114, 'Uuh': 116,
}
z_to_element = {v: k for k, v in element_to_z.items()}

### Utilities ###

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

def parse_stringfile(stringfile_path):
    """Parses data from string file.

    A string file is data presented as consecutive xyz data. The data could be
    three Cartesian coordinates for each atom, three atomic force vector
    components, or both coordinates and atomic forces in one line (referred to
    as extended xyz).
    
    Parameters
    ----------
    stringfile_path : :obj:`str`
        Path to string file.
    
    Returns
    -------
    :obj:`tuple` [:obj:`list`]
        Parsed atoms (as element symbols :obj:`str`), comments, and data as
        :obj:`float` from string file.
    """
    z, comments, data = [], [], []
    with open(stringfile_path, 'r') as f:
        for _, line in enumerate(f):
            line = line.strip()
            if not line:
                # Skips blank lines
                pass
            else:
                line_split = line.split()
                if len(line_split) == 1 \
                    and float(line_split[0]) % int(line_split[0]) == 0.0:
                    # Skips number of atoms line, adds comment line, and
                    # prepares next z and data item.
                    comment_line = next(f)
                    comments.append(comment_line.strip())
                    z.append([])
                    data.append([])
                else:
                    # Grabs z and data information.
                    z[-1].append(line_split[0])
                    data[-1].append([float(i) for i in line_split[1:]])
    return z, comments, data

def string_coords(z, R, decimal_points=9):
    """Puts atomic coordinates into a Python string. Typically used for 
    writing to an input file.
    
    Parameters
    atoms : :obj:`numpy.ndarray`
        A (n,) numpy array containing all ``n`` elements labled by their atomic 
        number.
    coords : :obj:`numpy.array`
        Contains atomic positions in a (n, 3) numpy array where the x, y, and z 
        Cartesian coordinates in Angstroms are given for the n atoms.
    decimal_points : :obj:`int`, optional
        Total number of decimal points to include. Defaults to 9.
    
    Returns
    -------
    :obj:`str`
        XYZ atomic coordinates as a string.
    """
    atom_coords_string = ''
    atom_index = 0
    while atom_index < len(z):
        atom_element = str(elements[z[atom_index]])
        format_string = '%.' + str(decimal_points) + 'f'
        coords_string = np.array2string(
            R[atom_index],
            suppress_small=True, separator='   ',
            formatter={'float_kind': lambda x: format_string % x}
        )[1:-1] + '\n'
        if len(atom_element) == 2:
            atom_spacing = '   '
        else:
            atom_spacing = '    '
        atom_coords_string += (atom_element + atom_spacing \
                               + coords_string).replace(' -', '-')
        atom_index += 1
    
    return atom_coords_string

def write_xyz(z, R, save_dir, file_name, comments=None, decimal_points=9):
    """Write XYZ file given atomic numbers and Cartesian coordinates.

    Parameters
    ----------
    z : :obj:`numpy.ndarray`
        A ``(n,)`` shape array of type :obj:`numpy.int32` containing atomic
        numbers of atoms in the structures in order as they appear.
    R : :obj:`numpy.ndarray`
        A :obj:`numpy.ndarray` with shape of ``(m, n, 3)`` where ``m`` is the
        number of structures and ``n`` is the number of atoms with three 
        Cartesian components.
    save_dir : :obj:`str`
        Path to directory to save XYZ file.
    file_name : :obj:`str`
        Name to save the file.
    comments : :obj:`list` [:obj:`str`], optional
        List of comments for each xyz structure.
    decimal_points : :obj:`int`, optional
        Total number of decimal points to include. Defaults to 9.
    """
    if R.ndim == 2:
        R = np.array([R])
    
    if save_dir[-1] != '/':
        save_dir += '/'
    
    num_atoms = len(z)
    with open(f'{save_dir}{file_name}.xyz', 'w') as f:
        for i in range(0, R.shape[0]):
            if comments is None:
                comment = '\n'
            else:
                if comments[i][:-2] != '\n':
                    comments[i] += '\n'
                comment = comments[i]
            f.writelines(
                [
                    f'{num_atoms}\n', comment,
                    string_coords(z, R[i], decimal_points=decimal_points)
                ]
            )

### Main ###

def main():

    parser = argparse.ArgumentParser(
        description='Partitions structure(s) for energy+gradient calculations.'
    )
    parser.add_argument(
        'lm_dir', metavar='lm_dir', type=str, nargs='?', default='.',
        help='Path to local minima directory from ABCluster.'
    )
    parser.add_argument(
        '--decimals', metavar='decimals', type=int, nargs='?', default=8,
        help='Number of decimal points to write in xyz files. Defaults to 8.'
    )

    args = parser.parse_args()

    print(f'\nABCluster Cleaner v{__version__}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)\n')

    # Ensures paths end in a '/'.
    lm_dir = args.lm_dir
    if lm_dir[-1] != '/':
        lm_dir += '/'
    
    # Gets all relevant files.
    abc_files = []
    for extension in file_extensions:
        abc_files.extend(get_files(lm_dir, extension, recursive=False))
    print(f'Found {len(abc_files)} files')
    
    # Cleans up directory
    removed_files = 0
    cleaned_files = 0
    print('Working on files ...')
    for abc_file in abc_files:
        extension = abc_file.split('.')[-1]

        # Remove files ending in '.cluster' or '.gjf'.
        if extension == 'cluster' or extension == 'gjf':
            os.remove(abc_file)
            removed_files += 1
        
        # Cleans and writes xyz files
        if extension == 'xyz':
            z, comments, data = parse_stringfile(abc_file)

            # Removes all coordinates for atom labels not in the element list.
            for structure in range(0, len(data)):
                indexes = []
                for i in range(0, len(z[structure])):
                    if z[structure][i] not in element_to_z.keys():
                        indexes.append(i)
                for i in sorted(indexes, reverse=True):
                    del z[structure][i]
                    del data[structure][i]
            
            # Checks if all structures have atoms in the same order.
            if len(set(tuple(i) for i in z)) == 1:
                z = z[0]
            else:
                print(abc_file)
                raise ValueError('Does not support XYZ files with multiple '
                                 'structures not with same atom order.')
            
            # Writes xyz file.
            abc_file_name = os.path.splitext(os.path.basename(abc_file))[0]
            z = np.array([element_to_z[i] for i in z])
            write_xyz(
                z, np.array(data), lm_dir, abc_file_name,
                comments=comments, decimal_points=args.decimals
            )
            cleaned_files += 1
    
    print(f'Removed {removed_files} files')
    print(f'Cleaned {cleaned_files} files')

if __name__ == "__main__":
    main()