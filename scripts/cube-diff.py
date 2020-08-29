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

import argparse
import re

def calc_density_diff(cube1_value, cube2_value):
    """Computes formated density difference.

    Parameters
    ----------
    cube1_value : str
        Single value of electron density.
    cube2_value : str
        Single value of electron density.
    
    Returns
    -------
    str
        cube1 - cube2 value.
    
    Examples
    --------
    >>> calc_density_diff('3.73097E-15', '2.43683E-15')
    '1.29414E-15'
    """
    
    diff = float(cube1_value) - float(cube2_value)

    diff = '{:0.5E}'.format(diff)
    return str(diff)


def cube_info(cube_path):
    """Parses info from ORCA generated cube file.

    Parameters
    ----------
    cube_path : str
        Path to cube file.
    
    Returns
    -------
    dict
        Cube info dictionary with keys:

        ``"atom_num"``
            Number of atoms in the system (`int`).
        ``"atoms"``
            Lines with atom information (`list` [`str`]).
        ``"dim1"``
            Number of grid points in x direction.
        ``"dim2"``
            Number of grid points in y direction.
        ``"dim3"``
            Number of grid points in z direction.
        ``"dim1_interval"``
            Spacing of grid points in x direction in Angstroms.
        ``"dim2_interval"``
            Spacing of grid points in y direction in Angstroms.
        ``"dim3_interval"``
            Spacing of grid points in z direction in Angstroms.
        ``"lines"``
            List of parsed lines for header (`list [`str`]).
        ``"min1"``
            Minimum x boundary in Angstroms (`float`).
        ``"min2"``
            Minimum y boundary in Angstroms (`float`).
        ``"min3"``
            Minimum z boundary in Angstroms (`float`).
    """

    cube_data = {'atoms': [], 'lines': []}
    with open(cube_path, 'r') as cube:

        line_num = 1
        for cube_line in cube:
            split = cube_line.split()
            if line_num < 3:
                pass
            elif line_num == 3:
                cube_data['atom_num'] = int(split[0])
                cube_data['min1'] = float(split[1])
                cube_data['min2'] = float(split[2])
                cube_data['min3'] = float(split[3])
            elif line_num == 4:
                cube_data['dim1'] = int(split[0])
                cube_data['dim1_interval'] = float(split[1])
            elif line_num == 5:
                cube_data['dim2'] = int(split[0])
                cube_data['dim2_interval'] = float(split[2])
            elif line_num == 6:
                cube_data['dim3'] = int(split[0])
                cube_data['dim3_interval'] = float(split[3])
            elif line_num > 6 and line_num <= 6 + cube_data['atom_num']:
                cube_data['atoms'].append(cube_line)
            else:
                break
            
            cube_data['lines'].append(cube_line)

            line_num += 1
    
    return cube_data

def write_header(cube1_path, cube2_path, save_path):
    """Writes density difference cube file header.

    Parameters
    ----------
    cube1_path : `str`
        Path to the positive cube file.
    cube2_path : `str`
        Path to the negative cube file.
    save_path : `str`
        Path to save file.
    
    Raises
    ------
    ValueError
        Raises if ``"min1"``, ``"min2"``, ``"min3"``, ``"dim1"``, ``"dim2"``,
        ``"dim3"``, ``"dim1_interval"``, ``"dim2_interval"``, or
        ``"dim3_interval"`` of the two cube files do not match.
    """
    
    cube1_name = cube1_path.split('/')[-1]
    cube2_name = cube2_path.split('/')[-1]
    
    header_lines = [
        'Electron density difference\n', f'({cube1_name}) - ({cube2_name})\n'
    ]

    cube1_data = cube_info(cube1_path)
    cube2_data = cube_info(cube2_path)

    # Checks if the the electron density specifications match.
    for key in cube1_data.keys():

        if key not in ['atoms', 'lines', 'atom_num']:
            if cube1_data[key] != cube2_data[key]:
                raise ValueError(f'{key} is not the same.')
    
    # TODO is there a case where we do a small fragment - large fragment?
    # If there is not, do a test here and raise ValueError.
    
    # Chooses header information from the cube file with the most atoms.
    # If doing electron density differences with fragments, one will have 
    # different numbers and should choose the largest system.
    if cube1_data['atom_num'] >= cube1_data['atom_num']:
        header_lines.extend(cube1_data['lines'][2:])
    else:
        header_lines.extend(cube2_data['lines'][2:])

    # Write the header.
    with open(save_path, 'w') as diff_file:
        diff_file.writelines(header_lines)


def skip_header(cube_file):
    """Moves cube file to first density line.

    Parameters
    ----------
    cube_file : `io.textiowrapper`
        Opened file variable.
    
    Returns
    -------
    `io.textiowrapper`
        File object with header info skipped and ready to read first density
        line.
    """

    # Skips two comment lines
    next(cube_file)
    next(cube_file)

    # Skips box and atom information
    atom_num, _, _, _ = cube_file.readline().split()
    for i in range(0, 3 + int(atom_num)):
        next(cube_file)

    return cube_file


def write_diff(cube1_path, cube2_path, save_path):
    """Compure and write electron density difference.

    Cube1 - cube2.

    Parameters
    ----------
    cube1_path : `str`
        Path to positive cube file.
    cube2_path : `str`
        path to negative cube file.
    save_path : `str`
        Path to save file.
    """

    
    # Prepares to write density data
    with open(save_path, 'a') as diff_file:
        with open(cube1_path, 'r') as cube1:
            cube1 = skip_header(cube1)
            with open(cube2_path, 'r') as cube2:
                cube2 = skip_header(cube2)

                # TODO Estimate percent done and print.
                # Could estimate number of lines by using dim numbers and 
                # how many numbers per line.

                # Writes density difference line by line
                for cube1_line, cube2_line, in zip(cube1, cube2):
                    cube1_line_split = cube1_line.split()
                    cube2_line_split = cube2_line.split()

                    item_index = 0
                    new_values = []
                    end_of_line = False

                    while item_index < len(cube1_line_split):
                        cube1_line_str = cube1_line_split[item_index]
                        cube2_line_str = cube2_line_split[item_index]


                        # Adds spaces
                        if cube1_line_str == '':
                            new_values.append('')
                        # Skips blank lines.
                        elif cube1_line_str.isspace():
                            pass
                        else:
                            
                            if item_index == (len(cube1_line_split) - 1):
                                end_of_line = True

                            density_diff_str = calc_density_diff(
                                cube1_line_str,
                                cube2_line_str
                            )

                            if end_of_line:
                                density_diff_str += '\n'
                            
                            new_values.append(density_diff_str)
                        
                        item_index += 1

                    # Formatting line
                    # TODO Fix formatting
                    diff_line = '   '.join(new_values)
                    diff_line = diff_line.replace(' -', '-')
                    
                    diff_file.write(diff_line)
                    
                

def main():
    
    parser = argparse.ArgumentParser(
        description='Creates density difference cube file (cube1 - cube2).'
    )
    parser.add_argument(
        'cube1', nargs=1, type=str,
        help='path to first cube file'
    )
    parser.add_argument(
        'cube2', nargs=1, type=str,
        help='path to second cube file'
    )
    parser.add_argument(
        '--save', metavar='save_name', nargs=1,
        default='electron-density-difference',
        help='name to save difference cube file, defaults to '
             '"electron-density-difference"'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', nargs=1, default='./',
        help='directory to save difference cube file'
    )
    args = parser.parse_args()

    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'

    save_path = f'{save_dir}{args.save}.cube'

    print('Writing header...')
    write_header(args.cube1[0], args.cube2[0], save_path)
    print('Writing density difference...')
    write_diff(args.cube1[0], args.cube2[0], save_path)

if __name__ == "__main__":
    main()