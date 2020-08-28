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

def calc_density_diff(cube1_str, cube2_str):
    
    cube1_value = float(cube1_str)
    cube2_value = float(cube2_str)
    diff = cube1_value - cube2_value

    diff = '{:0.5E}'.format(diff)
    return str(diff)

def write_diff(cube1_path, cube2_path, save_path):

    with open(save_path, 'w+') as diff_file:
        with open(cube1_path, 'r') as cube1:
            with open(cube2_path, 'r') as cube2:
                line_num = 1
                atom_num = 1

                for cube1_line, cube2_line, in zip(cube1, cube2):
                    cube1_line_split = cube1_line.split('   ')
                    cube2_line_split = cube2_line.split('   ')

                    # Gets actual atom number
                    if line_num == 3:
                        atom_num = int(cube1_line_split[1])

                    # Writes header
                    if line_num <= 6 + atom_num:
                        diff_file.write(cube1_line)
                    # Skip blank lines
                    elif len(cube1_line_split) == 1:
                        pass
                    # Writes volumetric data
                    else:
                        item_index = 0
                        new_values = []
                        end_of_line = False
                        while item_index < len(cube1_line_split):
                            cube1_line_str = cube1_line_split[item_index]
                            cube2_line_str = cube2_line_split[item_index]



                            if cube1_line_str == '':
                                new_values.append('')
                            else:
                                
                                if item_index == (len(cube1_line_split) - 1):
                                    end_of_line = True

                                    cube1_line_str = cube1_line_str[:-2]
                                    cube2_line_str = cube2_line_str[:-2]

                                density_diff_str = calc_density_diff(
                                    cube1_line_str,
                                    cube2_line_str
                                )

                                if end_of_line:
                                    density_diff_str += '\n'
                                
                                new_values.append(density_diff_str)
                            
                            item_index += 1

                        # Formatting line
                        diff_line = '   '.join(new_values)
                        diff_line = diff_line.replace(' -', '-')
                        
                        diff_file.write(diff_line)

                    line_num += 1
                    
                    #print(cube1_line)

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

    write_diff(args.cube1[0], args.cube2[0], save_path)

if __name__ == "__main__":
    main()