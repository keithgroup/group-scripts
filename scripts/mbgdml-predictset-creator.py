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
# mbGDML Predict Set Creator

Creates a mbGDML predict set.

## Requirements
- numpy
- mbgdml
"""

import os
import argparse
import json
import numpy as np
from mbgdml.data import mbGDMLPredictset



def main():

    parser = argparse.ArgumentParser(
        description='Creates mbGDML predict sets from data sets and models.'
    )
    parser.add_argument(
        'dataset', metavar='datasets', type=str, nargs='?',
        help='Path to mbGDML data set to predict energies and forces.'
    )
    parser.add_argument(
        '-m', '--models', metavar='models', nargs='+', default=[],
        help='Path to mbGDML data set to predict energies and forces.'
    )
    parser.add_argument(
        '--name', metavar='name', type=str, nargs='?', default='predictset',
        help='File name for the npz data set. Defaults to dataset.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save npz data set. Defaults to current directory.'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite npz data set.'
    )

    args = parser.parse_args()

    print('Making mbGDML predict sets')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)\n')

    
    # Ensures paths end in a '/'.
    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    
    # Checks to see if data set already exists.
    if os.path.isfile(f'{save_dir}{args.name}.npz') and not args.overwrite:
        print(f'{save_dir}{args.name}.npz already exists and overwrite is False.\n')
        raise FileExistsError

    # Creating the predict set.
    predictset = mbGDMLPredictset()

    print(f'Loading the {args.dataset} data set ... ', end='')
    predictset.load_dataset(args.dataset)
    print('done')

    print(f'Loading the {len(args.models)} models ... ', end='')
    predictset.load_models(args.models)
    print('done')

    print(f'Creating the predict set ...')
    predictset.save(args.name, predictset.predictset, save_dir)
    

    print(f'\nYour predict set is: {save_dir}{args.name}.npz')

if __name__ == "__main__":
    main()