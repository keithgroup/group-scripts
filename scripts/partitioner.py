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
# Partitioner

A partition, in this contex, is one or more molecules from a larger atomic
cluster. For example, a dimer partition of a tetramer would be one of the
possible combinations containing two molecules out of the possible four
molecules.

This script will take a single structure or trajectory, identify every possible
partition up to a specified size, then write xyz coordinates for each partition.
Also, a set of energy+gradient calculations using ORCA 4 can be prepared; these
calculations are pertinent for mbGDML data sets.

If multiple structures are desired, the atoms must be in the same order.

## Requirements
- mbgdml
"""

import os
import argparse
from mbgdml.utils import write_xyz
from mbgdml.partition import partition_structures
from mbgdml.data import mbGDMLDataset, structure
from mbgdml.calculate import partition_engrad


### Universal Information ###

calc_theory = 'MP2'
calc_basis = 'def2-TZVP'
calc_options = 'TightSCF'

def main():

    parser = argparse.ArgumentParser(
        description='Partitions structure(s) for energy+gradient calculations.'
    )
    parser.add_argument(
        'xyz_file', metavar='xyz_file', type=str, nargs='?',
        help='Path to xyz file with one or multiple structures to partition.'
    )
    parser.add_argument(
        '--dir_name', metavar='dir_name', type=str, nargs='?', default='partitions',
        help='Name of the top-level directory for all partitions.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save npz data set. Defaults to current directory.'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite npz data set.'
    )
    parser.add_argument(
        '-s', '--size', metavar='size', type=int, nargs='?', default=4,
        help='Maximum size of partitions.'
    )
    parser.add_argument(
        '--calcs', action='store_true',
        help='Write MP2/def2-TZVP ORCA 4 energy+gradient calculations.'
    )
    parser.add_argument(
        '--num_calcs', metavar='num_calcs', type=int, nargs='?', default=100,
        help='Number of structures to include in the energy+gradient calculations.'
    )
    parser.add_argument(
        '--cluster_name', metavar='cluster_name', type=str, nargs='?', default='cluster',
        help='Base name describing the original cluster.'
    )
    parser.add_argument(
        '--r_unit', metavar='r_unit', type=str, nargs='?', default='Angstrom',
        help='Units of distance for the structure(s): Angstrom or bohr.'
    )



    args = parser.parse_args()

    print(f'\nPartitioner v{__version__}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)\n')

    # Ensures paths end in a '/'.
    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    
    # Checks to see if partitions already exists.
    partitions_dir_path = f'{save_dir}{args.dir_name}/'
    if os.path.isdir(partitions_dir_path) and not args.overwrite:
        print(f'{partitions_dir_path} already exists and overwrite is False.\n')
        raise FileExistsError
    
    # Prepare to partition.
    print('Preparing partition directories ...')
    os.makedirs(partitions_dir_path, exist_ok=True)
    os.chdir(partitions_dir_path)

    # Parse file into mbGDML data set.
    print('Parsing XYZ file ... ')
    cluster_dataset = mbGDMLDataset()
    cluster_dataset.read_xyz(
        args.xyz_file, 'coords', r_unit=args.r_unit, energy_comments=False
    )
    print(f'Found {cluster_dataset.R.shape[0]} structure(s) '
          f'containing {cluster_dataset.z.shape[0]} atoms')
    
    # Partition all structures in the data set.
    print('Partitioning structure(s) ... ')
    cluster_partitions = partition_structures(
        cluster_dataset.z, cluster_dataset.R, args.size
    )
    print(f'Identified {len(cluster_partitions)} partitions')
    partition_sizes = []
    for size in range(1, args.size + 1):
        partition_sizes.append(
            len([i for i in cluster_partitions if len(i.split(',')) == size])
        )
        print(f'  {partition_sizes[-1]} {size}mer partitions')
    
    # Write all partitions and possibly energy+gradient calculations
    for partition in cluster_partitions:
        partition_size = len(partition.split(','))
        partition_name = f'{args.cluster_name}.mol{partition}'

        print(f'Working on {partition_name} partition')
        
        z = cluster_partitions[partition]['z']
        R = cluster_partitions[partition]['R']

        # XYZ file
        xyz_dir = f'{save_dir}{args.cluster_name}-{partition_size}.partitions/{partition_name}/'
        os.makedirs(xyz_dir, exist_ok=True)
        write_xyz(z, R, xyz_dir, partition_name)
        

        # Energy+Gradient calculation
        if args.calcs:
            # Adds '.first<num_calcs>' if the number of calcs is not the all
            # of the structures.
            if R.shape != R[:args.num_calcs].shape:
                partition_name += f'.first{args.num_calcs}'
            
            partition_dir = f'{xyz_dir}{partition_name}-calcs/'
            os.makedirs(partition_dir, exist_ok=True)
            calculation_name = f'{partition_name}-orca.engrad-mp2.def2tzvp.frozencore'

            partition_engrad(
                'orca',
                z,
                R[:args.num_calcs],
                partition_name,
                calculation_name,
                calculation_name,
                'MP2',
                'def2-TZVP',
                0,
                1,
                'smp',
                1,
                6,
                2,
                00,
                options='TightSCF',
                control_blocks=(
                    '%maxcore 8000\n\n'
                    '%scf\n    ConvForced true\nend'),
                calc_dir=partition_dir
            )

if __name__ == "__main__":
    main()