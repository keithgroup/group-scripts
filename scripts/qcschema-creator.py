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

__version__ = "0.0.1"

"""
# QCSchema Creator

Create QCSchema files (https://github.com/MolSSI/QCSchema) from output files.
This script uses cclib (https://github.com/cclib/cclib) to parse
data and implements custom parsers to retrive other information.

Unofficial modifications are made to QCSchema to meet the immediate needs of
our research. These modications are 
- Official QCSchemas do not support multiple configurations in a single JSON
file: one file per structure and calculation. This is inherently incompatible
with geometry optimizations, trajectories, or any other iterative procedure.
At of the time of writing (2021-01-03), there has been no concensus of
how to manage these files. Our immediate solution is to list each QCSchema
inside a list (i.e., [{},{},{}, ...]). Furthermore, only iterations that have
the same topology are supported (no checks are provided).

## Requirements
- cclib (>=1.7)
- numpy

## Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Custom parser for ORCA information such as integration grid, scf energy
  contributions (e.g., one-electron and two-electron energies).
- Debug option to raise errors instead of skipping over files.

### Changed
- Nest iterations into a list instead of having int labels.
- Standardized getting SCF, MP, and CC energies from cclib.
- Requires outfile path to initialize schema classes.

## [0.0.1] - 2021-01-03
### Added
- Support for output files containing multiple jobs or iterations
  (optimizations).
- ORCA Schema creator for single-point energies and gradients.

"""

import os
from functools import reduce
import argparse
from packaging import version
import cclib
import numpy as np
import json

# pylint: disable=no-member

### Hard coded information ###

_element_to_z = {
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
_z_to_element = {v: k for k, v in _element_to_z.items()}

# Program-specific keywords categorized into scf, moller-plesset, and coupled
# cluster methods as specified by QCSchema. All methods are in lowercase to 
# assist identification.
methods = {
    'orca': {
        'scf': [
            'rhf', 'rks', 'uhf', 'uks', 'rohf', 'roks', 'hf',
            'hfs', 'lda', 'lsd', 'vwn', 'vwn5', 'vwn3', 'pwlda', 'bp86', 'bp',
            'blyp', 'olyp', 'glyp', 'xlyp', 'pw91', 'mpwpw', 'mpwlyp', 'mpwlyp',
            'pbe', 'rpbe', 'revpbe', 'pwp', 'b1lyp', 'b3lyp', 'b3lyp/g',
            'o3lyp', 'x3lyp', 'b1p', 'b3p', 'b3pw', 'pw1pw', 'mpw1pw',
            'mpw1lyp', 'pbe0', 'pw6b95', 'bhandhylpl', 'tpss', 'tpssh', 'tpss0',
            'm06l', 'mo6', 'm062x', 'b97m-v', 'b97m-d3bj', 'scanfunc',
            'wb97', 'wb97x', 'wb97x-d3', 'wb97x-v', 'wb97x-d3bj', 'cam-b3lyp',
            'lc-blyp', 'b2plyp', 'b2plyp-d', 'b2plyp-d3', 'mpw2plpy',
            'mpw2plyp-d', 'b2gp-plyp', 'b2k-plyp', 'b2t-plyp', 'pwpb95',
            'dsd-blyp', 'dsd-pbep86', 'dsd-pbep95', 'wb2plyp', 'wb2gp-plyp'
        ],
        'moller-plesset': [
            'mp2', 'ri-mp2', 'scs-mp2', 'ri-scs-mp2', 'oo-ri-mp2',
            'oo-ri-scs-mp2', 'mp2-f12', 'mp2-f12', 'mp2-f12-ri', 'mp2-f12d-ri',
            'mp3', 'scs-mp3', 'dlpno-mp2', 'dlpno-scs-mp2', 'dlpno-mp2-f12',
            'dlpno-mp2-f12/d'
        ],
        'coupled cluster': [
            'ccsd', 'ccsd(t)', 'ccsd-f12', 'ccsd(t)-f12', 'ccsd-f12/ri',
            'ccsd(t)-f12d/ri', 'lpno-ccsd', 'dlpno-ccsd', 'dlpno-ccsd(t)',
            'dlpno-ccsd(t1)', 'dlpno-ccsd-f12', 'dlpno-ccsd-f12/d'
        ]
    }
}

# Program-specific, basis-set keywords in lowercase (to help with
# identification).
basis_sets = {
    'orca': [
        '3-21g', 'sto-3g', '3-21gsp', '4-22gsp', '6-31g', 'm6-31g', '6-311g',
        '6-31g*', '6-31g(d)', '6-31g**', '6-31g(d,p)', '6-31g(2d)',
        '6-31g(2df)', '6-31g(2d,p)', '6-31g(2d,2p)', '6-31g(2df,2dp)',
        '6-311g*', '6-311g(d)', '6-311g**', '6-311g(d,p)', '6-311g(2d)',
        '6-311g(2df)', '6-311g(2d,p)', '6-311g(2d,2p)', '6-311g(2df,2dp)',
        '6-311g(3df)','6-311g(3df,3pd)', '6-31+g', '6-31++g(d,p)',
        'def2-svp', 'def2-sv(p)', 'def2-tzvp', 'def2-tzvp(-f)', 'def2-tzvpp',
        'def2-qzvpp', 'sv', 'sv(p)', 'svp', 'tzv', 'tzv(p)', 'tzvp', 'tzvpp',
        'qzvp', 'qzvpp', 'ma-def2-svp', 'ma-def2-sv(p)', 'ma-def2-tavp',
        'ma-def2-tzvp(-f)', 'ma-def2-tzvpp', 'ma-def2-qzvpp', 'def2-svpd',
        'def2-tzvpd', 'def2-tzvppd', 'def2-qzvpd', 'def2-qzvppd',
        'dkh-def2-svp', 'zora-def2-svp', 'dkh-def2-sv(p)', 'zora-def2-sv(p)',
        'dkh-def2-tzvp', 'zora-def2-tzvp', 'dhk-def2-tzvp(-f)',
        'zora-def2-tzvp(-f)', 'dhk-def2-tzvpp', 'zora-def2-tzvpp',
        'dkh-def2-qzvpp', 'zora-def2-qzvpp', 'ma-dkh-def2-svp',
        'ma-zora-def2-svp', 'ma-dkh-def2-sv(p)', 'ma-zora-def2-sv(p)',
        'ma-dkh-def2-tzvp', 'ma-zora-def2-tzvp', 'ma-dhk-def2-tzvp(-f)',
        'ma-zora-def2-tzvp(-f)', 'ma-dhk-def2-tzvpp', 'ma-zora-def2-tzvpp',
        'ma-dkh-def2-qzvpp', 'ma-zora-def2-qzvpp', 'dkh-svp', 'zora-svp',
        'dkh-sv(p)', 'zora-sv(p)', 'dkh-tzvp', 'zora-tzvp', 'dhk-tzvp(-f)',
        'zora-tzvp(-f)', 'dhk-tzvpp', 'zora-tzvpp', 'dkh-qzvpp', 'zora-qzvpp',
        'ma-dkh-svp', 'ma-zora-svp', 'ma-dkh-sv(p)', 'ma-zora-sv(p)',
        'ma-dkh-tzvp', 'ma-zora-tzvp', 'ma-dhk-tzvp(-f)', 'ma-zora-tzvp(-f)',
        'ma-dhk-tzvpp', 'ma-zora-tzvpp', 'ma-dkh-qzvpp', 'ma-zora-qzvpp',
        'sarc-dkh-tzvp', 'sarc-dkh-tzvpp', 'sarc-zora-tzvp', 'sarc-zora-tzvpp',
        'sarc-dkh-svp', 'sarc-zora-svp', 'sarc2-dkh-qzv', 'sarc2-dkh-qzvp',
        'sarc2-zora-qzvp', 'sarc-zora-qzvp', 'pc-0' ,'pc-1', 'pc-2', 'pc-3',
        'pc-4', 'aug-pc-0' ,'aug-pc-1', 'aug-pc-2', 'aug-pc-3', 'aug-pc-4',
        'pcseg-0' ,'pcseg-1', 'pcseg-2', 'pcseg-3', 'pcseg-4', 'aug-pcseg-0',
        'aug-pcseg-1', 'aug-pcseg-2', 'aug-pcseg-3', 'aug-pcseg-4', 'pcsseg-0',
        'pcsseg-1', 'pcsseg-2', 'pcsseg-3', 'pcsseg-4', 'aug-pcsseg-0',
        'aug-pcsseg-1', 'aug-pcsseg-2', 'aug-pcsseg-3', 'aug-pcsseg-4', 'pcj-0',
        'pcj-1', 'pcj-2', 'pcj-3', 'pcj-4', 'aug-pcj-0', 'aug-pcj-1',
        'aug-pcj-2', 'aug-pcj-3', 'aug-pcj-4', 'sapporo-dzp-2012',
        'sapporo-tzp-2012', 'sapporo-qzp-2012', 'sapporo-dkh3-dzp-2012',
        'sapporo-dkh3-tzp-2012', 'sapporo-dkh3-qzp-2012', 'cc-pvdz', 'cc-pvtz',
        'cc-pvqz', 'cc-pv5z', 'cc-pv6z', 'aug-cc-pvdz', 'aug-cc-pvtz',
        'aug-cc-pvqz', 'aug-cc-pv5z', 'aug-cc-pv6z', 'cc-pcvdz', 'cc-pcvtz',
        'cc-pcvqz' 'cc-pcv5z', 'cc-pcv6z', 'aug-cc-pcvdz', 'aug-cc-pcvtz',
        'aug-cc-pcvqz' 'aug-cc-pcv5z', 'aug-cc-pcv6z', 'cc-pwcvdz', 'cc-pwcvtz',
        'cc-pwcvqz' 'cc-pwcv5z', 'aug-cc-pwcvdz', 'aug-cc-pwcvtz',
        'aug-cc-pwcvqz' 'aug-cc-pwcv5z', 'aug-cc-pwcvd(+d)z',
        'aug-cc-pwcvt(+d)z', 'aug-cc-pwcvq(+d)z' 'aug-cc-pwcv5(+d)z',
        'ano-pvdz', 'ano-pvtz', 'ano-pvqz', 'ano-pv5z', 'saug-ano-pvdz',
        'saug-ano-pvtz', 'saug-ano-pvqz', 'aug-ano-pvdz', 'aug-ano-pvtz',
        'aug-ano-pvqz', 'ano-rcc-full', 'ano-rcc-dzp', 'ano-rcc-tzp',
        'ano-rcc-qzp', 'd95', 'd95p', 'mini', 'minis', 'midi', 'minix',
        'wachters+f', 'partridge-1', 'partridge-2', 'partridge-3',
        'partridge-4', 'lanl2dz', 'lanl2tz', 'lanl2tz(f)', 'lanl08', 'lanl08(f)',
        'epr-ii', 'epr-iii', 'iglo-ii', 'iglo-iii', 'aug-cc-pvtz-j'
    ]
}










### Utility functions ###

def atoms_by_element(atom_list):
    """Converts a list of atoms identified by their atomic number to their
    elemental symbol in the same order.
    
    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.
    """
    return [_z_to_element[i] for i in atom_list]

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

def convert_forces(
    forces, e_units_calc, r_units_calc, e_units, r_units
):
    """Converts forces (or gradients) to specified units.

    Parameters
    ----------
    forces : :obj:`numpy.ndarray`
        An array with units of energy and distance matching `e_units_calc`
        and `r_units_calc`.
    e_units_calc : :obj:`str`
        Specifies package-specific energy units used in calculation. Available
        units are ``'eV'``, ``'hartree'``, ``'kcal/mol'``, and ``'kJ/mol'``.
    r_units_calc : :obj:`str`
        Specifies package-specific distance units used in calculation. Available
        units are ``'Angstrom'`` and ``'bohr'``.
    e_units : :obj:`str`
        Desired units of energy. Available units are ``'eV'``, ``'hartree'``,
        ``'kcal/mol'``, and ``'kJ/mol'``.
    r_units : obj:`str`
        Desired units of distance. Available units are ``'Angstrom'`` and
        ``'bohr'``.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Forces converted into the desired units.
    """
    #'ORCA': {'e_unit': 'hartree', 'r_unit': 'bohr'}
    if e_units not in ['eV', 'hartree', 'kcal/mol', 'kJ/mol']:
        raise ValueError(f'{e_units} is not an available energy unit.')
    if r_units not in ['Angstrom', 'bohr']:
        raise ValueError(f'{r_units} is not an available distance unit.')
    forces_conv = forces
    if e_units_calc != e_units:
        forces_conv = cclib.parser.utils.convertor(
            forces_conv, e_units_calc, e_units
        )
    if r_units_calc != r_units:
        forces_conv = cclib.parser.utils.convertor(
            forces_conv, r_units, r_units_calc
        )
    return forces_conv

def identify_package(outfile_path):
    """Identifies computational chemistry package.

    Only supported packaged should be included in `triggers`.

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    
    Returns
    -------
    :obj:`chemicalJSON`
        One of the supported ChemicalJSON classes.
    """
    with open(outfile_path, 'r') as f:
        for line in f:
            for parser, phrases, do_break in triggers:
                if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                    filetype = parser
                    if do_break:
                        return filetype

def read_json(json_path):
    """Read JSON file.
    
    Parameters
    ----------
    json_path : :obj:`str`
        Path to json file.
    
    Returns
    -------
    :obj:`dict`
        Contents of JSON file.
    """
    with open(json_path, 'r') as reader:
        json_dict = json.loads(reader.readline())
    
    return json_dict










### Manual parsing classes ###

class outfileParser:
    """Base class for parsing output files.

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    """

    def __init__(self, outfile_path):
        self.outfile_path = outfile_path
        self.file_name = '.'.join(outfile_path.split('/')[-1].split('.')[:-1])





class orcaParser(outfileParser):

    def __init__(self, outfile_path):
        super().__init__(outfile_path)
    
    def get_grid_info(self):
        """DFT integration grid information.

        All DFT calculations are performed with numerical integration, which
        means a integration grid must be specified. This is typically handeled
        with the Grid keyword. Furthermore, ORCA defaults to a multigrid
        approach where one grid is used for the SCF cycle, and a different
        (typically larger) grid is used for the final energy evaluation.

        After testing different keyword combinations, the Grid keywords are not
        always consistent with the final grid. Thus, we are going to directly
        parse the grid information from the output file instead of depending
        on the keywords.

        Returns
        -------
        :obj:`dict`
            Contains grid information using the following keys:

            ``'scf_grid_level'``
                :obj:`int` specifying ORCA default grid level (1 to 7) during
                the SCF cycle.
            ``'final_grid_level'``
                :obj:`int` specifying ORCA default grid level for the final
                evaluation.
        
        Notes
        -----
        Here are the main parameters specifying the ORCA default grid levels.

        +--------+---------------+----------+
        |  Grid  |  AngularGrid  |  IntAcc  |
        +========+===============+==========+
        |  1     |  Lebedev-50   |  4.67    |
        |  2     |  Lebedev-110  |  4.34    |
        |  3     |  Lebedev-194  |  4.34    |
        |  4     |  Lebedev-302  |  4.67    |
        |  5     |  Lebedev-434  |  5.01    |
        |  6     |  Lebedev-590  |  5.34    |
        |  7     |  Lebedev-770  |  5.67    |
        +--------+---------------+----------+
        """
        grid_info = {}
        lebedev_to_level = {
            '50': 1, '110': 2, '194': 3, '302': 4, '434': 5, '590': 6, '770': 7
        }

        with open(self.outfile_path, mode='r') as outfile:
            for line in outfile:
                # -------------------
                # DFT GRID GENERATION
                # -------------------
                # 
                # General Integration Accuracy     IntAcc      ...  4.340
                # Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
                # Angular Grid (max. acc.)         AngularGrid ... Lebedev-110
                if 'DFT GRID GENERATION' == line.strip():
                    while 'Angular Grid (max. acc.)' not in line:
                        line = next(outfile)
                    lebedev_num = line.strip().split('-')[-1]
                    grid_info['scf_grid_level'] = lebedev_to_level[lebedev_num]
                    continue
                
                if 'Setting up the final grid:' == line.strip():
                    while 'Angular Grid (max. acc.)' not in line:
                        line = next(outfile)
                    lebedev_num = line.strip().split('-')[-1]
                    grid_info['final_grid_level'] = lebedev_to_level[lebedev_num]
                    continue
        
        if 'final_grid_level' not in grid_info.keys():
            grid_info['final_grid_level'] = grid_info['scf_grid_level']
        return grid_info
    
    def _parse_other_scf(self, outfile, scf_info):
        """The nulear repulsion, one- and two-electron energy, and
        exchange-correlation energy after a SCF cycle.

        This is called directly after the ``'TOTAL SCF ENERGY'`` trigger, and 
        will terminate once the ``'SCF CONVERGENCE'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``scf_info`` here so any missing information (such as
        ``'scf_xc_energy'`` in MP2 calculations) is not an issue.

        Parameters
        ----------
        outfile
            An opened output file.
        scf_info : :obj:`dict`
            The scf_info dict that contains all QCSchema-supported SCF energies.
        
        Returns
        -------
        :obj:`dict`
            Available SCF energy components that could include the following
            keys.

            ``'scf_one_electron_energy'``
                The one-electron (core Hamiltonian) energy contribution to the
                total SCF energy.
            ``'scf_two_electron_energy'``
                The two-electron energy contribution to the total SCF energy.
            ``'nuclear_repulsion_energy'``
                The nuclear repulsion energy contribution to the total SCF
                energy.
            ``'scf_xc_energy'``
                The functional energy contribution to the total SCF energy.
        """
        for line in outfile:
            # Nuclear Repulsion  :  135.87324654 Eh    3697.29901 eV
            if 'Nuclear Repulsion' in line:
                if 'nuclear_repulsion_energy' not in scf_info.keys():
                    scf_info['nuclear_repulsion_energy'] = []
                scf_info['nuclear_repulsion_energy'].append(
                    float(line.split()[3])
                )

            # One Electron Energy: -674.26034691 Eh  -18347.55681 eV
            if 'One Electron Energy' in line:
                if 'scf_one_electron_energy' not in scf_info.keys():
                    scf_info['scf_one_electron_energy'] = []
                scf_info['scf_one_electron_energy'].append(
                    float(line.split()[3])
                )

            # Two Electron Energy:  245.90403408 Eh    6691.38895 eV
            if 'Two Electron Energy' in line:
                if 'scf_two_electron_energy' not in scf_info.keys():
                    scf_info['scf_two_electron_energy'] = []
                scf_info['scf_two_electron_energy'].append(
                    float(line.split()[3])
                )

            # E(XC)   :    -26.170406641397 Eh
            if 'E(XC)' in line:
                if 'scf_xc_energy' not in scf_info.keys():
                    scf_info['scf_xc_energy'] = []
                scf_info['scf_xc_energy'].append(
                    float(line.split()[2])
                )

            if 'SCF CONVERGENCE' == line.strip():
                break
        
        return scf_info
    
    def _parse_other_mp(self, outfile, mp_info):
        """Moller-Plesset calculation properties.

        This is called directly after the ``'ORCA  MP2 '`` trigger, and 
        will terminate once the ``'ORCA property calculations'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``mp_info`` here so any missing information is not an
        issue.

        Parameters
        ----------
        outfile
            An opened output file.
        mp_info : :obj:`dict`
            The mp_info dict that contains all QCSchema-supported mp energies.
        
        Returns
        -------
        :obj:`dict`
            Available MP energy components that could include the following
            keys.

            ``'mp2_correlation_energy'``
                The MP2 correlation energy.
        """
        for line in outfile:
            #  MP2 CORRELATION ENERGY   :     -3.132364939 Eh
            if 'MP2 CORRELATION ENERGY' in line:
                if 'mp2_correlation_energy' not in mp_info.keys():
                    mp_info['mp2_correlation_energy'] = []
                mp_info['mp2_correlation_energy'].append(
                    float(line.split()[4])
                )

            #      ***************************************
            #      *     ORCA property calculations      *
            #      ***************************************
            if '*     ORCA property calculations      *' == line.strip():
                break
        
        return mp_info
    
    def get_scf_info(self, iteration=-1):
        """Other scf energy components.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.

        Returns
        -------
        :obj:`dict`
            Available SCF energy components that could include the following
            keys.

            ``'scf_one_electron_energy'``
                The one-electron (core Hamiltonian) energy contribution to the
                total SCF energy.
            ``'scf_two_electron_energy'``
                The two-electron energy contribution to the total SCF energy.
            ``'nuclear_repulsion_energy'``
                The nuclear repulsion energy contribution to the total SCF
                energy.
            ``'scf_xc_energy'``
                The functional energy contribution to the total SCF energy.
        """
        
        if not hasattr(self, 'scf_info'):
            scf_info = {}
            with open(self.outfile_path, mode='r') as outfile:
                for line in outfile:

                    # ----------------
                    # TOTAL SCF ENERGY
                    # ----------------
                    if 'TOTAL SCF ENERGY' == line.strip():
                        # Iterative jobs will have multiple SCF calculations.
                        # Thus, we separated the parsing trigger and the actual
                        # parsing of the SCF energies. That way, the entire
                        # file is parsed instead of stopping at the first one.
                        scf_info = self._parse_other_scf(outfile, scf_info)

            self._scf_info = scf_info
        
        scf_info_iter = {}
        for k, v in self._scf_info.items():
            scf_info_iter[k] = v[iteration]
        return scf_info_iter
    
    def get_mp_info(self, iteration=-1):
        """Other MP energy components.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.

        Returns
        -------
        :obj:`dict`
            Available MP energy components that could include the following
            keys.

            ``'scf_one_electron_energy'``
                The one-electron (core Hamiltonian) energy contribution to the
                total SCF energy.
        """
        
        if not hasattr(self, 'mp_info'):
            mp_info = {}
            with open(self.outfile_path, mode='r') as outfile:
                for line in outfile:

                    # ----------------------------------------------------------
                    #                         ORCA  MP2 
                    # ----------------------------------------------------------
                    if 'ORCA  MP2' == line.strip():
                        # Iterative jobs will have multiple MP2 calculations.
                        # Thus, we separated the parsing trigger and the actual
                        # parsing of the MP2 energies. That way, the entire
                        # file is parsed instead of stopping at the first one.
                        mp_info = self._parse_other_mp(outfile, mp_info)

            self._mp_info = mp_info
        
        mp_info_iter = {}
        for k, v in self._mp_info.items():
            mp_info_iter[k] = v[iteration]
        return mp_info_iter











### QCSchema classes ###

class QCJSON:
    """Base quantum chemistry schema.

    Attributes
    ----------
    outfile : :obj:'cclib.parser.data.ccData_optdone_bool'
        Parsed data from output file.
    name : :obj:`str`
        Name identifying the QCschema. Defaults to ``'schema'`` until an output
        file is parsed by cclib.
    path : :obj:`str`
        Absolute path to the output file.
    multiple_iterations : :obj:`bool`
        If the QCSchema is from a calculation with multiple iterations.
    """

    def __init__(self):
        self.name = 'schema'

    def write(self, name, schema_dict, save_dir, prettify=True):
        """Writes QCSchema json.

        Parameters
        ----------
        name : :obj:`str`
            Name of file.
        save_dir : :obj:`str`
            Path to save directory.
        prettify : :obj:`bool`
            Indents JSON objects if True. If false the JSON file is only one
            line.
        """
        if save_dir[-1] != '/':
            save_dir += '/'
        if prettify:
            json_string = json.dumps(
                schema_dict, cls=cclib.io.cjsonwriter.JSONIndentEncoder,
                sort_keys=True, indent=4
            )
        else:
            json_string = json.dumps(
                schema_dict, cls=cclib.io.cjsonwriter.NumpyAwareJSONEncoder,
                sort_keys=True
            )
        with open(f'{save_dir}{name}.json', 'w') as f:
            f.write(json_string)
    
    def parse_output(self):
        """Parse output file using cclib.

        Parameters
        ----------
        outfile_path : :obj:`str`
            Path to computational chemistry output file.
        """
        filename_with_extension = os.path.basename(self.outfile_path)
        self.path = os.path.abspath(self.outfile_path)
        self.name = '.'.join(filename_with_extension.split('.')[:-1])
        self.outfile = cclib.io.ccread(self.outfile_path)
        if self.outfile.atomcoords.shape[0] > 1:
            self.multiple = True
        elif self.outfile.atomcoords.shape[0] == 1:
            self.multiple = False
    
    @property
    def schema(self):
        """Base QCSchema information.
        """
        if not hasattr(self, 'outfile'):
            raise AttributeError('No output file was parsed.')
        schema_dict = {
            'schema_name': 'qc_schema_output',
            'schema_version': 1,
            'qcschema_creator_version': __version__
        }
        return schema_dict
    
    def _get_scf_energy(self, iteration=-1):
        """The energy after SCF cycle.

        This is the energy for all DFT methods and HF energy for post-HF
        calculations.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            SCF energy in Hartree.
        """
        scfenergy = cclib.parser.utils.convertor(
            self.outfile.scfenergies[iteration], 'eV', 'hartree'
        )
        return scfenergy
    
    def _get_dispersion_energy(self, iteration=-1):
        """Dispersion energy corrections to DFT calculations.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            Dispersion energy correction in Hartree.
        """
        dispersion_energy = cclib.parser.utils.convertor(
            self.outfile.dispersionenergies[iteration], 'eV', 'hartree'
        )
        return dispersion_energy
    
    def _get_mp_energy(self, iteration=-1):
        """Total energies after Moller-Plesset corrections.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            Total energy from highest order MP corrections.
        """
        if self.outfile.mpenergies.ndim == 1:
            mpenergy_ev = self.outfile.mpenergies[iteration]
        elif self.outfile.mpenergies.ndim == 2:
            order = self.outfile.mpenergies.shape[1] + 1
            mpenergy_ev = self.outfile.mpenergies[iteration][order - 2]
        mpenergy_hartree = cclib.parser.utils.convertor(
            mpenergy_ev, 'eV', 'hartree'
        )
        return mpenergy_hartree
    
    def _get_cc_energy(self, iteration=-1):
        """Energy after all coupled cluster corrections.

        Only the highest order correction is included.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`float`
            CC energy in Hartree.
        """
        ccenergy = cclib.parser.utils.convertor(
            self.outfile.ccenergies[iteration], 'eV', 'hartree'
        )
        return ccenergy





class orcaJSON(QCJSON):
    """ORCA specific QCSchema information.

    Supported calculation types:

    * Single-point energies
    * Energy+Gradient
    * Optimizations

    Parameters
    ----------
    outfile_path : :obj:`str`
        Path to output file.

    Attributes
    ----------
    outfile_path : :obj:`str`
        Path to output file.
    method_type : :obj:`str`
        QCSchema method type of ``'scf'`` (e.g., DFT), ``'moller-plesset'``,
        or ``'coupled cluster'``.
    parser : :obj:`orcaParser`
        Manually parse information from output file.
    """

    def __init__(self, outfile_path):
        super().__init__()
        self.outfile_path = outfile_path
        self.parser = orcaParser(outfile_path)
    
    def get_topology(self, iteration=-1):
        """A full description of the overall molecule its geometry, fragments,
        and charges.

        Returned keys will be a top-level JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'molecule'``
                Atomic information about the system. Contains the following keys

                ``'geometry'``
                    (3 * nat, ) vector of XYZ coordinates [a0] of the atoms.
                ``'symbols'``
                    (nat, ) atom symbols in title case.
            ``'molecular_charge'``
                The overall charge of the molecule.
            ``'molecular_multiplicity'``
                The overall multiplicity of the molecule.
            ``'name'``
                The name of the molecule or system.
            ``'atomic_numbers'``
                (nat, ) atomic numbers, nuclear charge for atoms. Ghostedness
                should be indicated through ‘real’ field, not zeros here.
        """
        try:
            topology = {
                'molecule': {
                    'geometry' : self.outfile.atomcoords[iteration],
                    'symbols': atoms_by_element(self.outfile.atomnos.tolist())
                },
                'molecular_charge': self.outfile.charge,
                'molecular_multiplicity': self.outfile.mult,
                'name': self.name,
                'atomic_numbers': self.outfile.atomnos
            }
        except:
            error_out(self.path, 'Topology data not succussfully parsed.')
        return topology

    def get_model(self, iteration=-1):
        """Model chemistry used for the calculation.

        This will be placed under the ``'model'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'method'``
                They main quantum-chemical method used for the calculation. For
                DFT calculations this would be the functional.
            ``'basis'``
                The basis set used for the calculation.
        """
        model = {}
        _remove_keywords = []
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower in methods['orca']['scf']:
                self.method_type = 'scf'
                # Ensures dispersion method is included in 'properties' for
                # functionals that automatically include it.
                model['method'] = kw
                if kw_lower == 'wb97x-d3' or kw_lower == 'b2plyp-d3':
                    self.orca_keywords.append('d3')
                if  kw_lower == 'b97m-d3bj' or kw_lower == 'wb97x-d3bj':
                    self.orca_keywords.append('d3bj')
                _remove_keywords.append(kw)
                break
            elif kw_lower in methods['orca']['moller-plesset']:
                self.method_type = 'moller-plesset'
                if 'mp3' in kw_lower:
                    raise ValueError('MP3 is not supported.')
                model['method'] = kw
                _remove_keywords.append(kw)
                break
            elif kw_lower in methods['orca']['coupled cluster']:
                self.method_type = 'coupled cluster'
                model['method'] = kw
                _remove_keywords.append(kw)
                break
        
        # Need a separate keyword iteration for basis sets because
        # dispersion-corrected functionals messes up the indicies for pop.
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower in basis_sets['orca']:
                model['basis'] = kw
                _remove_keywords.append(kw)
        
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        return model
    
    def get_keywords(self, iteration=-1):
        """Package-specific job properties.

        This will be placed under the ``'keyword'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`
            Could contain the following keys:

            ``'dispersion'``
                Method of empirical dispersion included in job.
            ``'implicit_solvent'``
                Implicit (i.e., continuum) solvent model used in job.
            ``'solvent_name'``
                Name of the solvent (e.g., ``'water'``).
            ``'frozencore'``
                If the FrozenCore approximation was used. This defaults to
                ``True`` for MP and CC calculations in ORCA.
            ``'scf_convergence_tolerance'``
                Self-consistent field convergence tolerance keyword for ORCA.
        """
        keywords = {}
        _remove_keywords = []
        
        # Parses and categorizes calculation parameters.
        for kw in self.orca_keywords:
            kw_lower = kw.lower()

            # Empirical dispersion
            if kw_lower in ['d4', 'd3bj', 'd3', 'd3zero', 'd2']:
                if version.parse(self.outfile.metadata['package_version']) \
                   >= version.parse('4.0.0') and kw_lower == 'd3':
                    # In ORCA 4, D3 is D3BJ
                    keywords['dispersion'] = 'D3BJ'
                else:
                    keywords['dispersion'] = kw
                _remove_keywords.append(kw)
            
            # Implicit solvent models
            if 'cpcm' in kw_lower or 'c-pcm' in kw_lower:
                if 'smd' in self.outfile.metadata['input_file_contents']:
                    keywords['implicit_solvent'] = 'SMD'
                else:
                    keywords['implicit_solvent'] = 'CPCM'
                if '(' in kw_lower and ')' == kw_lower[-1]:
                    solvent_name = kw_lower[:-1].split('(')[-1]
                    keywords['solvent_name'] = solvent_name
                _remove_keywords.append(kw)
            
            if kw_lower == 'frozencore':
                keywords['frozencore'] = True
                _remove_keywords.append(kw)
            elif kw_lower == 'nofrozencore':
                keywords['frozencore'] = False
                _remove_keywords.append(kw)
            
            if kw_lower in ['sloppyscf', 'loosescf', 'normalscf', 'strongscf',
                            'tightscf', 'verytightscf', 'extremescf'] \
                or 'scfconv' in kw_lower:
                keywords['scf_convergence_tolerance'] = kw
                _remove_keywords.append(kw)
            
            # Removes grid information (will be manually parsed from outfile).
            if 'grid' in kw_lower:
                _remove_keywords.append(kw)
        
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        # Add default scf convergence if not already specified.
        if 'scf_convergence_tolerance' not in keywords.keys():
            keywords['scf_convergence_tolerance'] = 'NormalSCF'
        
        # Manually parses SCF and final grid.
        keywords = {
            **keywords, **self.parser.get_grid_info()
        }
        
        # Specifies if FrozenCore is used in MP or CC jobs.
        if self.method_type == 'moller-plesset':
            if version.parse(self.outfile.metadata['package_version']) \
                >= version.parse('4.0.0'):
                lower_list = [i.lower() for i in self.orca_keywords]
                if 'nofrozencore' in lower_list:
                    keywords['frozencore'] = False
                    self.orca_keywords.pop(lower_list.index('nofrozencore'))
                else:
                    if 'frozencore' not in keywords.keys():
                        keywords['frozencore'] = True

        # Add uncategorized calculation properties
        if len(self.orca_keywords) != 0:
            keywords['other'] = self.orca_keywords

        return keywords
    
    def get_properties(self, iteration=-1):
        """A list of valid quantum chemistry properties tracked by the schema.

        This will be placed under the ``'properties'`` JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'calcinfo_nbasis'``
                The number of basis functions for the computation.
            ``'calcinfo_nmo'``
                The number of molecular orbitals for the computation.
            ``'scf_total_energy'``
                The total electronic energy of the SCF stage of the calculation.
                This is represented as the sum of the … quantities.
            ``'scf_dispersion_correction_energy'``
                The dispersion correction appended to an underlying functional
                when a DFT-D method is requested.
            ``'scf_iterations'``
                The number of SCF iterations taken before convergence.
            ``'mp2_total_energy'``
                The total MP2 energy (MP2 correlation energy + HF energy).
        """
        properties = {}
        if self.method_type == 'scf':
            properties['scf_total_energy'] = self._get_scf_energy(
                iteration=iteration
            )
            properties['scf_dispersion_correction_energy'] = cclib.parser.utils.convertor(
                self.outfile.dispersionenergies[iteration], 'eV', 'hartree'
            )
            properties['scf_iterations'] = self.outfile.scfvalues[0].shape[0]
        elif self.method_type == 'moller-plesset':
            properties['scf_total_energy'] = self._get_scf_energy(
                iteration=iteration
            )
            properties['mp2_total_energy'] = self._get_mp_energy(
                iteration=iteration
            )
            properties = {
                **properties, **self.parser.get_mp_info(iteration=iteration)
            }
        elif self.method_type == 'coupled cluster':
            # TODO
            pass
        else:
            error_out(self.path, 'Unknown method type.')
        properties['calcinfo_nbasis'] = self.outfile.nbasis
        properties['calcinfo_nmo'] = self.outfile.nmo
        properties = {
            **properties, **self.parser.get_scf_info(iteration=iteration)
        }
        return properties
    
    def get_driver(self, iteration=-1):
        """The purpose of the calculation and its direct result.

        Returned keys will be a top-level JSON property.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        
        Returns
        -------
        :obj:`dict`

            ``'driver'``
                The purpose of the calculation. Could be ``'energy'``,
                ``'gradient'``, ``'optimization'``, or ``'frequency'``.
            ``'return_results'``
                The direct result of the driver calculation. For example, a
                energy+gradient calculation will return the gradient (the
                energies will be included in the ``'properties'`` property).
        """
        driver = {}
        _remove_keywords = []

        for kw in self.orca_keywords:
            kw_lower = kw.lower()

            # Gradients
            if kw_lower == 'engrad' or kw_lower == 'numgrad':
                self.calc_driver = 'gradient'
                driver['driver'] = self.calc_driver

                if self.outfile.grads.ndim == 3:
                    grads = self.outfile.grads[iteration]
                else:
                    raise ValueError('Please check gradient dimensions.')
                driver['return_result'] = convert_forces(
                    grads, 'hartree', 'bohr', 'hartree', 'Angstrom'
                )

                _remove_keywords.append(kw)
                break

            # Frequencies (not supported)
            elif kw_lower == 'freq' or kw_lower == 'numfreq':
                self.calc_driver = 'frequency'
                driver['driver'] = self.calc_driver

                _remove_keywords.append(kw)
                break

            # Optimizations
            elif kw_lower == 'opt' or kw_lower == 'copt' or kw_lower == 'zopt':
                self.calc_driver = 'optimization'
                driver['driver'] = self.calc_driver

                _remove_keywords.append(kw)
                break
            
            # Energies
            elif kw_lower == 'energy' or kw_lower == 'sp':
                # This driver will be captured by the len(driver) == 0 condition.
                _remove_keywords.append(kw)
                break
        
        # ORCA defaults to single-point energies if no keyword is present.
        if len(driver) == 0:
            self.calc_driver = 'energy'
            driver['driver'] = self.calc_driver
            if hasattr(self.outfile, 'ccenergies'):
                return_result = self.outfile.ccenergies[iteration]
            elif hasattr(self.outfile, 'mpenergies'):
                return_result = self._get_mp_energy(
                iteration=iteration
            )
            elif hasattr(self.outfile, 'scfenergies'):
                return_result = self._get_scf_energy(
                    iteration=iteration
                )
            driver['return_result'] = return_result
            pass

        # Removes triggered keywords so as not to be included in uncategorized
        # keywords.
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        return driver

    def get_provenance(self):
        """Information about the originating package of the calculation.

        Returns
        -------
        :obj:`dict`
            
            ``'creator'``
                Package name.
            ``'version'``
                Package version.
        """
        provenance = {
            'creator': 'ORCA',
            'version': self.outfile.metadata['package_version']
        }
        return provenance
    
    def get_json(self, debug=False):
        """QC schema of an ORCA output file.

        Calculations supported: single-point energies, gradients, optimizations.

        :type: :obj:`dict`
        """
        if not hasattr(self, '_schema'):

            all_jsons = []
            self.orca_keywords = self.outfile.metadata['keywords']
            for i in range(0, self.outfile.atomcoords.shape[0]):
                # Optimizations are iterative with only a single set of
                # keywords; this is different than consecutive jobs (e.g., 
                # energies) that have repeated keywords. So, we reinitialzie the
                # keywords for each optimization iteration.
                if hasattr(self, 'calc_driver') \
                   and self.calc_driver == 'optimization':
                    self.orca_keywords = self.outfile.metadata['keywords']
                try:
                    all_jsons.append(super().schema)
                    all_jsons[-1]['provenance'] = self.get_provenance()
                    all_jsons[-1] = {
                        **all_jsons[-1], **self.get_topology(iteration=i)
                    }
                    all_jsons[-1] = {
                        **all_jsons[-1], **self.get_driver(iteration=i)
                    }
                    all_jsons[-1]['model'] = self.get_model(iteration=i)
                    all_jsons[-1]['keywords'] = self.get_keywords(iteration=i)
                    all_jsons[-1]['properties'] = self.get_properties(iteration=i)
                    all_jsons[-1]['success'] = self.outfile.metadata['success']
                    if len(all_jsons) == 1:
                        self._schema = all_jsons[0]
                    else:
                        self._schema = all_jsons
                except:
                    if debug:
                        raise
                    else:
                        if self.path not in error_files:
                            error_out(self.path, 'Uncaught exceptions.')
                        self._schema = all_jsons
                        break
        return self._schema










### Runtime Functions ###

# Triggers to identify output files.
triggers = [
    (orcaJSON, ["O   R   C   A"], True)
]

def cclib_version_check():
    """Ensures cclib version is at least 1.7.

    Dispersion energies were not systematically included or parsed until version
    1.6.4. However, the __version__ property was incorrect until version 1.7.
    """
    cclib_version = cclib.__version__
    if version.parse(cclib_version) < version.parse('1.6.1'):
        raise ValueError(
            f'cclib version is {cclib_version}; need 1.7 or higher.'
        )

error_files = []

def error_out(outfile, error_message):
    """Controls what happens when an error is encountered.

    Prints error and does not write file.
    """
    global error_files
    error_files.append(outfile)
    print(f'\u001b[31;1mError: {error_message}\u001b[0m')










def main():
    
    parser = argparse.ArgumentParser(
        description='Creates QCSchema json files from output files. Files must '
                    'have "out" somewhere in the file name or as an extension.'
    )
    parser.add_argument(
        'outputs', metavar='outputs', type=str, nargs='?', default='.',
        help='Path to directory or specific computational chemistry output file.'
    )
    parser.add_argument(
        '--save_dir', metavar='save_dir', type=str, nargs='?', default='.',
        help='Path to save JSON files.'
    )
    parser.add_argument(
        '--remove', metavar='remove_files', nargs='+', default='',
        help='Ignore file paths matching the words in this string.'
    )
    parser.add_argument(
        '--name', metavar='name', nargs=1, default='data',
        help='Name of final combined file.'
    )
    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Recursively create schemas.'
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='Overwrite JSON files'
    )
    parser.add_argument(
        '--combine', action='store_true',
        help='Combine all JSON files into a '
        'single JSON file organized by directory structure.'
    )
    parser.add_argument(
        '-p', '--prettify', action='store_true',
        help='Prettify JSON files with indentation.'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Will not continue if an error is encountered'
    )
    args = parser.parse_args()
    print(f'QCSchema creator v{__version__}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)')
    print('Parsed data are converted to Hartrees and Angstroms\n')

    cclib_version_check()

    all_qcsjsons = []

    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    outputs = args.outputs
    if os.path.isfile(outputs):
        print(f'Making QCSchema for {outputs}')
        json_package = identify_package(outputs)
        out_json = json_package(outputs)
        out_json.get_json(debug=args.debug)  # Will trigger any errors before writing.
        if out_json.path not in error_files:
            all_qcsjsons.append(out_json)
    elif os.path.isdir(outputs):
        if outputs[-1] != '/':
            outputs += '/'
        
        if args.recursive:
            print(f'Looking for output files in {outputs}, recursively')
            all_outfiles = get_files(outputs, 'out', recursive=True)
        else:
            print(f'Looking for output files in {outputs}')
            all_outfiles = get_files(outputs, 'out', recursive=False)
        
        print(f'Found {len(all_outfiles)} output files')

        for outfile in all_outfiles:
            file_name = '.'.join(os.path.basename(outfile).split('.')[:-1])
            if save_dir == './':
                abs_path = os.path.dirname(os.path.abspath(outfile))
                if not args.overwrite \
                   and os.path.exists(f'{abs_path}/{file_name}.json'):
                    print(f'\n\u001b[36;1m{file_name}.json already exists.\u001b[0m')
                    continue
            else:
                if not args.overwrite \
                   and os.path.exists(f'{save_dir}/{file_name}.json'):
                    print(f'\n\u001b[36;1m{file_name}.json already exists.\u001b[0m')
                    continue
            print(f'\nMaking QCSchema for {file_name}')
            json_package = identify_package(outfile)
            out_json = json_package(outfile)
            out_json.parse_output()
            out_json.get_json(debug=args.debug)  # Will trigger any errors before writing.
            if out_json.path not in error_files:
                all_qcsjsons.append(out_json)
            
    else:
        raise ValueError(f'{outputs} is an unsupported type.')
    
    print('\nWriting all valid QCJSONs')
    for out_json in all_qcsjsons:
        if save_dir == './':
            abs_path = os.path.dirname(out_json.path)
            out_json.write(
                out_json.name, out_json.get_json(debug=args.debug), abs_path,
                prettify=args.prettify
            )
        else:
            out_json.write(
                out_json.name, out_json.get_json(debug=args.debug), save_dir,
                prettify=args.prettify
            )

    if args.combine:
        pass
    
    print(f'\n{len(error_files)} file(s) encountered errors and not written')
    for i in error_files:
        i_name = os.path.basename(i)
        print(f'\u001b[31;1m    {i_name}\u001b[0m')


if __name__ == "__main__":
    main()
