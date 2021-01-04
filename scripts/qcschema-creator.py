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
how to manage these files. Our immediate solution is to nest each iteration as
in its own QCSchema; each incrementally labeled an integer starting
from zero (e.g., '2'). Furthermore, only iterations that have the same topology
are supported (no checks are provided).

## Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.1] - 2021-01
### Added
- ORCA Schema creator for single-point energies and gradients.

"""

import os
from functools import reduce
import argparse
from packaging import version
import cclib
import numpy as np
import json

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

### QCSchema classes ###

class QCSchema:
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
    
    def parse_output(self, outfile_path):
        """Parse output file using cclib.

        Parameters
        ----------
        outfile_path : :obj:`str`
            Path to computational chemistry output file.
        """
        filename_with_extension = os.path.basename(outfile_path)
        self.path = os.path.abspath(outfile_path)
        self.name = '.'.join(filename_with_extension.split('.')[:-1])
        self.outfile = cclib.io.ccread(outfile_path)
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

class orcaSchema(QCSchema):
    """ORCA specific QCSchema information.

    Supported calculation types:

    * Single-point energies
    * Energy+Gradient
    """

    def __init__(self):
        super().__init__()
    
    def topology(self, iteration=-1):
        """Prepares basic topology information always available from cclib.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
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

    def model(self, iteration=-1):
        """Model chemistry used for the calculation.
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
    
    def keywords(self, iteration=-1):
        """Parses parameters for the calculation.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        """
        keywords = {
            'e_change': self.outfile.scftargets[iteration][0],
            'max_density_change': self.outfile.scftargets[iteration][1],
            'rms_density_change': self.outfile.scftargets[iteration][2]
        }
        _remove_keywords = []
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
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
            
            if kw_lower in ['normalscf', 'loosescf', 'sloppyscf', 'strongscf',
                            'tightscf', 'verytightscf', 'extremescf'] \
                or 'scfconv' in kw_lower:
                _remove_keywords.append(kw)
            
            if 'grid' in kw_lower:
                kw_lower_split = kw_lower.split('grid')
                if len(kw_lower_split) == 2 and int(kw_lower_split[-1]):
                    pass
        
        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
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
        keywords['other'] = self.orca_keywords
        return keywords
    
    def properties(self, iteration=-1):
        """Prepares basic properties information always available from cclib.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        """
        properties = {}
        if self.method_type == 'scf':
            properties['scf_total_energy'] = cclib.parser.utils.convertor(
                self.outfile.scfenergies[iteration], 'eV', 'hartree'
            )
            properties['scf_dispersion_correction_energy'] = cclib.parser.utils.convertor(
                self.outfile.dispersionenergies[iteration], 'eV', 'hartree'
            )
            properties['scf_iterations'] = self.outfile.scfvalues[0].shape[0]
        elif self.method_type == 'moller-plesset':
            if self.outfile.mpenergies.ndim == 1:
                mpenergy = self.outfile.mpenergies[iteration]
                scfenergy = self.outfile.scfenergies[iteration]
            elif self.outfile.mpenergies.ndim == 2 \
                 and self.outfile.mpenergies.shape[1] == 1:
                mpenergy = self.outfile.mpenergies[iteration][0]
                scfenergy = self.outfile.scfenergies[iteration]
            properties['scf_total_energy'] = cclib.parser.utils.convertor(
                scfenergy, 'eV', 'hartree'
            )
            properties['mp2_total_energy'] = cclib.parser.utils.convertor(
                mpenergy, 'eV', 'hartree'
            )
        elif self.method_type == 'coupled cluster':
            pass
        else:
            error_out(self.path, 'Unknown method type.')
        properties['calcinfo_nbasis'] = self.outfile.nbasis
        properties['calcinfo_nmo'] = self.outfile.nmo
        return properties
    
    def driver(self, iteration=-1):
        """The purpose of the calculation and its direct result.

        Parameters
        ----------
        iteration: :obj:`int`, optional
            Defaults to the last iteration.
        """
        driver = {}
        _remove_keywords = []
        for kw in self.orca_keywords:
            kw_lower = kw.lower()
            if kw_lower == 'engrad' or kw_lower == 'numgrad':
                driver['driver'] = 'gradient'
                if self.outfile.grads.ndim == 3:
                    grads = self.outfile.grads[iteration]
                else:
                    raise ValueError('Please check gradient dimensions.')
                driver['return_result'] = convert_forces(
                    grads, 'hartree', 'bohr', 'hartree', 'Angstrom'
                )
                _remove_keywords.append(kw)
                break
            elif kw_lower == 'freq' or kw_lower == 'numfreq':
                _remove_keywords.append(kw)
                break
            elif kw_lower == 'opt' or kw_lower == 'copt' or kw_lower == 'zopt':
                driver['driver'] = 'optimization'
                _remove_keywords.append(kw)
                break
            elif kw_lower == 'energy' or kw_lower == 'sp':
                driver['driver'] = 'energy'
                _remove_keywords.append(kw)
                break
        
        # ORCA defaults to single-point energies.
        if len(driver) == 0:
            driver['driver'] == 'energy'
            pass

        for kw in _remove_keywords:
            self.orca_keywords.remove(kw)
        
        return driver

    def provenance(self):
        """Program specification that performed the calculation.
        """
        provenance = {
            'creator': 'ORCA',
            'version': self.outfile.metadata['package_version']
        }
        return provenance
    
    @property
    def schema(self):
        """QC schema of an ORCA output file.

        Calculations supported: single-point energies.

        :type: :obj:`dict`
        """
        if not hasattr(self, '_schema'):
            schema_dicts = {}
            self.orca_keywords = self.outfile.metadata['keywords']
            for i in range(0, self.outfile.atomcoords.shape[0]):
                try:
                    schema_dicts[i] = super().schema
                    schema_dicts[i] = {
                        **schema_dicts[i], **self.topology(iteration=i)
                    }
                    schema_dicts[i] = {
                        **schema_dicts[i], **self.driver(iteration=i)
                    }
                    schema_dicts[i]['provenance'] = self.provenance()
                    schema_dicts[i]['model'] = self.model(iteration=i)
                    schema_dicts[i]['keywords'] = self.keywords(iteration=i)
                    schema_dicts[i]['properties'] = self.properties(iteration=i)
                    schema_dicts[i]['success'] = self.outfile.metadata['success']
                    if len(schema_dicts) == 1:
                        self._schema = schema_dicts[0]
                    else:
                        self._schema = schema_dicts
                except:
                    error_out(self.path, 'Uncaught exceptions.')
                    self._schema = {}
                    break
        return self._schema

### Runtime Functions ###

# Triggers to identify output files.
triggers = [
    (orcaSchema, ["O   R   C   A"], True)
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
        '--overwrite', action='store_true', help='Overwrite JSON files'
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
    args = parser.parse_args()
    print(f'QCSchema creator v{__version__}')
    print('Written by Alex M. Maldonado (@aalexmmaldonado)')
    print('Parsed data are converted to Hartrees and Angstroms\n')

    cclib_version_check()

    all_qcschemas = []

    save_dir = args.save_dir
    if save_dir[-1] != '/':
        save_dir += '/'
    outputs = args.outputs
    if os.path.isfile(outputs):
        print(f'Making QCSchema for {outputs}')
        schema = identify_package(outputs)
        out_schema = schema()
        out_schema.parse_output(outputs)
        out_schema.schema  # Will trigger any errors before writing.
        if out_schema.path not in error_files:
             all_qcschemas.append(out_schema)
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
            schema = identify_package(outfile)
            out_schema = schema()
            out_schema.parse_output(outfile)
            out_schema.schema  # Will trigger any errors before writing.
            if out_schema.path not in error_files:
                all_qcschemas.append(out_schema)
            
    else:
        raise ValueError(f'{outputs} is an unsupported type.')
    
    print('\nWriting all valid QCSchemas')
    for schema in all_qcschemas:
        if save_dir == './':
            abs_path = os.path.dirname(schema.path)
            schema.write(
                schema.name, schema.schema, abs_path, prettify=args.prettify
            )
        else:
            schema.write(
                schema.name, schema.schema, save_dir, prettify=args.prettify
            )

    if args.combine:
        pass
    
    print(f'\n{len(error_files)} file(s) encountered errors and not written')
    for i in error_files:
        i_name = os.path.basename(i)
        print(f'\u001b[31;1m    {i_name}\u001b[0m')


if __name__ == "__main__":
    main()
