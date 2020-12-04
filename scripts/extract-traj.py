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
import argparse
import numpy as np

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

def _string_coords(z, R):
    """Puts atomic coordinates into a Python string. Typically used for 
    writing to an input file.
    
    Parameters
    atoms : :obj:`numpy.ndarray`
        A (n,) numpy array containing all ``n`` elements labled by their atomic 
        number.
    coords : :obj:`numpy.array`
        Contains atomic positions in a (n, 3) numpy array where the x, y, and z 
        Cartesian coordinates in Angstroms are given for the n atoms.
    
    Returns
    -------
    :obj:`str`
        XYZ atomic coordinates as a string.
    """
    atom_coords_string = ''
    atom_index = 0
    while atom_index < len(z):
        atom_element = str(_z_to_element[z[atom_index]])
        coords_string = np.array2string(
            R[atom_index],
            suppress_small=True, separator='   ',
            formatter={'float_kind':'{:0.9f}'.format}
        )[1:-1] + '\n'
        atom_coords_string += (atom_element + '   ' \
                               + coords_string).replace(' -', '-')
        atom_index += 1
    
    return atom_coords_string

def _atomlabel_to_z(atom_labels):
    """Convert atom labels (e.g., O1, H1, H2) to their atomic number.

    Parameters
    ----------
    atom_labels : :obj:`numpy.ndarray`
        Strings of atom labels with one dimension.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Atomic numbers of atom labels.
    """
    if atom_labels.ndim != 1:
        raise ValueError('Array must have only one dimension.')
    z_list = []
    for atom_label in atom_labels:
        element = ''.join(filter(lambda x: not x.isdigit(), atom_label))
        z_list.append(_element_to_z[element])
    return np.array(z_list)
    

def get_traj(data):
    """
    Not all npz files will store the coordinates or atom specifications the
    same way. This function will combine the right data to get a full trajectory
    based on the package. If the "store-traj" script is modified to change the
    labels of data, this function will NOT work.

    Parameters
    ----------
    data : :obj:`dict`
        Dictionary of :obj:`numpy.ndarray` data from trajectory.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Atomic numbers of all atoms in system.
    :obj:`numpy.ndarray`
        Atomic Cartesian coordinates in the same order as the atomic numbers.
    """
    z = np.array([])
    R = None
    if data['package'][()].lower() == 'gamess':
        if data['n_qm'][()] != 0:
            z = np.concatenate((z, data['z_qm'].flatten()))
            if R is None:
                R = data['R_qm']
        if data['n_frag'] != 0:
            z_frag = _atomlabel_to_z(data['z_frag'])
            z = np.concatenate((z, z_frag.flatten()))
            if R is None:
                R = data['R_frag']
            else:
                R = np.concatenate((R, data['R_frag']), axis=1)
        if data['n_mm'] != 0:
            raise ValueError('MM atoms are not supported.')
    return z, R


def _write_traj(z, R, t, t_unit, E, E_unit, filename, save_dir):
    """

    Parameters
    ----------
    z : :obj:`numpy.ndarray`
        Atomic numbers of all atoms in system.
    R : :obj:`numpy.ndarray`
        Atomic Cartesian coordinates of the trajectory.
    t : :obj:`numpy.ndarray`
        Simulation time.
    t_unit : :obj:`numpy.ndarray`
        Units of time.
    E : :obj:`numpy.ndarray`
        Total energies of the snapshots.
    E_unit : :obj:`
    """
    if save_dir[-1] != '/':
        save_dir += '/'
    atom_num = z.shape[0]
    write_lines = []
    for i in range(R.shape[0]):
        write_lines.append(str(atom_num) + '\n')
        write_lines.append(f'time = {t[i]} {t_unit}; E_total = {E[i]} {E_unit}\n')
        write_lines.append(_string_coords(z, R[i]))
    with open(f'{save_dir}{filename}-traj.xyz', 'w') as f:
        f.writelines(write_lines)

def export_traj(z, R, data, filename, save_dir):
    """Handles writing and selecting timing and energies data.

    Parameters
    ----------
    z : :obj:`numpy.ndarray`
        Atomic numbers of all atoms in system.
    R : :obj:`numpy.ndarray`
        Atomic Cartesian coordinates of the trajectory.
    data : :obj:`dict`
        Dictionary of :obj:`numpy.ndarray` data from trajectory. Used to get
        energies and times for comments.
    filename : :obj:`str`
        Name of the npz file.
    save_dir : :obj:`str`
        Directory to save the trajectory.
    """

    if data['package'][()].lower() == 'gamess':
        E = data['E_total']
        E_unit = str(data['unit_E'][()])
        t = data['time']
        t_unit = str(data['unit_time'][()])
    _write_traj(z, R, t, t_unit, E, E_unit, filename, save_dir)


def main():
    
    parser = argparse.ArgumentParser(
        description=(
            'Extract data from a trajectory npz file.\n\n'
            'Data options are:\n'
            '* traj\n'
            '    * Will save a XYZ trajectory from the npz with energies as '
            'comments'
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'npz_traj', metavar='npz traj', type=str, nargs='?',
        help='Path to NumPy npz file with stored trajectory data.'
    )
    parser.add_argument(
        'traj_data', metavar='data', type=str, nargs='?',
        help='What data do you want to extract?'
    )
    args = parser.parse_args()

    npz_traj_path = args.npz_traj
    npz_traj_name = '.'.join(npz_traj_path.split('/')[-1].split('.')[:1])
    if not os.path.exists(npz_traj_path):
        raise ValueError(f'{npz_traj_path} does not exist.')
    data = dict(np.load(npz_traj_path))

    if args.traj_data == 'traj':
        z, R = get_traj(data)
        export_traj(z, R, data, npz_traj_name, '.')
    

if __name__ == "__main__":
    main()
