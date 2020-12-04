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

# pylint: disable=no-member

def assign_file(package, source):
    """Initializes package output class.

    Parameters
    ----------
    package : :obj:`str`
        Name of the package that generated the trajectory file.
    source : :obj:`str`
        Path to the trajectory file.
    
    Returns
    -------
    The class corresponding to the correct package.
    """
    if package.lower() == 'gamess':
        return GAMESS(source)
    else:
        raise ValueError(f'{package} is not supported.')

class OutFile:
    """General package class.
    """

    def __init__(self, package, source):
        self.package = package.lower()
        self.source = source
        self.file_name = '.'.join(source.split('/')[-1].split('.')[:1])
        self.helps = []
        self.data = {}
     
    def parse(self):
        """Parses trajectory file.
        """
        with open(self.source, mode='r') as trajfile:
            for trajline in trajfile:
                self.extract(trajfile, trajline)
    
    def save(self, save_dir):
        """Saves data dictionary to ``npz`` file.

        Parameters
        ----------
        save_dir : :obj:`str`
            Path to save directory.
        """
        if save_dir[-1] != '/':
            save_dir += '/'
        data_arrays = {}
        for key in self.data.keys():
            data_arrays[key] = np.array(self.data[key])
        helps = '\n'.join(self.helps)
        help_string = 'Description of arrays.\n----------------------\n' + helps
        data_arrays['help'] = np.array(help_string)
        np.savez_compressed(f'{save_dir}{self.file_name}.npz', **data_arrays)

        

class GAMESS(OutFile):
    """For GAMESS `.trj` files.
    """

    def __init__(self, source):
        super(GAMESS, self).__init__('GAMESS', source)
        self.helps.append('R_units: Units of distance')
        self.data['R_units'] = 'Angstrom'
        self.data['package'] = 'GAMESS'
    
    def extract(self, trajfile, line):
        """Extracts information from trajectory file line.
        """

        if '===== MD DATA PACKET =====' == line.strip():
            line = next(trajfile)

            # NAT=      11 NFRG=     288 NQMMM=       0
            if 'NAT' in line and 'NFRG' in line and 'NQMMM' in line:
                if 'n_qm' not in self.data.keys():
                    _, n_qm, _, n_frag, _, n_mm = line.split()
                    self.helps.append('n_qm: Number of atoms in QM region.')
                    self.data['n_qm'] = int(n_qm)
                    self.helps.append('n_frag: Number of fragments in simulation.')
                    self.data['n_frag'] = int(n_frag)
                    self.helps.append('n_mm: Number of MM atoms.')
                    self.data['n_mm'] = int(n_mm)
                    if int(n_mm) != 0:
                        print('NQMMM data parsing is not implemented.')
                line = next(trajfile)
            
            # TTOTAL=        0.00 FS    TOT. E=      -186333.219846 KCAL/MOL
            if 'TTOTAL' in line and 'TOT. E' in line:
                _, time, unit_time, _, _, E_total, unit_E = line.split()
                if 'time' not in self.data.keys():
                    self.helps.append('time: Elapsed simulation time.')
                    self.data['time'] = [float(time)]
                    self.helps.append('unit_time: Units of time.')
                    self.data['unit_time'] = unit_time
                    self.helps.append('E_total: Total energy of the system.')
                    self.data['E_total'] = [float(E_total)]
                    self.helps.append('unit_E: Units of energy.')
                    self.data['unit_E'] = unit_E
                else:
                    self.data['time'].append(float(time))
                    self.data['E_total'].append(float(E_total))
                line = next(trajfile)
            
            # POT. E=       -186333.219846 KCAL/MOL  BATHT=           300.000000
            if 'POT. E' in line and 'BATHT' in line:
                _, _, E_pot, _, _, T_bath = line.split()
                if 'E_pot' not in self.data.keys():
                    self.helps.append('E_pot: Potential energy of system.')
                    self.data['E_pot'] = [float(E_pot)]
                    self.helps.append('T_bath: Temperature bath set point.')
                    self.data['T_bath'] = float(T_bath)
                    self.helps.append('unit_T: Units of temperature.')
                    self.data['unit_T'] = 'K'
                else:
                    self.data['E_pot'].append(float(E_pot))
                line = next(trajfile)
            
            # KIN. E= 0.000000  TRANS KE= 0.000000  ROT KE= 0.000000 KCAL/MOL
            if 'KIN. E' in line and 'TRANS KE' in line and 'ROT KET' in line:
                _, _, E_kin, _, _, E_trans, _, _, E_rot, _ = line.split()
                if 'E_kin' not in self.data.keys():
                    self.helps.append('E_kin: Total kinetic energy of the system.')
                    self.data['E_kin'] = [float(E_kin)]
                    self.helps.append('E_trans: Translational kinetic energy.')
                    self.data['E_trans'] = [float(E_trans)]
                    self.helps.append('E_rot: Rotational kinetic energy.')
                    self.data['E_rot'] = [float(E_rot)]
                else:
                    self.data['E_kin'].append(float(E_kin))
                    self.data['E_trans'].append(float(E_trans))
                    self.data['E_rot'].append(float(E_rot))
                line = next(trajfile)

        # ----- QM PARTICLE COORDINATES FOR $DATA GROUP -----
        # C     6.0       -1.5146148541       -0.1669730010        0.3060154811
        # O     8.0       -1.9008047124       -0.1401553155       -0.8253740685
        # O     8.0       -0.9540157833       -0.4301183474        1.2998418078
        # B     5.0       -2.7901505687        1.9430321247        1.6931211193
        # H     1.0       -2.9334610667        0.8472872458        1.0641355653
        # H     1.0       -2.9277068935        1.6378289859        2.8445898128
        # H     1.0       -3.6719142898        2.8007307337        1.5110977164
        # O     8.0       -4.8023343241        1.0602827630       -0.5528402168
        # H     1.0       -4.2707220381        1.3296475494        0.1273999357
        # H     1.0       -5.1327052621        1.7821538295       -1.1056120755
        # H     1.0       -1.7806190780        2.5703408558        1.4177843584
        if '----- QM PARTICLE COORDINATES FOR $DATA GROUP -----' == line.strip():
            line = next(trajfile)

            if 'z_qm' not in self.data.keys():
                self.helps.append('z_qm: Atomic numbers of QM atoms.')
                self.data['z_qm'] = []
                self.helps.append('R_qm: Cartesian coordinates QM atoms.')
                self.data['R_qm'] = []

            self.data['R_qm'].append([])
            while '----- ' not in line:
                line_split = line.split()
                if len(self.data['R_qm']) == 1:
                    self.data['z_qm'].append(int(float(line_split[1])))
                self.data['R_qm'][-1].append([float(i) for i in line_split[2:]])
                line = next(trajfile)
        
        # ----- EFP PARTICLE COORDINATES FOR $EFRAG GROUP -----
        #  $EFRAG
        # COORD=CART POSITION=OPTIMIZE
        # FRAGNAME=TIP5P   !   1
        # O1               4.6692170327       -8.7877969884        4.1576974625
        # H2               4.0176908798       -9.2696333531        3.6482126651
        # H3               4.1548603695       -8.2734631691        4.7798914149
        # FRAGNAME=TIP5P   !   2
        # O1              -3.4493126255       -5.7822171821       -1.2823666153
        # H2              -2.9274425084       -6.5840925584       -1.2527890761
        # H3              -2.9373357511       -5.1825732769       -1.8250767000
        if '----- EFP PARTICLE COORDINATES FOR $EFRAG GROUP -----' == line.strip():
            for _ in range(2):
                line = next(trajfile)
            if 'COORD=CART POSITION=OPTIMIZE' != line.strip():
                raise ValueError(f'{line} is not supported.')
            line = next(trajfile)

            if 'z_frag' not in self.data.keys():
                self.helps.append('z_frag: Atom labels of fragments.')
                self.data['z_frag'] = []
                self.helps.append('R_frag: Cartesian coordinates all fragments.')
                self.data['R_frag'] = []
            
            self.data['R_frag'].append([])
            while '$END' != line.strip():
                if 'FRAGNAME=' in line and '!' in line:
                    line = next(trajfile)
                else:
                    line_split = line.split()
                    if len(self.data['R_frag']) == 1:
                        self.data['z_frag'].append(line_split[0])
                    self.data['R_frag'][-1].append([float(i) for i in line_split[1:]])
                    line = next(trajfile)
            line = next(trajfile)
        
        #     GRADIENT DATA (NOT USED BY RESTARTS)...
        # FRAGMENT #     1  TIP5P 
        #     EFCENT    8.7002553693  -16.6030919250    7.8688272886
        #      FORCE    0.0003305144   -0.0014810269    0.0004887864
        #       TORQ   -0.0018311880   -0.0001070612    0.0109241528
        # FRAGMENT #     2  TIP5P 
        #     EFCENT   -6.4089323372  -10.9481908902   -2.4775823330
        #      FORCE    0.0007830076    0.0047822137    0.0031648219
        #       TORQ    0.0008248002    0.0045364043    0.0020023457
        
        if 'GRADIENT DATA (NOT USED BY RESTARTS)...' == line.strip():
            line = next(trajfile)

            if 'frag_efcent' not in self.data.keys():
                self.helps.append('frag_efcent: Unknown.')
                self.data['frag_efcent'] = []
                self.helps.append('frag_force: Force vector on fragment.')
                self.data['frag_force'] = []
                self.helps.append('frag_torq: Torque on fragment.')
                self.data['frag_torq'] = []
            
            self.data['frag_efcent'].append([])
            self.data['frag_force'].append([])
            self.data['frag_torq'].append([])
            while '----- RESTART VELOCITIES FOR $MD GROUP -----' != line.strip():
                if 'FRAGMENT #' in line:
                    line = next(trajfile)
                else:
                    line_split = line.split()
                    if 'EFCENT' in line:
                        self.data['frag_efcent'][-1].append(
                            [float(i) for i in line_split[1:]]
                        )
                    elif 'FORCE' in line:
                        self.data['frag_force'][-1].append(
                            [float(i) for i in line_split[1:]]
                        )
                    elif 'TORQ' in line:
                        self.data['frag_torq'][-1].append(
                            [float(i) for i in line_split[1:]]
                        )
                    line = next(trajfile)
        
        if '----- RESTART VELOCITIES FOR $MD GROUP -----' == line.strip():
            pass
            

def main():
    
    parser = argparse.ArgumentParser(
        description=(
            'Parse and store MD output files into a NumPy npz file.\n\n'
            'Only GAMESS is implemented.'
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'package', metavar='package', type=str, nargs='?',
        help='What software package is the trajectory file from?'
    )
    parser.add_argument(
        'outfile', metavar='traj file', type=str, nargs='?',
        help='Path to trajectory file.'
    )
    args = parser.parse_args()

    package = args.package
    outfile_path = args.outfile
    if not os.path.exists(outfile_path):
        raise ValueError(f'{outfile_path} does not exist.')

    traj = assign_file(package, outfile_path)
    traj.parse()
    traj.save('.')

if __name__ == "__main__":
    main()
