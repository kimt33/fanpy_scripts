import os
import re
from glob import glob
from read_xyz import XYZReader
from input_g09 import GaussianInput
from subprocess import call
import subprocess

def run_fci_calcs(pattern: str, time=None, memory=None, outfile='outfile', units='AU'):
    """Run the calculations for the selected files/directories.

    Parameters
    ----------
    pattern : str
        Pattern for selecting the files.

    Notes
    -----
    Can only execute at the base directory.

    """
    cwd = os.getcwd()

    if time is not None and memory is not None:
        time = time.lower()
        if time[-1] == 'd':
            time = int(time[:-1]) * 24 * 60
        elif time[-1] == 'h':
            time = int(time[:-1]) * 60
        elif time[-1] == 'm':
            time = int(time[:-1])
        else:
            raise ValueError('Time must be given in minutes, hours, or days (e.g. 1440m, 24h, 1d).')
        if memory[-2:] not in ['mb', 'gb']:
            raise ValueError('Memory must be given as a MB or GB (e.g. 1024MB, 1GB)')
        submit_job = True
    elif memory is None:
        memory = "3gb"
    else:
        raise ValueError('You cannot provide only one of the time and memory.')

    for filename in glob(pattern):
        if os.path.commonpath([cwd, os.path.abspath(filename)]) != cwd:
            continue
        filename = os.path.abspath(filename)[len(cwd)+1:]

        _, system, basis, *wfn = filename.split(os.sep)
        if os.path.isdir(filename):
            os.chdir(filename)
        else:
            dirname, filename = os.path.split(filename)
            os.chdir(dirname)

        atoms, *_ = system.split('_')

        # get xyz
        with open(os.path.join('../', '..', '..', 'system.xyz2'), 'r') as f:
            xyz = f.readlines()
        xyz = [i.split() for i in xyz]

        # system info
        nelec = int(len(xyz))
        if basis == 'sto-6g':
            nbasis = nelec
        else:
            nbasis = 5 * nelec

        for spin in [1,3,5]:
            filename = 'fci_spin{0}.com'.format(spin)
            g09 = GaussianInput([i[0] for i in xyz], [i[1:] for i in xyz], method='cas({0},{1})'.format(nelec, nbasis), basis=basis, charge=0, spinmult=spin, route=f'Units={units} scf=tight geom=nocrowd', filename=filename, nosave=False)
            g09.write_input()

            if submit_job:
                #call(["sqsub", "-q", "NAP_9122", "-r", time, "-o", '{0}.out'.format(filename[:-4]),
                #      "--mpp={0}".format(memory), 'g16', filename])
                #with open('results.sh', 'w') as f:
                #    f.write('#!/bin/bash\n')
                #    f.write(f'g16 fci_spin{spin}.com\n')
                #subprocess.run(['sbatch', f'--time={time}', f'--output={outfile}', f'--mem={memory}',
                #                '--account=rmirandaquintana', './results.sh'])
                pass
            else:
                call(["g16", filename])
                #call(["sqsub", "-q", "serial", "-r", time, "-o", '{0}.out'.format(filename[:-4]),
                #      "--mpp={0}".format(memory), 'g09', filename])
                #call(["sqsub", "-q", "NAP_9122", "-r", time, "-o", '{0}.out'.format(filename[:-4]),
                #      "--mpp={0}".format(memory), 'g16', filename])
        if submit_job:
            with open('results.sh', 'w') as f:
                f.write('#!/bin/bash\n')
                f.write(f'g16 fci_spin1.com\n')
                f.write(f'g16 fci_spin3.com\n')
                f.write(f'g16 fci_spin5.com\n')
            subprocess.run(['sbatch', f'--time={time}', f'--output={outfile}', f'--mem={memory}',
                            '--account=rmirandaquintana', './results.sh'])

        # change directory
        os.chdir(cwd)
