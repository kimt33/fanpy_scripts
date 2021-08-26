import os
import re
import glob
import shutil
import subprocess
import numpy as np
import make_xyz
from make_com import make_com


def make_dirs(name: str, start_template: str, end_template: str, basis: str, num_steps=10, homedir='database'):
    """Make the directory, xyz, and gbs file for each step in the path between the two templates.

    Parameters
    ----------
    name : str
        Name of the calculations.
    start_template : str
        Name of the XYZ file that will be used as a template for the starting point.
    end_template : str
        Name of the XYZ file that will be used as a template for the end point.
    basis : str
        Name of the basis set.
    num_steps : int
        Number of steps in the path.

    Notes
    -----
    xyz file is stored as 'system.xyz2'.

    """
    start_template_index = int(start_template.split('_')[1])
    end_template_index = int(end_template.split('_')[1])

    dir_basename = os.path.join(homedir,
                                f'{name}_{start_template_index}{end_template_index}_{{}}_{basis}')
    for i, xyz in enumerate(make_xyz.xyz_from_templates(start_template, end_template, num_steps)):
        dirname = dir_basename.format(i)

        # if directory already exists
        if os.path.isdir(dirname):
            # find all others
            other_dirnames = glob.glob(dir_basename.format('*'))
            # check if xyz file matches any of the others
            has_matched = False
            for other_dirname in other_dirnames:
                xyzfile = os.path.join(other_dirname, 'system.xyz2')
                if os.path.isfile(xyzfile):
                    with open(xyzfile, 'r') as f:
                        if f.read() == xyz:
                            has_matched = True
                            break
            # if xyz file does not match any of the others
            else:
                # change the index (move to the end)
                i += len(other_dirnames)
                # update directory name
                dirname = dir_basename.format(i)

            if has_matched:
                continue

        # create directory
        os.mkdir(dirname)

        # make xyz
        xyzfile = os.path.join(dirname, 'system.xyz2')
        with open(xyzfile, 'w') as f:
            f.write(xyz)

        # make gbs
        gbsfile = os.path.join('basis', basis + '.gbs')
        # if basis set exists in basis directory
        if os.path.isfile(gbsfile):
            shutil.copyfile(gbsfile, os.path.join(dirname, basis + '.gbs'))
        else:
            print(f'Cannot find `.gbs` file that correspond to the given basis, {basis}')


def write_coms(pattern: str, memory='2gb', charge=0, multiplicity=1, units="AU"):
    """Write the Gaussian com files for the directories that match the given pattern.

    Parameters
    ----------
    pattern : str
        Pattern for selecting the directories.
    memory : str
        Amount of memory available to run the calculation.
        It should end in 'gb' or 'mb'
        Default is '2gb'
    charge : int
        Total charge of the molecule.
        Default is 0.
    multiplicity : int
        Spin (2n+1) of the molecule.
        Default is 0.

    Notes
    -----
    The directories are assumed to have the following format: `./system_templates_index_basis/`
    where

      - `system` is the chemical structure
      - `templates` is the indices for the templates used
      - `index` is the  index of the point in the path
      - `basis` is the basis set used

    Gaussian will only be used to run HF.

    """
    for parent in glob.glob(pattern):
        if not os.path.isdir(parent):
            continue
        *_, system_templates_index, basis = os.path.normpath(parent).split(os.path.sep)

        dirname = os.path.join(parent, 'hf')
        # make directory if it does not exist
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        # get xyz
        with open(os.path.join(parent, '..', 'system.xyz2'), 'r') as f:
            xyz_content = f.read()
        # get com content
        com_content = make_com(xyz_content, os.path.join("/home/kimt1/basis", basis), chkfile='hf_sp.chk', memory=memory,
                               title=f'HF/{basis} calculation for {system_templates_index}',
                               charge=charge, multiplicity=multiplicity, units=units)
        # make com file
        with open(os.path.join(dirname, 'hf_sp.com'), 'w')as f:
            f.write(com_content)


def make_orb_dirs(pattern: str, name: str):
    for parent in glob.glob(pattern):
        newdir = os.path.join(parent, name)
        if os.path.isdir(newdir):
            continue

        if not (os.path.isfile(os.path.join(parent, 'hf', 'oneint.npy')) and
                os.path.isfile(os.path.join(parent, 'hf', 'twoint.npy'))) and name != 'hf':
            continue

        os.mkdir(newdir)


def make_wfn_dirs(pattern: str, wfn_name: str, num_runs: int):
    """Make directories for running the wavefunction calculations.

    Parameters
    ----------
    pattern : str
        Pattern for the directories on which the new directories will be created.
        These directories must contain the files `oneint.npy` and `twoint.npy`.
    wfn_name : str
        Name of the wavefunction.
    num_runs : int
        Number of calculations that will be run.

    """
    for parent in glob.glob(pattern):
        if not os.path.isdir(parent):
            continue
        if not (os.path.isfile(os.path.join(parent, 'hf', 'oneint.npy')) and
                os.path.isfile(os.path.join(parent, 'hf', 'twoint.npy'))):
            continue

        newdir = os.path.join(parent, wfn_name)
        if not os.path.isdir(newdir):
            os.mkdir(newdir)

        for i in range(num_runs):
            try:
                os.mkdir(os.path.join(newdir, str(i)))
            except FileExistsError:
                pass


def write_wfn_py(pattern: str, nelec: int, wfn_type: str, optimize_orbs: bool=False,
                 pspace_exc=None, objective=None, solver=None, nproj=None,
                 ham_noise=None, wfn_noise=None,
                 solver_kwargs=None, wfn_kwargs=None,
                 load_orbs=None, load_ham=None, load_wfn=None, load_chk=None, load_prev=False,
                 memory=None, filename=None, ncores=1, exclude=None, old_fanpy=False):
    """Make a script for running calculations.

    Parameters
    ----------
    nelec : int
        Number of electrons.
    one_int_file : str
        Path to the one electron integrals (for restricted orbitals).
        One electron integrals should be stored as a numpy array of dimension (nspin/2, nspin/2).
    two_int_file : str
        Path to the two electron integrals (for restricted orbitals).
        Two electron integrals should be stored as a numpy array of dimension
        (nspin/2, nspin/2, nspin/2, nspin/2).
    wfn_type : str
        Type of wavefunction.
        One of `fci`, `doci`, `mps`, `determinant-ratio`, `ap1rog`, `apr2g`, `apig`, `apsetg`, and
        `apg`.
    optimize_orbs : bool
        If True, orbitals are optimized.
        If False, orbitals are not optimized.
        By default, orbitals are not optimized.
        Not compatible with solvers that require a gradient (everything except cma).
    pspace_exc : list of int
        Orders of excitations that will be used to build the projection space.
        Default is first and second order excitations of the HF ground state.
    objective : str
        Form of the Schrodinger equation that will be solved.
        Use `system` to solve the Schrodinger equation as a system of equations.
        Use `least_squares` to solve the Schrodinger equation as a squared sum of the system of
        equations.
        Use `variational` to solve the Schrodinger equation variationally.
        Must be one of `system`, `least_squares`, and `variational`.
        By default, the Schrodinger equation is solved as system of equations.
    solver : str
        Solver that will be used to solve the Schrodinger equation.
        Keyword `cma` uses Covariance Matrix Adaptation - Evolution Strategy (CMA-ES).
        Keyword `diag` results in diagonalizing the CI matrix.
        Keyword `minimize` uses the BFGS algorithm.
        Keyword `least_squares` uses the Trust Region Reflective Algorithm.
        Keyword `root` uses the MINPACK hybrd routine.
        Must be one of `cma`, `diag`, `least_squares`, or `root`.
        Must be compatible with the objective.
    ham_noise : float
        Scale of the noise to be applied to the Hamiltonian parameters.
        The noise is generated using a uniform distribution between -1 and 1.
        By default, no noise is added.
    wfn_noise : bool
        Scale of the noise to be applied to the wavefunction parameters.
        The noise is generated using a uniform distribution between -1 and 1.
        By default, no noise is added.
    load_orbs : str
        Numpy file of the orbital transformation matrix that will be applied to the initial
        Hamiltonian.
        If the initial Hamiltonian parameters are provided, the orbitals will be transformed
        afterwards.
    load_ham : str
        Numpy file of the Hamiltonian parameters that will overwrite the parameters of the initial
        Hamiltonian.
    load_wfn : str
        Numpy file of the wavefunction parameters that will overwrite the parameters of the initial
        wavefunction.
    load_chk : str
        Numpy file of the checkpoint for the optimization.
    memory : str
        Memory available to run the calculation.

    """
    cwd = os.getcwd()

    if optimize_orbs:
        optimize_orbs = ['--optimize_orbs']
    else:
        optimize_orbs = []

    if pspace_exc is None:
        pspace_exc = [1, 2, 3, 4]
    pspace_exc = [str(i) for i in pspace_exc]

    if objective is None:
        objective = 'variational'

    if solver is None:
        solver = 'cma'
    if solver == 'cma' and solver_kwargs is None:
        solver_kwargs = ("sigma0=0.01, options={'ftarget': None, 'timeout': np.inf, "
                         "'tolfun': 1e-11, 'verb_filenameprefix': 'outcmaes', 'verb_log': 0}")

    load_files = []
    if load_orbs:
        load_files += ['--load_orbs', load_orbs]
    if load_ham:
        load_files += ['--load_ham', load_ham]
    if load_wfn:
        load_files += ['--load_wfn', load_wfn]
    if load_chk:
        load_files += ['--load_chk', load_chk]

    if memory is None:
        memory = []
    else:
        memory = ['--memory', memory]

    kwargs = []
    if wfn_kwargs is not None:
        kwargs += ['--wfn_kwargs', wfn_kwargs]
    if solver_kwargs is not None:
        kwargs += ['--solver_kwargs', solver_kwargs]
    if ham_noise is not None:
        kwargs += ['--ham_noise', str(ham_noise)]
    if wfn_noise is not None:
        kwargs += ['--wfn_noise', str(wfn_noise)]

    for parent in glob.glob(pattern):
        if not os.path.isdir(parent):
            continue
        if exclude and any(i in os.path.abspath(parent) for i in exclude):
            #print(parent, 'x')
            continue

        os.chdir(parent)

        if filename is None:
            filename = 'calculate.py'

        # check if final directory is a number
        _, dirname = os.path.split(parent)
        
        if os.path.isfile('../hf/oneint.npy'):
            oneint = os.path.abspath('../hf/oneint.npy')
            twoint = os.path.abspath('../hf/twoint.npy')
            hf_energies = os.path.abspath('../hf/hf_energies.npy')
        else:
            oneint = os.path.abspath('../../hf/oneint.npy')
            twoint = os.path.abspath('../../hf/twoint.npy')
            hf_energies = os.path.abspath('../../hf/hf_energies.npy')

        nspin = np.load(oneint).shape[1] * 2
        nucnuc = np.load(hf_energies)[1]

        if load_prev:
            curr_sys = os.path.split(os.path.abspath("../../"))[1]
            re_curr_sys = re.search("(\w+_\w+_)(\d+)(_[\w-]+)", curr_sys)
            curr_index = int(re_curr_sys.group(2))
            prev_sys = re_curr_sys.group(1) + str(curr_index - 1) + re_curr_sys.group(3)
            prev_path = os.path.abspath(os.path.join("../../../", prev_sys, *parent.split(os.path.sep)[-3:]))
            load_ham = os.path.join(prev_path, "0", "hamiltonian.npy")
            load_wfn = os.path.join(prev_path, "0", "wavefunction.npy")
            kwargs += ["--load_ham", load_ham, "--load_wfn", load_wfn]

        save_ham = 'hamiltonian.npy'
        save_wfn = 'wavefunction.npy'
        save_chk = 'checkpoint.npy'

        #if load_chk:
        #    index = len(glob.glob('./results*.out'))
        #    for i in glob.glob('./*'):
        #        name, ext = os.path.splitext(i)
        #        if name[-1].isdigit():
        #            continue
        #        if 'checkpoint' in name:
        #            shutil.copy(i, name + str(index) + ext)
        #        else:
        #            shutil.move(i, name + str(index) + ext)
        #    #save_ham = 'hamiltonian{}.npy'.format(index+1)
        #    #save_wfn = 'wavefunction{}.npy'.format(index+1)
        #    #save_chk = 'checkpoint{}.npy'.format(index+1)

        #    if not glob.glob(load_chk) and load_chk != 'xxxx':
        #        os.chdir(cwd)
        #        continue

        print(parent)
        if nproj:
            pspace = ['--nproj', str(nproj)]
        else:
            pspace = ['--pspace', *pspace_exc]
        subprocess.run([#'wfns_make_script',
                        # 'python', '/home/kimt1/fanpy/wfns/scripts/make_script.py',
                        'fanpy_make_script' if old_fanpy else 'fanpy_make_fanci_script',
                        '--nelec', str(nelec), # '--nspin', str(nspin),
                        '--one_int_file', oneint, '--two_int_file', twoint,
                        '--nuc_repulsion', f'{nucnuc}', *optimize_orbs, '--wfn_type', wfn_type,
                        '--objective', objective,
                        '--solver', solver, *kwargs,
                        *load_files,
                        # '--save_ham', save_ham,
                        # '--save_wfn', save_wfn,
                        '--save_chk', save_chk,
                        '--filename', filename, *memory, 
                        *pspace,
                        # '--ncores', str(ncores),
                        ])

        os.chdir(cwd)


def run_calcs(pattern: str, time=None, memory=None, ncores=1, outfile='outfile', exclude=None, calc_range=None):
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
        memory = memory.upper()
        if memory[-2:] == 'GB':
            memory = str(int(float(memory[:-2]) * 1024)) + 'MB'
        elif memory[-2:] not in ['MB']:
            raise ValueError('Memory must be given as a MB or GB (e.g. 1024MB, 1GB)')
        print(memory)
        submit_job = True
    elif time is None and memory is None:
        submit_job = False
    else:
        raise ValueError('You cannot provide only one of the time and memory.')

    for filename in glob.glob(pattern):
        if os.path.commonpath([cwd, os.path.abspath(filename)]) != cwd:
            continue
        filename = os.path.abspath(filename)[len(cwd)+1:]
        print(filename)
        if exclude and any(i in filename for i in exclude):
            continue

        database, system, basis, *wfn = filename.split(os.sep)
        if os.path.isdir(filename):
            os.chdir(filename)
        else:
            dirname, filename = os.path.split(filename)
            os.chdir(dirname)

        if wfn[0] == 'hf' and os.path.splitext(filename)[1] == '.com':
            # write script (because sbatch only takes one command)
            if submit_job:
                with open('hf_sp.sh', 'w') as f:
                    f.write('#!/bin/bash\n')
                    f.write(f'g16 {filename}\n')
                command = ['hf_sp.sh']
            else:
                command = ['g16', filename]
        elif wfn[0] == 'hf' and os.path.splitext(filename)[1] == '.chk':
            command = ['formchk', filename]
        elif wfn[0] == 'hf' and os.path.splitext(filename)[1] == '.fchk':
            command = [os.environ.get('HORTONPYTHON'),
                       #'/home/kimt1/fanpy/wfns/wrapper/horton_gaussian_fchk.py',
                       '/home/kimt1/fanpy/fanpy/tools/wrapper/horton_gaussian_fchk.py',
                       'hf_energies.npy', 'oneint.npy', 'twoint.npy', 't_ab_mo.npy', 'fchk_file', filename]
        elif len(wfn) == 1:
            with open('results.sh', 'w') as f:
                f.write('#!/bin/bash\n')
                f.write('cd 50\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../51\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../52\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../53\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../54\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../55\n')
                f.write('python -u calculate.py > results.out\n')
                f.write('cd ../56\n')
                f.write('python -u calculate.py > results.out\n')
            subprocess.run(['chmod', 'u+x', 'results.sh'])
            command = ['./results.sh']
        elif len(wfn) == 2:
            if os.path.splitext(filename)[1] == '.py':
                with open('results.sh', 'w') as f:
                    f.write('#!/bin/bash\n')
                    f.write('cwd=$PWD\n')
                    if calc_range:
                        f.write(f'for i in {{{calc_range[0]}..{calc_range[1]}}}; do\n')
                    else:
                        f.write('for i in */; do\n')
                    f.write('    cd $i\n')
                    f.write('    python -u ../calculate.py > results.out\n')
                    f.write('    cd $cwd\n')
                    f.write('done\n')
            elif os.path.isfile('./calculate.py'):
                with open('results.sh', 'w') as f:
                    f.write('#!/bin/bash\n')
                    f.write('python -u ./calculate.py > results.out\n')
                    #f.write('python -u ./calculate.py\n')
            elif os.path.splitext(filename)[1] == '.com':
                with open('results.sh', 'w') as f:
                    f.write('#!/bin/bash\n')
                    f.write('g16 fci_spin1.com\n')
            else:
                os.chdir(cwd)
                continue
            subprocess.run(['chmod', 'u+x', 'results.sh'])
            command = ['./results.sh']

        #print(' '.join(['sbatch', f'--time={time}', f'--output={outfile}', f'--mem={memory}',
        #                 '--account=rrg-ayers-ab', *command]))
        print(os.getcwd())
        if submit_job:
            subprocess.run(['sbatch', f'--time={time}', f'--output={outfile}', f'--mem={memory}',
                            '--account=rmirandaquintana', f'--cpus-per-task={ncores}', *command])
        else:
            subprocess.run(command)

        # change directory
        os.chdir(cwd)
