import sys
sys.path.append("/home/kimt1/scripts_fanpy_db")

import os
import re
import json
import glob

import run_calc
import extract_energy
import edit_rbm


# RBM here is used as a multiplicative correction to an existing wavefunction
# given wavefunction
# |Psi> = \sum_i f(m_i) |m_i> where m_i is a slater determinent
# RBM multiplicative correction corresponds to
# |Psi_{rbm_corrected}> = \sum_i f(m_i) g(m_i) |m_i> where g(m_i) is the rbm overlap for given
# slater determinant

# system/job information
time = "7d"
memory = "8gb"  # total memory available
memory_for_wfn_cache = "6gb"  # memory allocated to storing wavefunction overlaps

# database structure:
# we will assume that the database/directories are structured as follows:
# db_directory/system/basis/wfn/repetitions/results
db_directory = "beh2"

# system information
nelec = 6
basis = "sto-6g"
units = "AU"
charge = 0
multiplicity = 1

# wavefunction to which rbm will multiplicatively correct
# we will assume that this wavefunciton has already been calculated
# i.e. there exists results for this wavefunction in the directory structure above
# from a given calculation result, there exists, ../../wfn/reptitions/results exists
wfn = "ap1rog"
wfn_dirname = "ap1rog"
# wfn_dirname is used to find the orbital optimized wavefunction results
# if orbital_optimization is False, then wfn_dirname = wfn (probably)
orbital_optimization = False

# rbm parameters
num_layers = 1  # number of layers in the rbm wavefunciton
orders = [1, 2, 3]  # orders of correlations accounted for in rbm wavefunction
orders_str = "1,2,3"

# initial guess parameters
wfn_noise = 0.0
ham_noise = 0.0  # not used if no orbital optimization

cwd = os.getcwd()

# 1. make Gaussian input
# pay special attention to units
#run_calc.write_coms(f'{db_directory}/*/{basis}/', memory='1gb', charge=charge, multiplicity=multiplicity, units=units)

# 2. run/submit Gaussian calculation
# If memory and runtime not given, runs in shell
#run_calc.run_calcs(f'{db_directory}/*/{basis}/hf/hf_sp.com')

# 3. make fchk
#run_calc.run_calcs(f'{db_directory}/*/{basis}/hf/hf_sp.chk')

# 4. make integrals (npy) using fchk and horton
#run_calc.run_calcs(f'{db_directory}/*/{basis}/hf/hf_sp.fchk')

# 5. make directories for wavefunction
run_calc.make_wfn_dirs(f'{db_directory}/*/{basis}/', f"{wfn_dirname}_rbm_{orders_str}_{num_layers}", 10, rep_dirname_prefix="energy_min") 

# 6. make rbm script (fanpy)
# define types of calculations to run in each directory label in "repetitions"
# i.e. following runs a energy minimization calculation in the directory energy_min
# key = name of the directory to be created/to store results
# value = 2-tuple of the objective and the optimization algorithm
# TODO: list supported objectives and algorithms

calc_type = {"*": ("energy", "minimize"), "projected": ("projected", "least_squares")}
for path_full in sorted(glob.glob(f'{db_directory}/*/{basis}/{wfn_dirname}_rbm_{orders_str}_{num_layers}'), reverse=True):
    print(path_full)
    if orbital_optimization:
        # replace energy objective to one_energy
        keys = []
        for key, (objective, algorithm) in calc_type.items():
            if "energy" in objective:
                keys.append(key)
        for key in keys:
            calc_type[key] = ("one_energy", *calc_type[key][1:])

    for dirpattern_results, (objective, solver) in calc_type.items():
        # we will use the "best" result from the wavefunction reference
        # here, we define "best" as one with lowest energy
        chk_energy = {}
        # NOTE: make sure to glob for files instead of directories
        # NOTE: we will assume that reference wavefunction stores its results in the same pattern as
        #       the current pattern
        for results_path in glob.glob(f'{path_full}/../{wfn_dirname}/{dirpattern_results}/results.out'):
            with open(results_path, 'r') as g:
                lines = g.readlines()
            # find/store last energy computed for given checkpoint
            energy = None
            for line in lines:
                # try and extract energy from current line
                temp_energy = extract_energy.extract_energy(line)
                if temp_energy:
                    energy = temp_energy
            if energy is not None:
                # NOTE: we are always assuming that the results_path will end in file (not directory)
                chk_energy[os.path.join(os.path.split(results_path)[0], 'checkpoint.npy')] = float(energy)
            else:
                print(f"Didn't get energy for directory, {results_path}")
        if not chk_energy:
            print(f"No energy in path {path_full}/../{wfn_dirname}/{dirpattern_results}/results.out")
            continue
        # checkpoint file with lowest energy
        checkpoint_file = os.path.join(cwd, min(chk_energy, key=lambda i: chk_energy[i]))

        for dirname_results in glob.glob(f'{path_full}/{dirpattern_results}'):
            if orbital_optimization:
                run_calc.write_wfn_py(
                    dirname_results, nelec, wfn, optimize_orbs=True,
                    objective=objective, solver=solver,
                    # FIXME: hardcoded to be full projection space
                    pspace_exc=list(range(1, nelec + 1)),
                    ham_noise=ham_noise, wfn_noise=wfn_noise,
                    memory=memory_for_wfn_cache,
                    load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint_file,
                    # use old fanpy
                    old_fanpy=True,
                )
            else:
                run_calc.write_wfn_py(
                    dirname_results, nelec, wfn, optimize_orbs=False,
                    objective=objective, solver=solver,
                    # FIXME: hardcoded full projection space
                    # projection space size
                    # 0 = largest possible projection space
                    # positive integer = exact size of projection space
                    # negative number = number of parameters * given number (as positive number)
                    nproj=0,
                    wfn_noise=wfn_noise,
                    memory=memory_for_wfn_cache,
                    load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint_file,
                    # use fast fanpy (uses Michelle's fanci)
                    old_fanpy=False,
                )

            # modify the file to incorporate rbm multiplicative correction
            edit_rbm.edit_file(
                f'{dirname_results}/calculate.py', rbm_kwargs='num_layers=1, orders=(1, 2)',
                old_fanpy=orbital_optimization, optimize_both=False, wfn_noise=wfn_noise,
            )

            # submit calculations
            # TODO: different ways of providing path results in different types of
            #       calculations/job submissions
            #run_calc.run_calcs(
            #    dirname_results, time=time, memory=memory, outfile='results.out'
            #)
