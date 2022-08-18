import sys
sys.path.append("/home/kimt1/scripts_fanpy_db")

import os
import re
import json
import glob

import run_calc
import extract_energy
import edit_cc


# CC here obtains the overlaps by explicitly evaluating the cumulent
# there are various flavors that ramon has developed

# supported wavefunctions:
# uncomment the you want
# "standardcc": CC by excitation orders
# "generalizedcc": generalized CC  # TODO: HOW TO DISTINGUISHh
# "basecc": arbitrary CC  # TODO: HOW TO DISTINGUISHh
# "senioritycc": CC by increasing seniority
# "pccd": PCCD
# "ap1rogsd": AP1roG + single and double excitation
# "ap1rogsd_spin": AP1roG + single and double excitations that preserves spin (alpha to alpha excitation and beta to beta excitation)
# "apsetgd": APsetG + double excitation
# "apsetgsd":  APsetG + single and double excitation
# "apg1rod": 1 reference orbital APG + double excitation
# "apg1rosd": 1 reference orbital APG + single and double excitation
# "ccsdsen0": CCS + seniority 0 double excitation
# "ccsdqsen0": CCS + seniority 0 double and quadruple excitation
# "ccsdtqsen0": CCST + seniority 0 double and quadruple excitation
# "ccsdtsen2qsen0": CCS + seniority 0 double and quadruple excitation + seniority 2 triple excitations
# "ccs": CCS
# "ccsd": CCSD
# "ccsdt": CCSDT
# "ccsdtq": CCSDTQ
wfn = "ccsdqsen0"

# wfn_dirname is used to find the orbital optimized wavefunction results
# if orbital_optimization is False, then wfn_dirname = wfn (probably)
wfn_dirname = wfn
# it is recommended to have orbital optimization off (False) if you have single excitations
orbital_optimization = False

# database structure:
# we will assume that the database/directories are structured as follows:
# db_directory/system/basis/wfn/repetitions/results.out
db_directory = "beh2"

# system information
nelec = 6
basis = "sto-6g"
units = "AU"
charge = 0
multiplicity = 1

# initial guess parameters
wfn_noise = 0.0
ham_noise = 0.0  # not used if no orbital optimization

cwd = os.getcwd()

# checkpoint (relative) path and the wavefunction
checkpoint_wfn = "ccsdsen0"
if checkpoint_wfn not in edit_cc.dict_str_cls:
    raise ValueError(f"Given checkpoint_wfn must be one of {list(edit_cc.dict_str_cls.keys())}")
# checkpoint is same wavefunction on different repetition
# checkpoint_path = f"../min_energy_0/checkpoint_{edit_cc.dict_str_cls[checkpoint_wfn].__name__}.npy"
# checkpoint is diff wavefucntion
checkpoint_path = f"../../ccsdsen0/*/checkpoint_{edit_cc.dict_str_cls[checkpoint_wfn].__name__}.npy"

# system/job information
time = "7d"
memory = "8gb"  # total memory available
memory_for_wfn_cache = "6gb"  # memory allocated to storing wavefunction overlaps

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
run_calc.make_wfn_dirs(f'{db_directory}/*/{basis}/', f"{wfn_dirname}", 10, rep_dirname_prefix="energy_min") 

# 6. make cc script (fanpy)
# define types of calculations to run in each directory label in "repetitions"
# i.e. following runs a energy minimization calculation in the directory energy_min
# key = name of the directory to be created/to store results
# value = 2-tuple of the objective and the optimization algorithm
# TODO: list supported objectives and algorithms

calc_type = {"energy_min_1": ("energy", "minimize"), "projected": ("projected", "least_squares")}
for path_full in sorted(glob.glob(f'{db_directory}/*/{basis}/{wfn_dirname}'), reverse=True):
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
        for dirname_results in glob.glob(f'{path_full}/{dirpattern_results}'):
            if orbital_optimization:
                run_calc.write_wfn_py(
                    dirname_results, nelec, wfn, optimize_orbs=True,
                    objective=objective, solver=solver,
                    # FIXME: hardcoded to be full projection space
                    pspace_exc=list(range(1, nelec + 1)),
                    ham_noise=ham_noise, wfn_noise=wfn_noise,
                    memory=memory_for_wfn_cache,
                    load_orbs=None, load_ham=None, load_wfn=None,
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
                    load_orbs=None, load_ham=None, load_wfn=None,
                    # use fast fanpy (uses Michelle's fanci)
                    old_fanpy=False,
                )

            # edit script to allow loading parameters of another CC wavefunction
            if checkpoint_path and checkpoint_wfn:
                edit_cc.edit_file(
                    f'{dirname_results}/calculate.py', checkpoint_wfn, checkpoint_path
                )

            # submit calculations
            # TODO: different ways of providing path results in different types of
            #       calculations/job submissions
            #run_calc.run_calcs(
            #    dirname_results, time=time, memory=memory, outfile='results.out'
            #)
