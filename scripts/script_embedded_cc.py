import sys
sys.path.append("/home/kimt1/scripts_fanpy_db")

import os
import re
import json
import glob

import run_calc
import extract_energy
import edit_cc


cwd = os.getcwd()

# CC here obtains the overlaps by explicitly evaluating the cumulent
# there are various flavors that ramon has developed

# supported wavefunctions:
# uncomment the you want
# "standardcc": CC by excitation orders
# "generalizedcc": generalized CC  # TODO: HOW TO DISTINGUISH
# "basecc": arbitrary CC  # TODO: HOW TO DISTINGUISH
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
wfns = ["apg1rosd", "ap1rogsd"]  # cc wavefunctions used in each subsystem in order

# NOTE: you do not actually provide number of electrons in each subsystem
#       because orbitals (and their occupations) partitioned "automatically"
# orbitals are localized, then are ordered by lowest to highest energy,
# electrons are assigned by aufbau principle
# NOTE: to have explicit control over subsystems (orbitals and number of electrons), you must edit
#       the produced scripts directly

# labels (index of subsystem) for each atom in system
atom_system_inds = [1, 0, 1, 1, 0, 0, 1, 1, 0 ,1]
# example:
# if [0, 1, 2, 2], 0th atom assigned to system 0
#                  1th atom assigned to system 1
#                  2th atom assigned to system 2
#                  3th atom assigned to system 2

# integral location
oneint_path = "../../hf_pm/oneint.npy"
twoint_path = "../../hf_pm/twoint.npy"
hf_energies_path = "../../hf_pm/hf_energies.npy"

# labels (index of ATOM) for each (spatial) orbital 
ao_inds_file = "../../hf_pm/ao_inds.npy"
# explicit orbital partition that WILL OVERRIDE ao_inds_file partition
ao_inds = None
#ao_inds = [9, 9, 0, 0, 0, 9, 9, 0, 0, 0]
# NOTE: if no localization scheme, then ao inds must be provided explicitly (list of ints)

# example:
# [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
# 0th-3th spatial orbital assigned to atom 0
# 4th-7th spatial orbital assigned to atom 1
# 8th-11th spatial orbital assigned to atom 2
# NOTE: it is not necessary that every atom "has" at least one orbital assigned to it
#       and by extension, it is not necessary that every system "has" at least one orbital assigned 
#       to it

# localization algorithm
# see https://pyscf.org/pyscf_api_docs/pyscf.lo.html
# "svd" - maximize overlap with minimal atomic orbital via svd (procrustes problem, quanmbo ish)
# "iao" - knizia's intrinsic atomic orbital
# "boys" - Boys' localization algorithm
# "pm" - Pipek Mezey algorithm
# "er" - Edmiston-Ruedenberg localization
# None - No localization
loc_type = "pm"
# FIXME: what happens if loc_type == none?

# NOTE: localization algorithms don't seem to be too reliable at the moment so orbital optimization 
#       is strongly recommended
# NOTE: if no localization, you MUST provide the ao_inds

# wfn_dirname is used to find the orbital optimized wavefunction results
# if orbital_optimization is False, then wfn_dirname = wfn (probably)
wfn_dirname = f"cc_embedded_{'_'.join(wfns)}"
# it is recommended to have orbital optimization off (False) if you have single excitations
orbital_optimization = True
if orbital_optimization:
    wfn_dirname += "_singleinterops"

# database structure:
# we will assume that the database/directories are structured as follows:
# db_directory/system/basis/wfn/repetitions/results.out
db_directory = "h10"
system_dirname = "h10_sheet*"

# system information
nelec = 10
basis = "sto-6g"
units = "Angstrom"
charge = 0
multiplicity = 1

# initial guess parameters
wfn_noise = 1e-3

# system/job information
time = "7d"
memory = "10gb"  # total memory available
memory_for_wfn_cache = "9gb"  # memory allocated to storing wavefunction overlaps

# 1. make Gaussian input
# pay special attention to units
#run_calc.write_coms(f'{db_directory}/*/{basis}/', memory='1gb', charge=charge, multiplicity=multiplicity, units=units)

# 2. run/submit Gaussian calculation
# If memory and runtime not given, runs in shell
#run_calc.run_calcs(f'{db_directory}/*/{basis}/hf/hf_sp.com')

# 3. make fchk
#run_calc.run_calcs(f'{db_directory}/*/{basis}/hf/hf_sp.chk')

# 4. make integrals (npy) using fchk and horton
#run_calc.run_calcs(f'{db_directory}/{system_dirname}/{basis}/hf/hf_sp.fchk')

# 5. localize orbitals
#run_calc.make_wfn_dirs(f'{db_directory}/{system_dirname}/{basis}/', f'hf_{loc_type}', 0)
#run_calc.run_calcs(f'{db_directory}/{system_dirname}/{basis}/hf_{loc_type}', arg=atom_system_inds, cubes=True, cube_dim=40)

# 6. make directories for wavefunction
run_calc.make_wfn_dirs(f'{db_directory}/{system_dirname}/{basis}/', f"{wfn_dirname}", 1, rep_dirname_prefix="energy_min_auto") 

# 7. make cc script (fanpy)
# define types of calculations to run in each directory label in "repetitions"
# i.e. following runs a energy minimization calculation in the directory energy_min
# key = name of the directory to be created/to store results
# value = 2-tuple of the objective and the optimization algorithm
# TODO: list supported objectives and algorithms

calc_type = {"energy_min_auto_0": ("energy", "minimize"), "projected": ("projected", "least_squares")}
for path_full in sorted(glob.glob(f'{db_directory}/{system_dirname}/{basis}/{wfn_dirname}'), reverse=True):
    print(path_full)
    for dirpattern_results, (objective, solver) in calc_type.items():
        for dirname_results in glob.glob(f'{path_full}/{dirpattern_results}'):
            for i,wfn in enumerate(wfns):
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
                    filename=f"template_{wfn}_{i}.py",
                    oneint_path=oneint_path, twoint_path=twoint_path, hf_energies_path=hf_energies_path,
                    load_orbs=None, load_ham=None, load_wfn=None,
                    # use fast fanpy (uses Michelle's fanci)
                    old_fanpy=False,
                )

            # create embedded wavefunction
            run_calc.write_embedding_wfn_py(
                dirname_results, 
                [f"template_{wfn}_{i}.py" for i, wfn in enumerate(wfns)],
                atom_system_inds,
                ao_inds_file,
                ao_inds=ao_inds,
                cc=True,
                loc_type=loc_type,
                orbital_optimization=orbital_optimization,
                wfn_noise=wfn_noise,
                filename='calculate.py',
            )

            # submit calculations
            # TODO: different ways of providing path results in different types of
            #       calculations/job submissions
            run_calc.run_calcs(
                dirname_results, time=time, memory=memory, outfile='results.out'
            )
