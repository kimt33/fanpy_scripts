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

# "ci_pairs": CI Pairs
# "cisd": CISD
# "fci": FCI
# "doci": DOCI
# "mps": MPS
# "determinant-ratio": quasiparticle determinant ratio
# "ap1rog": AP1roG
# "apr2g": APr2G
# "apig": APIG
# "apsetg": APsetG
# "apg": APG
# "network": Feed Forward  (Product Sum) Neural Network
# "rbm": RBM
wfn = "apg"

# order of fanpt
order = 1

# fanpt hyperparameter
# number of steps between lambda=0 and 1
nsteps = 100
# singular value tolerance for lstsq and pinv (see rtol in np.linalg.lstsq or np.linalg.pinv)
# NOTE: you have to be careful here because as the tolerance for optimizing the projected 
#       schrodinger equation loosens, this rtol must also looosen (b/c projected schrodinger
#       equation is no longer "exact")
singular_rtol = 1e-3
# solver to optimize projected schrodigner equation in between fanpt iterations
# to not optimize, use ""
# NOTE: for now, just use "least_squares" or ""
fanpt_solver = "least_squares"
#fanpt_solver = ""

# wfn_dirname is used to find the orbital optimized wavefunction results
# if orbital_optimization is False, then wfn_dirname = wfn (probably)
wfn_dirname = f"apg_fanpt_order{order}"
# orbital optimization is not supported
#orbital_optimization = False

# database structure:
# we will assume that the database/directories are structured as follows:
# db_directory/system/basis/wfn/repetitions/results.out
db_directory = "h8_stretch"

# system information
nelec = 8
basis = "sto-6g"
units = "Angstrom"
charge = 0
multiplicity = 1

# initial guess parameters
wfn_noise = 0.0
ham_noise = 0.0  # not used if no orbital optimization

cwd = os.getcwd()

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
# NOTE: optimizing the energy in between fanpt iterations was not found to be bad
#run_calc.make_wfn_dirs(f'{db_directory}/*/{basis}/', f"{wfn_dirname}", 10, rep_dirname_prefix="energy") 
run_calc.make_wfn_dirs(f'{db_directory}/*/{basis}/', f"{wfn_dirname}", 1, rep_dirname_prefix="projected") 

# 6. make cc script (fanpy)
# define types of calculations to run in each directory label in "repetitions"
# i.e. following runs a energy minimization calculation in the directory energy_min
# key = name of the directory to be created/to store results
# value = 2-tuple of the objective and the optimization algorithm
# TODO: list supported objectives and algorithms

# NOTE: optimizing the energy in between fanpt iterations was not found to be bad
# the optimization algorithm here will be used to optimize the wavefunction in between fanpt 
# iterations
# NOTE: to not have any optimization in between, set solver_kwargs=""
calc_type = {"projected_0": ("projected", "fanpt")}
for path_full in sorted(glob.glob(f'{db_directory}/*/{basis}/{wfn_dirname}'), reverse=True):
    print(path_full)
    for dirpattern_results, (objective, solver) in calc_type.items():
        for dirname_results in glob.glob(f'{path_full}/{dirpattern_results}'):
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
                fanpt=True,
                fanpt_kwargs=(
                    f"fill=fill, energy_active=True, resum=False, ref_sd=0, final_order={order}, "
                    f"lambda_i=0.0, lambda_f=1.0, steps={nsteps}, rcond={singular_rtol}"
                ),
# NOTE: you have to be careful here because as the tolerance for optimizing the projected 
#       schrodinger equation loosens, this rtol must also loosen (b/c projected schrodinger
#       equation is no longer "exact")
                fanpt_solver=fanpt_solver,
                fanpt_solver_kwargs = (
                    "'xtol':5.0e-4, 'ftol':1.0e-4, 'gtol':5.0e-4, 'max_nfev':fanci_wfn.nactive, "
                    "'verbose':2"
                ),
                # use fast fanpy (uses Michelle's fanci)
                old_fanpy=False,
            )

            # submit calculations
            # TODO: different ways of providing path results in different types of
            #       calculations/job submissions
            #run_calc.run_calcs(
            #    dirname_results, time=time, memory=memory, outfile='results.out'
            #)
