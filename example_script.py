import sys
sys.path.append("/home/kimt1/")

import os
import json
import run_calc
import run_fci_g09
import read_calc
import glob
import edit_calculate


cwd = os.getcwd()
basis = 'sto-6g'

## make com
#run_calc.write_coms('h8_octagon/*/{basis}/', memory='1gb', charge=0, multiplicity=1, units="AU")
#run_calc.write_coms('h8_stretch/*/{basis}/', memory='1gb', charge=0, multiplicity=1, units="Angstrom")
#
### run/submit com
#run_calc.run_calcs('h8_*/*/{basis}/hf/hf_sp.com')
#
## make fchk
#run_calc.run_calcs('h8_*/*/{basis}/hf/hf_sp.chk')
#
## make integrals (npy)
#run_calc.run_calcs('h8_*/*/{basis}/hf/hf_sp.fchk')
#
# make directories for wavefunction
#run_calc.make_wfn_dirs(f'h8*/*/{basis}/', 'oo-ap1rog', 10)
#run_calc.make_wfn_dirs(f'h8*/*/{basis}/', 'oo-apig', 10)
#run_calc.make_wfn_dirs(f'h8*/*/{basis}/', 'oo-doci', 10)
#run_calc.make_wfn_dirs(f'h8*/*/{basis}/', 'oo-apg', 10)
#run_calc.make_wfn_dirs(f'h8*/*/{basis}/', 'oo-pccd', 10)
#run_calc.make_wfn_dirs(f'h8*/{basis}/*/', 'fci', 10)
#run_calc.make_wfn_dirs(f'h8*/{basis}/*/', 'ccsd', 10)
#run_calc.make_wfn_dirs(f'h8*/{basis}/*/', 'ccsdt', 10)
#run_calc.make_wfn_dirs(f'h8*/{basis}/*/', 'ccsdtq', 10)

# make python script for optimizing wavefunction
# 0=cma, 1=bfgs, 2=energy+system+least_squares, 3=energy+system+trustregion
#calc_type = {0: ('one_energy', 'cma'), 1: ('one_energy', 'minimize'), 2: ('projected', 'least_squares'), 3: ('projected', 'trustregion')}
#calc_type = {6: ('one_energy', 'minimize')}
#for dir_name, (objective, solver) in calc_type.items():
#    run_calc.write_wfn_py(f'h8*/*/{basis}/oo-ap1rog/{dir_name}', 8, 'ap1rog', 
#                           optimize_orbs=True, pspace_exc=[1, 2,  3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                           ham_noise=1e-3, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
#    run_calc.write_wfn_py(f'h8*/*/sto-6g/oo-apig/{dir_name}', 8, 'apig', 
#                           optimize_orbs=True, pspace_exc=[1, 2,  3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
#    run_calc.write_wfn_py(f'h8*/*/sto-6g/oo-apg/{dir_name}', 8, 'apg', 
#                           optimize_orbs=True, pspace_exc=[1, 2,  3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
#    run_calc.write_wfn_py(f'h8*/*/sto-6g/oo-doci/{dir_name}', 8, 'doci', 
#                           optimize_orbs=True, pspace_exc=[1, 2,  3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
#    run_calc.write_wfn_py(f'h8*/*/sto-6g/oo-pccd/{dir_name}', 8, 'pccd', 
#                           optimize_orbs=True, pspace_exc=[1, 2,  3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
#    if dir_name in [2, 3]:
#        #edit_calculate.edit_file(f'h8*/*/sto-6g/oo-ap1rog/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#        #edit_calculate.edit_file(f'h8*/*/sto-6g/oo-apig/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#        #edit_calculate.edit_file(f'h8*/*/sto-6g/oo-apg/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#        edit_calculate.edit_file(f'h8*/*/sto-6g/oo-doci/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#        edit_calculate.edit_file(f'h8*/*/sto-6g/oo-pccd/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=True)
#
#    # ccsdtq, fci memory
#    run_calc.run_calcs(f'h8*/*/sto-6g/oo-ap1rog/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    run_calc.run_calcs(f'h8*/*/sto-6g/oo-apig/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    run_calc.run_calcs(f'h8*/*/sto-6g/oo-apg/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    run_calc.run_calcs(f'h8*/*/sto-6g/oo-doci/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    run_calc.run_calcs(f'h8*/*/sto-6g/oo-pccd/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#


#for basis in ['cc-pvdz', '6-31g**']:
#    run_calc.write_wfn_py(f'h8*/*/{basis}/ap1rog/', 8, 'ap1rog', 
#                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='projected', solver='trustregion',
#                           ham_noise=1e-3, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None)
#    run_calc.write_wfn_py(f'h8*/*/{basis}/apig/', 8, 'apig', 
#                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='projected', solver='trustregion',
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None)
#    run_calc.write_wfn_py(f'h8*/*/{basis}/apg/', 8, 'apg', 
#                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='projected', solver='trustregion',
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None)
#    run_calc.write_wfn_py(f'h8*/*/{basis}/doci/', 8, 'doci', 
#                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='projected', solver='trustregion',
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None)
#    run_calc.write_wfn_py(f'h8*/*/{basis}/pccd/', 8, 'pccd', 
#                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='projected', solver='trustregion',
#                           ham_noise=1e-1, wfn_noise=1e-2, memory='6gb',
#                           load_orbs=None, load_ham=None, load_wfn=None)

# run/submit calculations
#run_calc.run_calcs('h8*/*/sto-6g/ap1rog/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/apig/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/apg/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/doci/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/fci/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/pccd/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/ccsd/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/ccsdt/calculate.py', time='7d', memory='1gb', outfile='results.out')
#run_calc.run_calcs('h8*/*/sto-6g/ccsdtq/calculate.py', time='7d', memory='1gb', outfile='results.out')


# run calc using fanci
# 0=cma, 1=bfgs, 2=energy+system+least_squares, 3=energy+system+trustregion
# 5=cma with fanci
#calc_type = {6: ('energy', 'minimize')}
#for dir_name, (objective, solver) in calc_type.items():
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/ap1rog/{dir_name}', 8, 'ap1rog', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4],
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/apig/{dir_name}', 8, 'apig', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4],
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/apg/{dir_name}', 8, 'apg', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4],
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/doci/{dir_name}', 8, 'doci', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4, 5, 6, 7, 8],
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/fci/{dir_name}', 8, 'fci', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4, 5, 6, 7, 8],
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/ccsd/{dir_name}', 8, 'ccsd', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3] if dir_name >= 2 else [1, 2, 3, 4],
#    #                       wfn_kwargs = 'indices=None, refwfn=None, exop_combinations=None, refresh_exops=50000',
#    #                       wfn_noise=1e-4, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/pccd/{dir_name}', 8, 'pccd', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1,
#    #                       wfn_noise=1e-2, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/ccsdt/{dir_name}', 8, 'ccsdt', 
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4] if dir_name >= 2 else [1, 2, 3, 4],
#    #                       wfn_kwargs = 'indices=None, refwfn=None, exop_combinations=None, refresh_exops=50000',
#    #                       wfn_noise=1e-4, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #run_calc.write_wfn_py(f'h8*/*/{basis}/ccsdtq/{dir_name}', 8, 'ccsdtq',
#    #                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1, #pspace_exc=[1, 2, 3, 4, 5] if dir_name >= 2 else [1, 2, 3, 4],
#    #                       wfn_kwargs = 'indices=None, refwfn=None, exop_combinations=None, refresh_exops=50000',
#    #                       wfn_noise=1e-4, memory='6gb',
#    #                       load_orbs=None, load_ham=None, load_wfn=None)
#    #if dir_name in [7, 8]:
#    #    #edit_calculate.edit_file(f'h8*/*/{basis}/ap1rog/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    #edit_calculate.edit_file(f'h8*/*/{basis}/apig/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    #edit_calculate.edit_file(f'h8*/*/{basis}/apg/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/doci/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/fci/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/ccsd/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/ccsdt/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False, cc_chk=True)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/ccsdtq/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False, cc_chk=True)
#    #    edit_calculate.edit_file(f'h8*/*/{basis}/ccsdt/{dir_name}/calculate.py', truncate_projection=True, proj_seniority=False, cc_chk=True)
#    #edit_calculate.edit_file(f'h8*/*/{basis}/ccsdt/{dir_name}/calculate.py', cc_chk=True)
#    #edit_calculate.edit_file(f'h8*/*/{basis}/ccsdtq/{dir_name}/calculate.py', cc_chk=True)
#
#    #run_calc.run_calcs(f'h8*/*/{basis}/ap1rog/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/apig/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/apg/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/doci/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/fci/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/pccd/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/ccsd/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/ccsdt/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
#    #run_calc.run_calcs(f'h8*/*/{basis}/ccsdtq/{dir_name}/', time='7d', memory='7gb', outfile='results.out')




# RBM
# find smallest energy, edit
#for wfn in ['ap1rog', 'apig', 'doci', 'pccd', 'apg', 'oo-ap1rog', 'oo-apig', 'oo-doci', 'oo-pccd', 'oo-apg', 'ccsd', 'ccsdt', 'ccsdtq']:
#for wfn in ["ap1rogsd", "ap1rogsd_spin", "apsetgd", "apsetgsd", "apg1rod", "apg1rosd", "ccsdsen0", "ccsdqsen0", "ccsdtqsen0", "ccsdtsen2qsen0"]:
#    run_calc.make_wfn_dirs(f'h8*/*/{basis}/', f'{wfn}_rbm_1,2_1', 10)
#    for dirname in sorted(glob.glob(f'h8*/*/{basis}/{wfn}_rbm_1,2_1'), reverse=True):
#        if 'oo-' in wfn:
#            calc_type = {0: ('one_energy', 'minimize')}
#        else:
#            calc_type = {0: ('energy', 'minimize')}
#        for ind, (objective, solver) in calc_type.items():
#            # find smallest checkpoint with the energy
#            chk = {}
#            for i in glob.glob(f'{cwd}/{dirname}/../{wfn}/*/results.out'):
#                with open(i, 'r') as g:
#                    lines = g.readlines()
#                try:
#                    energy = None
#                    for line in lines:
#                        #if 'cma solver' in line:
#                        #    break
#                        temp_energy = edit_calculate.extract_energy(line)
#                        if temp_energy:
#                            energy = temp_energy
#                    if energy is not None:
#                        chk[i[:-11] + 'checkpoint.npy'] = float(energy)
#                    else:
#                        print("Didn't get energy")
#                except AttributeError:
#                    pass
#            checkpoint = min(chk, key=lambda i: chk[i])
#
#            if 'oo-' in wfn:
#                run_calc.write_wfn_py(f'{dirname}/{ind}', 8, wfn[3:], 
#                                       optimize_orbs=True, pspace_exc=[1, 2, 3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                                       ham_noise=0, wfn_noise=0, memory='6gb',
#                                       load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint, old_fanpy=True)
#            else:
#                run_calc.write_wfn_py(f'{dirname}/{ind}', 8, wfn, 
#                                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1,
#                                       wfn_noise=0, memory='6gb',
#                                       load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint)
#            edit_calculate.edit_file(
#                f'{dirname}/{ind}/calculate.py', truncate_projection=False, proj_seniority=False,
#                rbm=True, rbm_kwargs='num_layers=1, orders=(1, 2)',
#            )
#            run_calc.run_calcs(f'{dirname}/{ind}/', time='7d', memory='7gb', outfile='results.out')




## product-NN
## find smallest energy, edit
#for wfn in ['ap1rog', 'apig', 'doci', 'pccd', 'apg', 'oo-ap1rog', 'oo-apig', 'oo-doci', 'oo-pccd', 'oo-apg', 'ccsd', 'ccsdt', 'ccsdtq']:
#
#    run_calc.make_wfn_dirs(f'h8*/*/{basis}/', f'{wfn}_nn_1_2', 10)
#    for dirname in sorted(glob.glob(f'h8*/*/{basis}/{wfn}_nn_1_2'), reverse=True):
#        if 'oo-' in wfn:
#            calc_type = {0: ('one_energy', 'minimize')}
#        else:
#            calc_type = {0: ('energy', 'minimize')}
#        for ind, (objective, solver) in calc_type.items():
#            # find smallest checkpoint with the energy
#            chk = {}
#            for i in glob.glob(f'{cwd}/{dirname}/../{wfn}/*/results.out'):
#                with open(i, 'r') as g:
#                    lines = g.readlines()
#                try:
#                    energy = None
#                    for line in lines:
#                        #if 'cma solver' in line:
#                        #    break
#                        temp_energy = edit_calculate.extract_energy(line)
#                        if temp_energy:
#                            energy = temp_energy
#                    if energy is not None:
#                        chk[i[:-11] + 'checkpoint.npy'] = float(energy)
#                    else:
#                        print("Didn't get energy")
#                except AttributeError:
#                    pass
#            checkpoint = min(chk, key=lambda i: chk[i])
#
#            if 'oo-' in wfn:
#                run_calc.write_wfn_py(f'{dirname}/{ind}', 8, wfn[3:], 
#                                       optimize_orbs=True, pspace_exc=[1, 2, 3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
#                                       ham_noise=0, wfn_noise=0, memory='6gb',
#                                       load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint, old_fanpy=True)
#            else:
#                run_calc.write_wfn_py(f'{dirname}/{ind}', 8, wfn, 
#                                       optimize_orbs=False, objective=objective, solver=solver, nproj=-1,
#                                       wfn_noise=0, memory='6gb',
#                                       load_orbs=None, load_ham=None, load_wfn=None, load_chk=checkpoint)
#            edit_calculate.edit_file(
#                f'{dirname}/{ind}/calculate.py', truncate_projection=False, proj_seniority=False,
#                rbm_nn=True, rbm_nn_kwargs='num_layers=2',
#            )
#            run_calc.run_calcs(f'{dirname}/{ind}/', time='7d', memory='7gb', outfile='results.out')




#calc_type = {6: ('energy', 'minimize')}
##for wfn_type in ["ap1rogsd", "ap1rogsd_spin", "apsetgd", "apsetgsd", "apg1rod", "apg1rosd"]:
##for wfn_type in ["ccsdsen0", "ccsdqsen0", "ccsdtqsen0", "ccsdtsen2qsen0"]:
#for wfn_type in ["apsetgd", "apg1rod"]:
#    run_calc.make_wfn_dirs(f'h8*/*/{basis}/', f'{wfn_type}', 10)
#    for dir_name, (objective, solver) in calc_type.items():
#        run_calc.write_wfn_py(f'h8*/*/{basis}/{wfn_type}/{dir_name}', 8, wfn_type,
#                               optimize_orbs=False, objective=objective, solver=solver, nproj=-1,
#                               # wfn_kwargs = "ranks=None, indices=None, refwfn=None, exop_combinations=None, refresh_exops=None",
#                               wfn_kwargs = None,
#                               wfn_noise=1e-4, memory='6gb',
#                               load_orbs=None, load_ham=None, load_wfn=None)
#        run_calc.run_calcs(f'h8*/*/{basis}/{wfn_type}/{dir_name}/', time='7d', memory='7gb', outfile='results.out')

calc_type = {0: ('one_energy', 'minimize')}
for wfn_type in ["apsetgd", "apg1rod"]:
    run_calc.make_wfn_dirs(f'h8*/*/{basis}/', f'oo-{wfn_type}', 10)
    for dir_name, (objective, solver) in calc_type.items():
        run_calc.write_wfn_py(f'h8*/*/{basis}/oo-{wfn_type}/{dir_name}', 8, wfn_type,
                               optimize_orbs=True, pspace_exc=[1, 2, 3, 4, 5, 6, 7, 8], objective=objective, solver=solver,
                               wfn_kwargs = None,
                               wfn_noise=1e-4, memory='6gb',
                               load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
        run_calc.run_calcs(f'h8*/*/{basis}/oo-{wfn_type}/{dir_name}/', time='7d', memory='7gb', outfile='results.out')
