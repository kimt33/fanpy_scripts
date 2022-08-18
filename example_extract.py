from glob import glob
import re
import os
import numpy as np

def extract_energy(line):
    re_energy = re.search("^\s+\d+\s+\d+\s+([-\d\.\+e]+)\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+\d+:\d+\.?\d*", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^\(Mid Optimization\) Electronic [eE]nergy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Electronic FCI: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Final Electronic Energy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("Eigenvalue\s+([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)


def extract_hf_energy(line):
    re_energy = re.search("^Electronic HF: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

def extract_chk(line):
    re_chk = re.search(r"dirname, chk_point_file = os.path.split\(['\"](.+)['\"]\)", line)
    if re_chk:
        return re_chk.group(1)


systems = {
'h10_chain':'Angstrom', 
'h10_pyramid':'Angstrom', 
'h10_ring':'Angstrom', 
'h10_sheet':'Angstrom', 
}
wfns = {
        #'apg' : 'APG', 
        #'apg7': 'APG3.1', 'apg3.2': 'APG3.2', 'apg3.3': 'APG3.3', 'apg3.4': 'APG3.4', 'apg3.6': 'APG3.6',
        #'apg3.5': 'APG3.5', 'apg3.10': 'APG3.10', 'apg3.20': 'APG3.20', 'apg3.40': 'APG3.40', 'apg3.80': 'APG3.80',
        #'apg3.8': 'APG3.8', 'apg3.16': 'APG3.16', 'apg3.24': 'APG3.24', 'apg3.32': 'APG3.32', 'apg3.64': 'APG3.64',
        #'nn_1':'NN1', 'nn_2':'NN2', 'nn_3':'NN3', 'nn_4':'NN4', 'nn_5':'NN5',
        #'ap1rog': 'AP1roG', 'apig': 'APIG', 'doci': 'DOCI', 'pccd': 'PCCD',
        #'ccsd': 'CCSD', 'ccsdt': 'CCSDT', 'ccsdtq': 'CCSDTQ',
        #'oo-ap1rog': 'oo-AP1roG', 'oo-apig': 'oo-APIG', 'oo-doci': 'oo-DOCI', 'oo-pccd': 'oo-PCCD',
        #'ccsd': 'CCSD', 'ccsdt': 'CCSDT', 'ccsdtq': 'CCSDTQ',
'fci' : 'FCI',
'hf' : 'HF',
'ap1rog' : 'AP1roG',
'pccd' : 'PCCD',
'apig' : 'APIG',
'doci' : 'DOCI',
'apg' : 'APG',
'oo-ap1rog' : 'oo-AP1roG',
'oo-pccd' : 'oo-PCCD',
'oo-apig' : 'oo-APIG',
'oo-doci' : 'oo-DOCI',
'oo-apg' : 'oo-APG',
'ccsd' : 'CCSD',
'ccsdt' : 'CCSDT',
'ccsdtq' : 'CCSDTQ',

'ap1rogsd_spin' : 'AP1roGSD_spin',
'ap1rogsd' : 'AP1roGSD',
'apsetgd' : 'APsetGD',
'apsetgsd' : 'APsetGSD',
'apg1rod' : 'APG1roD',
'apg1rosd' : 'APG1roSD',

'ap1rog_rbm_1,2_1' : 'AP1roG_RBM_1,2_1',
'pccd_rbm_1,2_1' : 'PCCD_RBM_1,2_1',
'apig_rbm_1,2_1' : 'APIG_RBM_1,2_1',
'doci_rbm_1,2_1' : 'DOCI_RBM_1,2_1',
'apg_rbm_1,2_1' : 'APG_RBM_1,2_1',
'oo-ap1rog_rbm_1,2_1' : 'oo-AP1roG_RBM_1,2_1',
'oo-pccd_rbm_1,2_1' : 'oo-PCCD_RBM_1,2_1',
'oo-apig_rbm_1,2_1' : 'oo-APIG_RBM_1,2_1',
'oo-doci_rbm_1,2_1' : 'oo-DOCI_RBM_1,2_1',
'oo-apg_rbm_1,2_1' : 'oo-APG_RBM_1,2_1',
'ccsd_rbm_1,2_1' : 'CCSD_RBM_1,2_1',
'ccsdt_rbm_1,2_1' : 'CCSDT_RBM_1,2_1',
'ccsdtq_rbm_1,2_1' : 'CCSDTQ_RBM_1,2_1',

'ap1rogsd_spin_rbm_1,2_1' : 'AP1roGSD_spin_RBM_1,2_1',
'ap1rogsd_rbm_1,2_1' : 'AP1roGSD_RBM_1,2_1',
'apsetgd_rbm_1,2_1' : 'APsetGD_RBM_1,2_1',
'apsetgsd_rbm_1,2_1' : 'APsetGSD_RBM_1,2_1',
'apg1rod_rbm_1,2_1' : 'APG1roD_RBM_1,2_1',
'apg1rosd_rbm_1,2_1' : 'APG1roSD_RBM_1,2_1',

'ap1rog_nn_1_2' : 'AP1roG_NN_1_2',
'pccd_nn_1_2' : 'PCCD_NN_1_2',
'apig_nn_1_2' : 'APIG_NN_1_2',
'doci_nn_1_2' : 'DOCI_NN_1_2',
'apg_nn_1_2' : 'APG_NN_1_2',
'oo-ap1rog_nn_1_2' : 'oo-AP1roG_NN_1_2',
'oo-pccd_nn_1_2' : 'oo-PCCD_NN_1_2',
'oo-apig_nn_1_2' : 'oo-APIG_NN_1_2',
'oo-doci_nn_1_2' : 'oo-DOCI_NN_1_2',
'oo-apg_nn_1_2' : 'oo-APG_NN_1_2',
'ccsd_nn_1_2' : 'CCSD_NN_1_2',
'ccsdt_nn_1_2' : 'CCSDT_NN_1_2',
'ccsdtq_nn_1_2' : 'CCSDTQ_NN_1_2',

'ccsdsen0' : 'CCSD_sen0',
'ccsdqsen0' : 'CCSDQ_sen0',
'ccsdtqsen0' : 'CCSDTQ_sen0',
'ccsdtsen2qsen0' : 'CCSDT_sen2_Q_sen0',
        }
basis_sets = ['sto-6g']
for wfn in wfns:
    for basis in basis_sets:
        for system in systems:
            with open(f"data/{system}_{basis}_{wfn}_results.csv", "w") as f:
                f.write(
                    "File,Distance ({}),Nuclear Repulsion Energy (Hartree),{} Electronic Energy "
                    "(Hartree), {} Total Energy (Hartree),Counter\n".format(systems[system], wfns[wfn], wfns[wfn])
                )
                print(f"{system}_*/{basis}/{wfn}/*/results*.out")
                for i in glob(f"{system}_*/{basis}/{wfn}/*/results*.out"):
                    temp = os.path.split(i)
                    if os.path.isfile(
                       os.path.join(temp[0], temp[1].replace('results', 'calculate').replace('out', 'py'))
                    ):
                        calculate_file = os.path.join(temp[0], temp[1].replace('results', 'calculate').replace('out', 'py'))
                        with open(os.path.join(temp[0], temp[1].replace('results', 'calculate').replace('out', 'py')), 'r') as g:
                            if 'fixed' in g.read():
                                continue
                    elif os.path.isfile(
                        os.path.join(temp[0], '../', temp[1].replace('results', 'calculate').replace('out', 'py'))
                    ):
                        calculate_file = os.path.join(temp[0], '../', temp[1].replace('results', 'calculate').replace('out', 'py'))
                        with open(os.path.join(temp[0], '../', temp[1].replace('results', 'calculate').replace('out', 'py')), 'r') as g:
                            if 'fixed' in g.read():
                                continue
                    else:
                        continue
                    #print(i)
                    f.write(i.replace(',', '-'))
                    f.write(",")
                    #with open(os.path.join(os.path.split(i)[0], '../../../system.xyz2'), 'r') as g:
                    #    f.write(re.search("H ([-\d\.]+) ", g.readlines()[1]).group(1))
                    f.write(re.search(f'{system}_([\d\.]+)', i).group(1))
                    f.write(",")
                    nuc_energy = np.load(os.path.join(os.path.split(i)[0], '../../hf/hf_energies.npy'))[1]
                    f.write(str(nuc_energy))
                    f.write(",")
                    with open(i, 'r') as g:
                        lines = g.readlines()
                    try:
                        energy = None
                        counter = 0
                        for line in lines:
                            #if 'cma solver' in line:
                            #    break
                            temp_energy = extract_energy(line)
                            if temp_energy:
                                if energy and abs(float(energy) - float(temp_energy)) > 1e-6: 
                                    counter += 1
                                energy = temp_energy
                        if energy is not None:
                            f.write(energy)
                            f.write(f',{float(energy)+nuc_energy}')
                            f.write(f',{counter}')
                        else:
                            print(f"Didn't get energy: {i}")
                    except AttributeError:
                        pass

                    with open(calculate_file, 'r') as g:
                        lines = g.readlines()
                        counter = 0 
                        for line in lines:
                            chk = extract_chk(line)
                            if chk:
                                f.write(',')
                                f.write(chk)
                    f.write("\n")


for system in systems:
    for basis in basis_sets:
        with open(f"data/{system}_{basis}_fci_results.csv", "w") as f:
            f.write(
                "File,Distance ({}),Nuclear Repulsion Energy (Hartree),FCI Electronic Energy "
                "(Hartree),FCI Total Energy (Hartree)\n".format(systems[system])
            )
            for i in glob(f"{system}_*/{basis}/fci/*/fci_spin1.log"):
                f.write(i)
                f.write(",")
                #with open(os.path.join(os.path.split(i)[0], '../../../system.xyz2'), 'r') as g:
                #    f.write(re.search("H ([-\d\.]+) ", g.readlines()[1]).group(1))
                f.write(re.search(f'{system}_([\d\.]+)', i).group(1))
                f.write(",")
                nuc_energy = np.load(os.path.join(os.path.split(i)[0], '../../hf/hf_energies.npy'))[1]
                f.write(str(nuc_energy))
                f.write(",")
                with open(i, 'r') as g:
                    lines = g.readlines()
                try:
                    energy = None
                    for line in lines:
                        temp_energy = extract_energy(line)
                        if temp_energy:
                            energy = temp_energy
                    if energy is not None:
                        f.write(f'{float(energy)-nuc_energy}')
                        f.write(f',{float(energy)}')
                    else:
                        print(f"Didn't get energy: {i}")
                except AttributeError:
                    pass
                f.write("\n")


for system in systems:
    for basis in basis_sets:
        with open(f"data/{system}_{basis}_hf_results.csv", "w") as f:
            f.write(
                "File,Distance ({}),Nuclear Repulsion Energy (Hartree),HF Electronic Energy "
                "(Hartree),HF Total Energy (Hatree)\n".format(systems[system])
            )
            for i in glob(f"{system}_*/{basis}/hf/hf_sp.log"):
                #print(i)
                f.write(i)
                f.write(",")
                #with open(os.path.join(os.path.split(i)[0], '../../../system.xyz2'), 'r') as g:
                #    f.write(re.search("H ([-\d\.]+) ", g.readlines()[1]).group(1))
                f.write(re.search(f'{system}_([\d\.]+)', i).group(1))
                f.write(",")
                f.write(str(np.load(os.path.join(os.path.split(i)[0], 'hf_energies.npy'))[1]))
                f.write(",")
                f.write(str(np.load(os.path.join(os.path.split(i)[0], 'hf_energies.npy'))[0]))
                f.write(",")
                f.write(str(sum(np.load(os.path.join(os.path.split(i)[0], 'hf_energies.npy')))))
                f.write("\n")


