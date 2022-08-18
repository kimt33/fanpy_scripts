import re
import glob
import os
import numpy as np
import extract_energy
from fanpy.wfn.cc.ap1rog_generalized import AP1roGSDGeneralized
from fanpy.wfn.cc.ap1rog_spin import AP1roGSDSpin
from fanpy.wfn.cc.apg1ro_d import APG1roD
from fanpy.wfn.cc.apg1ro_sd import APG1roSD
from fanpy.wfn.cc.apset1rog_d import APset1roGD
from fanpy.wfn.cc.apset1rog_sd import APset1roGSD
from fanpy.wfn.cc.base import BaseCC
from fanpy.wfn.cc.ccsdq_sen0 import CCSDQsen0
from fanpy.wfn.cc.ccsd_sen0 import CCSDsen0
from fanpy.wfn.cc.ccsdtq_sen0 import CCSDTQsen0
from fanpy.wfn.cc.ccsdt_sen2_q_sen0 import CCSDTsen2Qsen0
from fanpy.wfn.cc.ccs import CCS
from fanpy.wfn.cc.generalized_cc import GeneralizedCC
from fanpy.wfn.cc.pccd_ap1rog import PCCD
from fanpy.wfn.cc.seniority_cc import SeniorityCC
from fanpy.wfn.cc.standard_cc import StandardCC


dict_str_cls = {
    "standardcc": StandardCC,
    "generalizedcc": GeneralizedCC,
    "basecc": BaseCC,
    "senioritycc": SeniorityCC,
    "pccd": PCCD,
    "ap1rogsd": AP1roGSDGeneralized,
    "ap1rogsd_spin": AP1roGSDSpin,
    "apset1rogd": APset1roGD,
    "apset1rogsd": APset1roGSD,
    "apg1rod": APG1roD,
    "apg1rosd": APG1roSD,
    "ccsdsen0": CCSDsen0,
    "ccsdqsen0": CCSDQsen0,
    "ccsdtqsen0": CCSDTQsen0,
    "ccsdtsen2qsen0": CCSDTsen2Qsen0,
    "ccs": CCS,
    "ccsd": StandardCC,
    "ccsdt": StandardCC,
    "ccsdtq": StandardCC,
}

dict_str_import = {
    "standardcc": "from fanpy.wfn.cc.standard_cc import StandardCC",
    "generalizedcc": "from fanpy.wfn.cc.generalized_cc import GeneralizedCC",
    "basecc": "from fanpy.wfn.cc.base import BaseCC",
    "senioritycc": "from fanpy.wfn.cc.seniority_cc import SeniorityCC",
    "pccd": "from fanpy.wfn.cc.pccd_ap1rog import PCCD",
    "ap1rogsd": "from fanpy.wfn.cc.ap1rog_generalized import AP1roGSDGeneralized",
    "ap1rogsd_spin": "from fanpy.wfn.cc.ap1rog_spin import AP1roGSDSpin",
    "apset1rogd": "from fanpy.wfn.cc.apset1rog_d import APset1roGD",
    "apset1rogsd": "from fanpy.wfn.cc.apset1rog_sd import APset1roGSD",
    "apg1rod": "from fanpy.wfn.cc.apg1ro_d import APG1roD",
    "apg1rosd": "from fanpy.wfn.cc.apg1ro_sd import APG1roSD",
    "ccsdqsen0": "from fanpy.wfn.cc.ccsdq_sen0 import CCSDQsen0",
    "ccsdsen0": "from fanpy.wfn.cc.ccsd_sen0 import CCSDsen0",
    "ccsdtqsen0": "from fanpy.wfn.cc.ccsdtq_sen0 import CCSDTQsen0",
    "ccsdtsen2qsen0": "from fanpy.wfn.cc.ccsdt_sen2_q_sen0 import CCSDTsen2Qsen0",
    "ccs": "from fanpy.wfn.cc.ccs import CCS",
}

def edit_file(filepattern, cc_wfn_str, cc_chk_pattern):
    """Edit CC script.

    1. Loads parameters from existing CC wavefunction into current (inside script that matches file 
    pattern CC wavefunction.

    Parameters
    ----------
    filepattern : str
        Glob pattern used to find scripts to edit
    cc_wfn_str : str
        CC wavefunction identifier that will load into CC wavefunction inside script.
    cc_chkpattern : str
        Relative path (wrt filepattern) to CC checkpoint file (for cc_wfn_str) from which parameters
        will be loaded.

    """
    if cc_wfn_str not in dict_str_cls:
        raise ValueError(f"Given cc wavefunction string must be one of {list(dict_str_cls.keys())}")

    cwd = os.getcwd()
    print(filepattern)
    for filename in glob.glob(filepattern):
        dirname = os.path.split(os.path.abspath(filename))[0]
        _, fcienergy_fanpy = extract_energy.extract_min_energy('{dirname}/../../fci/*/results.out')
        _, hfenergy_fanpy = extract_energy.extract_min_energy('{dirname}/../../fci/*/results.out')
        _, fcienergy_gauss = extract_energy.extract_min_energy('{dirname}/../../fci/*/*.log')
        _, hfenergy_gauss = extract_energy.extract_min_energy('{dirname}/../../fci/*/*.log')
        fcienergy = extract_energy.modmin(fcienergy_fanpy, fcienergy_gauss)
        hfenergy = extract_energy.modmin(hfenergy_fanpy, hfenergy_gauss)

        # FIXME: move out as standalone function
        # find checkpoint file from pattern with lowest energy
        chk_energy = {}
        # NOTE: make sure to glob for files instead of directories
        # NOTE: we will assume that reference wavefunction stores its results in the same pattern as
        #       the current pattern
        full_cc_chk_pattern = os.path.join(
            os.path.split(filepattern)[0],  # path to scripts
            os.path.split(cc_chk_pattern)[0],  # relative location of checkpoint file
            "results.out",  # energy output file
        )
        print(f"Glob pattern to find checkpoint files: {full_cc_chk_pattern}")
        for results_path in glob.glob(full_cc_chk_pattern):
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
                chk_energy[
                    os.path.join(
                        os.path.split(results_path)[0],
                        os.path.split(cc_chk_pattern)[1],
                    )
                ] = float(energy)
            else:
                print(f"Didn't get energy for directory, {results_path}")
        if not chk_energy:
            print(f"No energy in path {full_cc_chk_pattern}")
            continue
        # checkpoint file with lowest energy
        checkpoint_file = os.path.join(cwd, min(chk_energy, key=lambda i: chk_energy[i]))
       
        os.chdir(dirname)
        with open("calculate.py", 'r') as f:
            text = f.read()

        # define ranks for ccsd, ccsdt, ccsdtq
        if cc_wfn_str == "ccsd":
            ranks_str = ", ranks=[1, 2]"
        elif cc_wfn_str == "ccsdt":
            ranks_str = ", ranks=[1, 2, 3]"
        elif cc_wfn_str == "ccsdtq":
            ranks_str = ", ranks=[1, 2, 3, 4]"
        else:
            ranks_str = ""

        # load given wavefunction from checkpoint
        wfn_class_name = dict_str_cls[cc_wfn_str].__name__
        wfn_class_import = dict_str_import[cc_wfn_str]
        text = text.replace(
            "\nwfn =",
            f"""\n{wfn_class_import}
{cc_wfn_str} = {wfn_class_name}(nelec, nspin, params=None{ranks_str})
{cc_wfn_str}.assign_params(np.load("{checkpoint_file}"))
wfn ="""
        )
        # loads the given cc wavefunction into the one inside script
        if "wfn.assign_params(" in text:
            text = text.replace(
                "wfn.assign_params(",
                f"""wfn.assign_params({cc_wfn_str})
wfn.assign_params(""",
            )
        else:
            text = text.replace(
                "print('Wavefunction: ",
                f"""wfn.assign_params({cc_wfn_str})
print('Wavefunction: """,
            )

        with open('calculate.py', "w") as f:
            f.write(text)
        os.chdir(cwd)

