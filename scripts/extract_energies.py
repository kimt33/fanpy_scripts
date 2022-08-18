import sys
sys.path.append("/home/kimt1/scripts_fanpy_db")

import dataclasses
from glob import glob
import os
import re

import numpy as np
from extract_energy import extract_energy


# we will assume that the database/directories are structured as follows:
# db_directory/system/basis/wfn/repetitions/results.out

# create csv file for each calculation type in each system
# Columns: 
#    File
#    Alpha (units)
#    Nuclear Repulsion Energy (Hartree)
#    Wavefunction Electronic Energy(Hartree)
#    Wavefunction Total Energy (Hartree)
# NOTE: Since this is a csv, "," in the filename will be replaced with "-"

# systems
@dataclasses.dataclass
class System:
    """Description of system.

    Attributes
    ----------
    name : str
        Name of the system. This will be used to name the csv.
    unit : str
        Units used within the system. This will be used to label column in csv.
    alpha_name : str
        Name of the numerical metric of the system (e.g. bond distance, angle, etc). This is used to
        to label column in csv.
    system_pattern : str
        Glob pattern of the system directories. Used to find directories of the system.
    alpha_pattern : str
        Regex pattern to extract the numerical metric of the system (e.g. bond distance, angle, etc)
        This is used to extract numerical metric. Only extract first group.

    """
    name: str
    units: str
    alpha_name: str
    system_pattern: str
    alpha_pattern: str


systems = [
    System(
        name="h8_stretch", 
        units="Angstrom", 
        alpha_name="Distance",
        system_pattern="h8_stretch/h8_*_linear",
        alpha_pattern="h8_stretch/h8_([\d\.]+)_linear",
    )
]

# basis set
basis_sets = ["sto-6g"]

# calculation type
@dataclasses.dataclass
class CalculationType:
    """Description of calculation.

    Attributes
    ----------
    name : str
        Name of the calculation. This will be used to label column in csv.
    wfn_pattern : str
        Glob pattern for the wavefunction directories. Used to find directories of the wavefunction.
    repetition_pattern : str
        Glob pattern of the specific calculation type within wavefunction directories. We call this 
        level "repetitions". Used to find directories of desired results.out.
    output_pattern : str
        Glob pattern of the output file in repetition direcotires. Energies should be extracted
        from these output files.
    wfn_filename : str
        Wavefunction identifier to be used when naming the csv.
        Default value is wfn_pattern.
    identifier_filename : str
        Additional identifier of the calculation to be added when naming the csv.

    """

    name: str
    wfn_pattern: str
    repetition_pattern: str = "*"
    output_pattern: str = "results*.out"
    wfn_filename: str = "" 
    identifier_filename: str = ""
    # NOTE: default wfn_filename is wfn_pattern value b/c of __post_init__

    def __post_init__(self):
        if not self.wfn_filename:
            self.wfn_filename = self.wfn_pattern


calc_types = [
    # standard
    CalculationType(name="FCI", wfn_pattern="fci"),
    CalculationType(name="HF", wfn_pattern="fci", wfn_filename="hf"),
    # CalculationType(name="AP1roG", wfn_pattern="ap1rog"),
    # CalculationType(name="PCCD", wfn_pattern="pccd"),
    # CalculationType(name="APIG", wfn_pattern="apig"),
    # CalculationType(name="DOCI", wfn_pattern="doci"),
     CalculationType(name="APG", wfn_pattern="apg"),

    # orbital optimized
    # CalculationType(name="oo-AP1roG", wfn_pattern="oo-ap1rog"),
    # CalculationType(name="oo-PCCD", wfn_pattern="oo-pccd"),
    # CalculationType(name="oo-APIG", wfn_pattern="oo-apig"),
    # CalculationType(name="oo-DOCI", wfn_pattern="oo-doci"),
    # CalculationType(name="oo-APG", wfn_pattern="oo-apg"),

    # standard cc
    # CalculationType(name="CCSD", wfn_pattern="ccsd"),
    # CalculationType(name="CCSDT", wfn_pattern="ccsdt"),
    # CalculationType(name="CCSDTQ", wfn_pattern="ccsdtq"),

    # cc geminal
    # CalculationType(name="AP1roGSD_spin", wfn_pattern="ap1rogsd_spin"),
    # CalculationType(name="AP1roGSD", wfn_pattern="ap1rogsd"),
    # CalculationType(name="APsetGD", wfn_pattern="apsetgd"),
    # CalculationType(name="APsetGSD", wfn_pattern="apsetgsd"),
    # CalculationType(name="APG1roD", wfn_pattern="apg1rod"),
    # CalculationType(name="APG1roSD", wfn_pattern="apg1rosd"),
    # CalculationType(name="CCSD_sen0", wfn_pattern="ccsdsen0"),
    # CalculationType(name="CCSDQ_sen0", wfn_pattern="ccsdqsen0"),
    # CalculationType(name="CCSDTQ_sen0", wfn_pattern="ccsdtqsen0"),
    # CalculationType(name="CCSDT_sen2_Q_sen0", wfn_pattern="ccsdtsen2qsen0"),

    # rbm
    # CalculationType(name="AP1roG_RBM_1,2_1", wfn_pattern="ap1rog_rbm_1,2_1"),
    # CalculationType(name="PCCD_RBM_1,2_1", wfn_pattern="pccd_rbm_1,2_1"),
    # CalculationType(name="APIG_RBM_1,2_1", wfn_pattern="apig_rbm_1,2_1"),
    # CalculationType(name="DOCI_RBM_1,2_1", wfn_pattern="doci_rbm_1,2_1"),
    # CalculationType(name="APG_RBM_1,2_1", wfn_pattern="apg_rbm_1,2_1"),
    # CalculationType(name="oo-AP1roG_RBM_1,2_1", wfn_pattern="oo-ap1rog_rbm_1,2_1"),
    # CalculationType(name="oo-PCCD_RBM_1,2_1", wfn_pattern="oo-pccd_rbm_1,2_1"),
    # CalculationType(name="oo-APIG_RBM_1,2_1", wfn_pattern="oo-apig_rbm_1,2_1"),
    # CalculationType(name="oo-DOCI_RBM_1,2_1", wfn_pattern="oo-doci_rbm_1,2_1"),
    # CalculationType(name="oo-APG_RBM_1,2_1", wfn_pattern="oo-apg_rbm_1,2_1"),
    # CalculationType(name="CCSD_RBM_1,2_1", wfn_pattern="ccsd_rbm_1,2_1"),
    # CalculationType(name="CCSDT_RBM_1,2_1", wfn_pattern="ccsdt_rbm_1,2_1"),
    # CalculationType(name="CCSDTQ_RBM_1,2_1", wfn_pattern="ccsdtq_rbm_1,2_1"),
    # CalculationType(name="AP1roGSD_spin_RBM_1,2_1", wfn_pattern="ap1rogsd_spin_rbm_1,2_1"),
    # CalculationType(name="AP1roGSD_RBM_1,2_1", wfn_pattern="ap1rogsd_rbm_1,2_1"),
    # CalculationType(name="APsetGD_RBM_1,2_1", wfn_pattern="apsetgd_rbm_1,2_1"),
    # CalculationType(name="APsetGSD_RBM_1,2_1", wfn_pattern="apsetgsd_rbm_1,2_1"),
    # CalculationType(name="APG1roD_RBM_1,2_1", wfn_pattern="apg1rod_rbm_1,2_1"),
    # CalculationType(name="APG1roSD_RBM_1,2_1", wfn_pattern="apg1rosd_rbm_1,2_1"),

    # nn
    # CalculationType(name="AP1roG_NN_1_2", wfn_pattern="ap1rog_nn_1_2"),
    # CalculationType(name="PCCD_NN_1_2", wfn_pattern="pccd_nn_1_2"),
    # CalculationType(name="APIG_NN_1_2", wfn_pattern="apig_nn_1_2"),
    # CalculationType(name="DOCI_NN_1_2", wfn_pattern="doci_nn_1_2"),
    # CalculationType(name="APG_NN_1_2", wfn_pattern="apg_nn_1_2"),
    # CalculationType(name="oo-AP1roG_NN_1_2", wfn_pattern="oo-ap1rog_nn_1_2"),
    # CalculationType(name="oo-PCCD_NN_1_2", wfn_pattern="oo-pccd_nn_1_2"),
    # CalculationType(name="oo-APIG_NN_1_2", wfn_pattern="oo-apig_nn_1_2"),
    # CalculationType(name="oo-DOCI_NN_1_2", wfn_pattern="oo-doci_nn_1_2"),
    # CalculationType(name="oo-APG_NN_1_2", wfn_pattern="oo-apg_nn_1_2"),
    # CalculationType(name="CCSD_NN_1_2", wfn_pattern="ccsd_nn_1_2"),
    # CalculationType(name="CCSDT_NN_1_2", wfn_pattern="ccsdt_nn_1_2"),
    # CalculationType(name="CCSDTQ_NN_1_2", wfn_pattern="ccsdtq_nn_1_2"),

    # fanpt
    CalculationType(
        name="APG_FanPT_1",
        wfn_pattern="apg_fanpt_order1",
        repetition_pattern="projected_0",
        identifier_filename="100steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_2",
        wfn_pattern="apg_fanpt_order2",
        repetition_pattern="projected_0",
        identifier_filename="100steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_3",
        wfn_pattern="apg_fanpt_order3",
        repetition_pattern="projected_0",
        identifier_filename="100steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_4",
        wfn_pattern="apg_fanpt_order4",
        repetition_pattern="projected_0",
        identifier_filename="100steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_1",
        wfn_pattern="apg_fanpt_order1",
        repetition_pattern="projected_10_0",
        identifier_filename="10steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_2",
        wfn_pattern="apg_fanpt_order2",
        repetition_pattern="projected_10_0",
        identifier_filename="10steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_3",
        wfn_pattern="apg_fanpt_order3",
        repetition_pattern="projected_10_0",
        identifier_filename="10steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_4",
        wfn_pattern="apg_fanpt_order4",
        repetition_pattern="projected_10_0",
        identifier_filename="10steps_lstsq",
    ),
    CalculationType(
        name="APG_FanPT_1",
        wfn_pattern="apg_fanpt_order1",
        repetition_pattern="nosolve_10_0",
        identifier_filename="10steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_2",
        wfn_pattern="apg_fanpt_order2",
        repetition_pattern="nosolve_10_0",
        identifier_filename="10steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_3",
        wfn_pattern="apg_fanpt_order3",
        repetition_pattern="nosolve_10_0",
        identifier_filename="10steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_4",
        wfn_pattern="apg_fanpt_order4",
        repetition_pattern="nosolve_10_0",
        identifier_filename="10steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_1",
        wfn_pattern="apg_fanpt_order1",
        repetition_pattern="nosolve_100_0",
        identifier_filename="100steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_2",
        wfn_pattern="apg_fanpt_order2",
        repetition_pattern="nosolve_100_0",
        identifier_filename="100steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_3",
        wfn_pattern="apg_fanpt_order3",
        repetition_pattern="nosolve_100_0",
        identifier_filename="100steps_pure",
    ),
    CalculationType(
        name="APG_FanPT_4",
        wfn_pattern="apg_fanpt_order4",
        repetition_pattern="nosolve_100_0",
        identifier_filename="100steps_pure",
    ),
]


for system in systems:
    for basis in basis_sets:
        for calc_type in calc_types:
            csv_filename = "_".join(
                [
                    system.name,
                    basis,
                    calc_type.wfn_filename,
                    calc_type.identifier_filename,
                    "results.csv",
                ]
            )
            with open(f"data/{csv_filename}", "w") as f:
                f.write(
                    "File,{} ({}),Nuclear Repulsion Energy (Hartree),{} Electronic Energy "
                    "(Hartree), {} Total Energy (Hartree)\n".format(
                        system.alpha_name, system.units, calc_type.name, calc_type.name
                    )
                )
                glob_pattern = "{}/{}/{}/{}/{}".format(
                    system.system_pattern,
                    basis,
                    calc_type.wfn_pattern,
                    calc_type.repetition_pattern,
                    calc_type.output_pattern,
                )
                print(glob_pattern)
                for filepath in glob(glob_pattern):
                    #print(filepath)

                    # NOTE: replace , in filepath to -
                    f.write(filepath.replace(",", "-"))
                    f.write(",")
                    f.write(re.search(system.alpha_pattern, filepath).group(1))
                    f.write(",")
                    nuc_energy = np.load(
                        os.path.join(os.path.split(filepath)[0], "../../hf/hf_energies.npy")
                    )[1]
                    hf_energy = np.load(
                        os.path.join(os.path.split(filepath)[0], "../../hf/hf_energies.npy")
                    )[0]
                    f.write(str(nuc_energy))
                    f.write(",")

                    # extract energy from file
                    with open(filepath, "r") as g:
                        lines = g.readlines()
                    try:
                        energy = None
                        #counter = 0
                        for line in lines:
                            if calc_type.name == "HF":
                                temp_energy = hf_energy
                            else:
                                temp_energy = extract_energy(line)
                            if temp_energy:
                                # NOTE: if you want to check if energy is increasing
                                #if energy and abs(float(energy) - float(temp_energy)) > 1e-6:
                                #    counter += 1
                                energy = temp_energy
                        if energy is not None:
                            f.write(f"{energy},{float(energy)+nuc_energy}")
                            #f.write(f",{counter}")
                        else:
                            print(f"Didn't get energy: {filepath}")
                    except AttributeError:
                        pass
                    f.write("\n")
