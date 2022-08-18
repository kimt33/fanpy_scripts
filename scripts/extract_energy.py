import re
import glob
import os
import numpy as np

def extract_energy(line, cma=False):
    """Return energy inside given string if it exists.

    Supported formats:
    - fanpy output, optimized energy
    - cma optimization output
    - fanpy output, mid optimization energy
    - fanpy output, fci energy
    - fanpt output, after fanpt, after optimization
    - gaussian log, fci
    - gaussian log, correlated
    - gaussian log, hf

    """
    re_energy = re.search("^Final Electronic Energy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    if cma:
        re_energy = re.search("^\s+\d+\s+\d+\s+([-\d\.\+e]+)\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+", line)
        if re_energy:
            return re_energy.group(1)
        re_energy = re.search(
            "^\s+\d+\s+\d+\s+([-\d\.\+e]+)\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+\d+:\d+\.?\d*",
            line,
        )
        if re_energy:
            return re_energy.group(1)

    re_energy = re.search("^\(Mid Optimization\) Electronic Energy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Electronic FCI: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Energy after FanPT: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Energy after optimization: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search(r"Final Eigenvalue\s*([-\d\.\+D]+)", line)
    if re_energy:
        return re_energy.group(1).replace("D", "e")
 
    re_energy = re.search(r"E\(Corr\)=\s*([-\d\.\+D]+)", line)
    if re_energy:
        return re_energy.group(1).replace("D", "e")

    re_energy = re.search(r"E\(\w?HF\)=\s*([-\d\.\+D]+)", line)
    if re_energy:
        return re_energy.group(1).replace("D", "e")


def extract_hf(line):
    """Return energy inside string if it pertains to HF calculation inside fanpy output."""
    re_energy = re.search("^Electronic HF: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)


def get_results(filename):
    """Returns last energy computed from filename.

    Supported formats:
    - fanpy output, optimized energy
    - cma optimization output
    - fanpy output, mid optimization energy
    - fanpy output, fci energy

    """
    system = filename.split(os.path.sep)[0]
    with open(filename, 'r') as g:
        lines = g.readlines()
    energy = None
    for line in lines:
        temp_energy = extract_energy(line)
        if temp_energy:
            energy = temp_energy
    if energy is not None:
        return float(energy)

# TODO: consolidate hf energy with the others?
def get_results_hf(filename):
    """Return HF energy inside filename."""
    system = filename.split(os.path.sep)[0]
    with open(filename, 'r') as g:
        lines = g.readlines()
    energy = None
    for line in lines:
        temp_energy = extract_hf(line)
        if temp_energy:
            energy = temp_energy
    if energy is not None:
        return float(energy)


def modmin(a, b):
    """Returns min of a and b if it exists. If any of the two numbers is None, returns None."""
    if a is None and b is None:
        return None
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)


def extract_min_energy(filepattern):
    """Returns the minimum energy and the corresponding file that matches the given file pattern."""
    min_energy = None
    min_filename = None
    for filename in glob.glob(filepattern):
        energy = get_results(filename)
        if energy is None:
            continue
        # if min_energy is None and energy is not None
        if min_energy is None:
            min_energy = energy
            min_filename = filename
        # if min_energy is not None and energy < min_energy
        elif energy < min_energy:
            min_energy = energy
            min_filename = filename

    return min_filename, min_energy
