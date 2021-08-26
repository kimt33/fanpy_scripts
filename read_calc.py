import glob
import re
import os
import fnmatch
import numpy as np


def status(pattern: str):
    """Get the statuses of the calculations that match the given pattern.

    Pattern
    -------
    pattern : str
        Unix shell style wildcard pattern to find the calculations.

    """
    success = []
    opt_failed = []
    code_failed = []
    no_time = []
    running = []
    for filename in glob.glob(pattern):
        with open(filename, 'r') as f:
            results = f.read()

        if re.search('Optimization was successful', results):
            success.append(filename)
        elif re.search('Optimization was not successful: ', results):
            opt_failed.append(filename)
        elif re.search('Traceback (most recent call last):', results):
            code_failed.append(filename)
        elif re.search(r'slurmstepd: error: \*\*\* JOB .+ ON .+ CANCELLED AT .+ DUE TO TIME LIMIT',
                       results):
            no_time.append(filename)
        else:
            running.append(filename)

    return success, opt_failed, code_failed, no_time, running


def extract_results(pattern: str):
    """Get results from the calculations that match the given pattern.

    Pattern
    -------
    pattern : str
        Unix shell style wildcard pattern to find the calculations.

    """
    if pattern[-4:] in ['.out', '.log']:
        filenames, *_, no_time, running = status(pattern)
    elif pattern[-4:] == '.npy':
        filenames = glob.glob(pattern)
        running = []
        no_time = []

    cwd = os.getcwd()
    output = []
    for i, filename in enumerate(filenames + running + no_time):
        if i >= len(filenames):
            is_complete = False
        else:
            is_complete = True

        filename = os.path.abspath(filename)
        if os.path.commonpath([cwd, filename]) != cwd:
            continue
        filename = filename[len(cwd)+1:]
        print(filename)

        split_filename = filename.split(os.sep)
        dict_calc = {}
        if len(split_filename) in [5, 6]:
            with open(filename, 'r') as f:
                results = f.read()
            if results == '':
                continue

            _, system_basis, orbital, wfn, *index, filename = split_filename
            dict_calc['wfn'] = wfn
            dict_calc['index'] = index

            if filename == 'results.out':
                re_nuc = re.search(r'Nuclear-nuclear repulsion: (.+)', results)
                nuc_nuc = re_nuc.group(1)

                if is_complete:
                    re_energy = re.search(r'Final Electronic Energy: (.+)', results)
                    energy = re_energy.group(1)
                else:
                    lastline = re.split(r'\n', results)[-2]
                    if 'Iterat' in lastline:
                        continue
                    _, _, _, energy, _, sigma, *_ = re.split(r'\s+', lastline)
                    dict_calc['sigma'] = sigma
            elif filename == 'cc.log':
                re_nuc = re.search(r'\s+nuclear repulsion energy\s+(\d+\.\d+) Hartrees\.', results)
                nuc_nuc = re_nuc.group(1)
                energies = re.findall(r' E\(corr\)=\s+(-\d+\.\d+) ', results, flags=re.I)
                energy = energies[-1]
            elif filename == 'fci.log':
                print(system_basis)
                re_nuc = re.search(r'\s+nuclear repulsion energy\s+(\d+\.\d+) Hartrees\.', results)
                nuc_nuc = re_nuc.group(1)
                re_energy = re.search(r'EIGENVALUE\s+(-*\d+\.\d+)', results, flags=re.I)
                energy = re_energy.group(1)
        elif len(split_filename) == 4:
            energy, nuc_nuc = np.load(filename)

            _, system_basis, orbital, filename = split_filename
            dict_calc['wfn'] = 'hf'
        else:
            raise NotImplementedError(f'Unsupported file/directory: {filename}')
        system, basis = system_basis.rsplit('_', 1)
        dict_calc['system'] = system
        dict_calc['basis'] = basis
        dict_calc['orbital'] = orbital
        dict_calc['filename'] = filename
        dict_calc['energy'] = energy
        dict_calc['nuc_nuc'] = nuc_nuc

        output.append(dict_calc)

    return output


def select_results(results: dict, system: str, basis: str, orbital: str, wfn: str):
    """Select results that match the given options.

    Pattern
    -------
    results : dict
        Dictionary of the results from `extract_results`
    system : str
        Unix shell style wildcard pattern to select the systems.
    basis : str
        Basis set of the calculations.
    orbital : str
        Orbitals used for the calculation.
    wfn : str
        Wavefunction of the calculation.

    """
    output_x = []
    output_y = []
    output_error = []
    for result in results:
        if not (fnmatch.fnmatch(result['system'], system) and
                result['basis'] == basis and
                result['orbital'] == orbital and
                result['wfn'] == wfn):
            continue
        xval = result['system'].rsplit('_', 1)[1]
        output_x.append(xval)
        output_y.append(float(result['energy']) + float(result['nuc_nuc']))
        try:
            output_error.append(result['sigma'])
        except KeyError:
            output_error.append(0)
    return np.array(output_x, dtype=int), np.array(output_y, dtype=float), np.array(output_error, dtype=float)


def trim(x, y, keep='all'):
    new_x = []
    new_y = []
    for i in range(50):
        if np.all(x != i):
            continue
        # most frequent value
        y_vals = y[x == i]
        # round to 8th decimal
        y_vals = np.around(y_vals, 8)
        # remove repeated numbers
        y_vals_unique = np.unique(y_vals)
        # find frequency
        hist = np.array([np.sum(y_vals == i) for i in y_vals_unique])
        # if all are kept
        if keep == 'all':
            new_x.extend([i] * y_vals_unique.size)
            new_y.extend(y_vals_unique)
        # if only minimm is kept
        elif keep == 'min':
            new_x.append(i)
            new_y.append(min(y_vals_unique))
        # if only the most frequent is kept
        elif keep == 'frequent':
            new_x.append(i)
            new_y.append(y_vals_unique[np.argmax(hist)])

    new_x, new_y = np.array(new_x), np.array(new_y)
    return new_x, new_y
