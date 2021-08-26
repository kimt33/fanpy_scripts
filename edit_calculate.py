import re
import glob
import os
import numpy as np

def extract_energy(line):
    re_energy = re.search("^Final Electronic Energy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^\s+\d+\s+\d+\s+([-\d\.\+e]+)\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^\(Mid Optimization\) Electronic Energy: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)

    re_energy = re.search("^Electronic FCI: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)


def extract_hf(line):
    re_energy = re.search("^Electronic HF: ([-\d\.\+e]+)", line)
    if re_energy:
        return re_energy.group(1)


def get_results(filename):
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

def get_results_hf(filename):
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


from_text = "wfn.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext)))"
to_text = """from fanpy.wfn.geminal.apig import APIG
apig = APIG(nelec, nspin)
apig.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext))[:apig.nparams])
wfn.assign_params(apig)"""


def modmin(a, b):
    if a is None and b is None:
        return None
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)


def edit_file(filepattern, truncate_projection=True, hfupper=True, proj_seniority=False, energy_constraint=True, cma=False, apig_chk=False, doci_chk=False, apg_special=False, nn_chk=0, nn_num=0, apg_gradual=False, rbm=False, rbm_kwargs="num_layers=1, orders=(1, 2)", rbm_nn=False, rbm_nn_kwargs="num_layers=2", cc_chk=None):
    cwd = os.getcwd()
    print(filepattern)
    for filename in glob.glob(filepattern):
        dirname = os.path.split(os.path.abspath(filename))[0]
        fcienergy = 0
        hfenergy = 0
        for filename in glob.glob(f'{dirname}/../../fci/*/results.out'):
            fcienergy = modmin(get_results(filename), fcienergy)
            hfenergy = modmin(get_results_hf(filename), hfenergy)
        for filename in glob.glob(f'{dirname}/../../fci/*/*log'):
            fcienergy = modmin(get_results(filename), fcienergy)
            hfenergy = modmin(get_results_hf(filename), hfenergy)
       
        os.chdir(dirname)
        with open("calculate.py", 'r') as f:
            text = f.read()

        # truncate projection space
        if truncate_projection:
            from_text3 = "num_limit=None"
            to_text3 = "num_limit=wfn.nparams + ham.nparams"
            text = text.replace(from_text3, to_text3)

        if proj_seniority:
            text = text.replace('seniority=wfn.seniority', 'seniority=0')

        # system trust region
        if energy_constraint:
            text = text.replace('-np.inf', str(fcienergy))
            if hfupper:
                text = text.replace('np.inf', str(hfenergy))
            else:
                text = text.replace('np.inf', str(0))
           
        # apig checkpoint
        if apig_chk:
            from_text = "wfn.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext)))"
            to_text = """from fanpy.wfn.geminal.apig import APIG
    apig = APIG(nelec, nspin)
    apig.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext))[:apig.nparams])
    wfn.assign_params(apig)"""
            if apg_special:
                to_text +="""
    wfn.assign_params(wfn.params + 0.005 * 2 * (np.random.rand(*wfn.params.shape) - 0.5))"""
            text = text.replace(from_text, to_text)
           
        # doci checkpoint
        if doci_chk:
            from_text = "wfn.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext)))"
            to_text = """from fanpy.wfn.geminal.apig import APIG
    apig = APIG(nelec, nspin)
    apig.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext))[:apig.nparams])
    for sd in pspace:
        wfn.params[wfn.dict_sd_index[sd]] = apig.get_overlap(sd)
    """
            text = text.replace(from_text, to_text)
           
        # nn1 checkpoint
        if nn_chk:
            from_text = "wfn.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext)))"
            to_text = f"""wfn2 = NumpyNetwork(nelec, nspin, params=None, memory=None, num_layers={nn_chk})
    wfn2.assign_params(np.load(os.path.join(dirname, chk_point_file + '_wfn' + ext)))
    wfn2.normalize(pspace)
    wfn.output_scale = wfn2.output_scale
    if {nn_chk} > 1:
        wfn.assign_params(np.hstack([np.identity(wfn.nspin).flat] + [np.identity(wfn.nspin).flat / np.tanh(1)]*{nn_num-nn_chk-1} + [wfn2.params[:wfn2.nspin**2] / np.tanh(1), wfn2.params[wfn2.nspin**2:]]))
    else:
        wfn.assign_params(np.hstack([np.identity(wfn.nspin).flat, wfn2.params / np.tanh(1)]))"""
            text = text.replace(from_text, to_text)
            from_text = "# Normalize"
            to_text = """objective.assign_params(objective.params + 1e-3 * 2 * (np.random.rand(objective.params.size) - 0.5))
    # Normalize"""
            #text = text.replace(from_text, to_text)

            from_text = """param_selection = [(wfn, np.ones(wfn.nparams, dtype=bool)), (ham, np.ones(ham.nparams, dtype=bool))]"""
            to_text = """wfn_selection = np.zeros(wfn.nparams, dtype=bool)
    wfn_selection[:wfn.nspin**2] = True
    param_selection = [(wfn, wfn_selection), (ham, np.ones(ham.nparams, dtype=bool))]"""
            #text = text.replace(from_text, to_text)

        if apg_gradual:
            from_text = """# Normalize
    wfn.normalize(pspace)"""
            to_text = """for i in range(32, 15, -2):
        print('num_matching:',i)
        wfn = APG3(nelec, nspin, params=wfn.params, memory=None, tol=1e-8, num_matchings=i)
        wfn.normalize(pspace)
        objective = OneSidedEnergy(wfn, ham, param_selection=param_selection, tmpfile='checkpoint.npy',
                                   refwfn=pspace)
        print(objective.objective(objective.params))
        results = bfgs_minimize(objective, save_file='checkpoint.npy', method='BFGS', jac=objective.gradient,
                           options={'gtol': 1e-4, 'disp':True, 'maxiter':10})"""
            text = text.replace(from_text, to_text)

        # cma
        if cma:
            text = text.replace("# Solve", """from fanpy.solver.equation import minimize
    objective.print_energy = True
    ham.update_prev_params = True
    wfn.normalize(pspace)
    results = minimize(objective, save_file='checkpoint.npy', method='BFGS', jac=objective.gradient, options={'gtol': 1e-8, 'disp':True})
    objective.print_energy = False
    ham.update_prev_params = False
    wfn.normalize(pspace)
    # Solve""")
            text = text.replace("'tolfun': 1e-11", "'tolfun': 1e-5")
            #text = text.replace("'tolfun': 1e-11", "'tolfun': 1e-3, 'maxfevals':2000")
            text = text.replace("# Results", """from fanpy.solver.equation import minimize
    from fanpy.tools.math_tools import unitary_matrix
    objective.print_energy = True
    ham.update_prev_params = True
    ham._prev_unitary = ham._prev_unitary.dot(unitary_matrix(ham.params - ham._prev_params))
    ham._prev_params = ham.params.copy()
    ham.assign_params(ham.params)
    wfn.normalize(pspace)
    results = minimize(objective, save_file='checkpoint.npy', method='BFGS', jac=objective.gradient, options={'gtol': 1e-8, 'disp':True})
    # Results""")

        if rbm:
            nproj_txt = re.search(r'\s*(nproj = .+?)\n', text).group()
            if 'wfn' in nproj_txt:
                nproj_txt = nproj_txt.replace('wfn', 'wfn_rbm')
            else:
                nproj_txt = ''
            # find lowest energy checkpoint
            text = text.replace(
                "# Select parameters that will be optimized",
                f"""# Create product of wavefunction for rbm
wfn_orig = wfn
from fanpy.wfn.network.rbm import RestrictedBoltzmannMachine
wfn_rbm = RestrictedBoltzmannMachine(nelec, nspin, nbath=nspin, {rbm_kwargs})
wfn_rbm.assign_params(wfn_rbm.params + 0.001 * 2 * (np.random.rand(*wfn_rbm.params.shape) - 0.5))
from fanpy.wfn.composite.product import ProductWavefunction
wfn = ProductWavefunction([wfn_orig, wfn_rbm])
{nproj_txt}
# Select parameters that will be optimized""",
            )
            text = text.replace(
                    "(wfn, np.ones(wfn.nparams, dtype=bool))", "(wfn_rbm, np.ones(wfn_rbm.nparams, dtype=bool))"
                    )
        if rbm_nn:
            text = text.replace(
                "# Select parameters that will be optimized",
                f"""# Create product of wavefunction for rbm
wfn_orig = wfn
from fanpy.upgrades.numpy_network import NumpyNetwork
wfn_rbm = NumpyNetwork(nelec, nspin, {rbm_nn_kwargs})
params = []
scale = 1 / np.tanh(1)
if wfn_rbm.num_layers > 1:
    params.append(np.identity(wfn_rbm.nspin))
    for _ in range(wfn_rbm.num_layers - 2):
        params.append(np.identity(wfn_rbm.nspin) * scale)
output_weights = np.ones((wfn_rbm.nelec, wfn_rbm.nspin)) / nelec
if wfn_rbm.num_layers > 1:
    output_weights *= scale
params.append(output_weights)
wfn_rbm._template_params = params
wfn_rbm.assign_params(None)
wfn_rbm.assign_params(wfn_rbm.params + 0.001 * 2 * (np.random.rand(*wfn_rbm.params.shape) - 0.5))
from fanpy.wfn.composite.product import ProductWavefunction
wfn = ProductWavefunction([wfn_orig, wfn_rbm])

# Select parameters that will be optimized""",
            )
            text = text.replace(
                "(wfn, np.ones(wfn.nparams, dtype=bool))", "(wfn_rbm, np.ones(wfn_rbm.nparams, dtype=bool))"
            )

        if cc_chk:
            text = text.replace(
                "wfn = StandardCC(",
                f"""ccsd = StandardCC(nelec, nspin, params=None, ranks=[1, 2])
ccsd.assign_params(np.load("{cc_chk}"))
wfn = StandardCC("""
            )
            text = text.replace(
"wfn.assign_params(",
"""wfn.assign_params(ccsd)
wfn.assign_params(""",
            )

        with open('calculate.py', "w") as f:
            f.write(text)
        os.chdir(cwd)
