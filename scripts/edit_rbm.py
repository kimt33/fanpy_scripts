import re
import glob
import os
import numpy as np
import extract_energy


def edit_file(
        filepattern, rbm_kwargs="num_layers=1, orders=(1, 2)", wfn_noise=0.001, optimize_both=False,
        old_fanpy=False
    ):
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
       
        os.chdir(dirname)
        # read file to edit
        with open("calculate.py", 'r') as f:
            text = f.read()

        # add snippet to:
        #     instantiate RBM wavefunction
        #     create product wavefunction with original and RBM
        #     overwrite nproj if it is defined wrt wfn.nparams (to wfn_rbm.nparams)
        # FIXME: old fanpy does not use nproj
        if old_fanpy:
            # look for pspace instead
            nproj_txt = re.search(r'\s*(pspace = .+?)\nprint', text, flags=re.DOTALL).group(1)
        else:
            nproj_txt = re.search(r'\s*(nproj = .+?)\nprint', text, flags=re.DOTALL).group(1)
        if 'wfn' in nproj_txt:
            nproj_txt = nproj_txt.replace('wfn', 'wfn_rbm')
        else:
            nproj_txt = ''

        text = text.replace(
            "# Select parameters that will be optimized",
            f"""# Create product of wavefunction for rbm
from fanpy.wfn.network.rbm import RestrictedBoltzmannMachine
wfn_rbm = RestrictedBoltzmannMachine(nelec, nspin, nbath=nspin, {rbm_kwargs})
wfn_rbm.assign_params(wfn_rbm.params + {wfn_noise} * 2 * (np.random.rand(*wfn_rbm.params.shape) - 0.5))
from fanpy.wfn.composite.product import ProductWavefunction
wfn = ProductWavefunction([wfn_orig, wfn_rbm])
{nproj_txt}

# Select parameters that will be optimized""",
        )

        # replace all instances of wfn before rbm+product wavefunction to wfn_orig
        before_param_select_txt = re.search(
            r'Initialize wavefunction(.+)# Create product of wavefunction for rbm', text, flags=re.DOTALL
        ).group(1)
        text = text.replace(
            before_param_select_txt, 
            before_param_select_txt.replace("wfn.", "wfn_orig.").replace("wfn =", "wfn_orig ="),
        )

        # select wavefunction parameters to optimize
        if optimize_both:
            # optimize both original and rbm wavefunction parameters
            text = text.replace(
                "(wfn, np.ones(wfn.nparams, dtype=bool))", 
                "(wfn_orig, np.ones(wfn_orig.nparams, dtype=bool)), (wfn_rbm, np.ones(wfn_rbm.nparams, dtype=bool))",
            )
        else:
            # optimize only rbm wavefunction parameters
            text = text.replace(
                "(wfn, np.ones(wfn.nparams, dtype=bool))", 
                "(wfn_rbm, np.ones(wfn_rbm.nparams, dtype=bool))",
            )

        # overwrite
        with open('calculate.py', "w") as f:
            f.write(text)
        # reset directory location
        os.chdir(cwd)
