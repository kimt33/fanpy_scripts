Directory structure
===================
First level corresponds to some arbitrary groupings of systems. Think of it as the name of the
database.

Second level corresponds to the system. Think of it as molecular geometry.

Third level corresponds to the basis set.

Fourth level corresponds to the method. Methods with different hyperparameter settings are
considered different methods. For example, CCSD and CCSDT are considered to be different methods.
Different objectives such as projected vs variational may need to be differentiated here.

Fifth level corresponds to the a calculation for the given system at method/basis. Note that a given
calculation may be repeated multiple times as a result of checkpointing or restarting the
calculations, and each calculation may be stored as a separate directory.


Setup
=====

1. Load conda

`module load conda`

2. Activate virtual environment

`conda activate /home/kimt1/envs/fanpy`

3. Add commands to `.bashrc` if you'd like

```bash
  
    module load gaussian/16
    module load conda
    alias fanpy_activate='conda activate /home/kimt1/envs/fanpy'
```

Steps
=====

1. Make new database directory.

The directory on which you create this new database directory will be called base directory.


2. Make directory for each geometry.
   
Within each directory, make an ".xyz2" file which consists of the atom and its coordinates (space
separated) and different atoms separated by newline. Think of ".xyz" without the header. Though it
can be specified later on, *DEFAULT WILL ASSUME THAT IT IS IN ATOMIC UNITS*.

For example scripts, see scripts of the form `make_xxx_xyz.py`.


3. Make a directory for the basis set.
   
Within each directory, make/copy a ".gbs" file that you can obtain from
https://www.basissetexchange.org/ (Format: Gaussian). It contains the atomic basis set information
that will be used to run Gaussian calculations.

For example scripts, see scripts of the form `make_basis.py`.


4. Make a directory for HF. Within each directory, make a ".com" file to run HF.

This file (and the directory) can be generated using the python module `run_calc` in `/home/kimt1/`.
Use the function `write_coms`. Here's an example:
`run_calc.write_coms('paldus/*/sto-6g/', memory='1gb', charge=0, multiplicity=1, units="AU")`


5. Run HF calculations.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.com')`

Note that above runs Gaussian locally. To submit the job, specify the memory and time:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.com', time='1d', memory='3gb', outfile='results.out')`

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


6. Generate fchk from the chk files.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf_sp.chk')`

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


7. Generate one and two electron integrals as npy array.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf_sp.fchk')`

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


8. Make directories for the wavefunction/method.

You can use the `make_wfn_dirs` function in the `run_calc` module in `/home/kimt1/`. Here's an
example:
`run_calc.make_wfn_dirs('paldus/*/sto-6g/', 'ap1rog', 10)`

Note that you specify the number of directories that you create. This will become relevant a bit
 later.


9. Make scripts for running calculations. Specify the details of your calculation here.

You can use the `write_wfn_py` function in the `run_calc` module in `/home/kimt1/`. Here's an
example:
```python
    run_calc.write_wfn_py('paldus/*/sto-6g/ap1rog/', 4, 'ap1rog', 
                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='one_energy', solver='cma',
                           ham_noise=1e-3, wfn_noise=1e-2,
                           load_orbs=None, load_ham=None, load_wfn=None)
```

10. Run/submit the job

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/ap1rog/calculate.py', time='1d', memory='3gb', outfile='results.out')`
This will result in running the calculation *as many times as there are directories*. This means
that if 10 directories were created in step 8, the calculation will be run 10 times (in separate
directories) or until time runs out or until it gets killed.

To run each calculation separately, use the following commands:
```python
    run_calc.write_wfn_py('paldus/*/sto-6g/ap1rog/*/', 4, 'ap1rog', 
                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='one_energy', solver='cma',
                           ham_noise=1e-3, wfn_noise=1e-2,
                           load_orbs=None, load_ham=None, load_wfn=None)
    run_calc.run_calcs('paldus/*/sto-6g/ap1rog/*/calculate.py', time='1d', memory='3gb', outfile='results.out')
```

Using the same example as above, this will result in 10 jobs being submitted, each for running one
calculation.

To make a special script for a given directory, replace the `*` with the desired directory.

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


11. Rerun calculations

If your calculation crashed (e.g. due to insufficient resources), you can restart your calculation
from a checkpoint.

```python
    run_calc.write_wfn_py('paldus/*/sto-6g/ap1rog/0/', 4, 'ap1rog', 
                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='one_energy', solver='cma',
                           ham_noise=1e-3, wfn_noise=1e-2,
                           load_orbs=None, load_ham=None, load_wfn=None, load_chk='checkpoint.npy')
    run_calc.run_calcs('paldus/*/sto-6g/ap1rog/0/', time='1d', memory='3gb', outfile='results.out')
```

If you are restarting your calculation from the same directory, it might be a good idea to rename
checkpoint file and the output file so that they don't get overwritten.


12. Modify calculations

If you need to change the generated script, you can edit the generated script to fit your needs.
For example, if you want to start your calculation from a wavefunction that is related to your
wavefunction of interest (e.g. using optimized AP1roG to run APIG), then you must write some sort
of script to convert the parameters of one wavefunction into another.

To make the process a little bit easier, you can use the `edit_file` function in the
`edit_calculate` module in `/home/kimt1/`. Note that there's a lot of hard coding in this one, so
you'll like have to make changes yourself.

```python
    for dirname in sorted(glob('paldus/*/sto-6g/apig'), reverse=True):
        run_calc.write_wfn_py(f'{dirname}/0', 8, 'apig',
                              optimize_orbs=True, pspace_exc=[2, 4, 6, 8], objective='one_energy', solver='minimize',
                              ham_noise=1e-2, wfn_noise=1e-3, memory=None,
                              load_orbs=None, load_ham=None, load_wfn=None, load_chk='../../ap1rog/2/checkpoint.npy')
        # run_nn
        edit_calculate.edit_file(f'{dirname}/0/calculate.py', truncate_projection=False, proj_seniority=True, energy_constraint=False, cma=False, ap1rog_chk=True)
```

This module changes pretty frequently and isn't maintained that well. It'd be better for you to
simply skim it and make a module yourself. (All it does is replace the appropriate parts with
another).
