Directory structure
===================
Each "level" correspond to the nested directory. So second level describes a directory in the
directories in the first level.

Zeroth level correspond to the base directory on which all your database directories are located.

First level corresponds to some arbitrary groupings of systems. Think of it as the name of the
database.

Second level corresponds to the system. Think of it as molecular geometry.

Third level corresponds to the basis set.

Fourth level corresponds to the method. Methods with different hyperparameter settings are
considered different methods. For example, CCSD and CCSDT are considered to be different methods.
Different objectives such as projected vs variational may need to be differentiated here.

Fifth level corresponds to a calculation for the given system at method/basis. Note that a given
calculation may be repeated multiple times due to restarting the calculations or due to using
different hyperparameters for the calculations. Each calculation should be stored as a separate 
directory.

Example:
./paldus/p4_a20_alpha2.0/sto-6g/ap1rog/energy_min_1/
zeroth level = ./
  current directory is always considered the base directory
first level = paldus
  paldus is name for the class of model systems proposed by J. Paldus
second level = p4_a20_alpha2.5
  p4 is describes two H2 molecules being pushed together into a square
  a20 means that H2 molecules have bond length of 2 Bohr
  alpha2.5 means that the distance between the two H2 molecules is 2 Bohr
third level = sto-6g
fourth level = ap1rog
fifth level = energy_min_1
  1 means that this is the second calculation (zero-indexed) for minimizing the energy of 
  AP1roG/STO-6G for the system P4 with 2 Bohr bond length and 2.5 Bohr distance between the two H.

Note that the variations in fifth level can be moved into the fourth level. In this example, we
could have the following structure: ./paldus/p4_a20_alpha2.0/sto-6g/ap1rog_energy_min/1/

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
    alias cdw='cd /blue/rmirandaquintana/yourusername/'
    export HORTONPYTHON=/home/kimt1/envs/horton/bin/python
    export PYTHONPATH=$PYTHONPATH:/home/kimt1/

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

This file (and the directory) can be generated using the python module `/home/kimt1/run_calc.py`.
Use the function `write_coms`. Here's an example:
`run_calc.write_coms('paldus/*/sto-6g/', memory='1gb', charge=0, multiplicity=1, units="AU")`
You put this command in a python script

If you have nonstandard basis set directory name (i.e. not included in ~kimt1/basis/), then you need
to specify the name of the basis filename
`run_calc.write_coms('paldus/*/sto-6g/', memory='1gb', charge=0, multiplicity=1, units="AU", basisfile="basisfilename")`
Note that you do not include the '.gbs' in the basisfilename.

Take special care of the units.


5. Run HF calculations.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.com')`

Note that above runs Gaussian locally. To submit the job, specify the memory and time:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.com', time='1d', memory='3gb', outfile='results.out')`
You put this command in a python script

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


6. Generate fchk from the chk files.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.chk')`
You put this command in a python script

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


7. Generate one and two electron integrals as npy array.

You can use the `run_calcs` function in the `run_calc` module in `/home/kimt1/`. Here's an example:
`run_calc.run_calcs('paldus/*/sto-6g/hf/hf_sp.fchk')`
You put this command in a python script

*Note that `run_calcs` function assumes that you are in the base directory of the database.*


8. Make directories for the wavefunction/method.

You can use the `make_wfn_dirs` function in the `run_calc` module in `/home/kimt1/`. Here's an
example:
`run_calc.make_wfn_dirs('paldus/*/sto-6g/', 'ap1rog', 10)`
You put this command in a python script

Note that you specify the number of directories that you create. This will become relevant a bit
 later.


9. Make scripts for running calculations. Specify the details of your calculation here.

You can use the `write_wfn_py` function in the `run_calc` module in `/home/kimt1/`. Here's an
example:
```python
    run_calc.write_wfn_py('paldus/*/sto-6g/ap1rog/', 4, 'ap1rog', 
                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='one_energy', solver='cma',
                           ham_noise=1e-3, wfn_noise=1e-2,
                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
```
```python
    run_calc.write_wfn_py('paldus/*/sto-6g/ap1rog/0', 4, 'ap1rog', 
                           optimize_orbs=True, pspace_exc=[1, 2, 3, 4], objective='one_energy', solver='cma',
                           ham_noise=1e-3, wfn_noise=1e-2,
                           load_orbs=None, load_ham=None, load_wfn=None, old_fanpy=True)
```
You put this command in a python script

Here are the arguments and keyword arguments:
    nelec : int
        Number of electrons.
    wfn_type : str
        Type of wavefunction.
        One of `
        "ci_pairs", "cisd", "doci", "fci",
        "mps",
        "determinant-ratio",
        "ap1rog", "apr2g", "apig", "apsetg", "apg",
        "network", "rbm",
        "basecc", "standardcc", "generalizedcc", "senioritycc", "pccd", "ccsd", "ccsdt", "ccsdtq",
        "ap1rogsd", "ap1rogsd_spin", "apsetgd", "apsetgsd", "apg1rod", "apg1rosd",
        "ccsdsen0", "ccsdqsen0", "ccsdtqsen0", "ccsdtsen2qsen0".
    optimize_orbs : bool
        If True, orbitals are optimized.
        If False, orbitals are not optimized.
        By default, orbitals are not optimized.
        Not compatible with faster fanpy (i.e. `old_fanpy=False`)
    pspace_exc : list of int
        Orders of excitations that will be used to build the projection space.
        Default is first, second, third, and fourth order excitations of the HF ground state.
        Used for slower fanpy (i.e. `old_fanpy=True`)
    nproj : int
        Number of projection states that will be used.
        Default uses all possible projection states (i.e. Slater determinants.)
        Used for faster fanpy (i.e. `old_fanpy=False`)
        If 0, maximum possible projection space is used (all possible Slater determinants)
        If negative, then the projection space will have the size of `-nproj` times the number of
        parameters.
        If positive, then the projection space will have the size of `nproj`
    objective : str
        Form of the Schrodinger equation that will be solved.
        Use `projected` to solve the Schrodinger equation as a system of equations.
        Use `least_squares` to solve the Schrodinger equation as a squared sum of the system of
        equations.
        Use "one_energy" to minimize the energy projected on one sided.
        Use "variatioinal" to minimize the energy projected on both sided.
        Must be one of "system", "least_squares", "one_energy", and "variational".
        By default, the Schrodinger equation is solved as system of equations.
        "least_squares" is not supported in faster fanpy (i.e. `old_fanpy=False`)
    solver : str
        Solver that will be used to solve the Schrodinger equation.
        Keyword `cma` uses Covariance Matrix Adaptation - Evolution Strategy (CMA-ES).
        Keyword `diag` results in diagonalizing the CI matrix.
        Keyword `minimize` uses the BFGS algorithm.
        Keyword `least_squares` uses the Trust Region Reflective Algorithm.
        Keyword `root` uses the MINPACK hybrd routine.
        Must be one of `cma`, `diag`, `least_squares`, or `root`.
        Must be compatible with the objective.
    ham_noise : float
        Scale of the noise to be applied to the Hamiltonian parameters.
        The noise is generated using a uniform distribution between -1 and 1.
        By default, no noise is added.
    wfn_noise : bool
        Scale of the noise to be applied to the wavefunction parameters.
        The noise is generated using a uniform distribution between -1 and 1.
        By default, no noise is added.
    solver_kwargs : str
        String to be added as arguments and keyword arguments when running the solver.
        To be added after `solver(objective, `
        Default settings are used if not provided.
    wfn_kwargs : str
        String to be added as arguments and keyword arguments when instantiating the
        wavefunction.
        To be added after `Wavefunction(nelec, nspin, params=params, memory=memory, `
        Default settings are used if not provided.
    load_orbs : str
        Numpy file of the orbital transformation matrix that will be applied to the initial
        Hamiltonian.
        If the initial Hamiltonian parameters are provided, the orbitals will be transformed
        afterwards.
    load_ham : str
        Numpy file of the Hamiltonian parameters that will overwrite the parameters of the initial
        Hamiltonian.
    load_wfn : str
        Numpy file of the wavefunction parameters that will overwrite the parameters of the initial
        wavefunction.
    load_chk : str
        Base name for the numpy files of the checkpoint for the optimization.
        At the moment, multiple checkpoint files are created for each object with active parameters.
        They should have the name format of `checkpoint_{name}.npy`. Then, specify `checkpoint.npy`
        here.
    memory : str
        Memory available to run the calculation.
        e.g. '2gb'
        Default assumes no restrictions on memory usage
    filename : str
        Filename to save the generated script file.
        Default just prints it out to stdout.
    old_fanpy : bool
        Use old, slower (but probably more robust) fanpy.
        Default uses faster fanpy.
        Some features are not avaialble on new fanpy.

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
You put this command in a python script

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
You put this command in a python script

If you are restarting your calculation from the same directory, it might be a good idea to rename
checkpoint file and the output file so that they don't get overwritten.
You can just run the calculation from another directory also. 


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
You put this command in a python script

This module changes pretty frequently and isn't maintained that well. It'd be better for you to
simply skim it and make a module yourself. (All it does is replace the appropriate parts with
another).


Killing Jobs
============
To see details on the jobs you are running (e.g. directory in which the calculations are running), I
I use the following command:

```bash
    scontrol show jobid -dd | grep UserId=e.nash -B 1 -A 25
```

Then I find the jobid's of the jobs that I want to kill from the given information, and then use

```bash
    scancel jobid
```

to kill it.

This process is a little more streamline with the script `./scripts/scan_jobs.py`. Note that command
in line 57 is used to kill the given job. Make sure to put a condition (i.e. if statement) when
killing jobs to avoid killing all your jobs.

Example Use Case:
If I find out that a job for one system and wavefunction type will take too long to finish, I should
kill it. I will filter out the jobs that are associated with these calculations (from the directory
structure) and kill only these jobs.


Extracting Data
===============
First, make a directory at the "first level" calld "data".

To extract data, see the script `./example_extract.py`. Copy over or make your own version of this
script in the "first level"

Line 39 defines the name of the systems used and the correponding units. This corresponds to the
"second level" as described above. The script as it is now *will assume* that the name specified in
`systems` will be succeeded with exactly "\_number" (no backslash). where number corresponds to
value that identifies the system. If your naming structure is different, line 151 (and 201 and 237)
must be changed.

Line 45 defines the names of the methods used and the corresponding labels. This corresponds to the
"fourth level" as described above. This must match the names of the directories at this level.

Line 118 defines the names of the basis set used. This corresponds to the "third level" as described
above. This must match the names of the directories at this level.

Take a look at line 127 and 128 (and 196 and 231). This is how the results files are found. So if
you some other directory naming scheme, these lines need to be modified.

This script will generate a csv file for each combination of system, basis set, and method provided.
This csv file then can be read into an excel (or equivalent spreadsheet program) or be plotted
through matplotlib (or equivalent plotting module)


Common Issues
=============
- If you get an error something along the lines of "unable to generate `nproj=N` determinants", then
  that means that the number of projected states selected is greater than the number of states that
  can be generated by `FanCI`. You would need to decrease the number of projected states by
  specifying a smaller number in `nproj = {smaller_number}`. In the case you have as many states as
  there are parameters, you can change the projection state to the maximum number of states possible
  (i.e. `nproj = int(scipy.special.comb(nspin // 2, nelec - nelec // 2) * scipy.special.comb(nspin // 2, nelec // 2))`)
  or you can specify `nproj=0` in `write_wfn_py`

