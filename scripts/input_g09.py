import numpy as np
import os
import re

class GaussianInput:
    """Class for handling Gaussian inputs."""
    def __init__(self, atoms, coords, method='', basis='', route='', charge=0, spinmult=1, title='', filename='', template='default', nosave=False):
        """Initialize.

        Parameters
        ----------
        atoms : {list, tuple}
            Chemical symbols for each atom in the system.
            in proper capitalization
        coords
            array of coordinates (array of dim 3)
        NOTE: Defaults were arbitrary and depended only on the template I used to create it
        """
        if filename!='':
            self.filename=filename
        else:
            self.filename='gaussian.com'
        self.basename=self.filename[:-4] #NOTE: assume specified filename has three character extension
        #Link0
        self.chk=self.basename+'.chk'
        self.mem='1500MB'
        self.nosave=nosave
        #Route section
        if method!='':
            self.method=method
        else:
            self.method='uwb97xd'
        if basis!='':
            self.basis = basis
        else:
            self.basis = 'aug-cc-pvtz'
        #for ext in ['gbs','nwchem','davidh5']:
        #    if os.path.isfile(context.get_fn('basis/{0}.{1}'.format(self.basis,ext))):
        #        path = context.get_fn('basis/{0}.{1}'.format(self.basis,ext))
        #        self.basis = 'gen'
        #        break
        # this needs to be better
        if template == 'default':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry ' + route
        elif template == 'opt':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry opt=tight ' + route
        elif template == 'stable':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry stable=opt ' + route
        elif template == 'freq':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry freq ' + route
        elif template == 'pop':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry density=current pop=(chelpg,npa) IOp(6/80=1) ' + route
        route = route.split()

        self.keywords = {}
        self.routerest = []
        for i in route:
            if 'opt' == i.lower()[:3]:
                self.keywords['opt'] = i[4:]
            elif 'scf' == i.lower()[:3]:
                self.keywords['scf'] = i[4:]
            elif 'integral' == i.lower()[:8]:
                self.keywords['integral'] = i[9:]
            elif 'density' == i.lower()[:7]:
                self.keywords['density'] = i[8:]
            elif 'pop' == i.lower()[:3]:
                self.keywords['pop'] = i[4:]
            elif 'freq' == i.lower()[:4]:
                self.keywords['freq'] = i[5:]
            elif 'stable' == i.lower()[:6]:
                self.keywords['stable'] = i[7:]
            elif 'nosymmetry' == i.lower():
                self.keywords['nosymmetry'] = ''
            else:
                self.routerest.append(i)

        #title section
        chemformula = ''.join(atoms)
        self.title=title+' '+chemformula+' '+self.method+'/'+self.basis
        #molecule specification (default: 0 1)
        self.charge=charge #difference in charge from neutral
        self.spinmultiplicity=spinmult #1 for singlet, 2 for doublet, ...
        #atoms
        self.atoms=atoms
        self.coords=coords

    def write_input(self,filename=''):
        if filename == '':
            filename = self.filename
        with open(filename,'w') as fp:
            if self.nosave:
                fp.write('%oldchk='+self.chk+'\n')
                fp.write('%mem='+self.mem+'\n')
                fp.write('%nosave\n')
            else:
                fp.write('%chk='+self.chk+'\n')
                fp.write('%mem='+self.mem+'\n')
            # route
            fp.write('#p {0}/{1} '.format(self.method, self.basis))
            routerest = ''
            for keyword, val in self.keywords.items():
                if val == '':
                    routerest += keyword + ' '
                else:
                    routerest += keyword + '=' + val + ' '
            routerest += ' '.join(self.routerest)
            fp.write(routerest+'\n\n')
            # title
            fp.write(self.title+'\n\n')
            # mol spec
            fp.write(str(self.charge)+' '+str(self.spinmultiplicity)+'\n')
            # coord
            for atom,coord in zip(self.atoms,self.coords):
                line = atom + '{0:>17f}{1:>17f}{2:>17f}'.format(*[float(i) for i in coord])
                fp.write(line+'\n')
            fp.write('\n')
            # basis set
            if self.basis == 'gen':
                fp.write('''H     0   
S   8   1.00
    188.61445                 .00096385
     28.276596                .00749196
      6.4248300               .03759541
      1.8150410               .14339498
       .59106300              .34863630
       .21214900              .43829736
       .07989100              .16510661
       .02796200              .02102287
****
''')
                # fp.write(gbs_format())
            fp.write('\n\n\n')

    def scan_geometry_atom(self,atomindex,coord1,coord2,numsteps):
        """returns the generator for the input object where the coordinate of atom (specified by atomindex) changes incrementally following the path of coord1 to coord2

        """
        coord1=np.array(coord1)
        coord2=np.array(coord2)
        stepvec=(coord2-coord1)/(numsteps-1)
        coord_range=[coord1+i*stepvec for i in range(numsteps)]
        for coord in coord_range:
            self.update_coord_atom(coord.tolist(),atomindex)
            yield self
"""
a=g09_input(['H','O'],[[0.0,0.0,0.0],[2.4,0.0,0.0]])
a.write_input('testest.txt')
"""
