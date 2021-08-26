from smartreader import *
import re

class XYZReader:
    """reads xyz file

    """
    def __init__(self, inpfile):
        """ 

        Args:
            inpfile: string that describes the location of xyz file
        """
        self.filename = inpfile
        # FUCK TITLES
        #self.title=''
        self.num_atoms = []
        self.atoms = []
        self.coordinates = []
        self.extract_xyz()

    def extract_xyz(self):
        """ extracts the informations from self.filename

        """
        re_numatom = re.compile(r'^\s*(\d+)\s*$')
        #re_title = re.compile(r'^\s*(.*)\s*$')
        re_coords = re.compile(r'^\s*(\w+)\s+([\-\d\.eEdD]+)\s+([\-\d\.eEdD]+)\s+([\-\d\.eEdD]+)\s*$')
        item_numatom = Item(re_numatom, max_match=-1)
        #item_title = Item(re_title, header=item_numatom, max_match=1)
        item_coords = Item(re_coords, header=item_numatom, max_match=-1)
        with open(self.filename, 'r') as f:
            for line in f:
                #item_numatom(line)
                item_coords(line)
                #item_title(line)
        for i in item_coords.value:
            self.num_atoms.append(len(i))
            self.atoms.append(zip(*i)[0])
            self.coordinates.append([[float(k) for k in j] for j in zip(*zip(*i)[1:])])
