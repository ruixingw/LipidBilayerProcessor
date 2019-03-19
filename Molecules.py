import numpy as np
from typing import List, Callable


class Molecule(object):
    def __init__(self):
        self.atoms: List[Atom] = []

    @property
    def atom_names(self):
        atom_names: List[str] = []
        for atom in self.atoms:
            atom_names.extend(atom.atom_name)
        return atom_names

    @property
    def coordinates(self):
        coordinates: List[float] = []
        for atom in self.atoms:
            coordinates.extend(atom.coordinate)
        return coordinates

    def generate_coordinates(self) -> str:
        string: str = ""
        for atom in self.atoms:
            string += "{}     {:8.3f}     {:8.3f}     {:8.3f}\n".format(atom.atom_element, *atom.coordinate)
        return string


class Atom(Molecule):
    def __init__(self, atom_name: str, element: str, coordinate: List[float]):
        super(Atom, self).__init__()
        self.atoms: List[Atom] = [self]
        self.coordinate: List[float] = coordinate

        self.atom_name: str = atom_name
        self.atom_element = element

    def distance_to(self, another: "Atom") -> float:
        return np.linalg.norm(np.array(self.coordinate) - np.array(another.coordinate))


class Lipid(Molecule):

    def __init__(self, atom_names: List[str], coordinates: List[float]):
        super(Lipid, self).__init__()
        self.atoms: List[Atom] = []

        for i, atom in enumerate(atom_names):
            self.atoms.append(Atom(atom, atom[0], coordinates[3 * i:3 * i + 3]))

    def truncate(self, function: Callable[["Lipid"], "Lipid"]):
        return function(self)


class Water(Molecule):

    def __init__(self, coordinates: List[float]):
        super(Water, self).__init__()
        atom_names = ["OH2", "H1", "H2"]
        self.atoms: List[Atom] = []

        for i, atom in enumerate(atom_names):
            self.atoms.append(Atom(atom, atom[0], coordinates[3 * i:3 * i + 3]))


class Ion(Atom):

    def __init__(self, atom_name: str, coordinate: List[float]):
        super(Ion, self).__init__(atom_name, atom_name, coordinate)
        self.atom = self


class Cluster(object):

    def __init__(self, *args: Molecule):
        self.molecules: List[Molecule] = []
        for item in args:
            self.molecules.append(item)

    def generate_coordinates(self) -> str:
        string: str = ""
        for molecule in self.molecules:
            string += molecule.generate_coordinates()
        return string
