from typing import List
from Molecules import Lipid, Atom, Molecule
from LipidBilayer import LipidBilayer
import numpy as np


def trunc(this: Lipid) -> Lipid:
    carbon1 = this.atoms[37]
    carbon2 = this.atoms[84]
    nextCarbon1 = this.atoms[40]
    nextCarbon2 = this.atoms[87]
    this.atoms = this.atoms[0:41] + this.atoms[84:88]
    carbon1coord = np.array(carbon1.coordinate)
    carbon2coord = np.array(carbon2.coordinate)
    nextCarbon1coord = np.array(nextCarbon1.coordinate)
    nextCarbon2coord = np.array(nextCarbon2.coordinate)
    a = nextCarbon1coord - carbon1coord
    b = nextCarbon2coord - carbon2coord
    a = carbon1coord + a / np.linalg.norm(a) * 1.11
    b = carbon2coord + b / np.linalg.norm(b) * 1.11
    nextCarbon1.coordinate = list(a)
    nextCarbon2.coordinate = list(b)
    nextCarbon1.atom_name = "HT"
    nextCarbon1.atom_element = "H"
    nextCarbon2.atom_name = "HT"
    nextCarbon2.atom_element = "H"

    return this


def generateCoordinate(cluster: List[Molecule]) -> str:
    string = "#p blyp/def2svp opt freq scrf=(smd,solvent=water) em=gd3bj\n\ntitle\n\n0 1\n"
    for molecule in cluster:
        string += molecule.generate_coordinates()

    string += "\n\n"
    return string


if __name__ == "__main__":
    bilayer = LipidBilayer("../inputs/be.gro", 80, 131, 3168, 40, "Be")
    bilayer.read_file()
    bilayer.set_molecules()

    good_lipids = []
    for lipid in bilayer.lipids:
        lipid.truncate(trunc)
        good_lipid = True
        for atom in lipid.atoms:
            if not good_lipid:
                break
            for ion in bilayer.ions:
                if ion.distance_to(atom) < 3:
                    good_lipid = False
                    break
        if good_lipid:
            good_lipids.append(lipid)

    i = 0
    for lipid in good_lipids[0:5]:
        i += 1
        print(i)
        lipid.first_shell_waters = []
        for atom in lipid.atoms:
            for water in bilayer.waters:
                for water_atom in water.atoms:
                    if water_atom.distance_to(atom) < 3:
                        lipid.first_shell_waters.append(water)
                        break
        lipid.first_shell_waters = set(lipid.first_shell_waters)

    for i, lipid in enumerate(good_lipids[0:5]):
        with open("../outputs/" + str(i+1) + "-w.gjf", 'w') as f:
            f.write(generateCoordinate([lipid] + list(lipid.first_shell_waters)))
        with open("../outputs/" + str(i + 1) + ".gjf", 'w') as f:
            f.write(generateCoordinate([lipid]))



    # for water in bilayer.waters:
    #     for atom in water.atoms:
    #         if atom.atom_element == "H":
    #             todist.append(atom)
    #
    # for ion in bilayer.ions:
    #     ion.distance_to_lipids = []
    #     ion.distance_to_waters = []
    #     for lipid in bilayer.lipids:
    #         max_dist: float = 1.0e20
    #         for oxygen in lipid.atoms:
    #             if oxygen.atom_name in ["O11", "O12", "O13", "O14"]:
    #                 dist = ion.distance_to(oxygen)
    #                 if dist < max_dist:
    #                     max_dist = dist
    #         ion.distance_to_lipids.append([lipid, max_dist])
    #     for water in bilayer.waters:
    #         dist = ion.distance_to(water.atoms[0])
    #         ion.distance_to_waters.append([water, dist])
    #
    #     ion.distance_to_lipids = sorted(ion.distance_to_lipids, key=lambda x: x[1])
    #     ion.distance_to_waters = sorted(ion.distance_to_waters, key=lambda x: x[1])
    #
    # good_ion = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: []}
    # for ion in bilayer.ions:
    #     for i in range(0, 100):
    #         if ion.distance_to_lipids[i][1] > 2.5:
    #             break
    #     ion.bind_lipids = i
    #     good_ion[i].append(ion)
    #     for i in range(0, 100):
    #         if ion.distance_to_waters[i][1] > 3:
    #             break
    #     ion.bind_waters = i
    #     ion.bind_total = ion.bind_lipids + ion.bind_waters
