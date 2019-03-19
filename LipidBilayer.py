from typing import List
from Molecules import Lipid, Water, Ion


class LipidBilayer(object):

    def __init__(self, filename: str, num_lipids: int, natom_per_lipid: int, num_water: int, num_ions: int,
                 ion_name: str) -> None:
        self.filename: str = filename
        self.num_lipids: int = num_lipids
        self.natom_per_lipid: int = natom_per_lipid
        self.num_water: int = num_water
        self.num_ions: int = num_ions
        self.ion_name: str = ion_name
        self.residue_numbers: List[int] = []
        self.residue_names: List[str] = []
        self.atom_names: List[str] = []
        self.atom_numbers: List[int] = []
        self.coordinates: List[float] = []
        self.box: List[float] = []
        self.natoms: int = 0
        self.lipids: List[Lipid] = []
        self.waters: List[Water] = []
        self.ions: List[Ion] = []

        return

    def read_file(self) -> None:
        f = open(self.filename)
        next(f)  # skip title
        self.natoms = int(next(f))
        for i in range(1, self.natoms + 1):
            line = next(f)
            self.residue_numbers.append(int(line[0:5]))
            self.residue_names.append(line[5:10].strip())
            self.atom_names.append(line[10:15].strip())
            self.atom_numbers.append(int(line[15:20]))
            self.coordinates.append(round(float(line[20:28]) * 10, 2))
            self.coordinates.append(round(float(line[28:36]) * 10, 2))
            self.coordinates.append(round(float(line[36:44]) * 10, 2))
        tmp = next(f).split()
        self.box = list(map(lambda x: round(float(x) * 10, 2), tmp[0:3]))
        f.close()

        return

    def set_molecules(self) -> None:
        atom_names: List[str] = self.atom_names[0:self.natom_per_lipid]
        npl = self.natom_per_lipid
        for i in range(0, self.num_lipids):
            self.lipids.append(Lipid(atom_names, self.coordinates[3 * i * npl:3 * (i + 1) * npl]))

        water_indentation: int = 3 * self.num_lipids * self.natom_per_lipid

        for i in range(0, self.num_water):
            self.waters.append((Water(self.coordinates[9 * i + water_indentation:9 * (i + 1) + water_indentation])))

        ion_indentation: int = water_indentation + 9 * self.num_water

        for i in range(0, self.num_ions):
            self.ions.append(
                (Ion(self.ion_name, self.coordinates[3 * i + ion_indentation:3 * (i + 1) + ion_indentation])))

        return
