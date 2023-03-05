#! /usr/bin/env python
from pathlib import Path


def write_potcar(species: list[str], potcar="POTCAR"):
    """Writes the POTCAR file."""
    with open(potcar, "w") as potfile:
        for symbol in species:
            ppp = potcar_path_for_symbol(symbol)
            with open(ppp) as ppp_file:
                for line in ppp_file:
                    potfile.write(line)


# pot/Fe/POTCAR
def potcar_path_for_symbol(symbol: str):
    path = Path("pot")
    return str(path / symbol / "POTCAR")


def species_from_poscar(poscar="POSCAR") -> list[str]:
    with open(poscar) as fp:
        for _ in range(6):
            line = next(fp)
        return line.split()


if __name__ == "__main__":
    species = species_from_poscar()
    write_potcar(species)
