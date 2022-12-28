# [[file:../spdkit-python.note::fbe586a0][fbe586a0]]
from .spdkit import *

class Molecule(_Molecule):
    """ The main object featured in this library. This object
    represents a molecule with atoms and bonds.
    """
    def from_ase_atoms(ase_atoms):
        """Construct a molecule object from ase Atoms"""

        import ase

        assert isinstance(ase_atoms, ase.Atoms), "not an ase.Atoms object!"

        atoms = []
        for a in ase_atoms:
            new_atom = Atom.from_ase_atom(a)
            atoms.append(new_atom)

        # FIXME: lattice object
        # cell = Cell(tvs=ase_atoms.cell.tolist())
        # molecule = Molecule(atoms, cell=cell)
        return _Molecule.from_atoms(atoms)

class Atom(_Atom):
    def from_ase_atom(ase_atom):
        """Construct a molecule object from ase Atom object"""

        import ase
        assert isinstance(ase_atoms, ase.Atom), "not an ase.Atom object!"
        return _Atom(a.symbol, list(a.position))
# fbe586a0 ends here
