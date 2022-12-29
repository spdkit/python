# [[file:../spdkit-python.note::fbe586a0][fbe586a0]]
from .spdkit import *

class Molecule(_Molecule):
    """ The main object featured in this library. This object
    represents a molecule with atoms and bonds.
    """
    @classmethod
    def from_ase_atoms(cls, ase_atoms):
        """Construct a molecule object from ase Atoms"""

        import ase

        assert isinstance(ase_atoms, ase.Atoms), "not an ase.Atoms object!"

        atoms = []
        for a in ase_atoms:
            new_atom = Atom.from_ase_atom(a)
            atoms.append(new_atom)

        mol = cls(atoms)
        # lattice object
        if ase_atoms.pbc.all():
            lat = Lattice(ase_atoms.cell.tolist())
            mol.set_lattice(lat)
        return mol

    # @classmethod
    # def from_file(cls, path: str):
    #     return super(cls, cls).from_file(path)

    def to_ase_atoms(self):
        """Construct an ase `Atoms` object from `Molecule`"""
        import ase
        atoms = [ase.Atom(symbol=a.symbol(), position=a.position()) for _, a in self.atoms()]
        ase_atoms = ase.Atoms(atoms)
        lat = self.get_lattice()
        if lat:
            ase_atoms.set_pbc(True)
            ase_atoms.set_cell([lat.vector_a(), lat.vector_b(), lat.vector_c])
        return ase_atoms

class Atom(_Atom):
    """Atom is the smallest particle still characterizing a chemical element"""
    def from_ase_atom(ase_atom):
        """Construct a molecule object from ase Atom object"""

        import ase
        assert isinstance(ase_atom, ase.Atom), "not an ase.Atom object!"
        return _Atom(ase_atom.symbol, list(ase_atom.position))

    def to_ase_atom(self, ase_atom):
        """Construct an ase `Atom` object from native `Atom`"""
        Atom(ase_atom.symbol, ase_atom.position.tolist())

class Lattice(_Lattice):
    """Periodic 3D lattice"""
    pass
# fbe586a0 ends here
