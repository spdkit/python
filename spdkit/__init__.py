# [[file:../spdkit-python.note::fbe586a0][fbe586a0]]
from .spdkit import *


def from_ase_atoms(ase_atoms):
    """Construct a molecule object from ase Atoms"""

    import ase

    assert isinstance(ase_atoms, ase.Atoms), "not an ase.Atoms object!"

    atoms = [from_ase_atom(a) for a in ase_atoms]

    mol = Molecule.from_atoms(atoms)
    # lattice object
    if ase_atoms.pbc.all():
        lat = Lattice(ase_atoms.cell.tolist())
        mol.set_lattice(lat)
    return mol


def to_ase_atoms(mol: Molecule):
    """Construct an ase `Atoms` object from `Molecule`"""
    import ase

    atoms = [ase.Atom(symbol=a.symbol(), position=a.position()) for _, a in mol.atoms()]
    ase_atoms = ase.Atoms(atoms)
    lat = mol.get_lattice()
    if lat:
        ase_atoms.set_pbc(True)
        ase_atoms.set_cell([lat.vector_a(), lat.vector_b(), lat.vector_c()])
    return ase_atoms


def from_ase_atom(ase_atom) -> Atom:
    """Construct a molecule object from ase Atom object"""

    import ase

    assert isinstance(ase_atom, ase.Atom), "not an ase.Atom object!"
    return Atom(ase_atom.symbol, list(ase_atom.position))


def to_ase_atom(atom: Atom):
    """Construct an ase `Atom` object from native `Atom`"""
    import ase

    ase.Atom(symbol=atom.symbol(), position=atom.position())


def view_in_pymol(mol: Molecule):
    """View molecule object using pymol"""
    import subprocess, tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb") as f:
        molfile = f.name
        mol.to_file(molfile)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".py") as f:
            print(f"pymol.cmd.load('{molfile}')\n", file=f)
            print("pymol.cmd.show('sphere')", file=f)
            print("pymol.cmd.show('cell')", file=f)
            f.flush()
            subprocess.call(["pymol", "-q", f.name])
# fbe586a0 ends here
