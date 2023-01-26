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
# fbe586a0 ends here

# [[file:../spdkit-python.note::ec59e65f][ec59e65f]]
def view_in_pymol(mol: Molecule, rebond=False, format='pdb'):
    """View molecule object using pymol."""
    import subprocess, tempfile
    import time

    if rebond:
        # rebuild connectivity without periodic images
        lat = mol.unbuild_crystal()
        mol.rebond()
        if lat:
            mol.set_lattice(lat)

    with tempfile.NamedTemporaryFile(mode="w", suffix=f".{format}") as f:
        molfile = f.name
        title = mol.title
        lat = mol.get_lattice()
        mol.to_file(molfile)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".py") as f:
            print(f"pymol.cmd.load('{molfile}', '{title}')\n", file=f)
            # pymol doesnot recognize lattice data in mol2 format
            # set lattice/cell object
            if lat:
                a, b, c = lat.lengths()
                alpha, beta, gamma = lat.angles()
                print(
                    f"pymol.cmd.set_symmetry('all', {a}, {b}, {c}, {alpha}, {beta}, {gamma})",
                    file=f,
                )
            print("pymol.cmd.show('sphere')", file=f)
            print("pymol.cmd.show('stick')", file=f)
            print("pymol.cmd.show('cell')", file=f)
            print("pymol.cmd.label('all', 'ID')", file=f)
            # print("pymol.cmd.orient()", file=f)
            f.flush()
            # return subprocess.run(["pymol", "-J", f.name])
            # wait one second for pymol reading temp files
            p = subprocess.Popen(
                ["pymol", f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            time.sleep(2)
            return p


def view_traj_in_pymol(mols: list[Molecule], separated=True, format='pdb'):
    """View a list of molecule objects as trajectory using pymol."""
    import subprocess, tempfile, os
    import time

    with tempfile.TemporaryDirectory() as td:
        # create pymol script
        py = os.path.join(td, "script.py")
        fpy = open(py, "w")

        # create trajectory animation
        for i, m in enumerate(mols):
            i += 1
            molfile = os.path.join(td, f"{i}.{format}")
            title = m.title
            m.to_file(molfile)
            if separated:
                print(
                    f"pymol.cmd.load('{molfile}', object='traj', state=0)\n", file=fpy
                )
            else:
                print(f"pymol.cmd.load('{molfile}')\n", file=fpy)
            lat = m.get_lattice()
            # pymol doesnot recognize lattice data in mol2 format
            if lat:
                a, b, c = lat.lengths()
                alpha, beta, gamma = lat.angles()
                print(
                    f"pymol.cmd.set_symmetry('all', {a}, {b}, {c}, {alpha}, {beta}, {gamma})",
                    file=fpy,
                )

        print("pymol.cmd.set('movie_fps', 5)", file=fpy)
        print("pymol.cmd.show('sphere')", file=fpy)
        print("pymol.cmd.show('stick')", file=fpy)
        print("pymol.cmd.show('cell')", file=fpy)
        print("pymol.cmd.label('all', 'ID')", file=fpy)
        print("pymol.cmd.orient()", file=fpy)

        fpy.flush()
        # subprocess.run(["pymol", py])
        # time.sleep(50)
        p = subprocess.Popen(
            ["pymol", py], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # wait one second for pymol reading temp files
        time.sleep(2)
        return p
# ec59e65f ends here
