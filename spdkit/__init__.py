# [[file:../spdkit-python.note::fbe586a0][fbe586a0]]
from spdkit import *


def from_ase_atoms(ase_atoms):
    """Construct a Molecule object from ase Atoms"""

    import ase

    assert isinstance(ase_atoms, ase.Atoms), "not an ase.Atoms object!"

    atoms = [from_ase_atom(a) for a in ase_atoms]

    mol = Molecule.from_atoms(atoms)
    # lattice object
    if ase_atoms.pbc.all():
        lat = Lattice(ase_atoms.cell.tolist())
        mol.set_lattice(lat)
    return mol


def from_pmg_molecule(pmg_molecule):
    """Construct a Molecule object from pymatgen Molecule object"""

    import pymatgen

    assert isinstance(
        pmg_molecule, pymatgen.core.structure.Molecule
    ), "not a pymatgen.core.structure.Molecule object!"

    atoms = [from_pmg_site(s) for s in pmg_molecule]

    mol = Molecule.from_atoms(atoms)
    return mol


def from_pmg_structure(pmg_structure):
    """Construct a Molecule object from pymatgen Structure object"""

    import pymatgen

    assert isinstance(
        pmg_structure, pymatgen.core.structure.Structure
    ), "not a pymatgen.core.structure.Structure object!"

    atoms = [from_pmg_site(s) for s in pmg_structure]

    mol = Molecule.from_atoms(atoms)
    # lattice object
    lat = Lattice(pmg_struct.lattice.matrix.tolist())
    mol.set_lattice(lat)
    return mol


def to_ase_atoms(mol: Molecule):
    """Construct an ase `Atoms` object from `Molecule`"""
    import ase

    atoms = [to_ase_atom(a) for _, a in mol.atoms()]
    ase_atoms = ase.Atoms(atoms)
    lat = mol.get_lattice()
    if lat:
        ase_atoms.set_pbc(True)
        ase_atoms.set_cell([lat.vector_a, lat.vector_b, lat.vector_c])
    return ase_atoms


def to_pmg_molecule(mol: Molecule):
    """Construct a pymatgen `Molecule` object from `Molecule`"""

    import pymatgen

    sites = [to_pmg_site(a) for _, a in mol.atoms()]
    pmg_molecule = pymatgen.core.structure.Molecule.from_sites(sites)
    return pmg_molecule


def to_pmg_structure(mol: Molecule):
    """Construct a pymatgen `Structure` object from `Molecule`"""

    import pymatgen

    lat = mol.get_lattice()
    sites = [
        to_pmg_periodic_site(a, [lat.vector_a, lat.vector_b, lat.vector_c])
        for _, a in mol.atoms()
    ]
    pmg_structure = pymatgen.core.structure.Structure.from_sites(sites)
    return pmg_structure


def from_ase_atom(ase_atom) -> Atom:
    """Construct a Molecule object from ase Atom object"""

    import ase

    assert isinstance(ase_atom, ase.Atom), "not an ase.Atom object!"
    return Atom(ase_atom.symbol, list(ase_atom.position))


def from_pmg_site(pmg_site):
    """Construct an Atom object from pymatgen Site object"""

    import pymatgen

    assert isinstance(
        pmg_site, pymatgen.core.sites.Site
    ), "not a pymatgen.core.sites.Site object!"
    return Atom(pmg_site.specie.symbol, list(pmg_site.coords))


def to_ase_atom(atom: Atom):
    """Construct an ase `Atom` object from native `Atom`"""
    import ase

    return ase.Atom(symbol=atom.symbol, position=atom.position)


def to_pmg_site(atom: Atom):
    """Construct a pymatgein `Site` object from native `Atom`"""

    import pymatgen

    ele = pymatgen.core.periodic_table.Element(atom.symbol)
    pmg_site = pymatgen.core.sites.Site(ele, coords=atom.position)
    return pmg_site


def to_pmg_periodic_site(atom: Atom, lattice):
    """Construct a pymatgen `PeriodicSite` object from native `Atom`"""

    import pymatgen

    ele = pymatgen.core.periodic_table.Element(atom.symbol)
    lat = pymatgen.core.lattice.Lattice(lattice)
    pmg_site = pymatgen.core.sites.PeriodicSite(
        ele, coords=atom.position, lattice=lat, coords_are_cartesian=True
    )
    return pmg_site
# fbe586a0 ends here

# [[file:../spdkit-python.note::ec59e65f][ec59e65f]]
def view_in_pymol(mol: Molecule, rebond=False, format="pdb"):
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
        # for work with rpyc remote Molecule object, we do not call
        # to_file directly
        fmt = io.guess_format_from_path(molfile)
        s = mol.format_as(fmt)
        open(molfile, "w").write(s)

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


def view_traj_in_pymol(mols: list[Molecule], animated=True, format="mol2", sleep=5):
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
            # for work with rpyc remote Molecule object, we do not call to_file directly
            fmt = io.guess_format_from_path(molfile)
            s = m.format_as(fmt)
            open(molfile, "w").write(s)

            if animated:
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
        # wait a few seconds for pymol reading temp files
        time.sleep(sleep)
        return p
# ec59e65f ends here
