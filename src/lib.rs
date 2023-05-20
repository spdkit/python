// [[file:../spdkit-python.note::0cbf1e93][0cbf1e93]]
use pyo3::prelude::*;
use pyo3::types::PyType;
// 0cbf1e93 ends here

// [[file:../spdkit-python.note::787fe451][787fe451]]
mod analysis;
mod apps;
mod gosh;
mod gui;
mod io;
mod surface;
mod utils;

// mod htc;
// 787fe451 ends here

// [[file:../spdkit-python.note::bac0eb10][bac0eb10]]
mod common {
    pub use gut::prelude::*;
}
// bac0eb10 ends here

// [[file:../spdkit-python.note::95be1618][95be1618]]
use gchemol::Atom;

#[pyclass(name = "Atom")]
#[derive(Clone)]
/// Atom is the smallest particle still characterizing a chemical element.
#[pyo3(text_signature = "(symbol, position=[0, 0, 0])")]
// #[pyo3(signature = (position = [0, 0, 0]))]
pub struct PyAtom {
    inner: Atom,
}

#[pymethods]
impl PyAtom {
    #[new]
    /// Construct `Atom` object from `symbol` and `position`.
    #[args(position = "[0.0, 0.0, 0.0]")]
    pub fn new(symbol: String, position: [f64; 3]) -> Self {
        Self {
            inner: Atom::new(symbol, position),
        }
    }

    /// Return a copy of `Atom`.
    pub fn clone(&self) -> Self {
        let inner = self.inner.clone();
        Self { inner }
    }

    /// Return element symbol
    #[getter]
    pub fn symbol(&self) -> String {
        self.inner.symbol().to_string()
    }

    /// Set atom symbol.
    #[setter(symbol)]
    pub fn set_symbol(&mut self, symbol: String) {
        self.inner.set_symbol(symbol);
    }

    /// Return atom label
    #[getter]
    pub fn label(&self) -> String {
        self.inner.get_label().unwrap_or(self.inner.symbol()).to_owned()
    }

    /// Set atom label.
    #[setter(label)]
    pub fn set_label(&mut self, label: String) {
        self.inner.set_label(label);
    }

    /// Return atomic number
    #[getter]
    pub fn number(&self) -> usize {
        self.inner.number()
    }

    /// Get mass in atomic mass unit. Return None if atom is dummy.
    pub fn get_mass(&self) -> Option<f64> {
        self.inner.get_mass()
    }

    /// Return atom position in 3D Cartesian coordinates
    #[getter]
    pub fn position(&self) -> [f64; 3] {
        self.inner.position()
    }

    /// Change atom Cartesian position to `p` ([x, y, z]).
    #[setter]
    pub fn set_position(&mut self, p: [f64; 3]) {
        self.inner.set_position(p);
    }

    /// Return freezing mask array for Cartesian coordinates
    #[getter]
    pub fn freezing(&self) -> [bool; 3] {
        self.inner.freezing()
    }

    #[setter]
    pub fn set_freezing(&mut self, freeze: [bool; 3]) {
        self.inner.set_freezing(freeze);
    }

    /// Access Van der Waals radius of atom. Return None if no data available
    pub fn get_vdw_radius(&self) -> Option<f64> {
        self.inner.get_vdw_radius()
    }

    /// Access covalent radius of atom. Return None if no data available or atom is dummy.
    pub fn get_cov_radius(&self) -> Option<f64> {
        self.inner.get_cov_radius()
    }
}
// 95be1618 ends here

// [[file:../spdkit-python.note::c8807c91][c8807c91]]
use gchemol::Lattice;

/// Periodic 3D lattice
#[derive(Clone)]
#[pyclass(name = "Lattice", subclass)]
#[pyo3(text_signature = "(tvs)")]
pub struct PyLattice {
    inner: Lattice,
}

#[pymethods]
impl PyLattice {
    /// Construct Lattice from lattice matrix (3x3).
    #[new]
    pub fn new(tvs: [[f64; 3]; 3]) -> Self {
        let inner = Lattice::new(tvs);
        Self { inner }
    }

    /// Construct lattice from lattice parameters Unit cell angles in degrees, lengths in Angstrom
    #[staticmethod]
    pub fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let inner = Lattice::from_params(a, b, c, alpha, beta, gamma);
        Self { inner }
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_frac(p).into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    pub fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_cart(p).into()
    }

    /// Wrap point `p` in Cartesian coordinates into unit cell,
    /// obeying the periodic boundary conditions. Returns cartesian
    /// coordinates.
    pub fn wrap(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.wrap(p).into()
    }

    /// Wrap point `p` in fractional coordinates into unit cell,
    /// obeying the periodic boundary conditions. Returns fractional
    /// coordinates.
    pub fn wrap_frac(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.wrap_frac(p).into()
    }

    /// Return the shortest vector obeying the minimum image convention.
    pub fn apply_mic(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.apply_mic(p).into()
    }

    /// Returns Lattice vector a.
    #[getter]
    pub fn vector_a(&self) -> [f64; 3] {
        self.inner.vector_a().into()
    }

    /// Returns Lattice vector b.
    #[getter]
    pub fn vector_b(&self) -> [f64; 3] {
        self.inner.vector_b().into()
    }

    /// Returns Lattice vector c.
    #[getter]
    pub fn vector_c(&self) -> [f64; 3] {
        self.inner.vector_c().into()
    }

    /// Return the shortest distance between pi (point i) and the
    /// periodic images of pj (point j) under the minimum image
    /// convention
    pub fn distance(&self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        self.inner.distance(pi, pj).into()
    }

    /// Lattice length parameters: a, b, c.
    pub fn lengths(&self) -> [f64; 3] {
        self.inner.lengths()
    }

    /// Lattice angles parameters: alpha, beta, gamma in degree.
    pub fn angles(&self) -> [f64; 3] {
        self.inner.angles()
    }

    /// Return the volume of the unit cell the cache will be updated
    /// if necessary
    pub fn volume(&self) -> f64 {
        self.inner.volume()
    }

    /// Scale Lattice by a positive constant `v`
    pub fn scale_by(&mut self, v: f64) {
        self.inner.scale_by(v)
    }

    /// Scale Lattice in `a` direction by a positive constant `v`
    pub fn scale_by_a(&mut self, v: f64) {
        self.inner.scale_by_a(v)
    }

    /// Scale Lattice in `b` direction by a positive constant `v`
    pub fn scale_by_b(&mut self, v: f64) {
        self.inner.scale_by_b(v)
    }

    /// Scale Lattice in `c` direction by a positive constant `v`
    pub fn scale_by_c(&mut self, v: f64) {
        self.inner.scale_by_c(v)
    }
}
// c8807c91 ends here

// [[file:../spdkit-python.note::ef13f019][ef13f019]]
use gchemol::neighbors::{Neighbor, Neighborhood};

#[derive(FromPyObject)]
pub enum Selection {
    #[pyo3(transparent, annotation = "list")]
    List(Vec<usize>),
    #[pyo3(transparent, annotation = "str")]
    String(String),
}

impl Selection {
    fn try_into_list(self) -> Result<Vec<usize>> {
        let selection = match self {
            Selection::List(sel) => sel,
            Selection::String(s) => {
                let sel = gut::utils::parse_numbers_human_readable(&s)?;
                info!("parsed user selection: {sel:?}");
                sel
            }
        };
        Ok(selection)
    }
}

/// Represent a probe for searching nearby neighbors within distance
/// cutoff.
#[pyclass(name = "Neighborhood", subclass)]
pub struct PyNeighborProbe {
    inner: Neighborhood,
}

/// Helper struct for neighbors search result.
#[pyclass(name = "Neighbor", subclass)]
pub struct PyNeighbor {
    inner: Neighbor,
}

#[pymethods]
impl PyNeighbor {
    #[getter(node)]
    /// The node connected to the host point.
    pub fn node(&self) -> usize {
        self.inner.node
    }

    #[getter(distance)]
    ///The distance to the host point.
    pub fn distance(&self) -> f64 {
        self.inner.distance
    }

    #[getter(image)]
    /// Scaled displacment vector relative to origin cell if PBC enabled.
    pub fn image(&self) -> Option<[f64; 3]> {
        let image = self.inner.image?;
        Some(image.into())
    }
}

#[pymethods]
impl PyNeighborProbe {
    /// Return neighbors of a particle `p` within distance cutoff `r_cutoff`.
    ///
    /// # Parameters
    /// * p: the Cartesian position or probe
    /// * r_cutoff: the cutoff distance for searching neighbors
    pub fn search(&self, p: [f64; 3], r_cutoff: f64) -> Vec<PyNeighbor> {
        self.inner
            .search(p, r_cutoff)
            .map(|inner| PyNeighbor { inner })
            .collect()
    }

    /// Return an iterator of the nodes connected to the node `n`.
    pub fn neighbors(&self, n: usize, r_cutoff: f64) -> Vec<PyNeighbor> {
        self.inner.neighbors(n, r_cutoff).map(|inner| PyNeighbor { inner }).collect()
    }
}
// ef13f019 ends here

// [[file:../spdkit-python.note::969a9313][969a9313]]
use gchemol::prelude::*;
use gchemol::Molecule;
use gut::prelude::*;

/// The main object featured in this library. This object represents a
/// molecule, with atoms and bonds.
#[pyclass(name = "Molecule", subclass)]
#[pyo3(text_signature = "(title)")]
#[derive(Clone)]
pub struct PyMolecule {
    inner: Molecule,
}

#[pymethods]
impl PyMolecule {
    #[new]
    /// Construct an empty `Molecule` object named as `title`.
    pub fn new(title: String) -> Self {
        let inner = Molecule::new(&title);
        Self { inner }
    }

    /// Create a new copy of `Molecule.`
    pub fn clone(&self) -> Self {
        let inner = self.inner.clone();
        Self {inner}
    }

    /// Return the name of the molecule, which is typpically modified
    /// for safely storing in various chemical file formats.
    #[getter]
    pub fn title(&self) -> String {
        self.inner.title()
    }

    /// Set molecule's title to `title`.
    #[setter]
    pub fn set_title(&mut self, title: String) {
        self.inner.set_title(&title);
    }

    #[staticmethod]
    /// Build a molecule from atoms associated with serial numbers from 1.
    pub fn from_atoms(atoms: Vec<PyAtom>) -> Self {
        let inner = Molecule::from_atoms(atoms.into_iter().map(|m| m.inner));
        Self { inner }
    }

    /// Return its json representation of molecule object.
    pub fn to_json(&self) -> PyResult<String> {
        let json = gchemol::io::to_json(&self.inner)?;
        Ok(json)
    }

    /// Replace atom `i` with new `atom`.
    pub fn set_atom(&mut self, i: usize, atom: PyAtom) -> PyResult<()> {
        let a = self.inner.get_atom_mut(i).ok_or(format_err!("no atom {i}"))?;
        *a = atom.inner;
        Ok(())
    }

    /// Set label of atom `i`.
    pub fn set_atom_label(&mut self, i: usize, label: String) -> PyResult<()> {
        let a = self.inner.get_atom_mut(i).ok_or(format_err!("no atom {i}"))?;
        a.set_label(label);
        Ok(())
    }

    /// Get label of atom `i`.
    pub fn get_atom_label(&self, i: usize) -> PyResult<String> {
        let a = self.inner.get_atom(i).ok_or(format_err!("no atom {i}"))?;
        let s = a.get_label().unwrap_or_default();
        Ok(s.to_owned())
    }

    /// Set positions of atoms in sequential order.
    pub fn set_positions(&mut self, positions: Vec<[f64; 3]>) {
        self.inner.set_positions(positions);
    }

    /// Set fractional coordinates of atoms in sequence order.
    /// Panics if Molecule is aperiodic.
    pub fn set_scaled_positions(&mut self, frac_coords: Vec<[f64; 3]>) {
        self.inner.set_scaled_positions(frac_coords);
    }

    /// Return atom positions ordered by serial numbers.
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.inner.positions().collect()
    }

    /// Clean up molecule geometry using stress majorization algorithm.
    pub fn clean(&mut self) -> PyResult<()> {
        self.inner.clean()?;
        Ok(())
    }

    /// Renumber atoms consecutively from 1.
    pub fn renumber(&mut self) {
        self.inner.renumber();
    }

    /// Format molecule as string in specific molecular file format. Return
    /// error if cannot format molecule in `fmt`.
    pub fn format_as(&self, format: String) -> PyResult<String> {
        let s = self.inner.format_as(&format)?;
        Ok(s)
    }

    /// Write molecule to file with `path`. The molecule format will
    /// be determined based on file name extension.
    pub fn to_file(&self, path: String) -> PyResult<()> {
        self.inner.to_file(&path)?;
        Ok(())
    }

    /// Render molecule with template file from `path`. On success,
    /// return the formatted string.
    pub fn render_with(&self, path: String) -> PyResult<String> {
        let s = self.inner.render_with(path.as_ref())?;
        Ok(s)
    }

    /// Get the number of atoms.
    pub fn natoms(&self) -> PyResult<usize> {
        let n = self.inner.natoms();
        Ok(n)
    }

    /// Return chemical formula.
    pub fn formula(&self) -> PyResult<String> {
        Ok(self.inner.formula())
    }

    /// Unbuild current crystal structure leaving a nonperiodic
    /// structure. Return lattice if periodic, otherwise return None.
    pub fn unbuild_crystal(&mut self) -> Option<PyLattice> {
        let lat = self.inner.unbuild_crystal();
        lat.map(|inner| PyLattice {inner})
    }

    /// Create a Lattice from the minimal bounding box of the Molecule
    /// extended by a positive value of padding.
    ///
    /// NOTE: padding has to be large enough (> 0.5) to avoid self
    /// interaction with its periodic mirror.
    pub fn set_lattice_from_bounding_box(&mut self, padding: f64) -> PyResult<()> {
        self.inner.set_lattice_from_bounding_box(padding);
        Ok(())
    }

    /// Create a supercell version of new molecule.
    ///
    /// # Arguments
    ///
    /// * sa, sb, sc: an sequence of three scaling factors. E.g., [2,
    /// 1, 1] specifies that the supercell should have dimensions 2a x
    /// b x c
    pub fn supercell(&mut self, sa: usize, sb: usize, sc: usize) -> PyResult<Self> {
        if let Some(mol) = self.inner.supercell(sa, sb, sc) {
            let m = Self { inner: mol };
            Ok(m)
        } else {
            Err(pyo3::exceptions::PyException::new_err("not allowed for periodic structure"))
        }
    }

    /// Return the distance between atom i and atom j. For periodic
    /// structure, this method will return the distance under the
    /// minimum image convention.
    pub fn distance(&self, i: usize, j: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_distance(i, j)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}"))?;
        Ok(d)
    }

    /// Return the angle between atoms `i`, `j`, `k` in degrees,
    /// irrespective periodic images.
    pub fn angle(&self, i: usize, j: usize, k: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_angle(i, j, k)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}, {k}"))?
            .to_degrees();
        Ok(d)
    }

    /// Return the torsion angle between atoms `i`, `j`, `k`, `l` in
    /// degrees, irrespective periodic images.
    pub fn torsion(&self, i: usize, j: usize, k: usize, l: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_torsion(i, j, k, l)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}, {k}, {l}"))?
            .to_degrees();
        Ok(d)
    }

    /// Return molecule’s inertia matrix (3x3) in reference to molecule’s center of mass
    pub fn inertia_matrix(&self) -> PyResult<[[f64; 3]; 3]> {
        let im = self.inner.inertia_matrix();
        Ok(im)
    }

    /// Return the center of mass of molecule (COM).
    pub fn center_of_mass(&self) -> PyResult<[f64; 3]> {
        let com = self.inner.center_of_mass();
        Ok(com)
    }

    /// Return the center of geometry of molecule (COG).
    pub fn center_of_geometry(&self) -> PyResult<[f64; 3]> {
        let cog = self.inner.center_of_geometry();
        Ok(cog)
    }

    /// Set freezing flag to `freezed` for atoms in `selection`,
    /// useful for optimization or dynamic simulation.
    pub fn freeze_atoms(&mut self, selection: Selection, freezed: bool) -> PyResult<()> {
        for sn in selection.try_into_list()? {
            if let Some(a) = self.inner.get_atom_mut(sn) {
                a.set_freezing([freezed; 3]);
            } else {
                return Err(pyo3::exceptions::PyException::new_err("atom {sn} not found"));
            }
        }
        Ok(())
    }

    /// Get ONIOM layer of atom `sn`. If no layer information was set,
    /// will return `H`.
    pub fn get_oniom_layer(&self, sn: usize) -> PyResult<String> {
        use gchemol::io::formats::GaussianInputFile;

        if let Some(a) = self.inner.get_atom(sn) {
            let extra = GaussianInputFile::extra_atom_info(a);
            let l = extra.get_oniom_layer().unwrap_or("H");
            Ok(l.to_owned())
        } else {
            Err(pyo3::exceptions::PyException::new_err("atom {sn} not found"))
        }
    }

    /// Set ONIOM layer of atom `sn` to `layer` for Gaussian
    /// calculation. Possible layer includes `H`, `M`, `L`
    pub fn set_oniom_layer(&mut self, sn: usize, layer: &str) -> PyResult<()> {
        use gchemol::io::formats::GaussianInputFile;

        if let Some(a) = self.inner.get_atom_mut(sn) {
            let mut extra = GaussianInputFile::extra_atom_info(a);
            extra.set_oniom_layer(layer);
            extra.attach(a);
            Ok(())
        } else {
            Err(pyo3::exceptions::PyException::new_err("atom {sn} not found"))
        }
    }

    /// Find rings up to `nmax` atoms in `Molecule`.
    pub fn find_rings(&mut self, nmax: usize) -> PyResult<Vec<std::collections::HashSet<usize>>> {
        let rings = self.inner.find_rings(nmax);
        Ok(rings)
    }

    /// Return atom serial numbers.
    pub fn numbers(&self) -> Vec<usize> {
        self.inner.numbers().collect()
    }

    /// Return atom symbols
    pub fn symbols(&self) -> Vec<String> {
        self.inner.symbols().map(|x| x.to_owned()).collect()
    }

    /// Return atomic numbers
    pub fn atomic_numbers(&self) -> Vec<usize> {
        self.inner.atomic_numbers().collect()
    }

    /// Return an iterator over a tuple of atom serial number `n` and
    /// its associated `Atom` (n, Atom)
    pub fn atoms(&self) -> PyAtomsIter {
        let atoms: Vec<_> = self
            .inner
            .atoms()
            .map(|(i, atom)| (i, PyAtom { inner: atom.to_owned() }))
            .collect();
        PyAtomsIter {
            inner: atoms.into_iter(),
        }
    }

    /// Access the atom copy by atom serial number `n`. Return None if
    /// serial number `n` invalid.
    pub fn get_atom(&self, n: usize) -> Option<PyAtom> {
        let inner = self.inner.get_atom(n)?.clone();
        PyAtom { inner }.into()
    }
    
    /// Add atom a into molecule. If Atom numbered as a already exists in
    /// molecule, then the associated Atom will be updated with atom.
    pub fn add_atom(&mut self, n: usize, atom: PyAtom) {
        self.inner.add_atom(n, atom.inner)
    }
    
    /// Remove Atom a from Molecule. Return the removed Atom on success,
    /// and return None if Atom a does not exist.
    pub fn remove_atom(&mut self, n: usize) -> Option<PyAtom> {
        let inner = self.inner.remove_atom(n)?;
        PyAtom { inner }.into()
    }
    
    /// Remove atoms in `selection` from Molecule. Return the removed
    /// atoms on success, and return None if any atom does not exist.
    pub fn remove_atoms(&mut self, selection: Selection) -> Option<Vec<PyAtom>> {
        let mut removed = vec![];
        for n in selection.try_into_list().ok()? {
            let inner = self.inner.remove_atom(n)?;
            let a = PyAtom { inner };
            removed.push(a);
        }
        Some(removed)
    }

    /// Return a sub molecule induced by `atoms` in parent
    /// molecule. Return None if atom serial numbers are
    /// invalid. Return an empty Molecule if `atoms` empty.
    pub fn get_sub_molecule(&self, atoms: Selection) -> Option<Self> {
        let atoms = atoms.try_into_list().ok()?;
        let inner = self.inner.get_sub_molecule(&atoms)?;
        Self { inner }.into()
    }

    /// Set periodic lattice.
    pub fn set_lattice(&mut self, lat: PyLattice) {
        self.inner.set_lattice(lat.inner);
    }

    /// Get periodic lattice.
    pub fn get_lattice(&self) -> Option<PyLattice> {
        let lat = self.inner.get_lattice()?;
        PyLattice { inner: *lat }.into()
    }

    /// Return true if Molecule is a periodic structure.
    pub fn is_periodic(&self) -> bool {
        self.inner.is_periodic()
    }

    /// Return fractional coordinates relative to unit cell. Return
    /// None if not a periodic structure.
    pub fn get_scaled_positions(&self) -> Option<Vec<[f64; 3]>> {
        let scaled = self.inner.get_scaled_positions()?.collect_vec();
        scaled.into()
    }

    /// Get the number of bonds.
    pub fn nbonds(&self) -> usize {
        self.inner.nbonds()
    }
    
    /// Return the shortest distance counted in number of chemical
    /// bonds between two atoms. Return None if they are not
    /// connected.
    pub fn nbonds_between(&self, i: usize, j: usize) -> PyResult<Option<usize>> {
        let n = self.inner.nbonds_between(i, j);
        Ok(n)
    }
    
    #[pyo3(signature = (ignore_pbc=false, bond_tolerance=0.45))]
    #[pyo3(text_signature = "(ignore_pbc=False, bond_tolerance=0.45)")]
    /// Recalculates all bonds in molecule based on interatomic
    /// distances and covalent radii. For periodic system, the bonds
    /// are determined by applying miniumu image convention.
    pub fn rebond(&mut self, ignore_pbc: bool, bond_tolerance: f64) {
        std::env::set_var("GCHEMOL_REBOND_IGNORE_PBC", format!("{ignore_pbc}"));
        std::env::set_var("GCHEMOL_REBOND_BOND_TOLERANCE", format!("{bond_tolerance}"));
        self.inner.rebond();
    }
    
    /// Removes all existing bonds between atoms.
    pub fn unbound(&mut self) -> PyResult<()> {
        self.inner.unbound();
        Ok(())
    }
    
    /// Return a list of bonds in tuples (u, v, bond_order)
    pub fn bonds(&mut self) -> Vec<(usize, usize, f64)> {
        self.inner
            .bonds()
            .map(|(u, v, b)| if u > v { (u, v, b.order()) } else { (v, u, b.order()) })
            .collect()
    }
    
    /// Removes all bonds between atoms in `selection1` and `selection2`,
    /// in respect of pymol's `unbond` command.
    ///
    /// # Parameters
    /// * selection1: the first collection of atoms
    /// * selection2: the other collection of atoms
    ///
    /// # Reference
    /// * <https://pymolwiki.org/index.php/Unbond>
    #[pyo3(text_signature = "($self, selection1, selection2, /)")]
    pub fn unbond(&mut self, selection1: Selection, selection2: Selection) -> PyResult<()> {
        let atom_indices1 = selection1.try_into_list()?;
        let atom_indices2 = selection2.try_into_list()?;
        self.inner.unbond(&atom_indices1, &atom_indices2);
        Ok(())
    }
    
    /// Add a single bond between Atom `i` and Atom `j` into molecule.
    /// Panic if the specified atom a or b does not exist.
    ///
    /// # Parameters
    /// * kind: set bond kind using a string type. Possible values: dummy,
    ///   partial, single, aromatic, double, triple, quadruple
    #[pyo3(text_signature = "($self, i, j, /, kind = None)")]
    pub fn add_bond(&mut self, i: usize, j: usize, kind: Option<String>) {
        use gchemol::Bond;
    
        let bond = kind
            .map(|s| match s.as_str() {
                "dummy" => Bond::dummy(),
                "partial" => Bond::partial(),
                "single" => Bond::single(),
                "aromatic" => Bond::aromatic(),
                "double" => Bond::double(),
                "triple" => Bond::triple(),
                "quadruple" => Bond::quadruple(),
                _ => panic!("invalid bond kind {s}"),
            })
            .unwrap_or(Bond::single());
        self.inner.add_bond(i, j, bond)
    }
    
    /// Remove the bond between atom `i` and atom `j`.
    #[pyo3(text_signature = "($self, i, j, /)")]
    pub fn remove_bond(&mut self, i: usize, j: usize) {
        use gchemol::Bond;
    
        self.inner.remove_bond(i, j);
    }
    
    /// Return all directly bonded atoms in serial numbers with atom `n`.
    #[pyo3(text_signature = "($self, n)")]
    pub fn connected(&self, n: usize) -> Vec<usize> {
        self.inner.connected(n).collect()
    }

    /// Break molecule into multiple fragments based on its bonding
    /// connectivity.
    pub fn fragmented(&self) -> Vec<Self> {
        self.inner.fragmented().map(|inner| Self { inner }).collect()
    }
    
    /// Return all atoms that connected in the same fragment as atom
    /// `i`.
    pub fn connected_fragment_atoms(&self, i: usize) -> Vec<usize> {
        self.inner.connected_fragment_atoms(i).collect()
    }
    
    /// Reorder the atoms according to the ordering of keys. Keys define
    /// 1-to-1 mapping of atoms.
    ///
    /// # Parameters
    /// * keys: a list of numbers for sorting
    ///
    /// # NOTE
    /// * This method will cause serial numbers renumbered from 1.
    pub fn reorder(&mut self, keys: Vec<usize>) {
        self.inner.reorder(&keys);
    }
    
    /// Superimpose structure onto `mol_ref` which will be held fixed
    /// during alignment. Return superposition rmsd on done.
    ///
    /// # NOTE
    /// * The atoms must be in one-to-one mapping with atoms in `mol_ref`
    /// * The structure could be mirrored for better alignment.
    /// * Heavy atoms have more weights.
    #[pyo3(signature = (mol_ref, /, selection = None))]
    #[pyo3(text_signature = "($self, mol_ref, /, selection = None)")]
    pub fn superimpose_onto(&mut self, mol_ref: Self, selection: Option<Selection>) -> PyResult<f64> {
        let rmsd = if let Some(selected) = selection {
            let s = selected.try_into_list()?;
            self.inner.superimpose_onto(&mol_ref.inner, Some(&s))
        } else {
            self.inner.superimpose_onto(&mol_ref.inner, None)
        };
    
        Ok(rmsd)
    }

    pub fn educated_rebond(&mut self) {
        use educate::prelude::*;
        self.inner.educated_rebond();
    }
    
    pub fn educated_clean(&mut self) {
        use educate::prelude::*;
        self.inner.educated_clean();
    }
    
    pub fn educated_clean_selection(&mut self, selection: Selection) -> PyResult<()> {
        use educate::prelude::*;
        let selection = selection.try_into_list()?;
        self.inner.educated_clean_selection(&selection);
        Ok(())
    }
    
    /// Return a unique fingerprint of current molecule based on its bond
    /// graph. This fingerprint is independent of its 3d geometry or atom
    /// numbering.
    ///
    /// # NOTE
    ///   * This operation internally call `reorder_cannonically` method.
    pub fn fingerprint(&self) -> String {
        use spdkit::prelude::FingerPrintExt;
        self.inner.fingerprint()
    }
    
    /// Calculate disparity between `self` and `mol` using algorithm
    /// proposed by Lazauskas et al (DOI:10.1039/C6NR09072A)
    pub fn disparity_between(&self, mol: &PyMolecule) -> f64 {
        use spdkit::prelude::SimilarityExt;
        self.inner.disparity_between(&mol.inner)
    }
    
    /// This is an operation of reordering the atoms in a way that does not depend
    /// on where they were before. The bonding graph is important for this
    /// operation.
    pub fn reorder_cannonically(&mut self) -> (Vec<usize>, Vec<usize>) {
        use spdkit::prelude::FingerPrintExt;
        self.inner.reorder_cannonically()
    }
    
    /// Make `self` resemble `mol_ref` by applying rigid operations in
    /// permutation, translation or rotation, without changing inner
    /// 3D geometry. Equivalent atoms are recoginized based on
    /// connectivity. Return alignment rmsd on success.
    pub fn resemble_rigidly(&mut self, mol_ref: PyMolecule) -> Result<f64> {
        use spdkit::prelude::FingerPrintExt;
        self.inner.resemble_rigidly(&mol_ref.inner)
    }
    
    /// Convert `Molecule` to a graph object for distance geometry
    /// refinement.
    pub fn to_distance_geometry_graph(&self) -> DgGraph {
        let dg = self.inner.distance_geometry_graph();
        DgGraph { inner: dg }
    }

    /// Construct `Molecule` object from a file `path`. If the file
    /// contains multiple molecules, only the last one will be read.
    #[staticmethod]
    pub fn from_file(path: String) -> PyResult<Self> {
        let inner = Molecule::from_file(&path)?;
        // let instance = cls.call0()?;
        // let mol: Py<Self> = instance.extract()?;
        Ok(Self { inner })
    }
    
    /// Construct `Molecule` from string `src` in specific molecular file
    /// format `fmt`. Possible fmt includes `text/xyz`, `text/pxyz`,
    /// `vasp/input` etc.
    #[staticmethod]
    pub fn from_string(src: String, fmt: String) -> PyResult<Self> {
        let inner = Molecule::from_str(&src, &fmt)?;
        Ok(Self { inner })
    }
    
    #[staticmethod]
    /// Describe available backends for reading/writing molecule
    pub fn describe_available_backends() {
        gchemol::io::describe_backends();
    }

    /// Return selected atoms by expanding `n` chemical bonds away from
    /// the center atom `m`
    ///
    /// Note: the center atom m is put last in returned molecule.
    fn selection_by_expanding_bond(&self, m: usize, n: usize) -> Vec<usize> {
        self.inner.selection_by_expanding_bond(m, n)
    }
    
    /// Return selected atoms by cutoff distance `r` nearby central atom `n`
    fn selection_by_distance(&self, n: usize, r: f64) -> Vec<usize> {
        self.inner.selection_by_distance(n, r)
    }
    
    /// Return a `Neighborhood` struct for finding nearest neighbors.
    fn create_neighbor_probe(&self) -> PyNeighborProbe {
        let inner = self.inner.create_neighbor_probe();
        PyNeighborProbe { inner }
    }
    
    /// Return atoms with xyz coordinates freezed
    fn selection_freezed_atoms(&self) -> Vec<usize> {
        self.inner
            .atoms()
            .filter_map(|(i, a)| if a.is_fixed() { Some(i) } else { None })
            .collect()
    }
}
// 969a9313 ends here

// [[file:../spdkit-python.note::a31e85a4][a31e85a4]]
#[pyclass]
pub struct PyAtomsIter {
    inner: std::vec::IntoIter<(usize, PyAtom)>,
}

#[pymethods]
impl PyAtomsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(usize, PyAtom)> {
        slf.inner.next()
    }
}

#[pyclass]
pub struct PyMoleculeIter {
    iter: Box<dyn Iterator<Item = PyMolecule> + Send>,
}

#[pymethods]
impl PyMoleculeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyMolecule> {
        slf.iter.next()
    }
}
// a31e85a4 ends here

// [[file:../spdkit-python.note::df84f7ba][df84f7ba]]
use distances::prelude::*;

/// Represents a graph object for refinement of molecule structure
/// using distance geometry algorithm.
#[pyclass(mapping, module = "spdkit", subclass)]
#[derive(Clone)]
pub struct DgGraph {
    inner: distances::DgGraph,
}

#[pymethods]
impl DgGraph {
    /// Set weight for atom `i` and `j` for refinement of molecule structure.
    pub fn set_distance_weight(&mut self, i: usize, j: usize, w: f64) {
        self.inner.set_distance_weight(i, j, w);
    }

    /// Set distance bound (lower/upper) for atom `i` and `j` for
    /// refinement of molecule structure.
    pub fn set_distance_bound(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        self.inner.set_distance_bound(i, j, lb, ub);
    }

    /// Returns distance bound (lower/upper) for atom `i` and `j`.
    pub fn get_distance_bound(&self, i: usize, j: usize) -> (f64, f64) {
        self.inner.distance_bound(i, j)
    }

    /// Returns distance weight for atom `i` and `j`.
    pub fn get_distance_weight(&self, i: usize, j: usize) -> f64 {
        self.inner.distance_weight(i, j)
    }

    /// Constrain pairwise distances of atoms `selection` in Molecule `mol`.
    #[pyo3(text_signature = "($self, mol, selection)")]
    pub fn constrain_distances_within(&mut self, mol: PyMolecule, selection: Selection) -> PyResult<()> {
        let selection = selection.try_into_list()?;
        self.inner.constrain_distances_within(&mol.inner, &selection, 10.0);
        Ok(())
    }

    /// Refine molecule structure `mol` using distance geometry.
    #[pyo3(text_signature = "($self, mol)")]
    pub fn refine_molecule(&mut self, mol: &mut PyMolecule) {
        self.inner.refine_molecule(&mut mol.inner);
    }

    /// Refine molecule structure `mol` for atoms in selection only
    /// using distance geometry.
    #[pyo3(text_signature = "($self, mol, selection)")]
    pub fn refine_molecule_selection(&mut self, mol: &mut PyMolecule, selection: Selection) -> PyResult<()> {
        let selection = selection.try_into_list()?;
        self.inner.refine_molecule_selection(&mut mol.inner, &selection);
        Ok(())
    }
}
// df84f7ba ends here

// [[file:../spdkit-python.note::db9c632d][db9c632d]]
use distances::Interpolation;

/// Structure interpolation in distance space.
///
/// # Parameters
/// * r_cut: probe cutoff radius for nearest neighbours.
#[derive(Clone)]
#[pyclass(name = "Interpolation", subclass)]
#[pyo3(text_signature = "(r_cut=4.0)")]
pub struct PyInterpolation {
    r_cut: f64,
}

#[pymethods]
impl PyInterpolation {
    #[new]
    #[pyo3(signature = (r_cut = 4.0))]
    pub fn new(r_cut: f64) -> Self {
        Self { r_cut }
    }

    /// Interpolate a new Molecue at position `f` between `mol1` and
    /// `mol2`. `f` should be a positive number between 0 and 1.
    pub fn interpolate(&self, mol1: PyMolecule, mol2: PyMolecule, f: f64) -> PyMolecule {
        let inner = Interpolation::new(&mol1.inner, &mol2.inner)
            .with_probe_radius(self.r_cut)
            .interpolate(f);
        PyMolecule { inner }
    }

    /// Return a `Molecule` as the weighted average of `mol1` and `mol2` in
    /// distance space.
    ///
    /// # Parameters
    /// * w1, w2: the weights for `mol1` and `mol2`, respectively.
    pub fn weighted_average(&self, mol1: PyMolecule, mol2: PyMolecule, w1: f64, w2: f64) -> PyMolecule {
        let inner = Interpolation::new(&mol1.inner, &mol2.inner).weighted_average(w1, w2);
        PyMolecule { inner }
    }
}
// db9c632d ends here

// [[file:../spdkit-python.note::88853a11][88853a11]]
use distances::docs::dg_periodic_probe::ChemicalEnvironment;

#[pyclass(name = "ChemicalEnvironment", subclass)]
#[pyo3(text_signature = "(mol, rcut)")]
/// New ChemicalEnvironment for Molecule mol, without elastic connections.
pub struct PyChemicalEnvironment {
    mol: Molecule,
    inner: ChemicalEnvironment,
}

#[pymethods]
impl PyChemicalEnvironment {
    #[new]
    /// Define ChemicalEnvironment by probing local neighbors within
    /// distance cutoff r_cut
    pub fn probe(mol: PyMolecule, rcut: f64) -> Self {
        let mol = mol.inner;
        let inner = ChemicalEnvironment::probe(&mol, rcut);

        Self { mol, inner }
    }

    /// Reshape the structure of Molecule mol using probed coordination environment.
    pub fn reshape(&self, mol: &mut PyMolecule) -> usize {
        self.inner.reshape(&mut mol.inner)
    }

    /// Set parent `Molecule` to `mol`.  Molecule `mol` may have
    /// differnt 3D structure, but must share the same numbering
    /// system with its parent molecule.
    pub fn reparent(&mut self, mol: PyMolecule) {
        self.inner.reparent(mol.inner)
    }

    /// Update elastic constrains from `mol`, which maintaining original connectivity.
    pub fn update_constrains_from(&mut self, mol: &PyMolecule) -> Result<()> {
        self.inner.update_constrains_from(&mol.inner)?;
        Ok(())
    }

    /// Return distance constraint between node `u` and `v`.
    pub fn get_distance_constraint(&self, u: usize, v: usize) -> Option<f64> {
        self.inner.get_distance_constraint(u, v)
    }

    /// Constrain distance `u--v` to `d` when refinement.
    pub fn set_distance_constraint(&mut self, u: usize, v: usize, d: f64) {
        self.inner.set_distance_constraint(u, v, d)
    }

    /// Create central molecule from atom `i` with direct neighbors.
    pub fn create_central_molecule(&self, i: usize, mol_alt: Option<PyMolecule>) -> Option<PyMolecule> {
        let inner = self.inner.create_central_molecule(i, mol_alt.map(|m| m.inner).as_ref())?;
        let m = PyMolecule { inner };
        m.into()
    }

    /// Create a real auxiliary molecule with all periodic atoms.
    pub fn create_auxiliary_molecule(&self) -> Option<PyMolecule> {
        let inner = self.inner.create_auxiliary_molecule()?;
        PyMolecule { inner }.into()
    }

    /// Update immediate constraints of atom `i` using distances
    /// obtained from `mol_a`.
    ///
    /// # NOTE
    /// - `mol_b` must be a sub molecule of its parent, and has the same numbering system.
    /// - `mol_a` assumed in the same topology as `mol_b`
    /// - `mol_a  may have a different numbering system.
    pub fn update_constraints_for_center_from(&mut self, mol_a: PyMolecule, mol_b: PyMolecule, i: usize) -> Result<()> {
        self.inner.update_constraints_for_center_from(&mol_a.inner, &mol_b.inner, i)?;
        Ok(())
    }
}
// 88853a11 ends here

// [[file:../spdkit-python.note::fbe87af8][fbe87af8]]
#[pymodule]
#[pyo3(name = "spdkit")]
fn pyspdkit(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyMolecule>()?;
    m.add_class::<PyAtom>()?;
    m.add_class::<PyLattice>()?;
    m.add_class::<PyInterpolation>()?;

    // utils
    let s = utils::new(py, "utils")?;
    m.add_submodule(s)?;

    // surface
    let s = surface::new(py, "surface")?;
    m.add_submodule(s)?;

    // io
    let s = io::new(py, "io")?;
    m.add_submodule(s)?;

    // gosh, database
    let s = gosh::new(py, "gosh")?;
    m.add_submodule(s)?;

    // apps, applications
    let s = apps::new(py, "apps")?;
    m.add_submodule(s)?;

    // gui: visualization
    let s = gui::new(py, "gui")?;
    m.add_submodule(s)?;

    // analysis: bond valence, atom valence, general atom valence, ...
    let s = analysis::new(py, "analysis")?;
    m.add_submodule(s)?;

    // for ad-hoc experiments
    let dwim = PyModule::new(py, "dwim")?;
    dwim.add_class::<PyChemicalEnvironment>()?;
    m.add_submodule(dwim)?;

    Ok(())
}
// fbe87af8 ends here
