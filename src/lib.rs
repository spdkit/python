// [[file:../spdkit-python.note::0cbf1e93][0cbf1e93]]
use pyo3::prelude::*;
use pyo3::types::PyType;
// 0cbf1e93 ends here

// [[file:../spdkit-python.note::95be1618][95be1618]]
use gchemol::Atom;

#[pyclass(subclass)]
#[derive(Clone)]
pub struct _Atom {
    inner: Atom,
}

#[pymethods]
impl _Atom {
    #[new]
    /// Construct `Atom` object from `symbol` and `position`.
    fn new(symbol: String, position: [f64; 3]) -> Self {
        Self {
            inner: Atom::new(symbol, position),
        }
    }

    /// Return element symbol
    fn symbol(&self) -> String {
        self.inner.symbol().to_string()
    }

    /// Return atomic number
    fn number(&self) -> usize {
        self.inner.number()
    }

    /// Get mass in atomic mass unit. Return None if atom is dummy.
    fn get_mass(&self) -> Option<f64> {
        self.inner.get_mass()
    }

    /// Return atom position in 3D Cartesian coordinates
    fn position(&self) -> [f64; 3] {
        self.inner.position()
    }

    fn set_position(&mut self, p: [f64; 3]) {
        self.inner.set_position(p);
    }

    /// Return freezing mask array for Cartesian coordinates
    fn freezing(&self) -> [bool; 3] {
        self.inner.freezing()
    }

    /// Access Van der Waals radius of atom. Return None if no data available
    fn get_vdw_radius(&self) -> Option<f64> {
        self.inner.get_vdw_radius()
    }

    /// Access covalent radius of atom. Return None if no data available or atom is dummy.
    fn get_cov_radius(&self) -> Option<f64> {
        self.inner.get_cov_radius()
    }
}
// 95be1618 ends here

// [[file:../spdkit-python.note::c8807c91][c8807c91]]
use gchemol::Lattice;

/// Periodic 3D lattice
#[pyclass(subclass)]
#[derive(Clone)]
pub struct _Lattice {
    inner: Lattice,
}

#[pymethods]
impl _Lattice {
    /// Construct Lattice from lattice matrix (3x3).
    #[new]
    fn new(tvs: [[f64; 3]; 3]) -> Self {
        let inner = Lattice::new(tvs);
        Self { inner }
    }

    /// Construct lattice from lattice parameters Unit cell angles in degrees, lengths in Angstrom
    #[staticmethod]
    fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let inner = Lattice::from_params(a, b, c, alpha, beta, gamma);
        Self { inner }
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    fn to_frac(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_frac(p).into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_cart(p).into()
    }

    /// Returns Lattice vector a.
    fn vector_a(&self) -> [f64; 3] {
        self.inner.vector_a().into()
    }

    /// Returns Lattice vector b.
    fn vector_b(&self) -> [f64; 3] {
        self.inner.vector_b().into()
    }

    /// Returns Lattice vector c.
    fn vector_c(&self) -> [f64; 3] {
        self.inner.vector_c().into()
    }

    /// Return the shortest vector obeying the minimum image convention.
    fn apply_mic(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.apply_mic(p).into()
    }

    /// Return the shortest distance between pi (point i) and the
    /// periodic images of pj (point j) under the minimum image
    /// convention
    fn distance(&self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        self.inner.distance(pi, pj).into()
    }

    /// Lattice length parameters: a, b, c.
    fn lengths(&self) -> [f64; 3] {
        self.inner.lengths().into()
    }

    /// Return the volume of the unit cell the cache will be updated
    /// if necessary
    fn volume(&self) -> f64 {
        self.inner.volume()
    }

}
// c8807c91 ends here

// [[file:../spdkit-python.note::969a9313][969a9313]]
use gchemol::prelude::*;
use gchemol::Molecule;
use gut::prelude::*;

/// The main object featured in this library. This object represents a
/// molecule, with atoms and bonds.
#[pyclass(subclass)]
#[derive(Clone)]
pub struct _Molecule {
    inner: Molecule,
}

#[pymethods]
impl _Molecule {
    /// Construct `Molecule` object from a file `path`. If the file
    /// contains multiple molecules, only the last one will be read.
    #[staticmethod]
    fn from_file(path: String) -> PyResult<Self> {
        let inner = Molecule::from_file(&path)?;
        Ok(Self { inner })
    }

    #[new]
    fn from_atoms(atoms: Vec<_Atom>) -> Self {
        let inner = Molecule::from_atoms(atoms.into_iter().map(|m| m.inner));
        Self { inner }
    }

    /// Return its json representation of molecule object.
    fn to_json(&self) -> PyResult<String> {
        let json = gchemol::io::to_json(&self.inner)?;
        Ok(json)
    }

    /// Clean up molecule geometry using stress majorization algorithm.
    fn clean(&mut self) -> PyResult<()> {
        self.inner.clean()?;
        Ok(())
    }

    /// Renumber atoms consecutively from 1.
    fn renumber(&mut self) -> PyResult<()> {
        self.inner.renumber();
        Ok(())
    }

    /// Write molecule to file with `path`. The molecule format will
    /// be determined based on file name extension.
    fn to_file(&self, path: String) -> PyResult<()> {
        self.inner.to_file(&path)?;
        Ok(())
    }

    /// Render molecule with template file from `path`. On success,
    /// return the formatted string.
    fn render_with(&self, path: String) -> PyResult<String> {
        let s = self.inner.render_with(path.as_ref())?;
        Ok(s)
    }

    /// Get the number of atoms.
    fn natoms(&self) -> PyResult<usize> {
        let n = self.inner.natoms();
        Ok(n)
    }

    /// Return chemical formula.
    fn formula(&self) -> PyResult<String> {
        Ok(self.inner.formula())
    }

    /// Get the number of bonds.
    fn nbonds(&self) -> PyResult<usize> {
        let n = self.inner.nbonds();
        Ok(n)
    }

    /// Return the name of the molecule, which is typpically modified
    /// for safely storing in various chemical file formats.
    fn title(&self) -> PyResult<String> {
        let t = self.inner.title();
        Ok(t)
    }

    /// Return the shortest distance counted in number of chemical
    /// bonds between two atoms. Return None if they are not
    /// connected.
    fn nbonds_between(&self, i: usize, j: usize) -> PyResult<Option<usize>> {
        let n = self.inner.nbonds_between(i, j);
        Ok(n)
    }

    /// Unbuild current crystal structure leaving a nonperiodic structure
    fn unbuild_crystal(&mut self) -> PyResult<()> {
        self.inner.unbuild_crystal();
        Ok(())
    }

    /// Create a Lattice from the minimal bounding box of the Molecule
    /// extended by a positive value of padding.
    ///
    /// NOTE: padding has to be large enough (> 0.5) to avoid self
    /// interaction with its periodic mirror.
    fn set_lattice_from_bounding_box(&mut self, padding: f64) -> PyResult<()> {
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
    fn supercell(&mut self, sa: usize, sb: usize, sc: usize) -> PyResult<Self> {
        if let Some(mol) = self.inner.supercell(sa, sb, sc) {
            let m = Self { inner: mol };
            Ok(m)
        } else {
            Err(pyo3::exceptions::PyException::new_err("not allowed for periodic structure"))
        }
    }

    /// Recalculates all bonds in molecule based on interatomic
    /// distances and covalent radii. For periodic system, the bonds
    /// are determined by applying miniumu image convention.
    fn rebond(&mut self) -> PyResult<()> {
        self.inner.rebond();
        Ok(())
    }

    /// Removes all existing bonds between atoms.
    fn unbound(&mut self) -> PyResult<()> {
        self.inner.unbound();
        Ok(())
    }

    /// Removes all bonds between two selections to respect pymol's
    /// unbond command.
    ///
    /// # Parameters
    ///
    /// * atom_indices1: the first collection of atoms
    /// * atom_indices2: the other collection of atoms
    ///
    /// # Reference
    ///
    /// * <https://pymolwiki.org/index.php/Unbond>
    ///
    fn unbond(&mut self, atom_indices1: Vec<usize>, atom_indices2: Vec<usize>) -> PyResult<()> {
        self.inner.unbond(&atom_indices1, &atom_indices2);
        Ok(())
    }

    /// Return the distance between atom i and atom j. For periodic
    /// structure, this method will return the distance under the
    /// minimum image convention.
    fn distance(&self, i: usize, j: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_distance(i, j)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}"))?;
        Ok(d)
    }

    /// Return the angle between atoms `i`, `j`, `k` in degrees,
    /// irrespective periodic images.
    fn angle(&self, i: usize, j: usize, k: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_angle(i, j, k)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}, {k}"))?
            .to_degrees();
        Ok(d)
    }

    /// Return the torsion angle between atoms `i`, `j`, `k`, `l` in
    /// degrees, irrespective periodic images.
    fn torsion(&self, i: usize, j: usize, k: usize, l: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_torsion(i, j, k, l)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}, {k}, {l}"))?
            .to_degrees();
        Ok(d)
    }

    /// Return molecule’s inertia matrix (3x3) in reference to molecule’s center of mass
    fn inertia_matrix(&self) -> PyResult<[[f64; 3]; 3]> {
        let im = self.inner.inertia_matrix();
        Ok(im)
    }

    /// Return the center of mass of molecule (COM).
    fn center_of_mass(&self) -> PyResult<[f64; 3]> {
        let com = self.inner.center_of_mass();
        Ok(com)
    }

    /// Return the center of geometry of molecule (COG).
    fn center_of_geometry(&self) -> PyResult<[f64; 3]> {
        let cog = self.inner.center_of_geometry();
        Ok(cog)
    }

    /// Set freezing flag to `freezed` for atom `sn` when in
    /// optimization or dynamic simulation.
    fn freeze_atom(&mut self, sn: usize, freezed: bool) -> PyResult<()> {
        if let Some(a) = self.inner.get_atom_mut(sn) {
            a.set_freezing([freezed; 3]);
            Ok(())
        } else {
            Err(pyo3::exceptions::PyException::new_err("atom {sn} not found"))
        }
    }

    /// Get ONIOM layer of atom `sn`. If no layer information was set,
    /// will return `H`.
    fn get_oniom_layer(&self, sn: usize) -> PyResult<String> {
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
    fn set_oniom_layer(&mut self, sn: usize, layer: &str) -> PyResult<()> {
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
    fn find_rings(&mut self, nmax: usize) -> PyResult<Vec<std::collections::HashSet<usize>>> {
        let rings = self.inner.find_rings(nmax);
        Ok(rings)
    }

    /// Return atom serial numbers.
    fn numbers(&self) -> Vec<usize> {
        self.inner.numbers().collect()
    }

    /// Return an iterator over a tuple of atom serial number `n` and
    /// its associated `Atom` (n, Atom)
    fn atoms(&self) -> PyAtomsIter {
        let atoms: Vec<_> = self
            .inner
            .atoms()
            .map(|(i, atom)| (i, _Atom { inner: atom.to_owned() }))
            .collect();
        PyAtomsIter {
            inner: atoms.into_iter(),
        }
    }

    /// Return a sub molecule induced by `atoms` in parent
    /// molecule. Return None if atom serial numbers are
    /// invalid. Return an empty Molecule if `atoms` empty.
    fn get_sub_molecule(&self, atoms: Vec<usize>) -> Option<Self> {
        let inner = self.inner.get_sub_molecule(&atoms)?;
        Self { inner }.into()
    }

    /// Set periodic lattice.
    fn set_lattice(&mut self, lat: _Lattice) {
        self.inner.set_lattice(lat.inner);
    }

    /// Get periodic lattice.
    fn get_lattice(&self) -> Option<_Lattice> {
        let lat = self.inner.get_lattice()?;
        _Lattice { inner: *lat }.into()
    }

    /// Return true if Molecule is a periodic structure.
    fn is_periodic(&self) -> bool {
        self.inner.is_periodic()
    }

    /// Return fractional coordinates relative to unit cell. Return
    /// None if not a periodic structure.
    fn get_scaled_positions(&self) -> Option<Vec<[f64; 3]>> {
        let scaled = self.inner.get_scaled_positions()?.collect_vec();
        scaled.into()
    }

    // fn educated_rebond(&mut self) {
    //     use educate::prelude::*;
    //     self.inner.educated_rebond();
    // }
    // 
    // fn educated_clean(&mut self) {
    //     use educate::prelude::*;
    //     self.inner.educated_clean();
    // }
    // 
    // fn educated_clean_selection(&mut self, selection: Vec<usize>) {
    //     use educate::prelude::*;
    //     self.inner.educated_clean_selection(&selection);
    // }
    // 
    // /// Return unique fingerprint of current molecule
    // fn fingerprint(&self) -> String {
    //     use spdkit::prelude::FingerPrintExt;
    //     self.inner.fingerprint()
    // }
    // 
    // /// This is an operation of reordering the atoms in a way that does not depend
    // /// on where they were before. The bonding graph is important for this
    // /// operation.
    // fn reorder_cannonically(&mut self) -> Vec<usize> {
    //     use spdkit::prelude::FingerPrintExt;
    //     self.inner.reorder_cannonically()
    // }
    // 
    // /// Convert `Molecule` to a graph object for distance geometry
    // /// refinement.
    // fn to_distance_geometry_graph(&self) -> DgGraph {
    //     let dg = self.inner.distance_geometry_graph();
    //     DgGraph { inner: dg }
    // }
}
// 969a9313 ends here

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
    fn set_distance_weight(&mut self, i: usize, j: usize, w: f64) {
        self.inner.set_distance_weight(i, j, w);
    }

    /// Set distance bound (lower/upper) for atom `i` and `j` for
    /// refinement of molecule structure.
    fn set_distance_bound(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        self.inner.set_distance_bound(i, j, lb, ub);
    }

    /// Returns distance bound (lower/upper) for atom `i` and `j`.
    fn get_distance_bound(&self, i: usize, j: usize) -> (f64, f64) {
        self.inner.distance_bound(i, j)
    }

    /// Returns distance weight for atom `i` and `j`.
    fn get_distance_weight(&self, i: usize, j: usize) -> f64 {
        self.inner.distance_weight(i, j)
    }

    /// Refine molecule structure `mol` using distance geometry.
    fn refine_molecule(&mut self, mol: &mut _Molecule) {
        self.inner.refine_molecule(&mut mol.inner);
    }
}
// df84f7ba ends here

// [[file:../spdkit-python.note::c400da41][c400da41]]
#[pyfunction]
fn set_verbosity(level: u8) {
    let mut log = gut::cli::Verbosity::default();
    log.set_verbosity(level);
    log.setup_logger();
}
// c400da41 ends here

// [[file:../spdkit-python.note::a31e85a4][a31e85a4]]
#[pyclass]
struct PyAtomsIter {
    inner: std::vec::IntoIter<(usize, _Atom)>,
}

#[pymethods]
impl PyAtomsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(usize, _Atom)> {
        slf.inner.next()
    }
}

#[pyclass]
struct PyMoleculeIter {
    iter: Box<dyn Iterator<Item = _Molecule> + Send>,
}

#[pymethods]
impl PyMoleculeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<_Molecule> {
        slf.iter.next()
    }
}

#[pyfunction]
/// Read a list of `Molecule` from `path`. Returns an iterator over
/// `Molecule`, which allows reading a large file out of memory.
fn read(path: String) -> PyResult<PyMoleculeIter> {
    let mols = gchemol::io::read(path)?.map(|inner| _Molecule { inner });
    let mols = PyMoleculeIter { iter: Box::new(mols) };
    Ok(mols)
}
// a31e85a4 ends here

// [[file:../spdkit-python.note::fbe87af8][fbe87af8]]
#[pymodule]
#[pyo3(name = "spdkit")]
fn pyspdkit(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<_Molecule>()?;
    m.add_class::<_Atom>()?;
    m.add_class::<_Lattice>()?;
    m.add_function(wrap_pyfunction!(set_verbosity, m)?)?;

    let io = PyModule::new(py, "io")?;
    io.add_function(wrap_pyfunction!(read, io)?)?;
    m.add_submodule(io)?;

    Ok(())
}
// fbe87af8 ends here
