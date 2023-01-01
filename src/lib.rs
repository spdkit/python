// [[file:../spdkit-python.note::0cbf1e93][0cbf1e93]]
use pyo3::prelude::*;
use pyo3::types::PyType;
// 0cbf1e93 ends here

// [[file:../spdkit-python.note::95be1618][95be1618]]
use gchemol::Atom;

#[pyclass(name = "Atom")]
#[derive(Clone)]
/// Atom is the smallest particle still characterizing a chemical element.
#[pyo3(text_signature = "(symbol, position)")]
pub struct PyAtom {
    inner: Atom,
}

#[pymethods]
impl PyAtom {
    #[new]
    /// Construct `Atom` object from `symbol` and `position`.
    pub fn new(symbol: String, position: [f64; 3]) -> Self {
        Self {
            inner: Atom::new(symbol, position),
        }
    }

    /// Return element symbol
    pub fn symbol(&self) -> String {
        self.inner.symbol().to_string()
    }

    /// Return atomic number
    pub fn number(&self) -> usize {
        self.inner.number()
    }

    /// Get mass in atomic mass unit. Return None if atom is dummy.
    pub fn get_mass(&self) -> Option<f64> {
        self.inner.get_mass()
    }

    /// Return atom position in 3D Cartesian coordinates
    pub fn position(&self) -> [f64; 3] {
        self.inner.position()
    }

    /// Change atom Cartesian position to `p` ([x, y, z]).
    #[pyo3(text_signature = "($self, p, /)")]
    pub fn set_position(&mut self, p: [f64; 3]) {
        self.inner.set_position(p);
    }

    /// Set atom symbol.
    #[pyo3(text_signature = "($self, symbol, /)")]
    pub fn set_symbol(&mut self, symbol: String) {
        self.inner.set_symbol(symbol);
    }

    /// Return freezing mask array for Cartesian coordinates
    pub fn freezing(&self) -> [bool; 3] {
        self.inner.freezing()
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
    #[pyo3(text_signature = "($self, a, b, c, alpha, beta, gamma, /)")]
    pub fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let inner = Lattice::from_params(a, b, c, alpha, beta, gamma);
        Self { inner }
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    #[pyo3(text_signature = "($self, p, /)")]
    pub fn to_frac(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_frac(p).into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    #[pyo3(text_signature = "($self, p, /)")]
    pub fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.to_cart(p).into()
    }

    /// Returns Lattice vector a.
    pub fn vector_a(&self) -> [f64; 3] {
        self.inner.vector_a().into()
    }

    /// Returns Lattice vector b.
    pub fn vector_b(&self) -> [f64; 3] {
        self.inner.vector_b().into()
    }

    /// Returns Lattice vector c.
    pub fn vector_c(&self) -> [f64; 3] {
        self.inner.vector_c().into()
    }

    /// Return the shortest vector obeying the minimum image convention.
    pub fn apply_mic(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.apply_mic(p).into()
    }

    /// Return the shortest distance between pi (point i) and the
    /// periodic images of pj (point j) under the minimum image
    /// convention
    #[pyo3(text_signature = "($self, pi, pj, /)")]
    pub fn distance(&self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        self.inner.distance(pi, pj).into()
    }

    /// Lattice length parameters: a, b, c.
    pub fn lengths(&self) -> [f64; 3] {
        self.inner.lengths().into()
    }

    /// Return the volume of the unit cell the cache will be updated
    /// if necessary
    pub fn volume(&self) -> f64 {
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

    /// Return the name of the molecule, which is typpically modified
    /// for safely storing in various chemical file formats.
    pub fn title(&self) -> String {
        self.inner.title()
    }

    #[staticmethod]
    /// Build a molecule from atoms associated with serial numbers from 1.
    #[pyo3(text_signature = "($self, atoms, /)")]
    pub fn from_atoms(atoms: Vec<PyAtom>) -> Self {
        let inner = Molecule::from_atoms(atoms.into_iter().map(|m| m.inner));
        Self { inner }
    }

    /// Return its json representation of molecule object.
    pub fn to_json(&self) -> PyResult<String> {
        let json = gchemol::io::to_json(&self.inner)?;
        Ok(json)
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

    /// Write molecule to file with `path`. The molecule format will
    /// be determined based on file name extension.
    #[pyo3(text_signature = "($self, path, /)")]
    pub fn to_file(&self, path: String) -> PyResult<()> {
        self.inner.to_file(&path)?;
        Ok(())
    }

    /// Render molecule with template file from `path`. On success,
    /// return the formatted string.
    #[pyo3(text_signature = "($self, path, /)")]
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

    /// Unbuild current crystal structure leaving a nonperiodic structure
    pub fn unbuild_crystal(&mut self) -> PyResult<()> {
        self.inner.unbuild_crystal();
        Ok(())
    }

    /// Create a Lattice from the minimal bounding box of the Molecule
    /// extended by a positive value of padding.
    ///
    /// NOTE: padding has to be large enough (> 0.5) to avoid self
    /// interaction with its periodic mirror.
    #[pyo3(text_signature = "($self, padding, /)")]
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
    #[pyo3(text_signature = "($self, sa, sb, sc, /)")]
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
    #[pyo3(text_signature = "($self, i, j, /)")]
    pub fn distance(&self, i: usize, j: usize) -> PyResult<f64> {
        let d = self
            .inner
            .get_distance(i, j)
            .ok_or(format_err!("invalid serial numbers: {i}, {j}"))?;
        Ok(d)
    }

    /// Return the angle between atoms `i`, `j`, `k` in degrees,
    /// irrespective periodic images.
    #[pyo3(text_signature = "($self, i, j, k, /)")]
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
    #[pyo3(text_signature = "($self, i, j, k, l, /)")]
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

    /// Set freezing flag to `freezed` for atom `sn` when in
    /// optimization or dynamic simulation.
    #[pyo3(text_signature = "($self, sn, freezed, /)")]
    pub fn freeze_atom(&mut self, sn: usize, freezed: bool) -> PyResult<()> {
        if let Some(a) = self.inner.get_atom_mut(sn) {
            a.set_freezing([freezed; 3]);
            Ok(())
        } else {
            Err(pyo3::exceptions::PyException::new_err("atom {sn} not found"))
        }
    }

    /// Get ONIOM layer of atom `sn`. If no layer information was set,
    /// will return `H`.
    #[pyo3(text_signature = "($self, sn, /)")]
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
    #[pyo3(text_signature = "($self, sn, layer, /)")]
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
    #[pyo3(text_signature = "($self, namx, /)")]
    pub fn find_rings(&mut self, nmax: usize) -> PyResult<Vec<std::collections::HashSet<usize>>> {
        let rings = self.inner.find_rings(nmax);
        Ok(rings)
    }

    /// Return atom serial numbers.
    pub fn numbers(&self) -> Vec<usize> {
        self.inner.numbers().collect()
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
    #[pyo3(text_signature = "($self, n, /)")]
    pub fn get_atom(&self, n: usize) -> Option<PyAtom> {
        let inner = self.inner.get_atom(n)?.clone();
        PyAtom { inner }.into()
    }
    
    /// Add atom a into molecule. If Atom numbered as a already exists in
    /// molecule, then the associated Atom will be updated with atom.
    #[pyo3(text_signature = "($self, n, atom, /)")]
    pub fn add_atom(&mut self, n: usize, atom: PyAtom) {
        self.inner.add_atom(n, atom.inner)
    }
    
    /// Remove Atom a from Molecule. Return the removed Atom on success,
    /// and return None if Atom a does not exist.
    #[pyo3(text_signature = "($self, n, /)")]
    pub fn remove_atom(&mut self, n: usize) -> Option<PyAtom> {
        let inner = self.inner.remove_atom(n)?;
        PyAtom { inner }.into()
    }

    /// Return a sub molecule induced by `atoms` in parent
    /// molecule. Return None if atom serial numbers are
    /// invalid. Return an empty Molecule if `atoms` empty.
    #[pyo3(text_signature = "($self, atoms)")]
    pub fn get_sub_molecule(&self, atoms: Vec<usize>) -> Option<Self> {
        let inner = self.inner.get_sub_molecule(&atoms)?;
        Self { inner }.into()
    }

    /// Set periodic lattice.
    #[pyo3(text_signature = "($self, lat)")]
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
    
    /// Recalculates all bonds in molecule based on interatomic
    /// distances and covalent radii. For periodic system, the bonds
    /// are determined by applying miniumu image convention.
    pub fn rebond(&mut self) {
        self.inner.rebond();
    }
    
    /// Removes all existing bonds between atoms.
    pub fn unbound(&mut self) -> PyResult<()> {
        self.inner.unbound();
        Ok(())
    }
    
    /// Removes all bonds between `atom_indices1` and `atom_indices2`
    /// in respect of pymol's `unbond` command.
    ///
    /// # Parameters
    /// * atom_indices1: the first collection of atoms
    /// * atom_indices2: the other collection of atoms
    ///
    /// # Reference
    /// * <https://pymolwiki.org/index.php/Unbond>
    #[pyo3(text_signature = "($self, atom_indices1, atom_indices2, /)")]
    pub fn unbond(&mut self, atom_indices1: Vec<usize>, atom_indices2: Vec<usize>) -> PyResult<()> {
        self.inner.unbond(&atom_indices1, &atom_indices2);
        Ok(())
    }
    
    /// Add a single bond between Atom `i` and Atom `j` into molecule.
    /// Panic if the specified atom a or b does not exist
    #[pyo3(text_signature = "($self, i, j, /)")]
    pub fn add_bond(&mut self, i: usize, j: usize) {
        use gchemol::Bond;
    
        self.inner.add_bond(i, j, Bond::single())
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
    
    /// Reorder the atoms according to the ordering of keys. Keys define
    /// 1-to-1 mapping of atoms.
    ///
    /// # Parameters
    /// * keys: a list of numbers for sorting
    ///
    /// # NOTE
    /// * This method will cause serial numbers renumbered from 1.
    #[pyo3(text_signature = "($self, keys, /)")]
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
    #[pyo3(text_signature = "($self, mol_ref, /, selection = None)")]
    pub fn superimpose_onto(&mut self, mol_ref: Self, selection: Option<Vec<usize>>) -> f64 {
        use gchemol::geom::prelude::*;
        use gchemol::geom::Superimpose;
    
        let (positions_this, positions_prev, weights) = if let Some(selected) = selection {
            let this = selected.iter().map(|&i| self.get_atom(i).unwrap().position()).collect_vec();
            let prev = selected
                .iter()
                .map(|&i| mol_ref.get_atom(i).unwrap().position())
                .collect_vec();
            let weights = selected
                .iter()
                .map(|&i| self.get_atom(i).unwrap().get_mass().unwrap())
                .collect_vec();
            (this, prev, weights)
        } else {
            (
                self.inner.positions().collect_vec(),
                mol_ref.inner.positions().collect_vec(),
                self.inner.masses().collect_vec(),
            )
        };
        assert_eq!(positions_this.len(), positions_prev.len());
        assert_eq!(positions_this.len(), weights.len());
        let sp1 = Superimpose::new(&positions_this).onto(&positions_prev, weights.as_slice().into());
        let mut positions_this_mirrored = positions_this.clone();
        positions_this_mirrored.mirror_invert();
        let sp2 = Superimpose::new(&positions_this_mirrored).onto(&positions_prev, weights.as_slice().into());
        let positions_this_all = self.inner.positions().collect_vec();
        let (positions_new, rmsd) = if sp1.rmsd < sp2.rmsd {
            (sp1.apply(&positions_this_all), sp1.rmsd)
        } else {
            let mut positions_this_all_mirrored = positions_this_all.clone();
            positions_this_all_mirrored.mirror_invert();
            (sp2.apply(&positions_this_all_mirrored), sp2.rmsd)
        };
    
        self.inner.set_positions(positions_new);
        rmsd
    }

    pub fn educated_rebond(&mut self) {
        use educate::prelude::*;
        self.inner.educated_rebond();
    }
    
    pub fn educated_clean(&mut self) {
        use educate::prelude::*;
        self.inner.educated_clean();
    }
    
    pub fn educated_clean_selection(&mut self, selection: Vec<usize>) {
        use educate::prelude::*;
        self.inner.educated_clean_selection(&selection);
    }
    
    /// Return unique fingerprint of current molecule
    pub fn fingerprint(&self) -> String {
        use spdkit::prelude::FingerPrintExt;
        self.inner.fingerprint()
    }
    
    /// This is an operation of reordering the atoms in a way that does not depend
    /// on where they were before. The bonding graph is important for this
    /// operation.
    pub fn reorder_cannonically(&mut self) -> Vec<usize> {
        use spdkit::prelude::FingerPrintExt;
        self.inner.reorder_cannonically()
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
    #[pyo3(text_signature = "($self, path, /)")]
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
    #[pyo3(text_signature = "($self, src, fmt)")]
    pub fn from_string(src: String, fmt: String) -> PyResult<Self> {
        let inner = Molecule::from_str(&src, &fmt)?;
        Ok(Self { inner })
    }
    
    #[staticmethod]
    /// Describe available backends for reading/writing molecule
    pub fn describe_available_backends() {
        gchemol::io::describe_backends();
    }
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

    /// Refine molecule structure `mol` using distance geometry.
    pub fn refine_molecule(&mut self, mol: &mut PyMolecule) {
        self.inner.refine_molecule(&mut mol.inner);
    }
}
// df84f7ba ends here

// [[file:../spdkit-python.note::6144a4ef][6144a4ef]]
use gchemol_parser::TextViewer;

/// A simple line-based text viewer for quick peeking part of text
#[pyclass(name = "TextViwer", subclass)]
#[derive(Clone)]
pub struct PyTextViewer {
    inner: TextViewer,
}

#[pymethods]
impl PyTextViewer {
    /// Create a `TextViewer` from text string in `s`.
    #[staticmethod]
    #[pyo3(text_signature = "($self, s)")]
    pub fn from_string(s: String) -> Self {
        let inner = TextViewer::from_str(&s);
        Self { inner }
    }

    /// Return total number of lines.
    pub fn num_lines(&self) -> usize {
        self.inner.num_lines()
    }

    /// Get line number at cursor.
    pub fn current_line_num(&self) -> usize {
        self.inner.current_line_num()
    }

    /// Return all the text.
    pub fn text(&self) -> String {
        self.inner.text().to_owned()
    }

    /// Return the line at cursor.
    pub fn current_line(&self) -> String {
        self.inner.current_line().to_owned()
    }

    /// Move the cursor to line `n`, counting from line 1 at beginning
    /// of the text.
    #[pyo3(text_signature = "($self, n)")]
    pub fn goto_line(&mut self, n: usize) {
        self.inner.goto_line(n);
    }

    /// Move the cursor to the beginning of the previous line.
    pub fn goto_previous_line(&mut self) {
        self.inner.goto_previous_line();
    }

    /// Move the cursor to the beginning of the next line.
    pub fn goto_next_line(&mut self) {
        self.inner.goto_next_line();
    }

    /// Move the cursor to the beginning of the last line.
    pub fn goto_last_line(&mut self) {
        self.inner.goto_last_line();
    }

    /// Move the cursor to the beginning of the first line.
    pub fn goto_first_line(&mut self) {
        self.inner.goto_first_line();
    }

    /// Move the cursor to the line matching `pattern`. Regex pattern
    /// is allowed.
    #[pyo3(text_signature = "($self, pattern)")]
    pub fn search_forward(&mut self, pattern: String) -> PyResult<usize> {
        let n = self.inner.search_forward(&pattern)?;
        Ok(n)
    }

    /// Search backward from current point for `pattern`. Return
    /// current point after search.
    #[pyo3(text_signature = "($self, pattern)")]
    pub fn search_backward(&mut self, pattern: String) -> PyResult<usize> {
        let n = self.inner.search_backward(&pattern)?;
        Ok(n)
    }

    /// Peek line `n`.
    #[pyo3(text_signature = "($self, n)")]
    pub fn peek_line(&self, n: usize) -> String {
        self.inner.peek_line(n).to_owned()
    }

    /// Peek the text between line `n` and `m` (including line `m`)
    #[pyo3(text_signature = "($self, n, m)")]
    pub fn peek_lines(&self, n: usize, m: usize) -> String {
        self.inner.peek_lines(n, m).to_owned()
    }

    /// Select the next `n` lines from current point, including current line.
    #[pyo3(text_signature = "($self, n)")]
    pub fn selection(&self, n: usize) -> String {
        self.inner.selection(n).to_owned()
    }

    /// Select part of the string in next `n` lines (including
    /// currrent line), in a rectangular area surrounded by columns in
    /// `col_beg`--`col_end`.
    #[pyo3(text_signature = "($self, n, col_beg, col_end)")]
    pub fn column_selection(&self, n: usize, col_beg: usize, col_end: usize) -> String {
        self.inner.column_selection(n, col_beg, col_end)
    }
}
// 6144a4ef ends here

// [[file:../spdkit-python.note::7ff1511e][7ff1511e]]
use gchemol_parser::GrepReader;

/// Quick grep text by marking the line that matching a pattern,
/// suitable for very large text file.
#[pyclass(name = "GrepReader", subclass)]
pub struct PyGrepReader {
    inner: GrepReader,
}

#[pymethods]
impl PyGrepReader {
    /// Create a `GrepReader` from file in `path`.
    #[staticmethod]
    #[pyo3(text_signature = "($self, path)")]
    pub fn from_file(path: String) -> PyResult<Self> {
        let inner = GrepReader::try_from_path(path.as_ref())?;
        Ok(Self { inner })
    }

    /// Mark positions that matching pattern, so that we can seek these
    /// positions later. Return the number of marked positions.
    #[pyo3(text_signature = "($self, pattern, max_count = None)")]
    pub fn mark(&mut self, pattern: String, max_count: Option<usize>) -> PyResult<usize> {
        let n = self.inner.mark(&pattern, max_count)?;
        Ok(n)
    }

    /// Goto the start of the inner file.
    pub fn goto_start(&mut self) {
        self.inner.goto_start();
    }
    /// Goto the end of the inner file.
    pub fn goto_end(&mut self) {
        self.inner.goto_end();
    }

    /// Return the number of marked positions
    pub fn num_markers(&self) -> usize {
        self.inner.num_markers()
    }

    /// Goto the next position that marked. Return marker position on success.
    /// Return Err if already reached the last marker or other errors.
    pub fn goto_next_marker(&mut self) -> PyResult<u64> {
        let n = self.inner.goto_next_marker()?;
        Ok(n)
    }

    /// Goto the marked position in `marker_index`. Will panic if marker_index
    /// out of range.
    #[pyo3(text_signature = "($self, marker_index)")]
    pub fn goto_marker(&mut self, marker_index: isize) -> PyResult<u64> {
        let i = if marker_index < 0 {
            self.inner.num_markers() as isize + marker_index
        } else {
            marker_index
        };
        let n = self.inner.goto_marker(i as usize)?;
        Ok(n)
    }

    /// Return `n` lines in string on success from current
    /// position. Return error if reached EOF early.
    #[pyo3(text_signature = "($self, n)")]
    pub fn read_lines(&mut self, n: usize) -> PyResult<String> {
        let mut s = String::new();
        self.inner.read_lines(n, &mut s)?;
        Ok(s)
    }

    /// View next `n` lines like in a normal text viewer. This method
    /// will forward the cursor by `n` lines.
    #[pyo3(text_signature = "($self, n)")]
    pub fn view_lines(&mut self, n: usize) -> PyResult<PyTextViewer> {
        let inner = self.inner.view_lines(n)?;
        Ok(PyTextViewer { inner })
    }
}
// 7ff1511e ends here

// [[file:../spdkit-python.note::c400da41][c400da41]]
#[pyfunction]
#[pyo3(text_signature = "(level)")]
/// Enable logging by setting verbosity level to `level`.
pub fn set_verbosity(level: u8) {
    let mut log = gut::cli::Verbosity::default();
    log.set_verbosity(level);
    log.setup_logger();
}
// c400da41 ends here

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

// [[file:../spdkit-python.note::b8744829][b8744829]]
#[pyfunction]
/// Read a list of `Molecule` from `path`. Returns an iterator over
/// `Molecule`, which allows reading a large file out of memory.
#[pyo3(text_signature = "(path: String)")]
pub fn read(path: String) -> PyResult<PyMoleculeIter> {
    let mols = gchemol::io::read(path)?.map(|inner| PyMolecule { inner });
    let mols = PyMoleculeIter { iter: Box::new(mols) };
    Ok(mols)
}

#[pyfunction]
/// Write molecules into path. File format will be determined according to the path.
#[pyo3(text_signature = "(path: String, mols: List[Molecule])")]
pub fn write(path: String, mols: Vec<PyMolecule>) -> PyResult<()> {
    gchemol::io::write(path, mols.iter().map(|mol| &mol.inner))?;
    Ok(())
}
// b8744829 ends here

// [[file:../spdkit-python.note::fbe87af8][fbe87af8]]
#[pymodule]
#[pyo3(name = "spdkit")]
fn pyspdkit(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyMolecule>()?;
    m.add_class::<PyAtom>()?;
    m.add_class::<PyLattice>()?;
    m.add_function(wrap_pyfunction!(set_verbosity, m)?)?;

    let io = PyModule::new(py, "io")?;
    io.add_function(wrap_pyfunction!(read, io)?)?;
    io.add_function(wrap_pyfunction!(write, io)?)?;
    io.add_class::<PyTextViewer>()?;
    io.add_class::<PyGrepReader>()?;
    m.add_submodule(io)?;

    Ok(())
}
// fbe87af8 ends here
