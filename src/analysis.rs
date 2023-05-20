// [[file:../spdkit-python.note::106f6156][106f6156]]
use crate::common::*;
use pyo3::prelude::*;
// 106f6156 ends here

// [[file:../spdkit-python.note::cb1870b5][cb1870b5]]
use crate::PyMolecule;
use bvcalc::*;
use std::collections::HashMap;

/// Represent the evaluated atom properties such as atom valence,
/// general atom valence, coordination number etc.
#[pyclass(name = "BondValence")]
pub struct PyBondValence;

#[pymethods]
impl PyBondValence {
    /// Calculate atom valence, coordination number, and general atom
    /// valence of each atom in `mol`.
    #[staticmethod]
    #[pyo3(signature = (mol, /, distance_cutoff=4.0, bond_valence_cutoff=0.3))]
    #[pyo3(text_signature = "($self, mol, /, distance_cutoff=5.0, bond_valence_cutoff=0.3)")]
    fn evaluate(mol: PyMolecule, distance_cutoff: f64, bond_valence_cutoff: f64) -> HashMap<usize, Evaluated> {
        let opts = Options {
            distance_cutoff,
            bond_valence_cutoff,
            ..Default::default()
        };

        BondValence::evaluate(&mol.inner, &opts)
    }
}
// cb1870b5 ends here

// [[file:../spdkit-python.note::c84c0fc1][c84c0fc1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyBondValence>()?;
    // m.add_class::<mhm::PyMinimaHopping>()?;
    // m.add_function(wrap_pyfunction!(read, m)?)?;

    Ok(m)
}
// c84c0fc1 ends here
