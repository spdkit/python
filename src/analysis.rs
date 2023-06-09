// [[file:../spdkit-python.note::106f6156][106f6156]]
use crate::common::*;
use pyo3::prelude::*;
// 106f6156 ends here

// [[file:../spdkit-python.note::cb1870b5][cb1870b5]]
use crate::PyMolecule;
use bond_valence::*;
use std::collections::HashMap;

/// Represent the evaluated atom properties such as atom valence,
/// general atom valence, coordination number etc.
#[pyclass(name = "BondValenceModel")]
pub struct PyBondValenceModel {
    inner: BondValenceModel,
}

#[pymethods]
impl PyBondValenceModel {
    /// Calculate atom valence, coordination number, and general atom
    /// valence of each atom in `mol`.
    ///
    /// # Parameters
    ///   * mol: The molecule to evaluate.
    ///   * opts: optional bond valence parameters
    #[staticmethod]
    fn evaluate(mol: PyMolecule, opts: Option<Options>) -> HashMap<usize, Evaluated> {
        let opts = opts.unwrap_or_else(|| Options::default_from(&mol.inner));
        BondValenceModel::evaluate(&mol.inner, &opts)
    }

    /// Create chemical bonds from nearest neighbors based on bond
    /// valence values. The atom pair with larger bond valence will be
    /// set as chemical bond automatically determined by kmean
    /// clustering algorithm. Return the number of chemical bonds
    /// created.
    #[staticmethod]
    #[pyo3(signature = (mol, /, opts=None, ignore_pbc=false))]
    fn rebond(mol: &mut PyMolecule, opts: Option<Options>, ignore_pbc: bool) -> usize {
        let opts = opts.unwrap_or_else(|| Options::default_from(&mol.inner));
        let lat = if ignore_pbc { mol.unbuild_crystal() } else { None };
        let nbonds = BondValenceModel::rebond(&mut mol.inner, &opts);
        if ignore_pbc {
            if let Some(lat) = lat {
                mol.set_lattice(lat);
            }
        }
        nbonds
    }

    #[staticmethod]
    /// Return default bond valence options
    fn default_options(mol: PyMolecule) -> Options {
        Options::default_from(&mol.inner)
    }
}
// cb1870b5 ends here

// [[file:../spdkit-python.note::c84c0fc1][c84c0fc1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyBondValenceModel>()?;
    m.add_class::<Options>()?;
    // m.add_class::<mhm::PyMinimaHopping>()?;
    // m.add_function(wrap_pyfunction!(read, m)?)?;

    Ok(m)
}
// c84c0fc1 ends here
