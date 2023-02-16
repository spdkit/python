// [[file:../spdkit-python.note::*imports][imports:1]]
use crate::common::*;
use pyo3::prelude::*;
// imports:1 ends here

// [[file:../spdkit-python.note::45f23185][45f23185]]
use super::PyMolecule;

#[pyfunction]
#[pyo3(signature = (mol, /, r_cutoff=2.0, n_probes=500))]
#[pyo3(text_signature = "(mol, /, r_cutoff=2.0, n_probes=500)")]
/// Probe surface atoms approaching slab surface from top in z-axis
/// gradually
///
/// # Parameters
/// * mol: the molecule to be probed (slab model)
/// * r_cutoff: the cutoff radius for probing surface atoms
/// * n_probes: the number of random probe atoms
pub fn probe_surface_atoms(mol: &PyMolecule, r_cutoff: f64, n_probes: usize) -> Result<Vec<usize>> {
    let probed = spdkit_surface::probe::probe_surface_atoms(&mol.inner, r_cutoff, n_probes)?;
    Ok(probed)
}
// 45f23185 ends here

// [[file:../spdkit-python.note::4b881b59][4b881b59]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyTemplate>()?;
    m.add_function(wrap_pyfunction!(probe_surface_atoms, m)?)?;

    Ok(m)
}
// 4b881b59 ends here
