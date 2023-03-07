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

// [[file:../spdkit-python.note::3b89501f][3b89501f]]
#[pyfunction]
#[pyo3(signature = (mol, /, probe_element="C", maxcycle=1000))]
#[pyo3(text_signature = "(mol, /, probe_element='C', n_probes=1000)")]
pub fn probe_adsorption_sites(mol: &PyMolecule, probe_element: &str, maxcycle: usize) -> Result<Vec<PyMolecule>> {
    use ::surface::docs::surface::probe_adsorption_sites;
    let probed = probe_adsorption_sites(&mol.inner, probe_element, maxcycle)?;
    let probed = probed.into_iter().map(|inner| PyMolecule { inner }).collect();
    Ok(probed)
}
// 3b89501f ends here

// [[file:../spdkit-python.note::09bb2ce3][09bb2ce3]]
#[pyfunction]
/// Fragment molecule into connected parts by layer for periodic slab model.
pub fn fragment_atoms_by_layer(mol: &PyMolecule) -> Result<Vec<Vec<usize>>> {
    use spdkit_surface::fragment_atoms_by_layer;

    let probed = fragment_atoms_by_layer(&mol.inner)?.collect();
    Ok(probed)
}
// 09bb2ce3 ends here

// [[file:../spdkit-python.note::4b881b59][4b881b59]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyTemplate>()?;
    m.add_function(wrap_pyfunction!(probe_surface_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(probe_adsorption_sites, m)?)?;
    m.add_function(wrap_pyfunction!(fragment_atoms_by_layer, m)?)?;

    Ok(m)
}
// 4b881b59 ends here
