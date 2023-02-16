// [[file:../spdkit-python.note::1f520294][1f520294]]
use crate::common::*;
use pyo3::prelude::*;
// 1f520294 ends here

// [[file:../spdkit-python.note::4c0ea8d1][4c0ea8d1]]
#[pyfunction]
/// Enable logging by setting verbosity level to `level`.
pub fn set_verbosity(level: u8) {
    let mut log = gut::cli::Verbosity::default();
    log.set_verbosity(level);
    log.setup_logger();
}

/// Parse a list of numbers from a readable string `s`.
///
/// "2-5"   ==> [2, 3, 4, 5]
/// "1,3-5" ==> [1, 3, 4, 5]
/// "1 3,5" ==> [1, 3, 4, 5]
#[pyfunction]
pub fn parse_numbers_human_readable(s: String) -> Result<Vec<usize>> {
    use gut::utils::parse_numbers_human_readable;

    let selected = parse_numbers_human_readable(&s)?;
    Ok(selected)
}

#[pyfunction]
/// Show how to select atoms in jmol selection commands
pub fn jmol_selection_commands(selected: Vec<usize>) {
    let selected = selected.iter().map(|x| x.to_string()).join(",");
    println!("selectionhalo");
    println!("select none");
    println!("select atomno=[{selected}]");
}
// 4c0ea8d1 ends here

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

// [[file:../spdkit-python.note::0853f16f][0853f16f]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyTemplate>()?;
    m.add_function(wrap_pyfunction!(probe_surface_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(set_verbosity, m)?)?;
    m.add_function(wrap_pyfunction!(parse_numbers_human_readable, m)?)?;
    m.add_function(wrap_pyfunction!(jmol_selection_commands, m)?)?;

    Ok(m)
}
// 0853f16f ends here
