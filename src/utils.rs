// [[file:../spdkit-python.note::1f520294][1f520294]]
use crate::common::*;
use pyo3::prelude::*;
// 1f520294 ends here

// [[file:../spdkit-python.note::4c0ea8d1][4c0ea8d1]]
#[pyfunction]
/// Enable logging by setting verbosity level to `level`.
/// 0 => warn, 1 => info, 2 => debug, 3 or above => trace
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
/// Make an abbreviation of the long number list. Return the string
/// representation. For example: 1,2,3,6,7,8,9 ==> 1-3,6-9
pub fn abbreviate_numbers_human_readable(s: Vec<usize>) -> Result<String> {
    let selection = gut::utils::abbreviate_numbers_human_readable(&s)?;
    Ok(selection)
}

#[pyfunction]
/// Show how to select atoms using commands in jmol script console.
pub fn jmol_selection_commands(selected: Vec<usize>) {
    let selected = selected.iter().map(|x| x.to_string()).join(",");
    println!("selectionhalo");
    println!("select none");
    println!("select atomno=[{selected}]");
}

#[pyfunction]
/// Show how to select atoms using commands in pymol.
pub fn pymol_selection_commands(selected: Vec<usize>) {
    let selected = selected.iter().map(|x| x.to_string()).join("+");
    println!("select selection, id {selected}");
}

#[pyfunction]
/// Show codes to select atoms using GaussView
fn gaussview_selection_commands(selected: Vec<usize>) -> Result<()> {
    info!("Open Atom Selection dialog, and input the following code");
    println!("{}", gut::utils::abbreviate_numbers_human_readable(selected.as_slice())?);
    Ok(())
}
// 4c0ea8d1 ends here

// [[file:../spdkit-python.note::22b53bc5][22b53bc5]]
use super::PyMolecule;

#[pyfunction]
/// Rebond all molecules in parallel
fn par_rebond(mols: Vec<&PyCell<PyMolecule>>) -> Result<()> {
    use std::ops::DerefMut;

    // workaround for limitation in &mut PyMolecule
    let mut mols: Vec<_> = mols.into_iter().map(|cell| cell.try_borrow_mut().unwrap()).collect();
    let mut mols: Vec<_> = mols.iter_mut().map(|refr| refr.deref_mut()).collect();

    mols.as_mut_slice()
        .into_par_iter()
        .for_each(|mut mol| mol.rebond(false, None, None, None, None));

    Ok(())
}
// 22b53bc5 ends here

// [[file:../spdkit-python.note::0853f16f][0853f16f]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyTemplate>()?;
    m.add_function(wrap_pyfunction!(set_verbosity, m)?)?;
    m.add_function(wrap_pyfunction!(parse_numbers_human_readable, m)?)?;
    m.add_function(wrap_pyfunction!(abbreviate_numbers_human_readable, m)?)?;
    m.add_function(wrap_pyfunction!(jmol_selection_commands, m)?)?;
    m.add_function(wrap_pyfunction!(pymol_selection_commands, m)?)?;
    m.add_function(wrap_pyfunction!(gaussview_selection_commands, m)?)?;
    m.add_function(wrap_pyfunction!(par_rebond, m)?)?;

    Ok(m)
}
// 0853f16f ends here
