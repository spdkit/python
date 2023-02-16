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
// 4c0ea8d1 ends here

// [[file:../spdkit-python.note::0853f16f][0853f16f]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyTemplate>()?;
    m.add_function(wrap_pyfunction!(set_verbosity, m)?)?;
    m.add_function(wrap_pyfunction!(parse_numbers_human_readable, m)?)?;
    m.add_function(wrap_pyfunction!(abbreviate_numbers_human_readable, m)?)?;
    m.add_function(wrap_pyfunction!(jmol_selection_commands, m)?)?;

    Ok(m)
}
// 0853f16f ends here
