// [[file:../spdkit-python.note::d6f87faa][d6f87faa]]
use crate::common::*;
use pyo3::prelude::*;
// d6f87faa ends here

// [[file:../spdkit-python.note::97de14c4][97de14c4]]
use super::PyMolecule;
use gchemol::prelude::*;
use gchemol::Molecule;
use serde_json::Value;
use std::path::Path;
use zip::ZipArchive;

#[pyclass]
pub struct UserData {
    vars: Value,
    input_files: Vec<(String, String)>,
}

#[pymethods]
impl UserData {
    /// Construct `UserData` from a zip archive file.
    #[staticmethod]
    fn from_zip(path: &str) -> Result<Self> {
        use std::io::prelude::*;

        let zipfile = std::fs::File::open(path)?;
        let mut archive = ZipArchive::new(zipfile)?;
        let f = archive.by_name("script.txt")?;
        let mut vars: serde_json::Value = serde_json::from_reader(f)?;
        let input_file_names: Vec<String> = serde_json::from_value(vars["input_files"].take())?;
        let mut input_files = vec![];
        for fname in input_file_names {
            let mut s = String::new();
            let mut f = archive.by_name(&fname)?;
            f.read_to_string(&mut s)?;
            input_files.push((fname, s));
        }

        Ok(Self { vars, input_files })
    }

    /// Return variables in json represents (file in zip
    /// "script.txt").
    fn get_json_vars(&self) -> String {
        self.vars.to_string()
    }

    /// Return `Molecule` object extracted from file specified in
    /// `input_files` var.
    fn get_molecules(&self) -> Result<Vec<PyMolecule>> {
        let mut mols = vec![];

        for (fname, source) in self.input_files.iter() {
            // FIXME: guess format from file name
            let fmt = if fname.ends_with(".cif") {
                "text/cif"
            } else if fname.ends_with(".xyz") {
                "text/xyz"
            } else if fname.ends_with(".mol2") {
                "text/mol2"
            } else {
                unimplemented!("unkown {fname:?}");
            };

            let inner = Molecule::from_str(source, &fmt)?;
            mols.push(PyMolecule { inner });
        }

        Ok(mols)
    }
}
// 97de14c4 ends here

// [[file:../spdkit-python.note::016f6fa1][016f6fa1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<UserData>()?;
    // m.add_function(wrap_pyfunction!(optimize, m)?)?;
    Ok(m)
}
// 016f6fa1 ends here
