// [[file:../spdkit-python.note::5063691f][5063691f]]
use crate::common::*;
use pyo3::prelude::*;

use super::PyMolecule;
// 5063691f ends here

// [[file:../spdkit-python.note::c93dddfc][c93dddfc]]
use reqwest::blocking::Client;

/// A client for remote view using gchemol-view
#[pyclass(name = "GchemolViewClient")]
#[pyo3(text_signature = "(symbol, position=[0, 0, 0])")]
pub struct PyGchemolViewClient {
    client: Client,
    server: String,
}

#[pymethods]
impl PyGchemolViewClient {
    #[new]
    #[pyo3(signature = (port = 3039))]
    fn new(port: usize) -> Self {
        let server = format!("127.0.0.1:{port}");
        info!("connecting to {server:?}");
        let client = reqwest::blocking::Client::builder().build().expect("reqwest client");
        Self {
            client: client.into(),
            server,
        }
    }

    /// Load mols for remote view using gchemol-view services
    ///
    /// # NOTE
    ///   * currently only the first molecule in the list will be sent
    fn load(&mut self, mols: Vec<PyMolecule>) -> Result<()> {
        let server = &self.server;
        let uri = format!("http://{server}/view-molecule");
        if !mols.is_empty() {
            let mol = &mols[0];
            let resp = self.client.post(&uri).json(&mol.inner).send()?.text()?;
            println!("{resp:?}");
        } else {
            eprintln!("molecule list is empty.");
        }
        Ok(())
    }
}
// c93dddfc ends here

// [[file:../spdkit-python.note::18daa353][18daa353]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyGchemolViewClient>()?;
    // m.add_function(wrap_pyfunction!(read, m)?)?;

    Ok(m)
}
// 18daa353 ends here
