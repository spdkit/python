// [[file:../spdkit-python.note::be102058][be102058]]
use pyo3::prelude::*;

use gut::prelude::*;
// be102058 ends here

// [[file:../spdkit-python.note::6b44aa6c][6b44aa6c]]
use gosh_database::DbConnection;

#[pyclass(name = "DbConnection")]
#[derive(Clone)]
pub struct PyDbConnection {
    pub inner: DbConnection,
}

#[pymethods]
impl PyDbConnection {
    #[new]
    pub fn connect(path: String) -> Result<Self> {
        let inner = DbConnection::connect(&path)?;
        let r = Self { inner };
        Ok(r)
    }
}
// 6b44aa6c ends here

// [[file:../spdkit-python.note::83f2f6c1][83f2f6c1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyDbConnection>()?;
    Ok(m)
}
// 83f2f6c1 ends here
