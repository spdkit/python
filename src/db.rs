// [[file:../spdkit-python.note::be102058][be102058]]
use pyo3::prelude::*;

use gosh_database::CheckpointDb;
use gut::prelude::*;
// be102058 ends here

// [[file:../spdkit-python.note::5e90626e][5e90626e]]
use super::{Computed, PyComputed};

#[pyclass(name = "CheckpointDb")]
pub struct PyCheckpointDb {
    inner: CheckpointDb,
}

#[pymethods]
impl PyCheckpointDb {
    #[new]
    /// Construct Checkpoint from `path` to a file.
    pub fn new(path: String) -> Self {
        Self {
            inner: CheckpointDb::new(path),
        }
    }

    /// Load latest computed from checkpoint
    pub fn load_from_latest(&self) -> Result<PyComputed> {
        let inner = self.inner.load_from_latest::<Computed>()?;
        Ok(PyComputed { inner })
    }

    /// Load computed from checkpoint in slot `n`
    pub fn load_from_slot_n(&self, n: i32) -> Result<PyComputed> {
        let inner = self.inner.load_from_slot_n::<Computed>(n)?;
        Ok(PyComputed { inner })
    }

    /// Commit a checkpoint into database. Return true if committed, false otherwise.
    pub fn commit(&self, computed: super::PyComputed) -> Result<bool> {
        let r = self.inner.commit(&computed.inner)?;
        Ok(r)
    }

    pub fn list(&self) -> Result<bool> {
        let r = self.inner.list::<Computed>()?;
        Ok(r)
    }
}
// 5e90626e ends here

// [[file:../spdkit-python.note::83f2f6c1][83f2f6c1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    // m.add_class::<PyChemicalEnvironment>()?;
    Ok(m)
}
// 83f2f6c1 ends here
