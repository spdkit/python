// [[file:../spdkit-python.note::be102058][be102058]]
use pyo3::prelude::*;

use gut::prelude::*;

use super::PyMolecule;
// be102058 ends here

// [[file:../spdkit-python.note::6b44aa6c][6b44aa6c]]
use gosh::db::DbConnection;

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

// [[file:../spdkit-python.note::282503ba][282503ba]]
use gosh::db::prelude::Checkpoint;

#[pyclass(name = "Computed", subclass)]
#[derive(Clone)]
pub struct PyComputed {
    inner: Computed,
}

#[pymethods]
impl PyComputed {
    /// Get computed energy.
    pub fn get_energy(&self) -> Option<f64> {
        self.inner.get_energy()
    }

    /// Get computed forces.
    pub fn get_forces(&self) -> Option<Vec<[f64; 3]>> {
        self.inner.get_forces().map(|x| x.to_vec())
    }

    /// Get molecule structure.
    pub fn get_molecule(&self) -> Option<PyMolecule> {
        let mol = self.inner.get_molecule()?;
        Some(PyMolecule { inner: mol.to_owned() })
    }

    /// Converts computed results to string representation.
    pub fn to_string(&self) -> String {
        self.inner.to_string()
    }

    /// Converts computed results to json string.
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(&self.inner).unwrap()
    }

    // db related methods, which can only be put here for pyo3 limitations
    #[pyo3(text_signature = "($self, db)")]
    pub fn commit_checkpoint(&self, db: PyDbConnection) -> Result<()> {
        self.inner.commit_checkpoint(&db.inner)?;
        Ok(())
    }

    #[staticmethod]
    #[pyo3(text_signature = "(db, slot)")]
    /// Load computed from checkpoint `db` in slot `n`
    pub fn from_checkpoint_n(db: PyDbConnection, slot: i32) -> Result<Self> {
        let inner = Computed::from_checkpoint_n(&db.inner, slot)?;
        Ok(Self { inner })
    }

    #[staticmethod]
    #[pyo3(text_signature = "(db)")]
    /// Load latest computed from checkpoint `db`
    pub fn from_checkpoint(db: PyDbConnection) -> Result<Self> {
        Self::from_checkpoint_n(db, -1)
    }

    #[staticmethod]
    #[pyo3(text_signature = "(db)")]
    /// List available checkpoints in `db`.
    pub fn list_checkpoints(db: PyDbConnection) -> Result<()> {
        Computed::list_checkpoints(&db.inner)?;
        Ok(())
    }
}
// 282503ba ends here

// [[file:../spdkit-python.note::a6cc3f1c][a6cc3f1c]]
use gosh::model::{BlackBoxModel, ChemicalModel, Computed};

#[pyclass(name = "BlackBoxModel", subclass)]
#[pyo3(text_signature = "(dir)")]
pub struct PyBlackBoxModel {
    inner: BlackBoxModel,
}

#[pymethods]
impl PyBlackBoxModel {
    #[new]
    /// Construct BlackBoxModel model under directory context.
    pub fn from_dir(dir: String) -> Result<Self> {
        let inner = BlackBoxModel::from_dir(dir)?;
        let s = Self { inner };
        Ok(s)
    }

    /// Render `mol` to input string using this template.
    #[pyo3(text_signature = "($self, mol)")]
    pub fn render_input(&self, mol: PyMolecule) -> Result<String> {
        let r = self.inner.render_input(&mol.inner)?;
        Ok(r)
    }

    /// Compute `mol` using this model for properties.
    #[pyo3(text_signature = "($self, mol)")]
    pub fn compute(&mut self, mol: &PyMolecule) -> Result<PyComputed> {
        let inner = self.inner.compute(&mol.inner)?;
        Ok(PyComputed { inner })
    }

    /// Return the number of potentail evaluations
    pub fn number_of_evaluations(&self) -> usize {
        self.inner.number_of_evaluations()
    }

    /// Return bash script suitable for execution.
    #[pyo3(text_signature = "($self, mol)")]
    pub fn bash_script_for_execution(&mut self, mol: &PyMolecule) -> Result<String> {
        let s = self.inner.bash_script_for_execution(&mol.inner)?;
        Ok(s)
    }
}
// a6cc3f1c ends here

// [[file:../spdkit-python.note::83f2f6c1][83f2f6c1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyBlackBoxModel>()?;
    m.add_class::<PyComputed>()?;
    m.add_class::<PyDbConnection>()?;
    Ok(m)
}
// 83f2f6c1 ends here
