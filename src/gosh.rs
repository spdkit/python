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
use gosh::remote::JobHub;

#[pyclass(name = "BlackBoxModel", subclass)]
#[pyo3(text_signature = "(dir)")]
pub struct PyBlackBoxModel {
    inner: BlackBoxModel,
    jobhub: Option<JobHub>,
}

impl PyBlackBoxModel {
    fn get_jobhub(&mut self) -> Result<&mut JobHub> {
        let jobhub = self.jobhub.as_mut().ok_or(anyhow!("no scheduler for lazy computation"))?;
        Ok(jobhub)
    }
}

#[pymethods]
impl PyBlackBoxModel {
    #[new]
    /// Construct BlackBoxModel model under directory context.
    pub fn from_dir(dir: String) -> Result<Self> {
        let inner = BlackBoxModel::from_dir(dir)?;
        let s = Self { inner, jobhub: None };
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

    /// Set job execution scheduler from address in `scheduler_address`.
    #[pyo3(text_signature = "($self, scheduler_address)")]
    pub fn set_job_scheduler(&mut self, scheduler_address: String) {
        self.jobhub = JobHub::new(&scheduler_address).into();
    }

    /// Compute `mol` using this model for properties.
    #[pyo3(text_signature = "($self, mol)")]
    pub fn compute_lazily(&mut self, mol: &PyMolecule) -> Result<usize> {
        let s = self.bash_script_for_execution(&mol)?;
        let jobhub = self.get_jobhub()?;
        let n = jobhub.add_job(s);
        Ok(n)
    }

    /// Return computed result for job `jobid`
    #[pyo3(text_signature = "($self, jobid)")]
    pub fn get_computed_result(&mut self, jobid: usize) -> Result<PyComputed> {
        let jobhub = self.get_jobhub()?;
        let inner = jobhub.get_computed_from_job_out(jobid)?;
        let r = PyComputed { inner };
        Ok(r)
    }

    /// Compute all lazily scheduled jobs in parallel
    pub fn compute_scheduled(&mut self) -> Result<()> {
        let jobhub = self.get_jobhub()?;
        jobhub.run()?;
        Ok(())
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

// [[file:../spdkit-python.note::59752afa][59752afa]]
/// A job hub for parallel running of multiple jobs over remote
/// computational nodes
#[pyclass(name = "JobHub", subclass)]
#[pyo3(text_signature = "(scheduler_address)")]
pub struct PyJobHub {
    inner: JobHub,
}

#[pymethods]
impl PyJobHub {
    /// Create a job hub for background scheduler specified in
    /// `scheduler_address`.
    #[new]
    pub fn new(scheduler_address: String) -> Self {
        let inner = JobHub::new(&scheduler_address);
        Self { inner }
    }

    /// Return the number of parallel threads.
    #[staticmethod]
    pub fn num_threads() -> usize {
        JobHub::num_threads()
    }

    /// Add a new job into job hub for scheduling.
    #[pyo3(text_signature = "($self, cmd)")]
    pub fn add_job(&mut self, cmd: String) -> usize {
        self.inner.add_job(cmd)
    }

    /// Run all scheduled jobs with nodes in pool.
    pub fn run(&mut self) -> Result<()> {
        self.inner.run()?;
        Ok(())
    }

    /// Return raw output for job `jobid`.
    #[pyo3(text_signature = "($self, jobid)")]
    pub fn get_job_out(&mut self, jobid: usize) -> Result<String> {
        let s = self.inner.get_job_out(jobid)?;
        Ok(s)
    }
}
// 59752afa ends here

// [[file:../spdkit-python.note::83f2f6c1][83f2f6c1]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyBlackBoxModel>()?;
    m.add_class::<PyComputed>()?;
    m.add_class::<PyDbConnection>()?;
    m.add_class::<PyJobHub>()?;
    Ok(m)
}
// 83f2f6c1 ends here
