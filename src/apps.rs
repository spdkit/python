// [[file:../spdkit-python.note::a1bd7b8a][a1bd7b8a]]
use crate::common::*;
use pyo3::prelude::*;
// a1bd7b8a ends here

// [[file:../spdkit-python.note::825e4aec][825e4aec]]
use super::PyMolecule;
use rxview::cli::ReactionPreview;

/// A tool to interpolate reaction path with provided molecule
/// structures, similar to `reaction preview` tool as in Material
/// Studio.
#[pyclass(name = "ReactionPreview")]
pub struct PyReactionPreview {
    inner: ReactionPreview,
}

#[pymethods]
impl PyReactionPreview {
    #[new]
    fn new() -> Self {
        Self {
            inner: ReactionPreview::default(),
        }
    }

    #[pyo3(signature = (mol1, mol2, /, nimages=10))]
    #[pyo3(text_signature = "($self, mol1, mol2, /, nimages=10)")]
    /// Create reaction path with `nimages` molecules connecting `mol1` and `mol2`.
    fn create_reaction_path(&self, mol1: &PyMolecule, mol2: &PyMolecule, nimages: usize) -> Result<Vec<PyMolecule>> {
        let mols = self.inner.create_reaction_path(&mol1.inner, &mol2.inner, nimages)?;
        let mols = mols.into_iter().map(|inner| PyMolecule { inner }).collect();
        Ok(mols)
    }

    /// Define freezing atom list in human-readable scheme when create
    /// reaction path.
    ///
    /// Represents input scheme modifying which atoms are to be frozen
    /// during optimization. For examples, "3-6,17", will freeze atoms
    /// 3,4,5,6,17. If `atoms=auto`, we will deduce freezing atoms
    /// from reactant and product positions.
    #[pyo3(signature = (atoms="auto"))]
    #[pyo3(text_signature = "($self, atoms='auto')")]
    fn set_freezing(&mut self, atoms: &str) {
        self.inner.set_freezing(atoms);
    }

    /// Enable experimental interpolation algorithm for periodic system.
    #[pyo3(signature = (pbc=true))]
    #[pyo3(text_signature = "($self, pbc=True)")]
    pub fn set_pbc(&mut self, pbc: bool) {
        self.inner.set_pbc(pbc);
    }
}
// 825e4aec ends here

// [[file:../spdkit-python.note::55ce1aa4][55ce1aa4]]
mod mhm {
    use super::*;
    use minima_hopping::{MinimaHopping, MinimaHoppingOptions};
    use pyo3::types::PyDict;

    #[pyclass(name = "MinimaHopping", subclass)]
    pub struct PyMinimaHopping {
        options: MinimaHoppingOptions,
    }

    #[pymethods]
    impl PyMinimaHopping {
        #[staticmethod]
        #[pyo3(signature = (mol, bbm_dir, out, /, **options))]
        fn run(mol: &mut PyMolecule, bbm_dir: &str, out: &str, options: Option<&PyDict>) -> Result<()> {
            let mhm_options = if let Some(options) = options {
                pythonize::depythonize(options.as_ref())?
            } else {
                MinimaHoppingOptions::default()
            };
            mhm_options.validate()?;
            let mut mhm = MinimaHopping::new(mhm_options);
            mhm.run(&mut mol.inner, bbm_dir.as_ref(), out)?;
            Ok(())
        }
    }
}
// 55ce1aa4 ends here

// [[file:../spdkit-python.note::9d891d5f][9d891d5f]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyReactionPreview>()?;
    m.add_class::<mhm::PyMinimaHopping>()?;
    // m.add_function(wrap_pyfunction!(read, m)?)?;

    Ok(m)
}
// 9d891d5f ends here
