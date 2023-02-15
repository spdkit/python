// [[file:../spdkit-python.note::45ba73aa][45ba73aa]]
use crate::common::*;
use pyo3::prelude::*;
// 45ba73aa ends here

// [[file:../spdkit-python.note::377732f7][377732f7]]
use gchemol::io::Template;

#[pyclass(name = "Template")]
pub struct PyTemplate {
    inner: Template<'static>,
}

#[pymethods]
impl PyTemplate {
    /// Render template as string using vars from `json`.
    fn render_json(&self, json: &str) -> Result<String> {
        let s = self.inner.render_json(json)?;
        Ok(s)
    }

    /// Load template from file in `path`.
    fn from_file(&self, path: &str) -> Result<Self> {
        let inner = Template::try_from_path(path.as_ref())?;
        Ok(Self { inner })
    }
}
// 377732f7 ends here

// [[file:../spdkit-python.note::ca036abe][ca036abe]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyTemplate>()?;
    // m.add_function(wrap_pyfunction!(optimize, m)?)?;
    Ok(m)
}
// ca036abe ends here
