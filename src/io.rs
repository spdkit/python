// [[file:../spdkit-python.note::45ba73aa][45ba73aa]]
use crate::common::*;
use pyo3::prelude::*;
// 45ba73aa ends here

// [[file:../spdkit-python.note::4cfd5c0b][4cfd5c0b]]
use super::{PyMolecule, PyMoleculeIter};

#[pyfunction]
/// Read a list of `Molecule` from `path`. Returns an iterator over
/// `Molecule`, which allows reading a large file out of memory.
pub fn read(path: String) -> PyResult<PyMoleculeIter> {
    let mols = gchemol::io::read(path)?.map(|inner| PyMolecule { inner });
    let mols = PyMoleculeIter { iter: Box::new(mols) };
    Ok(mols)
}

/// Write molecules into path. File format will be determined according to the path.
#[pyfunction]
pub fn write(path: String, mols: Vec<PyMolecule>) -> PyResult<()> {
    gchemol::io::write(path, mols.iter().map(|mol| &mol.inner))?;
    Ok(())
}

/// Guess chemical file format from `path`
#[pyfunction]
pub fn guess_format_from_path(path: &str) -> Option<String> {
    let fmt = gchemol::io::guess_format_from_path(path.as_ref())?;
    Some(fmt)
}

#[pyfunction]
#[pyo3(signature = (pattern, /, root=".", recursive=true))]
#[pyo3(text_signature = "(pattern, root='.', recursive=True)")]
/// Recursively find all files in `root` dir with given file name
/// matching regex `pattern`
pub fn find_files(pattern: &str, root: &str, recursive: bool) -> Vec<String> {
    let root = root.as_ref();
    gchemol::io::find_files(pattern, root, recursive)
        .map(|p| p.to_string_lossy().to_string())
        .collect()
}
// 4cfd5c0b ends here

// [[file:../spdkit-python.note::377732f7][377732f7]]
use gchemol::io::Template;
use pyo3::types::PyDict;

#[pyclass(name = "Template")]
/// Template rendering using minijinja
pub struct PyTemplate {
    inner: Template<'static>,
}

#[pymethods]
impl PyTemplate {
    /// Render template using vars from keyword arguments
    #[pyo3(signature = (**keywords))]
    fn render(&self, keywords: Option<&PyDict>) -> Result<String> {
        let s = Python::with_gil(|py| -> PyResult<String> {
            let json = py.import("json")?;
            let s = json.call_method1("dumps", (keywords,))?;
            s.extract()
        })?;
        Ok(self.inner.render_json(&s)?)
    }

    /// Render template as string using vars from `json`.
    fn render_json(&self, json: &str) -> Result<String> {
        let s = self.inner.render_json(json)?;
        Ok(s)
    }

    /// Load template from file in `path`.
    #[staticmethod]
    fn from_file(path: &str) -> Result<Self> {
        let inner = Template::try_from_path(path.as_ref())?;
        Ok(Self { inner })
    }
}
// 377732f7 ends here

// [[file:../spdkit-python.note::c50e479a][c50e479a]]
use gchemol_parser::TextViewer;

/// A simple line-based text viewer for quick peeking part of text
#[pyclass(name = "TextViewer")]
#[derive(Clone)]
pub struct PyTextViewer {
    inner: TextViewer,
}

#[pymethods]
impl PyTextViewer {
    /// Create a `TextViewer` from text string in `s`.
    #[staticmethod]
    #[pyo3(text_signature = "($self, s)")]
    pub fn from_string(s: String) -> Self {
        let inner = TextViewer::from_str(&s);
        Self { inner }
    }

    /// Return total number of lines.
    pub fn num_lines(&self) -> usize {
        self.inner.num_lines()
    }

    /// Get line number at cursor.
    pub fn current_line_num(&self) -> usize {
        self.inner.current_line_num()
    }

    /// Return all the text.
    pub fn text(&self) -> String {
        self.inner.text().to_owned()
    }

    /// Return the line at cursor.
    pub fn current_line(&self) -> String {
        self.inner.current_line().to_owned()
    }

    /// Move the cursor to line `n`, counting from line 1 at beginning
    /// of the text.
    #[pyo3(text_signature = "($self, n)")]
    pub fn goto_line(&mut self, n: usize) {
        self.inner.goto_line(n);
    }

    /// Move the cursor to the beginning of the previous line.
    pub fn goto_previous_line(&mut self) {
        self.inner.goto_previous_line();
    }

    /// Move the cursor to the beginning of the next line.
    pub fn goto_next_line(&mut self) {
        self.inner.goto_next_line();
    }

    /// Move the cursor to the beginning of the last line.
    pub fn goto_last_line(&mut self) {
        self.inner.goto_last_line();
    }

    /// Move the cursor to the beginning of the first line.
    pub fn goto_first_line(&mut self) {
        self.inner.goto_first_line();
    }

    /// Move the cursor to the line matching `pattern`. Regex pattern
    /// is allowed.
    pub fn search_forward(&mut self, pattern: String) -> PyResult<usize> {
        let n = self.inner.search_forward(&pattern)?;
        Ok(n)
    }

    /// Search backward from current point for `pattern`. Return
    /// current point after search.
    pub fn search_backward(&mut self, pattern: String) -> PyResult<usize> {
        let n = self.inner.search_backward(&pattern)?;
        Ok(n)
    }

    /// Peek line `n`.
    pub fn peek_line(&self, n: usize) -> String {
        self.inner.peek_line(n).to_owned()
    }

    /// Peek the text between line `n` and `m` (including line `m`)
    pub fn peek_lines(&self, n: usize, m: usize) -> String {
        self.inner.peek_lines(n, m).to_owned()
    }

    /// Select the next `n` lines from current point, including current line.
    pub fn selection(&self, n: usize) -> String {
        self.inner.selection(n).to_owned()
    }

    /// Select part of the string in next `n` lines (including
    /// currrent line), in a rectangular area surrounded by columns in
    /// `col_beg`--`col_end`.
    pub fn column_selection(&self, n: usize, col_beg: usize, col_end: usize) -> String {
        self.inner.column_selection(n, col_beg, col_end)
    }
}
// c50e479a ends here

// [[file:../spdkit-python.note::048fa298][048fa298]]
use gchemol_parser::GrepReader;

/// Quick grep text by marking the line that matching a pattern,
/// suitable for very large text file.
#[pyclass(name = "GrepReader", subclass)]
#[pyo3(text_signature = "(path)")]
pub struct PyGrepReader {
    inner: GrepReader,
}

#[pymethods]
impl PyGrepReader {
    /// Create a `GrepReader` from file in `path`.
    #[new]
    pub fn from_file(path: String) -> PyResult<Self> {
        let inner = GrepReader::try_from_path(path.as_ref())?;
        Ok(Self { inner })
    }

    /// Mark positions that matching pattern, so that we can seek these
    /// positions later. Return the number of marked positions.
    #[pyo3(text_signature = "($self, pattern, max_count = None)")]
    pub fn mark(&mut self, pattern: String, max_count: Option<usize>) -> PyResult<usize> {
        let n = self.inner.mark(&pattern, max_count)?;
        Ok(n)
    }

    /// Goto the start of the inner file.
    pub fn goto_start(&mut self) {
        self.inner.goto_start();
    }
    /// Goto the end of the inner file.
    pub fn goto_end(&mut self) {
        self.inner.goto_end();
    }

    /// Return the number of marked positions
    pub fn num_markers(&self) -> usize {
        self.inner.num_markers()
    }

    /// Goto the next position that marked. Return marker position on success.
    /// Return Err if already reached the last marker or other errors.
    pub fn goto_next_marker(&mut self) -> PyResult<u64> {
        let n = self.inner.goto_next_marker()?;
        Ok(n)
    }

    /// Goto the marked position in `marker_index`. Will panic if marker_index
    /// out of range.
    pub fn goto_marker(&mut self, marker_index: isize) -> PyResult<u64> {
        let i = if marker_index < 0 {
            self.inner.num_markers() as isize + marker_index
        } else {
            marker_index
        };
        let n = self.inner.goto_marker(i as usize)?;
        Ok(n)
    }

    /// Return `n` lines in string on success from current
    /// position. Return error if reached EOF early.
    pub fn read_lines(&mut self, n: usize) -> PyResult<String> {
        let mut s = String::new();
        self.inner.read_lines(n, &mut s)?;
        Ok(s)
    }

    /// View next `n` lines like in a normal text viewer. This method
    /// will forward the cursor by `n` lines.
    pub fn view_lines(&mut self, n: usize) -> PyResult<PyTextViewer> {
        let inner = self.inner.view_lines(n)?;
        Ok(PyTextViewer { inner })
    }

    /// Return text from current position to the next marker or file
    /// end. It method will forward the cursor to the next marker.
    pub fn read_until_next_marker(&mut self) -> PyResult<String> {
        let mut s = String::new();
        self.inner.read_until_next_marker(&mut s)?;
        Ok(s)
    }

    /// View all lines until next marker like in a normal text viewer.
    /// It method will forward the cursor to the next marker.
    pub fn view_until_next_marker(&mut self) -> PyResult<PyTextViewer> {
        let inner = self.inner.view_until_next_marker()?;
        Ok(PyTextViewer { inner })
    }
}
// 048fa298 ends here

// [[file:../spdkit-python.note::d5cde32b][d5cde32b]]
fn plot_3d(z: Vec<Vec<f64>>, x: Vec<f64>, y: Vec<f64>) {
    use plotly::layout::Axis;
    use plotly::*;

    let mut plot = Plot::new();

    let trace1 = Surface::new(z).x(x).y(y).cauto(false);
    plot.add_trace(trace1);
    plot.show();
}

#[pyclass]
#[derive(Debug, Default)]
struct OptimizationTrajactory {
    energy: Vec<f64>,
    fmax: Vec<f64>,
}

#[pymethods]
impl OptimizationTrajactory {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an image data in `energy` for total energy and `fmax`
    /// for max force.
    pub fn append_data(&mut self, energy: f64, fmax: f64) {
        self.energy.push(energy);
        self.fmax.push(fmax);
    }

    /// Plot using plotly and write html file in `path`.
    pub fn plot(&self, path: String) {
        use plotly::common::{AxisSide, Mode};
        use plotly::layout::Axis;
        use plotly::{Layout, Plot, Scatter};

        let n = self.energy.len();
        let x = (0..n).collect_vec();
        let y = self.energy.clone();
        let trace1 = Scatter::new(x.clone(), y).name("energy").mode(Mode::Lines).y_axis("y1");
        let y = self.fmax.clone();
        let trace2 = Scatter::new(x.clone(), y).name("fmax").mode(Mode::Lines).y_axis("y2");
        let mut plot = Plot::new();
        plot.add_trace(trace1);
        plot.add_trace(trace2);
        let layout = Layout::new()
            .y_axis2(Axis::new().side(AxisSide::Right).overlaying("y"))
            .title("<b>Optimization trajectory</b>".into());
        plot.set_layout(layout);

        plot.use_local_plotly();
        plot.write_html(&path);

        info!("plot was wrote to {path}");
    }

    /// Return string for pretty printing table.
    pub fn table(&self) -> String {
        use tabled::TableIteratorExt;
        use tabled::{Style, Table, Tabled};

        #[derive(Tabled)]
        struct OptHist {
            step: usize,
            energy: f64,
            fmax: f64,
        }

        let hist = self.energy.iter().zip(&self.fmax).enumerate().map(|(i, (e, f))| OptHist {
            step: i,
            energy: *e,
            fmax: *f,
        });

        Table::new(hist).with(Style::markdown()).to_string()
    }
}
// d5cde32b ends here

// [[file:../spdkit-python.note::ca036abe][ca036abe]]
pub fn new<'p>(py: Python<'p>, name: &str) -> PyResult<&'p PyModule> {
    let m = PyModule::new(py, name)?;
    m.add_class::<PyTemplate>()?;
    m.add_class::<PyTextViewer>()?;
    m.add_class::<PyGrepReader>()?;
    m.add_class::<OptimizationTrajactory>()?;
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    m.add_function(wrap_pyfunction!(find_files, m)?)?;
    m.add_function(wrap_pyfunction!(guess_format_from_path, m)?)?;

    Ok(m)
}
// ca036abe ends here
