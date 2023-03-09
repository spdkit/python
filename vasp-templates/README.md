This directory contains templates using [minijinja](https://docs.rs/minijinja/latest/minijinja/syntax/index.html) for VASP
calculations refactored from [vaspkit](https://vaspkit.com) generated INCAR files. You can
play with MiniJinja online in the [browser playground](https://mitsuhiko.github.io/minijinja-playground/).

Directory layout

-   sp: single point static calculation
-   opt: structure optimization
-   md: molecule dynamics
-   freq: frequency calculation

Example usage in spdkit-python scripting

    from spdkit import io, gosh, Molecule
    
    # initialize minijinja template from file
    t = io.Template.from_file("/path/to/vasp/sp/INCAR.jinja")
    # change VASP computational keywords dynamically
    s = t.render(ENCUT=500, ISPIN=2)
    # write it to INCAR file in BBM directory
    open("/path/to/vasp/sp/INCAR", "w").write(s)
    
    # construct Molecule object from file
    m = Molecule.from_file("path/to/mol-file")
    # construct BBM from the directory
    vasp = gosh.BlackBoxModel("/path/to/vasp/sp")
    # call vasp bbm compute molecule `m`
    r = vasp.compute(m)
    # get energy for further process
    energy = r.get_energy()

