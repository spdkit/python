#! /usr/bin/env bash

## 0. load VASP related environment variables
source ~/apps/env.rc
apps load vasp

## 1. Prepare input files required for VASP calculation
# Create POSCAR from standard input (stdin), formatted with user template
# (input.tera)
cat > POSCAR

# prepare POTCAR from species in POSCAR
# ppp file layout: pot/{symbol}/POTCAR
if [ ! -f POTCAR ]; then
    apps load python
    # NOTE: this script is also in scripts/ dir
    generate-potcar-from-poscar.py
fi

## Prepare other input files
# copy important files into the .tmp* scratch directory, which will be
# automatically removed if job finished successful.
cp "$BBM_TPL_DIR/INCAR" .
cp "$BBM_TPL_DIR/KPOINTS" .

## 2. How to run vasp
# PLEASE CHANGE
# submit vasp, ignoring stdout produced from vasp exe.
hostname > "$BBM_JOB_DIR"/vasp.log
vasp630 >> "$BBM_JOB_DIR"/vasp.log

## 3. Post-processes
# VASP completed. Save intermediate structures to job starting dir
cp OUTCAR CONTCAR OSZICAR "$BBM_JOB_DIR/" 2>/dev/null

## 4. extract energy and forces from OUTCAR to stdout
gosh-adaptor vasp OUTCAR
