# [[file:spdkit-python.note::91bf6351][91bf6351]]
build-python:
	rm -f /scratch/cargo/structure-predication/wheels/*.whl
	maturin build --release
	pip install --force-reinstall /scratch/cargo/structure-predication/wheels/*.whl

build-docs:
	cd docs && make singlehtml
# 91bf6351 ends here
