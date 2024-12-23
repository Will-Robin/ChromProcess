.PHONY: tutorials docs format type_check

docs: docs/index.html
	cd docs && make

tutorials:
	jupyter nbconvert --execute --to notebook --inplace Tutorials/ChromProcess_Introduction_1.ipynb
	jupyter nbconvert --execute --to notebook --inplace Tutorials/ChromProcess_Introduction_2.ipynb
	jupyter nbconvert --execute --to notebook --inplace Tutorials/ChromProcess_Introduction_3_deconvolution.ipynb
	jupyter nbconvert --execute --to notebook --inplace Tutorials/ChromProcess_Tutorial_5_Background_subtraction.ipynb

format:
	ruff format src/ChromProcess

type_check:
	mypy src/ChromProcess

lint:
	ruff check src/ChromProcess
