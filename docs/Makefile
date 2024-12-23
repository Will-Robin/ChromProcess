CHROMPROCESS=../src/ChromProcess
CHROMPROCESS_MOD=../src/ChromProcess/*.py

CLASSES=$(CHROMPROCESS)/Classes

LOADING=$(CHROMPROCESS)/Loading

WRITERS=$(CHROMPROCESS)/Writers

PROCESSING=$(CHROMPROCESS)/Processing

UTILS=$(CHROMPROCESS)/Utils

CHROMPROCESS_PACKAGE=$(LOADING) $(WRITERS) $(PROCESSING) $(UTILS) $(CHROMPROCESS_MOD)

all: index.html

index.html: $(CHROMPROCESS_PACKAGE)
	pdoc ../ChromProcess --docformat "numpy" -o ./

clean:
	rm *.html
	rm ChromProcess/*.html
	rm search.js
