# PBAlgebras --- build the paper, run the checks, regenerate the tables.

PYTHON ?= python3
PAPER  := paper/main.pdf

.PHONY: all paper test verify quick clean distclean help

help:
	@echo "make paper    build paper/main.pdf"
	@echo "make test     run the verification test suite"
	@echo "make verify   recompute critical constants, regenerate the paper's table"
	@echo "make quick    a CI-sized verify"
	@echo "make all      verify + paper"
	@echo "make clean    remove LaTeX build products"

all: verify paper

paper: $(PAPER)

$(PAPER): paper/main.tex paper/pbalgebras.sty paper/refs.bib \
          $(wildcard paper/sections/*.tex) paper/generated/constants.tex
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

test:
	cd verify && $(PYTHON) -m pytest tests -q

verify:
	cd verify && $(PYTHON) experiments/run_all.py

quick:
	cd verify && $(PYTHON) experiments/run_all.py --quick

t2:
	cd verify && $(PYTHON) experiments/t2_exact.py

clean:
	cd paper && latexmk -C >/dev/null 2>&1 || true
	rm -rf verify/**/__pycache__ verify/.pytest_cache

distclean: clean
	rm -f $(PAPER)
