# Ultratherm
Used for design of DNA/RNA temperature sensitive structures, including DNA thermometers, RNA thermometers, and DNA:DNA or RNA:RNA thermosensitive duplexes.

Design is structure-agnostic. Scoring is performed either with NUPACK or Vienna.
NUPACK will not perform pseudoknot structure prediction.

## Installing
Use the environment.yml to make an appropriately-named conda env (prefer mamba).
*Technical aside: Vienna on bioconda is nonfunctional so the environment.yml specifies the PyPi version.*

```
mamba env create -f environment.yml
conda activate ultratherm
```

Install nupack from local using pip3 with the env active. It should appear in the env list.

```
python3 -m pip install -U nupack -f ./nupack-4.0.1.8/package
```

As of NUPACK 4.0 Python >3.7 is required (a standard installation will use this).

**Operation is not guaranteed with any other version of Python!**

## Development
pdoc is used to generate documentation.

Install pdoc to the *same* conda environment, make sure main() and test() are commented out in main.py, and run:
```
pdoc --docformat google -o docs *.py
```
