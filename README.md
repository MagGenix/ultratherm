# Ultratherm
Used for design of DNA/RNA temperature sensitive structures, including DNA thermometers, RNA thermometers, and DNA:RNA thermosensitive heteroduplexes.

*Note: so far, nupack and vienna folding of non-pseudoknotted DNA and RNA thermosensitive structures has been implemented. Homo/heteroduplex folding is to come.*

Design is structure-agnostic. Scoring is performed either with NUPACK or Vienna.
NUPACK will not perform pseudoknot structure prediction and cannot fold heteroduplices.

## Installing
Use the environment.yml to make an appropriately-named conda env (prefer mamba).
*Technical aside: Vienna on bioconda is nonfunctional so the environment.yml specifies the PyPi version.*

```
mamba env create -f environment.yml
conda activate ultratherm
```

Install nupack from local using pip3 with the env active. It should appear in the env list.

```
pip3 install -U nupack -f ./nupack-4.0.1.8/package
```

As of RNAlib2.6.2 Python 3.9 is required (a standard installation will use this).

**Operation is not guaranteed with any other version of Python!**

## Development
pdoc is used to generate documentation.

Install pdoc to the *same* conda environment, make sure main() and test() are commented out in main.py, and run:
```
pdoc --docformat google -o docs *.py
```
