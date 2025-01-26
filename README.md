# Ultratherm
Used for design of DNA/RNA temperature sensitive structures, including DNA thermometers, RNA thermometers, and DNA:DNA or RNA:RNA thermosensitive duplexes.

Design is structure-agnostic. Scoring is performed with Vienna by default.
NUPACK may also be used for scoring, though it is not included with a default installation of Ultratherm (users must supply their own copy of NUPACK).
Note that NUPACK will not perform pseudoknot structure prediction.

## Installing
Use the environment.yml to make an appropriately-named conda env (prefer mamba).
*Technical aside: Vienna on bioconda is nonfunctional so the environment.yml specifies the PyPi version.*

```
mamba env create -f environment.yml
conda activate ultratherm
```

### Adding NUPACK to your Ultratherm environment (optional)
As a reminder, NUPACK is not included by the default installation above.
If you intend to use NUPACK for scoring, you must supply your own files for NUPACK as follows:

Install nupack from local using pip3 with the env active. It should appear in the env list.

```
python3 -m pip install -U nupack -f ./nupack-4.0.1.8/package
```

### Python version requirements
As of NUPACK 4.0 and ViennaRNA 2.6.4 Python >=3.8 is required (a standard installation will use 3.10).
*Note: NUPACK 4.x is not distributed with Windows wheels, so using NUPACK as an Ultratherm scorer on Windows is not guaranteed. You may consider using WSL or compiling NUPACK from your source.*

**Operation is not guaranteed with any other version of Python!**

## Development
pdoc is used to generate documentation.

Install pdoc to the *same* conda environment, make sure main() and test() are commented out in main.py, and run:
```
pdoc --docformat google -o docs *.py
```

## Contributors
Sreenivas Eadara developed this software while at MagGenix, Inc. from 2023 to February 2025.

## Acknowledgements

Thanks to "ViennaRNA Package 2.0" authors Ronny Lorenz, Stephan H Bernhart, Christian Höner zu Siederdissen, Hakim Tafer, Christoph Flamm, Peter F Stadler & Ivo L Hofacker.<br>
Source: (https://github.com/ViennaRNA/ViennaRNA)<br>
License: (https://github.com/ViennaRNA/ViennaRNA/blob/master/license.txt)
