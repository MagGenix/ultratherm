# Ultratherm
Used for design of DNA/RNA temperature sensitive structures, including DNA thermometers, RNA thermometers, and DNA:RNA thermosensitive heteroduplexes.

*Note: so far, nupack folding of non-pseudoknotted DNA and RNA thermosensitive structures has been implemented.*
*Vienna implementation to come.*

Design is structure-agnostic. Scoring is performed either with NUPACK or Vienna.
NUPACK will not perform pseudoknot structure prediction and cannot fold heteroduplices.

## Installing
Use the environment.yml to make an appropriately-named conda env (prefer mamba).

```
mamba env create -f environment.yml
conda activate ultratherm
```

Install nupack from local using pip3 with the env active. It should appear in the env list.

```
pip3 install -U nupack -f ./nupack-4.0.1.8/package
```

Install the ViennaRNA package using instructions provided by ViennaRNA to get RNAlib2.6.2 (or newer).
Note that this requires Python 3.9, so the environment has been updated to use Python3.9.

In my case (Darwin x86), I downloaded and installed precompiled binaries provided at:
https://www.tbi.univie.ac.at/RNA/

This completed with the following note:
>In case you have installed the Perl or Python bindings for RNAlib, you might want to adjust your PERL5LIB and/or PYTHONPATH, respectively. You can do this by adding the following lines to your ~/.bash_profile, or ~/.bashrc: 

```
export PYTHONPATH=/usr/local/lib/python3.9/site-packages:${PYTHONPATH}
```