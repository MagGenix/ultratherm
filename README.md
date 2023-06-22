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
```

Install nupack from local using pip3 with the env active. It should appear in the env list.

```
pip3 install -U nupack -f ./nupack-4.0.1.8/package
```