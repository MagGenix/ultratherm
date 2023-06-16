# Ultratherm
Used for design of DNA/RNA temperature sensitive structures, including DNA thermometers, RNA thermometers, and DNA:RNA thermosensitive heteroduplexes.

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

Install nuad now that both nupack and vienna are installed, with the environment active.

```
pip3 install nuad
```