# Ultratherm Data Specification
Sreenivas Eadara

## Purpose
Human-readable sequence formats are required to allow users to view the results of ultratherm.<br>
However, ultratherm (and perhaps other programs) may require additional data associated with each nucleotide of each sequence, which are required in order for a human-readable sequence format to be used as an input file. This specification is designed to utilize the FASTQ standard format to store sequences and 5 bits of data per nucleotide. The protocol may be expanded in the future.

## Specification
The FASTQ format allows phred_quality scores, represented as a string of ASCII characters between ! and ~, to be stored alongside nucleotide sequences, such that each nucleotide has an associated quality (from 0 to 93 inclusive).


In this specification, only a subset of the FASTQ characters (characters between ! and ~) are to be used. Only alphanumeric characters are permitted. 5 consecutive bits are encoded per character in b1, b2, b3, b4, b5. The sequence of bits is as follows:

`[10XXXXX]`

Where `XXXXX` are the 5 least significant bits.

In ultratherm, this will correspond to the following (0: false, 1: true)<br>
(note that ASCII is 7-bit. When stored in a byte the MSB is always 0).

| b5              | b4        | b3           | b2        | b1        |
|-----------------|-----------|--------------|-----------|-----------|
| 0 (False)       | 0 (False) | 0 (False)    | 0 (False) | 0 (False) |
| 1 (True)        | 1 (True)  | 1 (True)     | 1 (True)  | 1 (True)  |
| score_this_nucl | is_rna    | score_region | no_indel  | no_mod    |

## PacBio PHRED compatibility

It is possible to map each byte of data (i.e. each ASCII character) to a PHRED score. The ASCII ‘@’, which is the first character in this specification, corresponds to ASCII 64 and PHRED score 31. In order to obtain a PHRED score from a list of bits b1-b5, one just needs to convert the list of bits to a 5-bit integer and add 31, obtaining the PHRED score.

## Multiplex Sequences

This specification supports the addition of a specific gap character such that multiple sequences (i.e. hybrids) can be stored within the same record. This is accomplished with the ‘&’ character as follows:

```
@16 0.9314516162927949
GAAGCAGACAAGAGCAGUCCUAGUCCUGAUUGCCAAAAGGCGAAUUAUUA&GGACUAGCUCUU
+
___XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX&LLLLLLLLLLLL
```

This can be accomplished only because & is a permissible PHRED score (5) and such characters are permitted in BioPython as well as other parsers. The ‘&’ character was selected because ViennaRNA permits multiple strand entry using the ‘&’ character and this increases the convenience for manual review of ultratherm sequences.
