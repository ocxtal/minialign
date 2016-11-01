
# minialign

Minialign is a little bit fast and moderately accurate nucleotide sequence alignment tool designed for PacBio and Nanopore long reads. It is built on three key algorithms, minimizer-based index of the [minimap](https://github.com/lh3/minimap) overlapper, array-based seed chaining, and SIMD-parallel Smith-Waterman-Gotoh extension.

## Getting started

C99 compiler (gcc / clang / icc) is required to build the program.

```
$ make && make install	# PREFIX=/usr/local by default
$ minialign <reference.fa> <reads.fq> > out.sam		# read-to-ref alignment
```

Reference sequence index can be stored in separate. Using prebuilt index saves around a minute per run for human (~3G) genome.

```
$ minialign -d index.mai <reference.fa>	# build index
$ minialign -l index.mai <reads.fq> > out.sam	# mapping on prebuilt index
```

Frequently used options are: scoring parameters, minimum length/score cut-offs, and number of threads.

```
$ minialign -a1 -b2 -p2 -q1		# match, mismatch, gap-open and gap-extend
$ minialign -m400	# set minimum score threshold to 400
$ minialign -r0.8	# set report threshold at 0.8 x highest score for every read
$ minialign -t10	# minialign is now 10x faster!!!
```

## Benchmarks

All the benchmarks were took on Intel i5-6260U (Skylake, 2C4T, 2.8GHz, 4MBL3) with 32GB (DDR4, 2133MHz) RAM.

### Speed

|                      Time (sec.)                     |  minialign  |   DALIGNER  |   BWA-MEM   |
|:----------------------------------------------------:|:-----------:|:-----------:|:-----------:|
| E.coli (MG1655) x100 simulated read (460Mb) to ref.  |        16.7 |        39.5 |        6272 |
| S.cerevisiae (sacCer3) x100 sim. (1.2Gb) to ref.     |        43.0 |           - |       10869 |
| D.melanogaster (dm6) x20 sim. (2.75Gb) to ref.       |         139 |           - |       31924 |
| Human (hg38) x3 sim. (9.2Gb) to ref.                 |        1571 |           - |           - |

Notes: PBSIM (PacBio long-read simulator), [modified version based on 1.0.3](https://github.com/ocxtal/pbsim/tree/nfree) not to generate reads containing N's, was used to generate read sets. Parameters (len-mean, len-SD, acc-mean, acc-SD) were fixed at (20k, 2k, 0.88, 0.07) in all the samples. Minialign and DALIGNER were run with default parameters except for the multithreading options, `-t4` and `-T4` respectively. BWA-MEM was run with `-t4 -A1 -B2 -O2 -E1 -L0`, where scoring (mismatch and gap-open) parameters adjusted based on the presets of `-xpacbio`. Index construction (minialign and BWA-MEM) and format conversion time (DALIGNER: fasta -> DB, las -> sam) are excluded from measurements. Peak RAM usage was around 12GB in human read-to-ref mapping with four threads.

### Read-lendth vs. sensitivity trend

![length-sensitivity plot](https://github.com/ocxtal/minialig/blob/master/pic/len_sens.png)

Notes: Sensitivity is defined as: the number of reads whose originating location is correctly identified (including secondary mappings) / the total number of reads. Reads are generated from hg38 without ALT / random contigs using PBSIM with the same parameters as the speed benchmark. Reads were mapped onto the reference with ALT / random contigs included. Minialign was run with the same parameters as in the speed benchmark except for the minimum mapped region length `-m`, set to the half of the mean read length.

### Speed vs. sensitivity trend

![speed-sensitivity plot](https://github.com/ocxtal/minialig/blob/master/pic/spd_sens.png)

Notes: k-mer length and window size parameters were altered in the benchmark. Details of the correspondences between point and parameters were found in scripts directory.

## Algorithm overview

### Minimizer-based index structure

Indexing routines: minimizer calculation, hash table construction, and hash table retrieval are roughly diverted from the original minimap program except that the position of sequence direction flag in hash table is moved from the least significant bit to the int32_t sign bit. See descriptions in the minimap [repository](https://github.com/lh3/minimap) and [paper]() for the details of the invertible hash function used to generate minimizers.

### Seed chaining

Collected seeds (minimizers) were first sorted by its (rpos - 2*qpos) value, resulting in lining up along with 15-degree leaned lines from the diagonal. Then (rpos - qpos/2) values are evaluated from head to tail on the sorted elements, and chained when the current seed is inside a 30-degree-angled parallelogram window at the right-bottom direction of the previous one. Chaining is iteratively performed on the remaining seed array and terminated when no seed is left in it. Finally, collected chains are sorted by their rough path lengths: (re - rs + qe - qs).

### Smith-Waterman-Gotoh extension

The second head seed of each chain is extended upward (3' on the reference side) then downward from the maximum score position found. If the resulted path is shorter than the seed chain span, similar extension is repeatedly performed on the next one (three times at the maximum). Each extension is carried out by the [GABA library](https://github.com/ocxtal/libgaba), which implements the adaptive-banded semi-global Smith-Waterman-Gotoh algorithm with difference recurrences.

## Notes, issues and limitations

* SDUST masking is removed from the original minimap implementation.
* Repetitive seed-hit region detection is also removed.
* Large gap open penalty (> 5) and large X-drop penalty (> 64) are disallowed due to the limitation of the GABA library.
* Index file format is incompatible with of the minimap.

## Updates

* 2016/11/2 (0.3.0) First release of 0.3 series, with a better chaining algorithm.
* 2016/10/5 (0.2.1) Last tagged commit of the version 0.2 series.
* 2016/9/13 (0.2.0) First tagged commit (unstable).

## Gallery

#### *Fast and Accurate* logo

![metcha hayaiyo](https://github.com/ocxtal/minialig/blob/master/pic/hayai.pdf)

#### Intel nuc, my main development machine

![he is also powerful](https://github.com/ocxtal/minialig/blob/master/pic/nuc.jpg)

## Copyright and license

The original source codes of the minimap program were developed by Heng Li and licensed under MIT, modified by Hajime Suzuki. The other codes, libgaba and ptask, were added by Hajime Suzuki. The whole repository except for the pictures in the gallery section (contents of pic directory) is licensed under MIT, following that of the original repository.
