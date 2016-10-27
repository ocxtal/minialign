
# minialign

Minialign is a fast and accurate nucleotide sequence alignment tool designed for PacBio and Nanopore long reads. It is built on three key algorithms, minimizer-based index of the [minimap](https://github.com/lh3/minimap) overlapper, array-based seed chaining, and SIMD-parallel Smith-Waterman-Gotoh extension.

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

All the benchmarks were took on Intel i5-6260U (Skylake, 2C4T, 2.8GHz, 4MBL3) with 32GB (DDR4, 2133MHz) memory.

### Speed

|                 Time (sec.)                  |  minialign  |   DALIGNER  |   BWA-MEM   |
|:--------------------------------------------:|:-----------:|:-----------:|:-----------:|
| E.coli (MG1655) x100 simulated read to ref.  |       18.3s |       39.5s |       6272s |
| E.coli (MG1655) x2000 simulated read to ref. |        477s |           - |         10h |
| C.serevisiae x100 simulated read to ref.     |       84.2s |           - |      10869s |
| D.melanogaster x20 simulated read to ref.    |        654s |           - |      31924s |
| Human (hg39) x20 simulated read to ref.      |           - |           - |           - |

Notes: PBSIM (PacBio long-read simulator) was used to generate read sets. Parameter sets (len-mean, len-SD, acc-mean, acc-SD) were set to (20k, 2k, 0.88, 0.07) in both samples. Minialign was run with default parameters except `-t4`, and BWA-MEM was run with `-t4 -A1 -B2 -O2 -E1 -L0`. Index construction time (minialign and BWA-MEM) and format conversion time (DALIGNER: fasta -> DB, las -> sam) are excluded from results.

### Read-lendth vs. sensitivity trend


Notes: Sensitivity is defined as: number of reads whose original location is correctly detected / total reads. Reads are generated from hg39 using PBSIM with the same parameters as the speed benchmarks. ALT/random contigs were excluded from reference sequences in read generation and included in mapping.

### Speed vs. sensitivity trend


Notes:

## Algorithm overview

### Minimizer-based index structure

Indexing routines: minimizer calculation, hash table construction, and hash table retrieval are roughly diverted from the original minimap program. The position of direction flag in hash table is moved from LSb in the original to the sign bit of int32_t in this code. See descriptions in the minimap [repository](https://github.com/lh3/minimap) and [paper]() for the details of the invertible hash function used to generate minimizers.

### Seed chaining

Collected seeds (minimizers) were first sorted by its (rpos - 2*qpos) value, resulting in lining up along with 15-degree leaned lines from the diagonal. Then (rpos - qpos/2) values are evaluated from head to tail, and chained if the current seed is inside the 30-degreed parallelogram window at the right-bottom direction of the previous one. Chaining is iteratively performed on remaining seed array and terminated when no seed is left in it. Finally, collected chains are sorted by their rough path lengths: (re - rs + qe - qs).

### Smith-Waterman-Gotoh extension

The head seed of each chain is extended upward on the reference side then downward from the maximum score position found. If the resulted path is shorter than the seed span, similar extension is repeatedly performed on the next one (three times at the maximum). Each extension is carried out by the GABA library, which implements adaptive-banded semi-global Smith-Waterman-Gotoh algorithm with difference recurrence.

## Notes, issues and limitations

* SDUST masking is removed from the original minimap implementation.
* Repetitive seed hit detection is also removed.
* Large gap open penalty (> 5) is disallowed due to the limitation of the GABA library.

## Gallery

* Intel nuc

* Fast and accurate

## License

MIT (following the original license)
