# minialign

Minialign is a fork of the [minimap](https://github.com/lh3/minimap) long-read overlapper, with Smith-Waterman alignment calculation and sam output.

## Algorithm

See `algorithm overview' of the minimap repository for the detail of the seed gathering and chaining. See [libgaba](https://github.com/ocxtal/libgaba) for the detail of the adaptive-banded alignment algorithm. Note that the overall sensitivity is slightly low due to the minimizer-based indexing.

## Usage

C99 compiler (gcc / clang) is required.

```
$ make
$ ./minialign <reference.fa> <reads.fq> > out.sam
```

## License

MIT (following the original license)
