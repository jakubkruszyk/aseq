# Aseq
Library for arithmetic sequences detection in set of numbers (vectors) and selection of those sequences that will cover all given vectors.

## Usage of binary
```shell
aseq <path_to_vectors> <filter> [options]
```

where possible filters are:
- `greedy` - filter will used sequnces that covers the most of remaining vectors
- `hybrid <cutoff>` - filter works like greedy version, but during each iteration sequences with ratio of covered vectors to sequence length below cutoff value are discarded. That implies that cutoff value should be from range <0:1>

## Program outputs:
All sequences `aseq` prints contains indexes of vectors in input file. For example, consider file test.txt:
```txt
111
222
333
444
```
Sequence [0, 1, 3] translates to [111, 222, 444]. Additionaly software prints "de-mapped" version where indexes are replaced with actual vectors.
