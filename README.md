# Project 1 - Sequence Alignment
## Author: Parker Grosjean
## BMI 203 UCSF

![BuildStatus](https://github.com/pgrosjean/Project1/workflows/HW1/badge.svg?event=push)

## This Repository

This repo holds all of the information needed to run the Smith-Waterman and Needleman-Wunsch alignment for local and global alignment respectively.
This repository is set up with three main sections locations for the (1) Docs (2) Align Module and (3) the pytest functionalities.

### Docs
To view the API docs for the align module and specifically the classes contained within align.algs please see the [API doc](https://github.com/pgrosjean/Project1/blob/main/docs/build/html/align.html) by running "open index.rst" in the location of the file pointed to by the API doc link.

### Align Module
The [Align Module](https://github.com/pgrosjean/Project1/tree/main/align) holds the algs submodule that contains both the base class PairwiseAligner and the inhereted classes for SmithWaterman and NeedlemanWunsch. See the API docs for more in depth descriptions of how to use this module.

### Pytest Functionalities
All of the necessary fasta files for testing are in [this folder](https://github.com/pgrosjean/Project1/tree/main/files_test) and the actual unit tests are implemented [here](https://github.com/pgrosjean/Project1/blob/main/test/test_align.py).

### Testing
Testing is as simple as running
```
python -m pytest
```
from the root directory of this project.
This testing corresponds to the badge at the top of this readme.
