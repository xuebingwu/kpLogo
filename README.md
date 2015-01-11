# PKA: positional kmer analysis

PKA is tool for identifying well-positioned short motifs from a set of aligned sequences. 

Developed by Xuebing Wu (wuxbl@wi.mit.edu) at Whitehead Institute.

## Installation
1. Download the source code from github: https://github.com/xuebingwu/PKA
2. Unzip, enter the subfolder 'src', then type 'make':

```sh
cd src
make
``` 
3. The compiled programs are in the folder 'bin'.
4. To recompile, type 'make clean' then 'make'.

## Test examples

Go to the folder 'bin', and run the following examples:

1. unweighted sequences in fasta, default background (relative to other positions)

```sh
./PKA ../test/input/input-fixed-length.fa -o ../test/output/unweighted-default
```

2. weighted sequences 

```sh
 ./PKA ../test/input/input.weighted.txt -weighted 
```

Please refer to the manual for other examples and explanation of the output files.


## Project home page

Source code, compiled binaries (Linux), documentation, and examples can be found on github page:

https://github.com/xuebingwu/xtools


## Cite

This software is currently unpublished, so please just cite the homepage (thanks!).

## Copyright

Copyright (c) 2014 Xuebing Wu. See LICENSE.txt for further details.
