# PKA: positional kmer analysis

PKA is tool for identifying well-positioned short motifs from a set of aligned sequences.

Developed by Xuebing Wu (wuxbl@wi.mit.edu) at Whitehead Institute.

## Web server

PKA is also available as an online tool (under-construction): http://iona.wi.mit.edu/xuebing/PKA/

## Installation
Download the source code from github: https://github.com/xuebingwu/PKA/releases/

Unzip, enter the subfolder 'src', then type 'make':

```sh
cd src
make
``` 

The compiled programs are in the folder 'bin'.

To recompile, type 'make clean' then 'make'.

## Test examples

Go to the folder 'bin', and run the following examples:

Unweighted sequences in fasta, default background (relative to other positions)

```sh
./PKA ../test/input/input-fixed-length.fa -o ../test/output/unweighted-default
```

Weighted sequences 

```sh
 ./PKA ../test/input/input.weighted.txt -weighted 
```

Please refer to the manual for other examples and explanation of the output files.


## Project home page

Source code, compiled binaries (Linux), documentation, and examples can be found on github page:

https://github.com/xuebingwu/PKA


## Cite

This software is currently unpublished, so please just cite the homepage (thanks!).

## Copyright

Copyright (c) 2014 Xuebing Wu. See LICENSE.txt for further details.

PKA reuses the following library/code: 

1. Boost library (version 1.57.0): http://www.boost.org/  

2. The ushuffle code written by Minghui Jiang: http://digital.cs.usu.edu/~mjiang/ushuffle/


