# kpLogo: k-mer probability logo

kpLogo is tool for identifying well-positioned short motifs from a set of aligned sequences.

Developed by Xuebing Wu (wuxbl@wi.mit.edu) at Whitehead Institute.

## Web server

kpLogo is also available as an online tool: http://kpLogo.wi.mit.edu

## Installation
Download the source code from github: https://github.com/xuebingwu/kpLogo/releases/

Unzip, enter the subfolder 'src', then type 'make':

```sh
cd src
make
``` 

The compiled programs are in the folder 'bin'.

To recompile, type 'make clean' then 'make'.

## Project home page

Source code, compiled binaries (Linux), documentation, and examples can be found on github page:

https://github.com/xuebingwu/kpLogo


## Cite

Xuebing Wu and David Bartel (2017) kpLogo: positional k-mer analysis reveals hidden specificity in biological sequences, submitted.

## Copyright

Copyright (c) 2014 Xuebing Wu. See LICENSE.txt for further details.

PKA reuses the following library/code: 

1. Boost library (version 1.57.0): http://www.boost.org/  

2. The ushuffle code written by Minghui Jiang: http://digital.cs.usu.edu/~mjiang/ushuffle/

