# xtools

[Xuebing](http://www.mit.edu/~wuxbl)'s tools for computational biology (C/C++)


## Prerequisite Software/Package

1. C++ library Boost, only tested [Boost version 1.57.0](http://sourceforge.net/projects/boost/files/boost/1.57.0/)
2. [R](http://www.r-project.org/) is required for making figures and calculating FDR from raw p values.

## Installation
1. Download the source code from xtools page on github: https://github.com/xuebingwu/xtools
2. Unzip and then enter the subfolder 'include':

```sh
unzip xtools-master.zip
``` 
3. Download the C++ library boost (tested version 1.57.0) from http://www.boost.org/ or http://sourceforge.net/projects/boost/files/boost/1.57.0/
2. Unzip and move the subfolder 'boost_1_57_0/boost' to the folder 'xtools-master/include'.
3. Go to folder 'xtools-master/src', type 'make' to compile the source code:
```sh
cd xtools-master/src
make 
``` 
4. The compiled programs are in the folder 'bin'.
5. To recompile, type 'make clean' then 'make'.


## Project home page

Source code, compiled binaries (Linux), documentation, and examples can be found on github page:

https://github.com/xuebingwu/xtools


## Cite

This software is currently unpublished, so please just cite the homepage (thanks!).

## Copyright

Copyright (c) 2014 Xuebing Wu. See LICENSE.txt for further details.
