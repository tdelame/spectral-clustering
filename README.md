# Spectral Clustering

## What is it
I needed to experiment some spectral clustering algorithms for my research 
topics. I deciced to let the code available as it may help other people to
understand spectral clustering to implement their own algorithms. 

This code implements different variations of the spectral clustering method. It
is based on the 2007 tutorial of Ulrike von Luxburg (to see a pdf of this 
tutorial, follow this link http://goo.gl/W2sRH1). I haven't got the time to 
profile, optimize and document the code yet. If you find any mistake, have some
questions or have ideas from improvements (there might be a lot to do:-}), feel
free to contact me. 

Also, if you have some knowledge of the underlying data, you can improve some 
key steps of the implementation. For example, if you know that no cluster should
have parts of two different connected component of the similarity graph, you 
should set the number of clusters to be at least the multiplicity of the eigen 
value 0 (this multiplicity is the number of connected components). Another 
example: if some cluster is way too small for your applications, you should 
consider to merge this cluster with adjacent one and use another eigen value to
generate another cluster. 

## Configuration
The code has only been tested on GNU/Linux platforms, compiled with gcc 5.2. to
5.3. Three libaries are used:
- eigen (latest version)
- nanoflann (one header included in the sources)
- spectral (a header only library included in the sources).

It is worth noting that the nanoflann library, while being advertised as a 
header only simple version of the flann library, also solved many issues I had 
with flann (radius search missed many results, for no known reasons). This is 
why I recommend the use of nanoflann everywhere the flann library is used.

## Installation
First check out this repository
```
git clone https://github.com/tdelame/spectral-clustering.git
```

Then edit the file _Makefile_ to define some variables for the compilation. 
Select a build type by typing one of the two following lines: 
```
make release
make debug
```

Once the build type is active, you can build the command line application by
typing ``make``.
    
## Usage
Before running the binary, you need to set some environnement variables. They
are written in the file 'SourceMe.sh'. Also, this script ensure the creation of
output directory. Once you set those environnement variables to the correct 
values for your system, simply type:
```
source SourceMe.sh
```

## Licensing
This code is release under the MIT Licence. Please see the file called LICENCE.

## Contacts
Thomas Delame, tdelame@gmail.com
