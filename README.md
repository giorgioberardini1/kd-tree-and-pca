# Kd-tree - PCA

  

A Kd-tree is a k-dimensional binary search tree. Given a dataset D, that is a set of a k-dimensional points, a kd-tree can index these points allowing to answer fast to the following question:



  

*given a point Q and a radius r, which are the points of D that are less than r from Q?*

  
  

It has been proved that when k is big enough an indexing technique is not very efficient.

  

A solution is to reduce the dataset dimension considering only the principal component, this is also known as *Principal Component Analysis* (PCA).


The algorithm used to implement the PCA is the NIPALS algorithm. 

In order to obtain an improvement in performance in terms of time, we used the following optimization techniques both in C  and  Assembly environment ([ILP](https://en.wikipedia.org/wiki/Instruction-level_parallelism) - [SIMD](https://en.wikipedia.org/wiki/SIMD)):


- Caching
- Loop interchange
- Code vectorization
- Loop unrolling


 

  

## Getting Started

  

  

### Prerequisites

  

Software neeeded:

  

- NASM (www.nasm.us)

- GCC (gcc.gnu.org)

  

for example on Ubuntu you can type:

```

sudo apt-get install nasm

  

sudo apt-get install gcc

```

  

It might be necessary to run:

  

sudo apt-get install lib32gcc-4.8-dev

sudo apt-get install libc6-dev-i386

  
  
  

## Running the tests

  

  

If you want to compile and run use:

  

    ./runkdtreepca32 D [-pca <h>] [-kdtree [-rq <r>]]

  

for 64 bit version use:

  

    ./runkdtreepca64 D [-pca <h>] [-kdtree [-rq <r>]]

  
  

- D: dataset name (without '.ds')

- pca <<h>h>: number of pca component

- kdtree: index the tree

- rq <<r>r>: range query with radius r

Note 1: if you are running in pca mode you're going to decompose the original dataset into two small datasets. They will be stored in two files named 'datasetname.U' and 'datasetname.V' both bytes typed row major.

  
  

Note2: if you are running in range query mode you're going to generate a list of query answers that are stored in one file named 'datasetname.qa' bytes typed row major.

  
  

### About the dataset

  

The dataset need to be bytes typed row major in this format:

```

  

numberOfColumns numberOfRows Point1FirstComponent Point1SecondComponent and so on..

  

```

  

you can find 2 datasets in the GitHub folder of the following dimension:

  

- 8000x128 (test1.ds)

- 5000x36 (test2.ds)

  
  
  
  

## Performance and Results

  

We got pretty good results in term of performance, we ranked 1st over 19 groups for the fastest implementation.

  
