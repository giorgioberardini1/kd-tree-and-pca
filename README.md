# Kd-tree and PCA 

A Kd-tree is a k-dimensional binary search tree. Given a dataset D, that is a set of a k-dimensional points, build a kd-tree means to index a dataset allowing to answer fast to the following question: 

  *given a point Q and a radius r, what are the points of D which are less than r from Q?*


It has been proved that when k is big enough an indexing technique is not very efficient.

A solution is to reduce the dataset dimension considering only the principal component, this is also known as *Principal Component Analysis* (PCA).   


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


### Generate the executable
  
 
 Use the following command for compiling: 

```
nasm -f elf32 kdtreepca32.nasm && gcc -O0 -m32 -msse kdtreepca32.o kdtreepca32c.c -o kdtreepca32c && ./kdtreepca32c
```

  
or, if you want to compile and run just type: 

  

```
./runkdtreepca32
```

for the 64 bit version use: 

```
./runkdtreepca64
```

## Running the tests

  

If you want to run the program use:

    ./runkdtreepca32 D [-pca <h>] [-kdtree [-rq <r>]]


 - D: dataset name (without '.ds')
 - pca <<h>h>: number of pca component 
 - kdtree: index the tree
 - rq <<r>r>: range query with radius r 
 
  
Note 1: if you are running in pca mode you're going to decompose the original dataset in two small datasets. They will be stored in two files named 'datasetname.U' and 'datasetname.V' each one bytes typed row major. 

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

 We got pretty good results in term of performance, we've got 1st place (team #17) over 19 groups for the fastest implementation. 

![Ciao](https://imagizer.imageshack.com/img923/8595/wjp7Kt.png)

  


