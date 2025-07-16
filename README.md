# rmatrix linear algebra package
A _very_ basic linear algebra package meant for use with the ROOT library. I wrote this because I disliked ROOT's lin algebra package. 
Right now, it supports addition of matrices, multplication with vectors, element-wise access, and, most importantly, solving of linear systems. 

## Compile for use with ROOT
I am using this package for development with other projects on WSL, so I set up the CMakeLists.txt file to compile these files as ROOT libraries. To use them, run the following commands (on a bash terminal).

Fist, create an empty directory, and clone the repo: 

    cd /path/to/local/clone
    git clone https://github.com/sethcarl2000/rmatrix

Then, compile the files: 

    mkdir build
    cmake -B build -S .
    cmake --build build
    cmake --install build

Then, to make sure you can use the RMatrix class in any root macro, just invoke: 

    setenv LD_LIBRARY_PATH="/path/to/local/clone/build:${LD_LIBRARY_PATH}"

Or, if you're in a Conda environment (as I am), activate the environment (which has ROOT installed as a package!) and run: 

    conda env config vars set LD_LIBRARY_PATH="/path/to/local/clone/build:${LD_LIBRARY_PATH}"

And it will ask to you deactivate/reactivate your env to have the change take effect. after this, you should be able to create and use an RMatrix object just like any other ROOT class. 
