# rmatrix linear algebra package
A _very_ basic linear algebra package meant for use with the ROOT library. I wrote this because I disliked ROOT's lin algebra package. 
Right now, it supports addition of matrices, multplication with vectors, element-wise access, and, most importantly, solving of linear systems. 

## Basic usage
Suppose you have a 3x3 system of equations, something like  Ax = b. To use this package to find x = (A^-1)b, use: 

    //create the matrix (3 rows, 3 cols)
    RMatrix A( 3, 3,  { 1.0,  0.2, -0.3, 
                        0.2, -2.1,  3.0, 
                        1.1,  0.0,  1.5 }); 

    //Create the b vector: 
    vector<double> b = {-9.1, -2.0, -0.3 }; 

    //Solve the system: 
    auto x = A.Solve(b); 

    printf("Answer (x): {%+.4f,%+.4f,%.4f}\n", x[0], x[1], x[2]); 

    //check our answer: 
    auto b_test = A * x; 

    for (int i=0; i<b.size(); i++) printf("Error: %+.10f\n", b_test.at(i) - b.at(i)); 

This should return: 

    Answer (x): {-9.0268,+9.2636,6.4197}
    Error: +0.0000000000
    Error: +0.0000000000
    Error: -0.0000000000

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
