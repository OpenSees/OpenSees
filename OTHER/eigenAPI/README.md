This folder contains header-only sources meant to interface OpenSees with eigen in an effective way. 
We define macros for vectors, matrices and tensors, as well as implementations of commont operations
that are useful for implementing new elements and materials. 

It uses a git submodule to get the source for eigen 

After a fresh git clone of the OpenSees repo, one must execute

    git submodule init
    git submodule update 

To get the  eigen sources in the right folder