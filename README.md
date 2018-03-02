Hashing Source Code
===================

Source code used in the hashing experiments of the thesis ``On the Analysis of
Two Randomized Algorithms: Multi-Pivot Quicksort and Efficient Hash Functions''.

Dependencies
------------

Needs g++, cmake in version >= 2.8, libboost-random and libpapi for performance
measurements. (Enabled by default, can be changed in CMakeLists.txt.)

How to build
------------
Use the following commands on the top-level directory of the project.

>> mkdir build; cd build
>> cmake ..
>> make

For a release build with optimization flags, use 
>> cmake -DCMAKE_BUILD_TYPE=Release

After successful compilation, the executable is located at
build/src/hashingtest.

Example Calls
-------------

Example calls can be found in the directory ``examples''.

Written November 28, 2014 by Martin Aumüller.
