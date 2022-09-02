# sk-test

An attempt in translating the Fortran code in dftb to C++ for reading the inputs from slater koster files, e.g. B-B.skf

helloworld.cpp
- c++ Module feature testing

memory_tracker.h
- C++ code for tracking memory allocation

SlakoCont_init in 
- Tried usinh this-> pointer for constucting init
- Later followed Modern c++ guide for constructor
- Still has operator error due to TSlater creating an instance

Trying to convert from dftbplus:
- sk.F90
- orbital.F90
- pairrepulsive.F90
- parser.F90(readSKFiles)
- slakocont.F90
- linkedlist.F90
- linkedlisti1.F90

Trying to follow styles from
- [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html#thread_local)
- [Header files (C++) | Microsoft Docs](https://docs.microsoft.com/en-us/cpp/cpp/header-files-cpp?view=msvc-170#include-guards)
