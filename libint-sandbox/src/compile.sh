#!/usr/bin/bash

# NOTE: libint is was compiled following the documentation. To summarize, I ran the following steps:
#       1. Install eigen
#           a. download the source code
#           b. mkdir build; cmake ..; cmake --build . --target install
#       2. Install boost
#           a. On Ubuntu: sudo apt install build-essential libboost-system-dev libboost-thread-dev libboost-program-options-dev libboost-test-dev libboost-filesystem-dev
#       3. Compile/Install libint:
#           a. download latest release
#           b. mkdir build; cmake ..; 
#           c. cmake --build . --target check
#           d. cmake --build . --target install
#               - This will install the headers and library in /usr/local/include and /usr/local/lib, respectively

RUN_LEVEL="-O0"

if [[ $# -gt 1 ]]; then
    echo "Invalid number of arguments. You must pass only a single argument and it must be either 'debug' or 'release'"
    exit 1
elif [[ $# -eq 0 ]]; then
    echo "WARNING: No run level specified. Running in DEBUG by default"
elif [[ $# -eq 1 ]]; then
    if [ "$1" = "debug" ]; then
        RUN_LEVEL="-O0"
    elif [ "$1" = "release" ]; then
        RUN_LEVEL="-O3"
    else
        echo "Invalid argument: '$1'. You must pass only a single argument and it must be either 'debug' or 'release'"
        exit 1
    fi
fi
    
set -x

g++ -std=c++23                      \
    -I/usr/local/include            \
    -I/usr/local/include/eigen3     \
    -L/usr/local/lib                \
    ${RUN_LEVEL}                    \
    -Wall -Wextra -Wpedantic        \
    main.cpp                        \
    -l:libint2.a

