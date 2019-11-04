# Monte Carlo methods for radiative transfer

This project is an exploration of Monte Carlo (MC) methods for radiative
transfer. It is developed as part of a graduate course on radiative transfer.
A report on the development as well as the solutions to the simulation exercises
can be found [here](https://simonpf.github.io/llrte).

## Dependencies

- CMake >= 3
- GCC >= 7.4.0
- To generate plots: Python 3, matplotlib, numpy

## Installation

To build the executable for the exercises:

````
$ git clone https://github.com/simonpf/llrte
$ cd llrte
$ mkdir build
$ cd build
$ cmake ..
$ make
````

To run the simulations and generate plots:

````
$ make exercise_1
$ make exercise_2
````
