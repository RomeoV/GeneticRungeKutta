# GeneticRungeKutta

This is just a toy project of mine, where I implemented a [genetic algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm)
to optimze a [Runge Kutta Tableau](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods),
which represents a (or rather **the**) method to solve ordinary differential equations numerically.
Currently, only explicit schemes are supported and it's probably going to stay that way. Maybe one day I'll write an implicit solver.

This tool also makes use of C++-17 parallel algorithms using Intel's TBB library, so in theory this program is highly paralellized without the use of OpenMP or something similar.
In practice I haven't seen great performance differences, but that might be due to my machine.

To compile this, you need a pretty recent compiler, [Intel's TBB library](https://github.com/intel/tbb) and C++17 support.
