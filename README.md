## triangulator

A simple 2D Delaunay triangulation library written in C++.

It uses the [Bowyer-Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm) as well as [Ruppert's algorithm](https://en.wikipedia.org/wiki/Ruppert%27s_algorithm) to create quality meshes from 2D polygons.

The implementation is partly based on [PyDelaunay2D](https://github.com/jmespadero/pyDelaunay2D) by Jose M. Espadre, a didactic implementation of the Bowyer-Watson-algorithm written in Python.

The library uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), a free C++ template library for linear algebra.

### Building

The library can be added to your project as a submodule or by placing the repos contents in your project tree. It includes a `CMakeLists.txt` that defines the target `triangulator` so you can include it in your CMake build by using `add_subdirectory`.

### Usage

An example application using `triangulator` can be found [here](https://github.com/cp3-ws1920/applications).

TOOD: Add some documentation.
