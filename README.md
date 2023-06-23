# Fortran bindings to read GF3D produced Subset files

Only really tested with `gfortran`. There maybe portability issues with `INQUIRE`
routine.

Organization of the code could be better. I'm not quite sure whether the going
down the `submodule` route is the right thing. I have to write an extra
interface for every function in the main module. Either way, that can wait.
Right now, I'm just keeping the convention with module names so that
`src/<dir>/<dir>.F90` has submodules `src/<dir>/<dir>_<sub>.F90`. I don't think
this is very clean, but an ok way to get around the `*.mod` having unique names
convention.

If you want to build a shared library, please use the `-fPIC` flag
```bash
export CFLAGS="-fPIC"
export FFLAGS="-fPIC"
```

Whether you use your own built parallel HDF5 or a cluster non-mpi HDF5, both
should work fine! Just make sure that the compilers are the same so that there
are no issues.

Make sure the `HDF5_ROOT` variable is set.
```bash
export HDF5_ROOT=/path/to/hdf5
```

`cmake` will automatically find the include dirs and libraries related to the
`HDF5` installation.

Then, you can build the package using

```bash
cd /path/to/gf3df/
cmake -S . -B build
cmake --build build
```

For using the modules in your own package add `/path/to/gf3df/build/include` to
your compilation and in fortran you should be able to just `use gf3d`.


## Building with HDF5

```bash
export HDF5_ROOT="/home/lsawade/ph5py-testing/hdf5/build/phdf5"

# For now add some debugging flags
export FFLAGS="-Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fba
cktrace"

rm -rf build && \
    cmake -S . -B build && \
    cmake --build build && \
    ./build/bin/gf3d single_element.h5
```

Just building and running

```bash
cmake --build build && ./bin/gf3d single_element.h5
```

