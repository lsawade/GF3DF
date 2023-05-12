# Fortran bindings to read GF3D produced Subset files

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
rm -rf build && \
    cmake -S . -B build && \
    cmake --build build && \
    ./build/bin/gf3d single_element.h5
```

Just building and running

```bash
cmake --build build && ./bin/gf3d single_element.h5
```

