

Make sure the `HDF5_ROOT` variable is set.
```bash
export HDF5_ROOT=/path/to/hdf5
```

`cmake` will automatically find the include dirs and libraries related to the
`HDF5` installation.

Then, you can build the package using

```bash
cmake -S . -B build
cmake --build build
```