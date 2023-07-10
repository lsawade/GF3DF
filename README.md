# GF3DF

Fortran extraction tool for [Specfem](https://github.com/SPECFEM/specfem3d_globe)
and [GF3D](https://github.com/lsawade/GF3D) generated Green function database
subset files. The main intent of this repo to include 3D Green functions in the
[GCMT project](https://www.globalcmt.org)'s CMT inversion routine.


## Quickstart

```bash
git clone
cd /path/to/gf3df/
cmake -S . -B build
cmake --build build
```

## Dependencies

- **Fortran compiler** -- Tested only gfortran so far; most recently GNU Fortran
  13.1.0. Have not tried any other setups.
- **CMake** -- currently I'm running with Cmake 3.26.4, and I haven't really tested
  any other cmake setups. The package is not really using any "crazy" CMake
  functionalities, so I'm not too worried about the whole thing.
- **HDF5** -- I have tested with a few different HDF5 packages, but since the
  function that I'm using are pretty basic reading functions it has been working
  with all of them (1.10, 1.12, 1.14).

Random tests that worked as well:
- Using a parallel hdf5 implementation


## Shared library

If you want to build a shared library, please use the `-fPIC` flag
```bash
export CFLAGS="-fPIC"
export FFLAGS="-fPIC"
```

## HDF5

Whether you use your own built parallel HDF5 or a cluster non-mpi HDF5, both
should work fine! Just make sure that the compilers are the same so that there
are no issues.

Make sure the `HDF5_ROOT` variable is set.
```bash
export HDF5_ROOT=/path/to/hdf5
```

`cmake` will automatically find the include dirs and libraries related to the
`HDF5` installation.



For using the modules in your own package add `/path/to/gf3df/build/include` to
your compilation and in fortran you should be able to just `use gf3d`.


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



---

## Notes on code structure

Organization of the code could be better. I'm not quite sure whether the going
down the `submodule` route is the right thing, I did however start changing the
library into something that is a set of `submodules` with interfaces -- if
needed -- in the main module. I'll try to rearrange the software a little more.


## Known Issues

- **CMake does sometimes not find HDF5 correctly**. On my Mac, I have
  had the issue that `find_package` in `CMake` failed or messed up the
  `HDF5_INCLUDE_DIRS` path. Then it is really the easiest to some variables via
  command line. First the normal paths
  ```bash
  export HDF5_ROOT=<path/to/hdf5>
  export HDF5_LIBRARIES="${HDF5_ROOT}/lib"
  export HDF5_INCLUDE_DIRS="${HDF5_ROOT}/include"
  export PATH="${HDF5_ROOT}/bin:${PATH}"
  ```
  and the an easy way to collect all the flags is to simply use the HDF5 fortran
  wrapper as a reference:
  ```bash
  export FFLAGS=${FFLAGS} $(echo "${FFLAGS}" $(h5fc -show | cut -w -f '2-60'))
  ```
  Another way would be to just set `export FC=h5fc`.

- Only really tested with `gfortran`. There maybe portability issues with `INQUIRE`
routine.