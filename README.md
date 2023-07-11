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

## Usage in your package

An example of an external program that uses GF3D is given in
`./external_program`.

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

## Optimization flags

To get a speed increase of 3x, please use the following environment variable
before making.

```bash
export FFLAGS="-O3"
```

executables are larger do to unrolled loops and other global optimizations.
I tried using both `-O4` and `-O5`, as well. There is no advantage in `-O4` and
the aggressive option `-O5` is running even slower. Since extraction time for
17 stations is at 0.65 seconds for a 4h seismogram. I will leave it at that for
now. I could to some advanced profiling later on.

## Usage in your own module

For using the modules in your own package add `/path/to/gf3df/build/include` to
your compilation and in fortran you should be able to just `use gf3d`. An
example use is shown in `./app/gf3d-get-sdp-demo`


## Shared library

If you want to build a shared library, please use the `-fPIC` flag
```bash
export CFLAGS="-fPIC"
export FFLAGS="-fPIC"
```

## HDF5

Whether you use your own built parallel HDF5 or a cluster non-mpi HDF5, both
should work fine! Just make sure that the compilers are the same so that there
are no issues. In general, `cmake` will automatically find the include dirs and libraries related to the `HDF5` installation. However, sometimes it's not
quite the best at setting the link paths. If that's the case simply set.

```bash
export HDF5_ROOT=/path/to/hdf5
```

If that does not work, see below in the `Known Issues section`.


## Debugging setep

```bash
export HDF5_ROOT="/home/lsawade/ph5py-testing/hdf5/build/phdf5"

# For now add some debugging flags
export FFLAGS="-Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fba
cktrace"

# Remove build completely and rebuild.
rm -rf build && \
    cmake -S . -B build && \
    cmake --build build && \
    ./build/bin/gf3d-get-sdp-demo single_element.h5 3

```

---

## Notes on code structure

Organization of the code could be better. I'm not quite sure whether the going
down the `submodule` route is the right thing, I did however start changing the
library into something that is a set of `submodules` with interfaces -- if
needed -- in the main module. I'll try to rearrange the software a little more.


## Known Issues

- **CMake does sometimes not find HDF5 correctly**. I had a lot of trouble
  getting it right on my Mac. And eventually, I created a fix in the CMake file.
  Should you have trouble. Instead of updating the `./CMakeLists.txt`, I would
  start with setting all the `HDF5` paths like so

  ```bash
  export HDF5_ROOT=<path/to/hdf5>
  export HDF5_LIBRARIES="${HDF5_ROOT}/lib"
  export HDF5_INCLUDE_DIRS="${HDF5_ROOT}/include"
  export PATH="${HDF5_ROOT}/bin:${PATH}"
  ```
  Then try to compiles it again. If that should also fail, you could manually
  link HDF5 using the HDF5 fortran compiler wrapper as a reference:
  ```bash
  export FFLAGS=${FFLAGS} $(echo "${FFLAGS}" $(h5fc -show | cut -w -f '2-60'))
  ```
  If you end up doing this, make sure to remove all lines with
  `find_package(HDF5 ...)`

- Only really tested with `gfortran`. There maybe portability issues with `INQUIRE`
routine.