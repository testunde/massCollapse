# massCollapse
Simulates the collapse off mass particles by their own gravitational 
fields.

## Build Instructions

## Dependencies
OpenCV v4 (for older versions the source code may be adopted)

OpenMP (only for gcc compiler and if parallelisation desired)

OpenCL v1.2 (check that you installed the right OpenCL version depending on your graphics card!)

### Configure
Use a out-of-tree build to not pollute your checkout:
```
mkdir build
cd build
```

Depending if you want to use a debug or release build, either run

```
# debugging
cmake -GNinja -DCMAKE_BUILD_TYPE=Debug -DSANITIZE_ADDRESS=On ..
```

or

```
# profiling
cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
```

or

```
# public release
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
```

### Actual Build
```
ninja
```
Usually, you do NOT need to re-run cmake.

### Run
```
./masscollapse
```

## Developer Notes
Do `export OMP_NUM_THREADS=X` before execution to run program on X threads, if OpenMP is available and compiled with gcc.

## Auto-Format
```
cd build
ninja format
```
