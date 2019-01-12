# massCollapse
Simulates the collapse off mass particles by their own gravitational 
fields.

## Build Instructions

## Dependencies
OpenCV v4 (can be disabled in src/global.h; for older versions the source code may be adopted)

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
-

## Auto-Format
```
cd build
ninja format
```
