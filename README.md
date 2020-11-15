# audacity-noise-reduction
2-pass noise reduction, pulled out of Audacity

# Usage

## Driver
The following commands will build and run the driver, i.e. the main application.
```bash
$ ./build_driver.sh
$ ./build/noisereduction_driver --help
```

In the directory `test` there are some examples to test the driver.  
They can be run with the corresponding shell script `run.sh`.  
For example:
```bash
$ cd test/separate-track-for-profiling/mono
$ ./run.sh
```

## Test
The following command will build and run automated tests.
```bash
$ ./test.sh
```

## Clean
The following command will clean the build directory.
```bash
$ ./clean.sh
```

# Dependencies
- CMake
- libsndfile

## Linux
```bash
$ sudo apt-get install cmake
$ sudo apt-get install libsndfile-dev
```

## MacOS
```bash
$ brew install cmake
$ brew install libsndfile
```
