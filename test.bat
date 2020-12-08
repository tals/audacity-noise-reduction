mkdir build
cd build
cmake ..
cmake --build . --target test_runner
cd Debug
mkdir temp
.\test_runner.exe
PAUSE
