mkdir build
cd build
cmake ..
cmake --build . --target noisereduction_driver
cd ..
copy /y "C:\Program Files\Mega-Nerd\libsndfile\bin\libsndfile-1.dll" ".\build\Debug\libsndfile-1.dll"
PAUSE
