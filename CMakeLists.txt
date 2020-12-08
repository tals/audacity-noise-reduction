cmake_minimum_required(VERSION 3.1)

if(NOT WIN32)
	# do not use clang on Windows,
	# will use MSVC anyway with default Visual Studio generator, even if clang is specified
	set(CMAKE_C_COMPILER clang)
	set(CMAKE_CXX_COMPILER clang++)
endif()

if(NOT WIN32)
	# mmap not supported for Windows
	set(USE_MMAP TRUE)
endif()

project(AudacityNoiseReduction VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF) # disable compiler-specific extensions

if(MSVC)
	# for MSVC, only warning level 3, is much too verbose with -Wall
	add_compile_options(-W3)
else()
	add_compile_options(-Wall)
endif()

################################################################################
# Create core library
set(CORE_SRC_FILES
    noisereduction/NoiseReduction.cpp
    noisereduction/InputTrack.cpp
    noisereduction/OutputTrack.cpp
    noisereduction/TrackUtils.cpp
    noisereduction/RealFFTf.cpp
    noisereduction/Extensions.cpp
)
if(USE_MMAP)
	list(APPEND CORE_SRC_FILES noisereduction/SndMmap.cpp)
endif()

# Add headers as well so they are listed in Visual Studio
set(CORE_HEADER_FILES
	noisereduction/NoiseReduction.h
	noisereduction/InputTrack.h
	noisereduction/OutputTrack.h
	noisereduction/TrackUtils.h
	noisereduction/SndContext.h
	noisereduction/RealFFTf.h
	noisereduction/Types.h
	noisereduction/MemoryX.h
	noisereduction/cxxopts.hpp
	noisereduction/loguru.hpp	
)
if(USE_MMAP)
	list(APPEND CORE_HEADER_FILES noisereduction/SndMmap.h)
endif()
	
add_library(noisereduction_core STATIC ${CORE_SRC_FILES} ${CORE_HEADER_FILES})

if(USE_MMAP)
	target_compile_definitions(noisereduction_core PUBLIC USE_MMAP)
endif()

find_package(Threads REQUIRED)

if(WIN32)
	# default installation dir for libsndfile on Windows
	set(LIBSNDFILE_DIR "C:\\Program Files\\Mega-Nerd\\libsndfile")
endif()

if(WIN32)
	# on Windows, need full name of library and path hint to find it
	find_library(LIBSNDFILE_LIBRARY libsndfile-1 REQUIRED
		HINTS "${LIBSNDFILE_DIR}\\lib"
	)
else()
	find_library(LIBSNDFILE_LIBRARY sndfile REQUIRED)
endif()

set(CORE_LIBRARIES
    Threads::Threads
    ${LIBSNDFILE_LIBRARY}
)

if(NOT WIN32)
	# do not need dl library on Windows
	find_library(DL_LIBRARY dl REQUIRED)
	list(APPEND CORE_LIBRARIES ${DL_LIBRARY})
endif()

if(WIN32)
	# must add include dir on Windows
    target_include_directories(noisereduction_core PUBLIC "${LIBSNDFILE_DIR}\\include")
endif()

target_link_libraries(noisereduction_core PUBLIC ${CORE_LIBRARIES})

################################################################################
# Create driver executable
add_executable(noisereduction_driver noisereduction/main.cpp)

target_link_libraries(noisereduction_driver PRIVATE noisereduction_core)

################################################################################
# Create test executable
set(TEST_SRC_FILES
    noisereduction/tests/test_main.cpp
    noisereduction/tests/test_audacity.cpp
)
if(USE_MMAP)
    list(APPEND TEST_SRC_FILES noisereduction/tests/test_mmap_snd.cpp)
endif()

# Add headers as well so they are listed in Visual Studio
set(TEST_HEADER_FILES
	noisereduction/tests/catch.hpp
)

add_executable(test_runner ${TEST_SRC_FILES} ${TEST_HEADER_FILES})

# Define samples directory relative to executable
if(MSVC)
	# two levels up because MSVC puts exe under build/Debug
	target_compile_definitions(test_runner PRIVATE SAMPLES_DIR="../../samples")
else()
	target_compile_definitions(test_runner PRIVATE SAMPLES_DIR="../samples")
endif()

target_link_libraries(test_runner PRIVATE noisereduction_core)
