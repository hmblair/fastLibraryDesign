
cmake_minimum_required(VERSION 3.14)

PROJECT(fastLibraryDesign)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

add_executable(fastLibraryDesign main.cpp)
target_link_libraries(fastLibraryDesign argparse)
