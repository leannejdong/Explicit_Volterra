cmake_minimum_required(VERSION 3.12)
project(Explicit_Volterra)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(Explicit_Volterra main.cpp volterra.cpp)

target_compile_options(
        Explicit_Volterra PRIVATE -Wall -Wextra -pedantic -Werror
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g -O1")