project(kima)
# Just for Cmake to stop complaining
cmake_minimum_required(VERSION 3.6)

add_executable(kima)

# Define source files
set(source_dir ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(kima_sources
    ${source_dir}/distributions/Fixed.cpp
    ${source_dir}/data.cpp
    # ${source_dir}/RVConditionalPrior.cpp
    # ${source_dir}/RVmodel.cpp
    ${source_dir}/system.cpp
    ${source_dir}/planet.cpp
    ${source_dir}/main.cpp
)
target_sources(kima PUBLIC ${kima_sources})
# Include .h from kima
target_include_directories(kima PUBLIC ${source_dir})

# Set libs folder for all external code
set(third_party_dir ${CMAKE_CURRENT_SOURCE_DIR}/third-party)
set(libs_dir ${CMAKE_CURRENT_SOURCE_DIR}/libs)

# ----------------------- EIGEN ---------------------------------
# Include Eigen
target_include_directories(kima PUBLIC ${third_party_dir}/eigen)

# ----------------------- DNEST4 --------------------------------
# Compile Dnest4 lib - Building everything here to not mess with DNest4 files
add_library(dnest4 STATIC)

# Define DNest4 sources and add include dirs
set(dnest4_dir ${third_party_dir}/DNest4/code)
set(dnest4_sources
    ${dnest4_dir}/DNest4.h
    ${dnest4_dir}/Barrier.cpp
    ${dnest4_dir}/CommandLineOptions.cpp
    ${dnest4_dir}/Level.cpp
    ${dnest4_dir}/LikelihoodType.cpp
    ${dnest4_dir}/Options.cpp
    ${dnest4_dir}/RNG.cpp
    ${dnest4_dir}/Utils.cpp
    ${dnest4_dir}/Distributions/Cauchy.cpp
    ${dnest4_dir}/Distributions/ContinuousDistribution.cpp
    ${dnest4_dir}/Distributions/Exponential.cpp
    ${dnest4_dir}/Distributions/Gaussian.cpp
    ${dnest4_dir}/Distributions/Jeffreys.cpp
    ${dnest4_dir}/Distributions/Kumaraswamy.cpp
    ${dnest4_dir}/Distributions/Laplace.cpp
    ${dnest4_dir}/Distributions/LogUniform.cpp
    ${dnest4_dir}/Distributions/Pareto.cpp
    ${dnest4_dir}/Distributions/Rayleigh.cpp
    ${dnest4_dir}/Distributions/Triangular.cpp
    ${dnest4_dir}/Distributions/Uniform.cpp
    ${dnest4_dir}/RJObject/ConditionalPriors/BasicCircular.cpp
    ${dnest4_dir}/RJObject/ConditionalPriors/ClassicMassInf.cpp
    ${dnest4_dir}/RJObject/ConditionalPriors/ClassicMassInf1D.cpp
    ${dnest4_dir}/RJObject/ConditionalPriors/ConditionalPrior.cpp
)
target_sources(dnest4 PRIVATE ${dnest4_sources})
target_include_directories(dnest4 PUBLIC ${dnest4_dir})

# Link DNest4 lib. Pthread comes from the PUBLIC interface
target_link_libraries(kima PUBLIC dnest4)

# ------------------------ spdlog lib ---------------------------------
add_subdirectory(${third_party_dir}/spdlog-1.4.2)
target_link_libraries(kima PRIVATE spdlog)

# ---------------- OS-specific definitions and libs -------------------
# If we are on windows, handle the lack of unistd and getopt files from DNest4
if (WIN32)
    set(libs_sources
        ${libs_dir}/unistd/getopt.c
    )
    target_sources(dnest4 PRIVATE ${libs_sources})
    target_include_directories(dnest4 PRIVATE ${libs_dir}/unistd)
endif()

# OS-specific compile options and definitions (can't target release and debug specifically yet). 
if (WIN32)
    # WINDOWS: Define _USE_MATH_DEFINES to have M_PI and M_PI_2 definitions
    target_compile_definitions(kima PUBLIC _USE_MATH_DEFINES)
    
    # Debug flags
    # target_compile_options(dnest4 PUBLIC -Wall)
endif(WIN32)
if(UNIX)
    # UNIX: Add pthread for threading
    target_link_libraries(dnest4 PUBLIC pthread)
    #target_compile_options(kima PRIVATE -pthread)

    # Debug flags
    # target_compile_options(dnest4 PUBLIC -Wall -Wextra -Wpedantic -Wconversion)
endif(UNIX)

# Common definitions for all OSs
target_compile_definitions(kima PUBLIC EIGEN_MPL2_ONLY)
target_compile_features(kima PUBLIC cxx_std_11 c_std_11)