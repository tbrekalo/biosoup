cmake_minimum_required(VERSION 3.11)
project(biosoup VERSION 0.10.1
                LANGUAGES CXX
                DESCRIPTION "Biosoup is a c++ collection of header only data structures used for storage and logging in bioinformatics tools.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(biosoup_main_project ON)
endif ()

option(biosouppy "Generate python bindings for nucleic_acid" OFF)
option(biosouppy_ext_dir "Used by setup.py")

option(biosoup_install "Generate install target" OFF)
option(biosoup_build_tests "Build unit tests" OFF)

include(FetchContent)
if (biosoup_build_tests)
  find_package(GTest 1.10.0 QUIET)
  if (NOT GTest_FOUND)

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

if (biosouppy)
  find_package(pybind11 2.9.1 QUIET)
  if (NOT pybind11_FOUND)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.9.1)

    FetchContent_GetProperties(pybind11)
    if (NOT pybind11_POPULATED)
      FetchContent_Populate(pybind11) 
      add_subdirectory(
        ${pybind11_SOURCE_DIR}
        ${pybind11_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()
endif ()

add_library(biosoup INTERFACE)
add_library(${PROJECT_NAME}::biosoup ALIAS biosoup)

target_include_directories(biosoup INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

if (biosouppy)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

  pybind11_add_module(biosouppy
    ${PROJECT_SOURCE_DIR}/pybind/main.cc)
  target_link_libraries(biosouppy PRIVATE biosoup)
  add_library(${PROJECT_NAME}::biosouppy ALIAS biosouppy)
  
  if (biosouppy_ext_dir)
    install(TARGETS biosouppy
            LIBRARY DESTINATION ${biosouppy_ext_dir})
  endif()
endif ()

if (biosoup_install)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)

  install(
    TARGETS biosoup
    EXPORT ${PROJECT_NAME}Targets)
  install(
    DIRECTORY include/biosoup
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif ()

if (biosoup_build_tests)
  add_executable(biosoup_test
    test/nucleic_acid_test.cpp
    test/overlap_test.cpp
    test/progress_bar_test.cpp
    test/sequence_test.cpp
    test/timer_test.cpp)

  target_link_libraries(biosoup_test
    biosoup
    GTest::Main)
endif ()
