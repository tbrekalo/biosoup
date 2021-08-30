find_package(benchmark 1.5.6 QUIET)
if (NOT benchmark_FOUND)
  include(FetchContent)
  
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark
    GIT_TAG v1.5.6
  )

  if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
    FetchContent_MakeAvailable(benchmark)
  else()
    FetchContent_Populate(benchmark)
    add_subdirectory(
      benchmark_SOURCE_DIR
      benchmark_BINARY_DIR
      EXCLUDE_FROM_ALL
    )
  endif()
endif()

set(${PROJECT_NAME}_benchmark_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/nucleic_acid_benchmark.cpp
)

add_executable(${PROJECT_NAME}_benchmark ${${PROJECT_NAME}_benchmark_SOURCES})
target_link_libraries(${PROJECT_NAME}_benchmark biosoup benchmark::benchmark)
