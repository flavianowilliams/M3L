cmake_minimum_required(VERSION 3.30.1)

project(python_fortran LANGUAGES C Fortran)

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()

# imprimindo informações importantes
message(STATUS "Diretório do código fonte: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Diretório de arquivos binários: ${CMAKE_CURRENT_BINARY_DIR}")

# carrengando o pacote numpy 
find_package(Python COMPONENTS NumPy REQUIRED)

message(STATUS "Definindo diretório de headers do pacote f2py: ${F2PY_INCLUDE_DIR}")

# definindo o diretório do include fortranobject.c 
execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Definindo diretório de headers do numpy: ${NUMPY_INCLUDE_DIR}")

# definindo diretório dos includes do NumPy
execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

#flags de compilação
set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" CACHE STRING "" FORCE)

# configurando ambiente CMAKE 
include(cmake/config.cmake)

set(filename "fortran")
#
foreach(VALUE IN LISTS filename)
  function_module(${VALUE})
endforeach()

#set(filename "modulo1")
#
#function_module(filename)

