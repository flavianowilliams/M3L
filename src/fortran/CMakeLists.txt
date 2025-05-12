cmake_minimum_required(VERSION 3.30.1)

project(python_fortran LANGUAGES C Fortran)

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()

# configurando ambiente de instalação 
include(cmake/config.cmake)

# imprimindo informações importantes
message(STATUS "Diretório do código fonte: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Diretório de arquivos binários: ${CMAKE_CURRENT_BINARY_DIR}")

# setando instruções para o compilador fortran 
set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall" CACHE STRING "" FORCE)

# carrengando o pacote numpy 
find_package(Python COMPONENTS NumPy REQUIRED)

# definindo o diretório do include fortranobject.c 
execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Definindo diretório de headers do pacote f2py: ${F2PY_INCLUDE_DIR}")

# definindo diretório dos includes do NumPy
execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Definindo diretório de headers do numpy: ${NUMPY_INCLUDE_DIR}")

# gerando arquivos module.c f2pywrappers2.f90 
#set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${FORTRAN_FILES}")

add_custom_command(
  OUTPUT ${PROJECT_NAME}module.c ${PROJECT_NAME}-f2pywrappers2.f90
  DEPENDS ${FORTRAN_FILES}
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py "${FORTRAN_FILES}" -m ${PROJECT_NAME} --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} 
)

message(STATUS "Dependências: ${PROJECT_NAME}module.c e ${PROJECT_NAME}-f2pywrappers2.f90 em ${CMAKE_CURRENT_BINARY_DIR}")

# gerando biblioteca principal
add_library(
  ${PROJECT_NAME} MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}module.c"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-f2pywrappers2.f90"
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  "${FORTRAN_FILES}"
)

message(STATUS "Definindo biblioteca ${PROJECT_NAME}")

# linkando bibliotecas do python na biblioteca principal
target_link_libraries(${PROJECT_NAME} PRIVATE Python::NumPy)

message(STATUS "Linkando biblioteca Python.h em ${PROJECT_NAME}")

## adicionando diretorios de headers no target principal
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC
  ${F2PY_INCLUDE_DIR}
  ${NUMPY_INCLUDE_DIR}
  ${Python_INCLUDE_DIRS}
)

message(STATUS "Adicionando diretório de headers em ${PROJECT_NAME}")

# definindo nome do target principal
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

message(STATUS "Alterando o sufixo para ${PROJECT_NAME}")

# definindo diretorio de instalação 
install(TARGETS ${PROJECT_NAME} DESTINATION .)

message(STATUS "Definindo local de instalação do target ${PROJECT_NAME}")
