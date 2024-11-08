cmake_minimum_required(VERSION 3.17.2...3.29)

project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

set(filename "md_fortran")

message("Carregando interpretador python e módulos Module e Numpy!!!!")

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

message("OK!")
message("Definindo diretório do arquivo fortranobject.c como pré-requisito para compilação do módulo.")

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

message("OK!")
message("Gerando arquivos ${filename}module.c e ${filename}-f2pwrappers2.f90 como pré-requisito para compilação do módulo.")

add_custom_command(
  OUTPUT ${filename}module.c ${filename}-f2pywrappers2.f90
  DEPENDS ${filename}.f90
  VERBATIM
  COMMAND "${PYTHON_EXECUTABLE}" -m numpy.f2py "${CMAKE_CURRENT_SOURCE_DIR}/${filename}.f90" -m ${filename} --lower
)

message("OK!")
message("Criando a biblioteca * em Python::Module")

python_add_library(
  ${filename} MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/${filename}module.c"
  "${CMAKE_CURRENT_BINARY_DIR}/${filename}-f2pywrappers2.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/${filename}.f90" 
  WITH_SOABI
)

target_link_libraries(${filename} PRIVATE fortranobject)

message("OK!")
message("Gerando as regras de instalação para os alvos pré-definidos.")

install(TARGETS ${filename} DESTINATION .)
