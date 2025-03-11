# definindo versão mínima do cmake 
cmake_minimum_required(VERSION 3.30.1)

# definindo nome do projeto
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

# definindo nome do target principal
set(filename "md_fortran")

#
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

#
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# adicionando biblioteca fortranobject
add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")

# linkando biblioteca NumPy no target fortranobject
target_link_libraries(fortranobject PUBLIC Python::NumPy)

# incluindo diretório dos headers em fortranobject
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")

# configurando propriedades no target fortranobject
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

# setando instruções para o compilador fortran 
set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall" CACHE STRING "" FORCE)

#aplicando comando customizado para gerar as saídas md_fortranmodule.c e md_fortran-f2pywrappers2.f90
add_custom_command(
  OUTPUT ${filename}module.c ${filename}-f2pywrappers2.f90
  DEPENDS ${filename}.f90
  VERBATIM
  COMMAND "${PYTHON_EXECUTABLE}" -m numpy.f2py "${CMAKE_CURRENT_SOURCE_DIR}/${filename}.f90" -m ${filename} --lower
)

python_add_library(
  ${filename} MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/${filename}module.c"
  "${CMAKE_CURRENT_BINARY_DIR}/${filename}-f2pywrappers2.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/${filename}.f90" 
  WITH_SOABI
)

# linkando a biblioteca fortranobject no target filename
target_link_libraries(${filename} PRIVATE fortranobject)

# instalando target filename na pasta local -- necessário o comando make install
install(TARGETS ${filename} DESTINATION .)
