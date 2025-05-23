macro(function_module VAR)

  message("")
  message("***************************************************************************************************")

  set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${VAR}_module.f90")

  set(modulename ${VAR})

  message(STATUS "Módulo: ${modulename}")

  message("--------------------------------------------------------------------------------------------------")

  message(STATUS "Source: ${fortran_src_file}")

  message(STATUS "Dependências: ${modulename}module.c e ${modulename}-f2pywrappers2.f90 em ${CMAKE_CURRENT_BINARY_DIR}")

  #gerando arquivos module.c e f2pywrappers2 necessários para f2py
  add_custom_command(
    OUTPUT ${modulename}module.c ${modulename}-f2pywrappers2.f90
    DEPENDS ${fortran_src_file}
    COMMAND ${Python_EXECUTABLE} -m numpy.f2py "${fortran_src_file}" -m ${modulename} --lower
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} 
  )

  # gerando biblioteca principal

  add_library(
    ${modulename} MODULE
    "${CMAKE_CURRENT_BINARY_DIR}/${modulename}module.c"
    "${CMAKE_CURRENT_BINARY_DIR}/${modulename}-f2pywrappers2.f90"
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    "${fortran_src_file}"
  )

  message(STATUS "Definindo biblioteca ${modulename}")

  # linkando bibliotecas do python na biblioteca principal
  target_link_libraries(${modulename} PRIVATE Python::NumPy)

  message(STATUS "Linkando biblioteca Python.h em ${modulename}")

  ## adicionando diretorios de headers no target principal
  target_include_directories(
    ${modulename}
    PUBLIC
    ${F2PY_INCLUDE_DIR}
    ${NUMPY_INCLUDE_DIR}
    ${Python_INCLUDE_DIRS}
  )

  message(STATUS "Adicionando diretório de headers em ${modulename}")

  # definindo nome do target principal
  set_target_properties(${modulename} PROPERTIES PREFIX "")

  message(STATUS "Alterando o sufixo para ${modulename}")

  # definindo diretorio de instalação 
  install(TARGETS ${modulename} DESTINATION .)

  message(STATUS "Definindo local de instalação do target ${modulename}")

  message("***************************************************************************************************")

endmacro()

