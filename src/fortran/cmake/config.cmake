# definindo modulos
#set(FILES_LIST "fortran_module.f90")
#set(FORTRAN_FILES "")
#foreach(VALUE IN LISTS FILES_LIST)
#  set(ITEM "${CMAKE_CURRENT_SOURCE_DIR}/${VALUE}")
#  list(APPEND FORTRAN_FILES ${ITEM})
#endforeach()

#set(FILES_LIST "ff_module.f90;fortran_module.f90")
#set(FORTRAN_FILES "")
#foreach(VALUE IN LISTS FILES_LIST)
#  set(ITEM "${CMAKE_CURRENT_SOURCE_DIR}/${VALUE}")
#  list(APPEND FORTRAN_FILES ${ITEM})
#endforeach()

set(FORTRAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/fortran_module.f90")
# definindo funcoes
function(teste)
  message("tudo blz!")
endfunction()
