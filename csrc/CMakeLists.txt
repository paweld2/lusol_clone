cmake_minimum_required(VERSION 3.8)

INCLUDE(FindPkgConfig)
include(GNUInstallDirs)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


PROJECT(lusol)

ENABLE_TESTING()

CONFIGURE_FILE("lusol.pc.in" "lusol.pc" @ONLY)

set(LUSOL_HEADER_FILES commonlib.h myblas.h lusol.h lusolio.h lusolmain.h myblas.h)
include_directories ( . )

ADD_LIBRARY(lusol SHARED lusol.c lusolio.c mmio.c myblas.c commonlib.c hbio.c ${LUSOL_HEADER_FILES})
target_link_libraries(lusol ${BLAS_LIBRARIES})
target_link_libraries(lusol ${LAPACK_LIBRARIES})
target_link_libraries(lusol dl)
target_link_libraries(lusol m)

ADD_EXECUTABLE(lusolmain lusolmain.c ${LUSOL_HEADER_FILES})
target_link_libraries(lusolmain lusol)

INSTALL(TARGETS lusol lusolmain
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(FILES ${LUSOL_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lusol)
INSTALL(FILES "${CMAKE_BINARY_DIR}/csrc/lusol.pc" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")