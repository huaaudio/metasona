# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(
  metasona
  VERSION 0.1.0
  DESCRIPTION "Source-informed psychoacoustic metrics library"
  LANGUAGES C
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(MS_BUILD_SHARED "Build a shared native library" ON)
option(MS_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(MS_WARNINGS_AS_ERRORS "Treat project warnings as errors" ON)
set(MS_ABI_SOVERSION 1 CACHE STRING "Native shared-library ABI version")

set(MS_LIBRARY_TYPE STATIC)
if(MS_BUILD_SHARED)
  set(MS_LIBRARY_TYPE SHARED)
endif()

# Keep pinned upstream code unchanged and isolate its warning policy from
# project code. Prefix and hide its symbols in both static and shared builds.
add_library(ms_pocketfft OBJECT src/third_party/pocketfft/pocketfft.c)
add_library(ms_psychoacoustics OBJECT src/roughness_dw.c src/tonality_aures1985.c)
target_include_directories(ms_psychoacoustics PRIVATE src/third_party/pocketfft src include)
target_compile_definitions(ms_psychoacoustics PRIVATE
  MS_BUILDING_LIBRARY=1 TONALITY_AURES1985_STATIC=1
  roughness_dw=ms_roughness_kernel tonality_aures1985=ms_tonality_kernel)
set_target_properties(ms_psychoacoustics PROPERTIES
  C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO
  POSITION_INDEPENDENT_CODE YES C_VISIBILITY_PRESET hidden)
set_target_properties(ms_pocketfft PROPERTIES
  C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO
  POSITION_INDEPENDENT_CODE YES C_VISIBILITY_PRESET hidden)
if(MSVC)
  target_compile_options(ms_pocketfft PRIVATE /fp:precise /Brepro)
  target_compile_options(ms_psychoacoustics PRIVATE /fp:precise /Brepro)
else()
  target_compile_options(ms_pocketfft PRIVATE -ffp-contract=off)
  target_compile_options(ms_psychoacoustics PRIVATE -ffp-contract=off)
endif()

add_library(metasona ${MS_LIBRARY_TYPE}
  src/metasona.c
  src/ms_loudness.c
  src/ms_signal_metrics.c
  src/ms_sharpness.c
  $<TARGET_OBJECTS:ms_pocketfft>
  $<TARGET_OBJECTS:ms_psychoacoustics>
)
target_compile_definitions(metasona PRIVATE
  TONALITY_AURES1985_STATIC=1
  roughness_dw=ms_roughness_kernel tonality_aures1985=ms_tonality_kernel)
foreach(symbol make_cfft_plan destroy_cfft_plan cfft_backward cfft_forward cfft_length
               make_rfft_plan destroy_rfft_plan rfft_backward rfft_forward rfft_length)
  target_compile_definitions(ms_pocketfft PRIVATE "${symbol}=ms_pocketfft_${symbol}")
  target_compile_definitions(ms_psychoacoustics PRIVATE "${symbol}=ms_pocketfft_${symbol}")
  target_compile_definitions(metasona PRIVATE "${symbol}=ms_pocketfft_${symbol}")
endforeach()
add_library(MetaSona::metasona ALIAS metasona)

set_target_properties(metasona PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED YES
  C_EXTENSIONS NO
  POSITION_INDEPENDENT_CODE YES
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
  WINDOWS_EXPORT_ALL_SYMBOLS OFF
)
if(NOT SKBUILD)
  set_target_properties(metasona PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${MS_ABI_SOVERSION}
  )
endif()
target_compile_definitions(metasona PRIVATE MS_BUILDING_LIBRARY=1)
if(NOT MS_BUILD_SHARED)
  target_compile_definitions(metasona PUBLIC MS_STATIC=1)
  target_compile_definitions(ms_psychoacoustics PRIVATE MS_STATIC=1)
endif()
target_include_directories(metasona
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(MSVC)
  target_compile_options(metasona PRIVATE /W4 /permissive- /sdl /fp:precise /Brepro)
  target_link_options(metasona PRIVATE /Brepro)
  if(MS_WARNINGS_AS_ERRORS)
    target_compile_options(metasona PRIVATE /WX)
  endif()
else()
  target_compile_options(metasona PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes
    -Wsign-conversion -Wmissing-prototypes -Wdouble-promotion -Wformat=2
    -ffp-contract=off
  )
  if(MS_WARNINGS_AS_ERRORS)
    target_compile_options(metasona PRIVATE -Werror)
  endif()
  target_link_libraries(metasona PRIVATE m)
endif()

if(MS_ENABLE_SANITIZERS)
  if(MSVC)
    target_compile_options(ms_pocketfft PRIVATE /fsanitize=address)
    target_compile_options(ms_psychoacoustics PRIVATE /fsanitize=address)
    target_compile_options(metasona PRIVATE /fsanitize=address)
    target_link_options(metasona PRIVATE /INFERASANLIBS)
  elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(ms_pocketfft PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_compile_options(ms_psychoacoustics PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_compile_options(metasona PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(metasona PRIVATE -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "Sanitizers are not configured for ${CMAKE_C_COMPILER_ID}")
  endif()
endif()

if(SKBUILD)
  set(MS_RUNTIME_INSTALL_DIR "metasona/_native")
  set(MS_LIBRARY_INSTALL_DIR "metasona/_native")
else()
  set(MS_RUNTIME_INSTALL_DIR "${CMAKE_INSTALL_BINDIR}")
  set(MS_LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
endif()

install(TARGETS metasona
  EXPORT MetaSonaTargets
  RUNTIME DESTINATION "${MS_RUNTIME_INSTALL_DIR}" COMPONENT Runtime
  LIBRARY DESTINATION "${MS_LIBRARY_INSTALL_DIR}" COMPONENT Runtime
    NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Development
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development)

configure_package_config_file(
  cmake/MetaSonaConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MetaSonaConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MetaSona
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MetaSonaConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(EXPORT MetaSonaTargets
  NAMESPACE MetaSona::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MetaSona
  COMPONENT Development
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MetaSonaConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MetaSonaConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MetaSona
  COMPONENT Development
)
install(FILES LICENSE NOTICE THIRD_PARTY.md
  DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/MetaSona"
  COMPONENT Runtime
)
install(FILES LICENSE NOTICE THIRD_PARTY.md
  DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/MetaSona"
  COMPONENT Development
)

install(DIRECTORY LICENSES/
  DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/MetaSona/LICENSES"
  COMPONENT Runtime
)
install(DIRECTORY LICENSES/
  DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/MetaSona/LICENSES"
  COMPONENT Development
)
