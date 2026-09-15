# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)
project(metasona_consumer_example LANGUAGES C)

find_package(MetaSona 0.1 CONFIG REQUIRED)

add_executable(metasona_c_example metasona_example.c)
target_link_libraries(metasona_c_example PRIVATE MetaSona::metasona)
target_compile_features(metasona_c_example PRIVATE c_std_11)
