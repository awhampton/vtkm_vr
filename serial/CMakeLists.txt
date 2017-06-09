cmake_minimum_required(VERSION 3.6)

set(VTKm_DIR "/home/users/hampton2/vtk-m/build/lib")
set(VTK_DIR "/home/users/hampton2/VTK-build")

find_package(VTKm REQUIRED QUIET
  )

find_package(VTK REQUIRED)

add_executable(vr vr.cxx)

target_include_directories(vr PRIVATE ${VTKm_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})
target_link_libraries(vr PRIVATE ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
target_compile_options(vr PRIVATE ${VTKm_COMPILE_OPTIONS})
