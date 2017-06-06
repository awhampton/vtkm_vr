cmake_minimum_required(VERSION 3.6)

set(VTKm_DIR "/home/users/hampton2/vtk-m/build_tbb/lib")
set(VTK_DIR "/home/users/hampton2/VTK-build")

find_package(VTKm REQUIRED QUIET
  OPTIONAL_COMPONENTS Serial CUDA TBB
  )

find_package(VTK REQUIRED)

#if(VTKm_CUDA_FOUND)
#  set(old_nvcc_flags ${CUDA_NVCC_FLAGS})
#  set(old_cxx_flags ${CMAKE_CXX_FLAGS})
#  vtkm_setup_nvcc_flags( old_nvcc_flags old_cxx_flags)
#  vtkm_disable_troublesome_thrust_warnings()

  # Cuda compiles do not respect target_include_directories
#  cuda_include_directories(${VTKm_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})
#  cuda_add_executable(vr vr.cu)

#  set(CUDA_NVCC_FLAGS ${old_nvcc_flags})
#  set(CMAKE_CXX_FLAGS ${old_cxx_flags})
#else()
  add_executable(vr vr.cxx)
#endif()

target_include_directories(vr PRIVATE ${VTKm_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})
target_link_libraries(vr PRIVATE ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
target_compile_options(vr PRIVATE ${VTKm_COMPILE_OPTIONS})
