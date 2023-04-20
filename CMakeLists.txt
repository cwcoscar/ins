cmake_minimum_required(VERSION 3.0.2)
project(ins)

add_definitions(-std=c++14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  message_generation
  tf
  sensor_msgs
  geometry_msgs
  ublox_msgs
)

## Generate messages in the 'msg' folder
add_message_files(DIRECTORY msg)

## Generate added messages and services with any dependencies listed here
generate_messages(
  DEPENDENCIES
  std_msgs
)

catkin_package(
 CATKIN_DEPENDS 
 roscpp 
 rospy 
 std_msgs 
 message_runtime
 sensor_msgs
)

###########
## Build ##
###########
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

add_executable(main src/main.cpp src/ins.cpp)

target_link_libraries(main ${catkin_LIBRARIES})

add_dependencies(main
  ${catkin_EXPORTED_TARGETS}
)

#############
## Install ##
#############

install(TARGETS
  main
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )

