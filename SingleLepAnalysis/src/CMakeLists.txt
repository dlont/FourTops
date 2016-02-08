execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE CMAKE_CXX_FLAGS_TEMP)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE CMAKE_LD_FLAGS_TEMP)
string(STRIP ${CMAKE_CXX_FLAGS_TEMP} CMAKE_CXX_FLAGS)
string(STRIP ${CMAKE_LD_FLAGS_TEMP} CMAKE_LD_FLAGS)

set(CMAKE_CXX_COMPILER "g++")

message( "CMAKE_C_FLAGS: ${CMAKE_CXX_FLAGS}" )
message( "CMAKE_LD_FLAGS: ${CMAKE_LD_FLAGS}" )

# Specify TopBrussels include and lib path
set(TopBrussels_DIR "/user/dlontkov/CMSSW_7_6_3/src/TopBrussels/")

#set(TopBrussels_DIR "" CACHE "TopBrussels installation paths")
if(IS_DIRECTORY ${TopBrussels_DIR})
        message(STATUS "Found dir: ${TopBrussels_DIR}")
        set(TopBrussels_LIBRARY_DIR ${TopBrussels_DIR})
        set(TopBrussels_INCLUDE_DIR ${TopBrussels_DIR})
else()
        message(SEND_ERROR "ERROR: ${TopBrussels_DIR} does not exist")
endif()

# include TopBrussels header search path
include_directories(${TopBrussels_INCLUDE_DIR})
# link against TopBrussels
#link_directories(${TopBrussels_LIBRARY_DIR}/)

# include TopBrussels header search path
include_directories(${BSMFourTops_SOURCE_DIR}/SingleLepAnalysis/interface)

add_library(FourTopSingleLep SHARED CutsTable.cc EventBDT.cc HadronicTopReco.cc Trigger.cc Zpeak.cc)
# link against root
target_link_libraries(FourTopSingleLep PRIVATE ${CMAKE_LD_FLAGS})
install(TARGETS FourTopSingleLep DESTINATION ~/lib)

