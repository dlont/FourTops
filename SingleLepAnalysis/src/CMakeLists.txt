execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE CMAKE_CXX_FLAGS_TEMP)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE CMAKE_LD_FLAGS_TEMP)
string(STRIP ${CMAKE_CXX_FLAGS_TEMP} CMAKE_CXX_FLAGS)
string(STRIP ${CMAKE_LD_FLAGS_TEMP} CMAKE_LD_FLAGS)

set(CMAKE_CXX_COMPILER "g++")

message( "CMAKE_C_FLAGS: ${CMAKE_CXX_FLAGS}" )
message( "CMAKE_LD_FLAGS: ${CMAKE_LD_FLAGS}" )

# include TopBrussels header search path
include_directories(${TopBrussels_SOURCE_DIR})
include_directories(${TopBrussels_SOURCE_DIR}/BSMFourTops)

# add TopTreeAna76 TopTreeAnaContent76 lib search path
link_directories(~/lib)

# include TopBrussels header search path
#include_directories(${BSMFourTops_SOURCE_DIR}/SingleLepAnalysis/interface)
#message(STATUS "lib SingleLepton dir: ${BSMFourTops_SOURCE_DIR}")

# define output target
add_library(FourTopSingleLep SHARED CutsTable.cc EventBDT.cc HadronicTopReco.cc Trigger.cc Zpeak.cc)

# link against root and TopTreeAna76 TopTreeAnaContent76
target_link_libraries(FourTopSingleLep PUBLIC TopTreeAna76 TopTreeAnaContent76 PRIVATE ${CMAKE_LD_FLAGS})
install(TARGETS FourTopSingleLep DESTINATION ~/lib)

