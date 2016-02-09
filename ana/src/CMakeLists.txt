execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE CMAKE_CXX_FLAGS_TEMP)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE CMAKE_LD_FLAGS_TEMP)
string(STRIP ${CMAKE_CXX_FLAGS_TEMP} CMAKE_CXX_FLAGS)
string(STRIP ${CMAKE_LD_FLAGS_TEMP} CMAKE_LD_FLAGS)

set(CMAKE_CXX_COMPILER "g++")

# include TopBrussels SingleLep header search path
include_directories(${TopBrussels_SOURCE_DIR})				# TopAnalysisBase/interface TopTreeProducer/interface
include_directories(${TopBrussels_SOURCE_DIR}/BSMFourTops)		# SingleLep/interface/
include_directories(${TopBrussels_SOURCE_DIR}/BSMFourTops/ana/inc)

# add FourTopSingleLep lib search path
link_directories(${TopBrussels_SOURCE_DIR}/BSMFourTops/SingleLepAnalysis/src)

# add TopTreeAna76 TopTreeAnaContent76 lib search path
link_directories(~/lib)

# define output target
add_executable(SLMACROLocal FourTop_EventSelection_SingleLepton_Run2_Local.cc)

# link against root and FourTopSingleLep TopTreeAna76 TopTreeAnaContent76
target_link_libraries(SLMACROLocal PUBLIC FourTopSingleLep TopTreeAna76 TopTreeAnaContent76 PRIVATE ${CMAKE_LD_FLAGS})
install(TARGETS SLMACROLocal DESTINATION ~/CMSSW_7_6_3/src/TopBrussels/BSMFourTops/)
