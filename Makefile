ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
UNAME = $(shell uname -s)
ifeq ($(UNAME), Darwin)
DllSuf        = dylib
endif
ifeq ($(UNAME), Linux)
DllSuf        = so
endif

OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
ROOTLIBS_NoTMVA = $(shell root-config --libs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam
ROOTGLIBS       = $(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
ROOTLIBS       += -L$(ROOFITSYS)/lib

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES -I..
CXX           = g++
CXXFLAGS	  = -g3 -gdwarf-2 -O -Wall -fPIC $(DEFINES)  -I./TMVA/include
ifeq ($(UNAME), Darwin)
CXXFLAGS      += -I/opt/local/include
endif
LD		      = g++
LDFLAGS		  = -g3 -gdwarf-2 -O -Wall -fPIC
ifeq ($(UNAME), Darwin)
SOFLAGS       = -dynamiclib
endif
ifeq ($(UNAME), Linux)
SOFLAGS       = -shared
endif

CXXFLAGS	  += $(ROOTCFLAGS)
LIBS		  = -I./TMVA/include -L./TMVA/lib $(ROOTLIBS) -lEG -I.. -L. -L libTopTreeAna76 -L TopTreeAnaContent76
LIBS_NoTMVA   = $(ROOTLIBS_NoTMVA) -lEG -I.. -L. -L libTopTreeAna76 -L TopTreeAnaContent76
ifeq ($(UNAME), Darwin)
LIBS          += -I/opt/local/include
endif
GLIBS		  = $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES       = $(wildcard SingleLepAnalysis/src/*cc)
HEADERS       = $(wildcard SingleLepAnalysis/interface/*.h)
OBJECTS		  = $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		  = $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	  = $(SOURCES:.$(SrcSuf)=.$(DllSuf))

all:  libFourTopSingleLep.$(DllSuf)
	cp libFourTopSingleLep.$(DllSuf) ~/lib/

SLMACRO: ana/src/FourTop_EventSelection_SingleLepton_Run2_New.cc libFourTopSingleLep.$(DllSuf)
	@echo "compiling grid macro"
	g++ -g3 -gdwarf-2 -std=c++11 -I.. -I. -L ~/lib -l TopTreeAnaContent76 -l TopTreeAna76 -l FourTopSingleLep -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` ana/src/FourTop_EventSelection_SingleLepton_Run2_New.cc -o SLMACRO

SLMACROLocal: ana/src/FourTop_EventSelection_SingleLepton_Run2_Local.cc libFourTopSingleLep.$(DllSuf)
	@echo "compiling local macro"
	g++ -g3 -gdwarf-2 -std=c++11 -I.. -I. -L ~/lib -l TopTreeAnaContent76 -l TopTreeAna76 -l FourTopSingleLep -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` ana/src/FourTop_EventSelection_SingleLepton_Run2_Local.cc -o SLMACROLocal

SLMACROEL: ana/src/FourTop_EventSelection_SingleLepton_Run2_New.cc libFourTopSingleLep.$(DllSuf)
	@echo "compiling grid electron macro"
	g++ -g3 -gdwarf-2 -std=c++11 -I.. -I. -L ~/lib -l TopTreeAnaContent76 -l TopTreeAna76 -l FourTopSingleLep -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` ana/src/FourTop_EventSelection_SingleLepton_Run2_New.cc -o SLMACROEL

SLMACROLocalEl: ana/src/FourTop_EventSelection_SingleLepton_Run2_Local.cc libFourTopSingleLep.$(DllSuf)
	@echo "compiling local electron macro"
	g++ -g3 -gdwarf-2 -std=c++11 -I.. -I. -L ~/lib -l TopTreeAnaContent76 -l TopTreeAna76 -l FourTopSingleLep -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` ana/src/FourTop_EventSelection_SingleLepton_Run2_Local.cc -o SLMACROLocalEl

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(DEPENDS) macros/*.exe *Dict.* *.$(DllSuf) core 

.SUFFIXES: .$(SrcSuf) .C .o .$(DllSuf)

###

Dict.$(SrcSuf): $(HEADERSDIC) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERSDIC) ./LinkDef.h

libFourTopSingleLep.$(DllSuf): $(OBJECTS)
	@echo "Building libFourTopSingleLep..."
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

# specific stuff for btag eff analysis ONLY

ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO 

macros/%.exe: macros/%.cc $(HEADERS) libFourTopSingleLep.$(DllSuf)
	$(LD) -lFourTopSingleLep $(LIBS) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@

SOURCES_MACROS = $(wildcard macros/*.cc)

macros: $(SOURCES_MACROS:.cc=.exe)
