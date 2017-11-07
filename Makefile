PROJNAME := a.out

# include debug information: 1=yes, 0=no
DEBUG :=  

DEPEND := .dependencies

BINDIR := bin
INCDIR := src
SRCDIR := src
OBJDIR := obj

ADDSRCDIR := 
ADDINCDIR := 

# CXX      := $(shell which g++)
CXX := g++
LIBS     :=  -lm
CXXFLAGS :=  -I$(INCDIR)c -g
#-I$(ADDINCDIR)

ITPPFLAGS := `itpp-config --cflags`
ITPPLIBS  := `itpp-config --static --libs`
LIBS      += $(ITPPLIBS)
CXXFLAGS  += $(ITPPFLAGS)

ifdef DEBUG
#SUFFIX  := .dbg
SUFFIX   :=
#CXXFLAGS += -g -Wall
CXXFLAGS += -Wall
else
SUFFIX   :=
CXXFLAGS += -O2
endif

OBJSUF := .o$(SUFFIX)

SRC    := $(wildcard $(SRCDIR)/*.cpp) 
ADDSRC := $(wildcard $(ADDSRCDIR)/*.cpp)
OBJ    := $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o$(SUFFIX)) $(ADDSRC:$(ADDSRCDIR)/%.cpp=$(OBJDIR)/%.o$(SUFFIX)) 
BIN    := $(BINDIR)/$(PROJNAME)$(SUFFIX)

.phony: default clean tags bin depend

default: depend bin 

clean:
	@echo remove all objects
	@rm -f $(OBJDIR)/*

tags:
	@echo update tag table
	@ctags inc/*.h src/*.cpp

bin:    $(OBJ)
	@echo
	@echo 'creating binary "$(BIN)"'
	@$(CXX) -o $(BIN) $(OBJ) $(LIBS)
	@echo '... done'
	@echo

depend:
	@echo
	@echo 'checking dependencies'
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $(SRC) $(ADDSRC)                  \
         | sed '\''s@\(.*\)\.o[ :]@$(OBJDIR)/\1.o$(SUFFIX):@g'\''               \
         >$(DEPEND)'
	@echo

$(OBJDIR)/%.o$(SUFFIX): $(SRCDIR)/%.cpp
	@echo 'compiling object file "$@" ...'
	@$(CXX) -c -o $@ $(CXXFLAGS) $<

$(OBJDIR)/%.o$(SUFFIX): $(ADDSRCDIR)/%.cpp
	@echo 'compiling object file "$@" ...'
	@$(CXX) -c -o $@ $(CXXFLAGS) $<

-include $(DEPEND)

