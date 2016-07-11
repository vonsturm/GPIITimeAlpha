###################################################################
# This Makefile was created using the bat-project script
# for project TimeAlpha
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -O2 -Wall -fPIC -Wno-deprecated
LD           = /opt/exp_software/gerda/common/sw/binutils/linux-scientific-6-x86_64/2.25/bin/ld -m elf_x86_64
LDFLAGS      = -g -O2 

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

# BAT
CXXFLAGS += $(shell bat-config --cflags)
LIBS := $(shell bat-config --libs)

# CUBA
CXXFLAGS += -I$(CUBA_BASE_DIR)/include
LIBS += -L$(CUBA_BASE_DIR)/lib -lcuba

# ROOT
CXXFLAGS += $(shell root-config --cflags)
LIBS += $(shell root-config --libs) -lMinuit

# GERDA-ADA
CXXFLAGS += -I$(GERDA_BASE_DIR)/include/gerda-ada/
LIBS += -L$(GERDA_BASE_DIR)/lib -lgerda-ada-core -lgerda-ada-calib-ged -lgerda-ada-dataprod -lgerda-ada-evtviewer -lgerda-ada-monitoring -lgerda-ada-psd-base -lgerda-ada-stats

# GELATIO
CXXFLAGS += -I$(GERDA_BASE_DIR)/include/gelatio
LIBS += -L$(GERDA_BASE_DIR)/lib -lGELATIODecoders -lGELATIOManagement -lGELATIOModules -lGELATIOUtilities


# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \
	runTimeAlpha.cxx TimeAlpha.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
MYPROGS     = runTimeAlpha

GARBAGE      = $(CXXOBJS) *.o *~ link.d $(MYPROGS)

# targets
all : runTimeAlpha

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm -f $(GARBAGE)

runTimeAlpha : $(CXXOBJS)
	$(CXX) $(LDFLAGS) $(CXXOBJS) $(LIBS) -o runTimeAlpha

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS)
	@echo c++ objs  : $(CXXOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
