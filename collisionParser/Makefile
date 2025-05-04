# Detect compiler and flags used for ROOT
CC              = $(shell root-config --cxx)  # Get the C++ compiler used for ROOT
CCFLAGS         = $(shell root-config --cflags) -O2  # Use ROOT’s compiler flags

INCLUDE         = -I$(shell root-config --incdir)  # Explicitly add -I
ROOTLIBS        = $(shell root-config --libs)
LIBDIRS         = -L$(shell root-config --libdir)  # Explicitly add -L

collisionParser:	collisionParser.cpp
		$(CC) $(CCFLAGS) $(INCLUDE) $(LIBDIRS) $(ROOTLIBS) \
		-o collisionParser collisionParser.cpp

clean:
		rm -f collisionParser
