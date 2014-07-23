## *******************************************************************
##   (C) Copyright 2014 Leiden Institute of Advanced Computer Science
##   Universiteit Leiden
##   All Rights Reserved
## *******************************************************************
## Extractor (library)
## *******************************************************************
## FILE INFORMATION:
##   File:     Makefile
##   Author:   Jonathan K. Vis
##   Revision: 2.01a
##   Date:     2014/07/22
## *******************************************************************
## DESCRIPTION:
##  Build the Extractor library for python.
## *******************************************************************

SOURCES=extractor.cc
TARGET=_extractor.so
DEBUG=debug.cc

CXX=g++
CFLAGS=-c -fpic -Wall -Wextra -O3 -D__debug__
LDFLAGS=-Wall -O3 -shared

SWIG=swig
SWIGFLAGS=-c++ -python
INCLUDES=-I/usr/include/python2.7

WRAPPER=$(SOURCES:.cc=.py) $(SOURCES:.cc=)_wrap.cxx
OBJECTS=$(SOURCES:.cc=.o) $(WRAPPER:.cxx=.o)

.PHONY: all clean debug rebuild tarball

all: $(TARGET)

$(TARGET): $(OBJECTS) $(WRAPPER)
	$(CXX) $(LDFLAGS) $(filter-out $(SOURCES:.cc=.py),$(OBJECTS)) -o $@
	chmod -x $(TARGET)

debug: $(SOURCES:.cc=.o) $(DEBUG)
	$(CXX) $(LDFLAGS:-shared=) $(SOURCES:.cc=.o) $(DEBUG) -o $@

%_wrap.cxx: %.i
	$(SWIG) $(SWIGFLAGS) $(SOURCES:.cc=.i)

%.py: %_wrap.cxx
	@

%.o: %.cc
	$(CXX) $(CFLAGS) -o $@ $<

%.o: %.cxx
	$(CXX) $(CFLAGS:-Wextra=) $(INCLUDES) -o $@ $<

clean:
	rm -f $(DEBUG:.cc=.o) $(filter-out $(SOURCES:.cc=.py),$(OBJECTS)) $(WRAPPER) $(TARGET) debug

rebuild: clean all

tarball:
	tar -czf $(SOURCES:.cc=.tgz) $(SOURCES) $(SOURCES:.cc=.i) $(SOURCES:.cc=.h) $(DEBUG) Makefile

