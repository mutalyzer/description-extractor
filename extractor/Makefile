## *******************************************************************
##   (C) Copyright 2013 Leiden Institute of Advanced Computer Science
##   Universiteit Leiden
##   All Rights Reserved
## *******************************************************************
## Extractor (library)
## *******************************************************************
## FILE INFORMATION:
##   File:     Makefile
##   Author:   Jonathan K. Vis
##   Revision: 1.02a
##   Date:     2013/07/24
## *******************************************************************
## DESCRIPTION:
##  Build the Extractor library for python.
## *******************************************************************

SOURCES=extractor.cc
TARGET=_extractor.so
DEBUG=debug.cc

CXX=g++
CFLAGS=-c -fpic -Wall -Wextra -O3
LDFLAGS=-shared -Wall -O3

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

%.py: %.i
	$(SWIG) $(SWIGFLAGS) $(SOURCES:.cc=.i)

%.o: %.cc
	$(CXX) $(CFLAGS) -o $@ $<

%.o: %.cxx
	$(CXX) $(CFLAGS:-Wextra=) $(INCLUDES) -o $@ $<

clean:
	rm -f $(filter-out $(SOURCES:.cc=.py),$(OBJECTS)) $(WRAPPER) $(TARGET) debug

rebuild: clean all

tarball:
	tar -czf $(SOURCES:.cc=.tgz) $(SOURCES) $(SOURCES:.cc=.i) $(SOURCES:.cc=.h) $(DEBUG) Makefile

