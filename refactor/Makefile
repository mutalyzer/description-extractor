CC=g++
CFLAGS=-c -Wall -Wextra -O3 -D__debug__
LDFLAGS=-Wall -Wextra -O3
SOURCES=debug.cc extractor.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=debug

.PHONY: all clean rebuild

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

rebuild: clean all

