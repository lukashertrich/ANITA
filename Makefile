CXX = g++ 
INC_DIR = ./include
CXXFLAGS = -c -Wall -std=c++11 -I$(INC_DIR)
LDFLAGS = -L ./lib
LDLIBS = 
SOURCES = main.cpp
OBJECTS = $(SOURCES: .cpp = .o)
EXECUTABLE = earthmodel

.PHONY: all

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.p:
	$(CC) $(CFLAGS) $< -o $@