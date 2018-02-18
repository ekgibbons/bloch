TARGET = bloch_cardiac

INSTALL_PATH = ~/python_utils/simulations

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES   = $(wildcard $(SRCDIR)/*.c)
SOURCES_CPP = $(wildcard $(SRCDIR)/*.cpp)
INCUDES   = $(wildcard $(SRCDIR)/*.h)
INCUDES_HPP = $(wildcard $(SRCDIR)/*.hpp)
OBJECTS_C    = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
OBJECTS_CPP  = $(SOURCES_CPP:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# include base path
BASE_INCLUDE = /home/mirl/egibbons/.conda/envs/ekg/include
PYBIND_INCLUDE = /home/mirl/egibbons/.conda/envs/ekg/include/python3.5m
INCLUDE_FLAGS = -I$(BASE_INCLUDE) -I$(PYBIND_INCLUDE)


# library base path
LIBRARY_PATH = -L/home/mirl/egibbons/.conda/envs/ekg/lib

# libraries to link
LD_LIBS = -larmadillo

CPPFLAGS = -O3 -Wall -Wconversion -fPIC
CXXFLAGS = -std=c++11

# CFLAGS = 
# LDFLAGS = -shared -Wl,--export-dynamic 
# LDFLAGS += -Wl,--unresolved-symbols=report-all
LDFLAGS = -shared

# .PHONY: all
# all: $(TARGET) install # run

WARNINGS=-Wall
PYBIND_INCLUDES=`python3 -m pybind11 --includes`
FLAGS=-O3 -shared -std=c++11 -fPIC
EXTENSION=`python3-config --extension-suffix`

$(TARGET): $(OBJECTS_C) $(OBJECTS_CPP)
	$(CXX) $(OBJECTS_C) $(OBJECTS_CPP) $(LD_LIBS) $(LIBRARY_PATH) $(LDFLAGS) -o $(TARGET).so
	mv $(TARGET).so python/

$(OBJECTS_C): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CPPFLAGS) $(INCLUDE_FLAGS) -c -o $@ $< 

$(OBJECTS_CPP): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

.PHONY: install
install:
	cp python/$(TARGET).so $(INSTALL_PATH)

.PHONY: clean
clean:
	rm -rf *~ *.o *.so
	rm -rf $(OBJECTS_C)
	rm -rf $(OBJECTS_CPP)
	rm -rf $(TARGET).so
	rm -rf python/$(TARGET).so
	rm -rf python/*.pdf python/*.jpg python/*.png

