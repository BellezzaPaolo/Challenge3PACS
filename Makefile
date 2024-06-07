# Specify which compiler
CXX=mpicxx -std=c++20

# Initialize the directory PACS_ROOT and add all the compiler flags
PACS_ROOT?=../pacs-examples/Examples
CXXFLAG= -fopenmp -pipe -DMPI -I$(PACS_ROOT)/include

# Installation of the muParser
# Usually could be disable since I'm trouble with installation
MUPARSER_LIBDIR=-L$(PACS_ROOT)/lib
LIBRARY_NAME=muparser
DYNAMIC_LIBFILE=lib$(LIBRARY_NAME).so
LDLIBS+=-L. -l$(LIBRARY_NAME) $(MUPARSER_LIBDIR) -lmuparser -L$(PACS_ROOT)/lib -lpacs
CXXFLAGS+=$(MUPARSER_INCLUDE)


EXEC=main
SRC=main.cpp problem.cpp
OBJECTS=$(SRC:.cpp=.o)
HEADERS=problem.hpp

# Make all the object files linking the headers and finally make the executable
$(EXEC) : $(OBJECTS)
	$(CXX) $(CXXFLAG) $^ -o $@ $(LDLIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAG) -c $< -o $@


# Rule to remove all the object file and the executable
clean:
	rm -f *.o $(EXEC)

# Rule for the optimization of the exxecutable
optimize: CXXFLAGS += -O3
optimize: $(EXEC)
