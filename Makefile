#define the compiler
CXX=g++ -std=c++20 
CXXFLAG= -fopenmp 
#-lmuparser

#define the executablee
EXEC=main

#define the headers, the objects, the source files to get the executable
SRC=main.cpp problem.cpp
OBJECTS=$(SRC:.cpp=.o)
HEADERS=problem.hpp

# Link the executable
$(EXEC) : $(OBJECTS)
	$(CXX) $(CXXFLAG) $^ -o $@

# Compile the source files into object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAG) -c $< -o $@


# Clean rule to remove object files and the executable
clean:
	rm -f *.o $(EXEC)

# Define rule to optimize the build
optimize: CXXFLAGS += -O3
optimize: $(EXEC)