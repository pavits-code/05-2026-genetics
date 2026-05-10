CXX ?= g++
CXXFLAGS ?= -Wall -Werror -pedantic -g --std=c++17 -Wno-sign-compare -Wno-comment
CXXFLAGS += -fsanitize=address

# Run regression test
test: genetics.exe
	./genetics.exe

# Compile the main executable
genetics.exe: genetics.cpp
	$(CXX) $(CXXFLAGS) genetics.cpp -o genetics.exe

# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump