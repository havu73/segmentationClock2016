#include <cerrno> // Needed for errno, EEXIST
#include <cstdio> // Needed for fopen, fclose, fseek, ftell, rewind
#include <sys/stat.h> // Needed for mkdir
#include <unistd.h> // Needed for read, write, close

#include "io.hpp" // Function declarations
#include "sim.hpp" // Needed for anterior_time
#include "structs.hpp"
#include "main.hpp"

using namespace std;

extern terminal* term; // Declared in init.cpp

