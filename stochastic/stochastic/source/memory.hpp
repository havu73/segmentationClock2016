#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <cstdlib> // Needed for size_t

using namespace std;

void* mallocate(size_t);
void mfree(void*);
#if defined(MEMTRACK)
	void print_heap_usage();
#endif

#endif
