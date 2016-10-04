#include <unistd.h> // Needed for getpid

#include "init.hpp" // Function declarations
#include "io.hpp"

using namespace std;

terminal* term = NULL; // The global terminal struct

/* init_terminal creates and initializes a new terminal struct
	parameters:
	returns: nothing
	notes:
	todo:
*/
void init_terminal () {
	if (term != NULL) {
		delete term;
	}
	term = new terminal();
}

/* free_terminal frees the terminal from memory and resets the terminal text color to its default value
	parameters:
	returns: nothing
	notes:
	todo:
*/
void free_terminal () {
	cout << term->reset;
	delete term;
}

/* init_verbosity sets the verbose stream to /dev/null if verbose mode is not enabled
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void init_verbosity (input_params& ip) {
	if (!ip.verbose) {
		term->set_verbose_streambuf(ip.null_stream->rdbuf());
	}
}

