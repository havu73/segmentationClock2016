#include <cerrno> // Needed for errno, EEXIST
#include <cstdio> // Needed for fopen, fclose, fseek, ftell, rewind
#include <sys/stat.h> // Needed for mkdir
#include <unistd.h> // Needed for read, write, close
#include <stdio.h>
#include <string.h>

#include "io.hpp" // Function declarations
#include "main.hpp"
#include "macros.hpp"

using namespace std;

extern terminal* term; // Declared in init.cpp

/* copy string from value to field.
 * Used in accept_input_params (init.cpp) to store users' input files names into corresponding field in ip
 * Ex: store_filename(&(ip.passed_file_name), value)
 * Params: field: pointers to an array of characters, most likely pointers to a field in ip
 * 		   value: an array of char (from users' arguments)
 */
void store_filename (char** field, const char* value) {
	strcpy(*field, value);
}

/* open_file opens the file with the given name and stores it in the given output file stream
	parameters:
		file_pointer: a pointer to the output file stream to open the file with
		file_name: the path and name of the file to open
		append: if true, the file will appended to, otherwise any existing data will be overwritten
	returns: nothing
	notes:
	todo:
*/
void open_file (ofstream* file_pointer, char* file_name, bool append) {
	try {
		if (append) {
			cout << term->blue << "Opening " << term->reset << file_name << " . . . ";
			file_pointer->open(file_name, fstream::app);
		} else {
			cout << term->blue << "Creating " << term->reset << file_name << " . . . ";
			file_pointer->open(file_name, fstream::out);
		}
	} catch (ofstream::failure) {
		cout << term->red << "Couldn't write to " << file_name << "!" << term->reset << endl;
		exit(EXIT_FILE_WRITE_ERROR);
	}
	term->done();
}
