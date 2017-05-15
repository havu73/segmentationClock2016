
#include <cerrno> // Needed for errno, EEXIST
#include <cstdio> // Needed for fopen, fclose, fseek, ftell, rewind
#include <sys/stat.h> // Needed for mkdir
#include <unistd.h> // Needed for read, write, close
#include <stdio.h>
#include <string.h>

#include "io.hpp" // Function declarations
#include "main.hpp"
#include "macros.hpp"
#include "memory.hpp"

using namespace std;

extern terminal* term; // Declared in init.cpp
char* mutant_dir_name [] = {"wildtype"};
char* state_file_name [] = {"mHer1", "mHer7", "mDelta", "pHer1", 
	"pHer7", "pDelta", "dimerH1H1", "dimerH1H7", "dimerH7H7", "geneH1", "geneH1N", "geneH1PH11", "geneH7",
	"geneH7N", "geneH7PH11", "geneDelta", "geneDeltaPH11"};
	
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


/* Read parameters from pipe, and store into pr structure. Used when sres sends parameters into simulation. 
 * Params: pr: structure declared in main (main.cpp) to store parameters
 * 		   ip: input params created and modified in main (main.cpp), containing file descriptors to get input (ip.pipe)in)
 * Notes: in sres, the first number we pipe in it the total number of rates, so the first number we check here
 * is also the number of rates, if this number is not NUM_PARAMS (macros.hpp), program will be frozen. 
 * in sres, the second number piped in is the number of param set (1), so the second number we read from pipe is the number of param sets,
 * stored in ip.num_sets
 * Next, we read the parameters.
 */
void read_pipe (parameters& pr, input_params& ip) {
	// Read how many rates per set will be piped in
	int num_pars = 0;
	read_pipe_int(ip.pipe_in, &num_pars);
	if (num_pars != NUM_PARAMS) {
		cout << term->red << "An incorrect number of rates will be piped in! This simulation requires " << NUM_PARAMS << " rates per set but the sampler is sending " << num_pars << " per set." << term->reset << endl;
		exit(EXIT_INPUT_ERROR);
	}

	// Read how many sets will be piped in
	int num_sets = 0;
	read_pipe_int(ip.pipe_in, &num_sets);
	if (num_sets < ip.num_sets) {
		cout << term->red << "The number of num_sets provided by pipe is " << num_sets << ", but you specified " << ip.num_sets << " sets. "<< term->reset << endl;
		exit(EXIT_INPUT_ERROR);
	}
	
	// Read every set and store into pr.data
	for (int i = 0; i < ip.num_sets; i++) {
		if (read(ip.pipe_in, pr.data[i], sizeof(double) * NUM_PARAMS) == -1){
			term->failed_pipe_read();
			exit(EXIT_PIPE_READ_ERROR);
		}
	}
}

/* read int from file descriptor fd and store that number into int pointer address
 * Helper function of read_pipe (io.cpp)
 * params: fd: file descriptor, most likely ip.pipe_in
 * 		
 */
void read_pipe_int (int fd, int* address) {
	if (read(fd, address, sizeof(int)) == -1) {
		term->failed_pipe_read();
		exit(EXIT_PIPE_READ_ERROR);
	}
}

/* Put the file's content into input_data* ifd.buffer
 * Called in read_experiment_data, read_sim_data (init.cpp) to read users' file input
 * Params: input_data* structure declared in main(main.cpp)
 */
void read_file (input_data* ifd) {
	cout << term->blue << "Reading file " << term->reset << ifd->filename << " . . . ";
	// Open the file for reading
	FILE* file = fopen(ifd->filename, "r");
	if (file == NULL) {
		cout << term->red << "Couldn't open " << ifd->filename << "!" << term->reset << endl;
		exit(EXIT_FILE_READ_ERROR);
	}
	
	// Seek to the end of the file, grab its size, and then rewind
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	ifd->size = size;
	rewind(file);
	// Allocate enough memory to contain the whole file
	ifd->buffer = (char*)mallocate(sizeof(char) * size + 1);
	// Copy the file's contents into the buffer
	long result = fread(ifd->buffer, 1, size, file);
	if (result != size) {
		cout << term->red << "Couldn't read from " << ifd->filename << term->reset << endl;
		exit(EXIT_FILE_READ_ERROR);
	}
	ifd->buffer[size] = '\0';
	
	// Close the file
	if (fclose(file) != 0) {
		cout << term->red << "Couldn't close " << ifd->filename << term->reset << endl;
		exit(EXIT_FILE_READ_ERROR);
	}
	term->done();
}

/* parse_param_line reads a line in the given parameter sets buffer and stores it in the given array of doubles
	parameters:
		params: the array of doubles to store the parameters in
		buffer_line: the buffer with the line to read
		index_buffer: the index of the buffer to start from
	returns: true if a line was found, false if the end of the file was reached without finding a valid line
	notes:
		The buffer should contain one parameter set per line, each set containing comma-separated floating point parameters.
		Blank lines and lines starting with # will be ignored.
		Each line must contain the correct number of parameters or the program will exit.
		index_buffer is a reference, allowing this function to store where it finished parsing.
	todo:
*/
bool parse_param_line (parameters& pr, int j, char* buffer_line, int& index_buffer) {
	static const char* usage_message = "There was an error reading the given parameter sets file.";
	int index_params = 0; // Current index in params
	int index_digits = index_buffer; // Index of the start of the digits to read
	int i = index_buffer; // Current index in buffer_line
	int line_start = i; // The start of the line, used to tell whether or not a line is empty
	for (; not_EOL(buffer_line[i]); i++) {
		if (buffer_line[i] == '#') { // Skip any lines starting with #
			for(; not_EOL(buffer_line[i]); i++);
			line_start = i + 1; 
			index_digits = i + 1;
			i++;
		} else if (buffer_line[i] == ',' ) { // Indicates the end of the digits to read
			if (sscanf(buffer_line + index_digits, "%lf", &pr.data[j][index_params++]) < 1) { // Convert the string of digits to a double when storing it in params
				usage(usage_message);
			}
			index_digits = i + 1;
		}
	}
	index_buffer = i + 1;
	if (i - line_start > 0) { // This line has content
		if (sscanf(buffer_line + index_digits, "%lf", &pr.data[j][index_params++]) < 1) {
			usage(usage_message);
		}
		if (index_params != NUM_PARAMS) {
			cout << term->red << "The given parameter sets file contains sets with an incorrect number of rates! This simulation requires " << NUM_PARAMS << " per set but at least one line contains " << index_params << " per set." << term->reset << endl;
			exit(EXIT_INPUT_ERROR);
		}
		return true;
	} else if (buffer_line[index_buffer] != '\0' && buffer_line[index_buffer] != '\n') { // There are more lines to try to parse
		return parse_param_line(pr, j, buffer_line, index_buffer);
	} else { // The end of the buffer was found
		return false;
	}
}

/* parse_ranges_file reads the given buffer and stores every range found in the given ranges array
	parameters:
		ranges: the array of pairs in which to store the lower and upper bounds of each range
		buffer: the buffer with the ranges to read
	returns: nothing
	notes:
		The buffer should contain one range per line, starting the name of the parameter followed by the bracked enclosed lower and then upper bound optionally followed by comments.
		e.g. 'msh1 [30, 65] comment'
		The name of the parameter is so humans can conveniently read the file and has no semantic value to this parser.
		Blank lines and lines starting with # will be ignored. Anything after the upper bound is ignored.
	todo:
*/
void parse_ranges_file (pair <double, double> ranges[], char* buffer) {
	int i = 0;
	int rate = 0;
	for (; buffer[i] != '\0'; i++) {
		// Ignore lines starting with #
		while (buffer[i] == '#') {
			while (buffer[i] != '\n' && buffer[i] != '\0') {i++;}
			i++;			
		}
		
		// Ignore whitespace before the opening bracket
		while (buffer[i] != '[' && buffer[i] != '\0') {i++;}
		if (buffer[i] == '\0') {break;}
		i++;
		
		// Read the bounds
		ranges[rate].first = atof(buffer + i);
		while (buffer[i] != ',') {i++;}
		i++;
		ranges[rate].second = atof(buffer + i);
		if (ranges[rate].first < 0 || ranges[rate].second < 0) { // If the ranges are invalid then set them to 0
			ranges[rate].first = 0;
			ranges[rate].second = 0;
		}
		rate++;
		// Skip any comments until the end of the line
		while (buffer[i] != '\n' && buffer[i] != '\0') {i++;}
	}
}

/* not_EOL returns whether or not a given character is the end of a line or file (i.e. '\n' or '\0', respectively)
	parameters:
		c: the character to check
	returns: true if c is the end of a line or file, false otherwise
	notes:
		When reading input file strings, use this instead of a straight newline check to avoid EOF (end of file) issues.
	todo:
*/
bool not_EOL (char c) {
	return (c != '\n' && c != '\0' && c != EOF);
}

/*
 * Create directory set_setIndex_scnORfib inside simulation. this directory contains concentrations files. 
 * Called in simulate_all_params(sim.cpp)
 * params:
 * 		set_index: index of the current param set (found in simulate_all_params)
 * 		ip
 */
void create_set_directory (int set_index, input_params& ip){
	char * dir_name = new char[40];
	
	cout << term->blue << "Creating set_" << set_index << " directory if necessary . . . " << term->reset << endl;
	if (ip.has_out_dir){
		sprintf(dir_name, "%s/set_%d", ip.out_dir, set_index);
	}
	else{
		sprintf(dir_name, "set_%d", set_index);
	}
	if (mkdir(dir_name, 0775) != 0 && errno != EEXIST){
		cout << term->red << "Couldn't create '" << dir_name << "' directory!" << term->reset << endl;
		delete [] dir_name;
		exit(EXIT_FILE_WRITE_ERROR);
	}
	
	delete [] dir_name;
}


void print_concentrations(input_params& ip, int set_index, int mutant_index, con_levels& cons){
	char* file_name = new char [50];
	
	// create directory name
	if (ip.has_out_dir){
		sprintf(file_name, "%s/set_%d/%s.txt", ip.out_dir, set_index, mutant_dir_name[mutant_index]);
	}
	else{
		sprintf(file_name, "set_%d/%s.txt", set_index, mutant_dir_name[mutant_index]);
	}
	
	ofstream file_cons;
	open_file(&file_cons, file_name, true);
	int num_record = cons.num_sim_steps / ip.record_gran;
	for (int i = 0 ; i < num_record; i++){
		double time = i * ip.record_gran * ip.step_size;
		file_cons << time;
		for (int j = 0; j < NUM_STATES; j++){
			int index = i * ip.record_gran;
			double con = cons.sim_data[j][index];
			file_cons << "," << con;
		}
		file_cons << "\n" ;
	}
	delete [] file_name;
}

void write_pipe(double* score, input_params& ip){
	if (write(ip.pipe_out, score, sizeof(double)) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
	// Close the pipe
	if (close(ip.pipe_out) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
}
