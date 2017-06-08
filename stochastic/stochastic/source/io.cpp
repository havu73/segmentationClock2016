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

void open_file_quiet (ofstream* file_pointer, char* file_name){
	try{
		file_pointer->open(file_name, fstream::app);
	}catch (ofstream::failure){
		cout << term->red << "Couldn't write to " << file_name << "!" << term->reset << endl;
		exit(EXIT_FILE_WRITE_ERROR);
	}
}

/* Read parameters from pipe, and store into pr structure. Used when sres sends parameters into simulation. 
 * Params: pr: structure declared in main (main.cpp) to store parameters
 * 		   ip: input params created and modified in main (main.cpp), containing file descriptors to get input (ip.pipe)in)
 * Notes: in sres, the first number we pipe in it the total number of rates, so the first number we check here
 * is also the number of rates, if this number is not NUM_RATES (macros.hpp), program will be frozen. 
 * in sres, the second number piped in is the number of param set (1), so the second number we read from pipe is the number of param sets,
 * stored in ip.num_sets
 * Next, we read the parameters.
 */
void read_pipe (parameters& pr, input_params& ip) {
	// Read how many rates per set will be piped in
	int num_pars = 0;
	read_pipe_int(ip.pipe_in, &num_pars);
	if (num_pars != NUM_RATES) {
		cout << term->red << "An incorrect number of rates will be piped in! This simulation requires " << NUM_RATES << " rates per set but the sampler is sending " << num_pars << " per set." << term->reset << endl;
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
		if (read(ip.pipe_in, pr.data[i], sizeof(double) * NUM_RATES) == -1){
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
		if (index_params != NUM_RATES) {
			cout << term->red << "The given parameter sets file contains sets with an incorrect number of rates! This simulation requires " << NUM_RATES << " per set but at least one line contains " << index_params << " per set." << term->reset << endl;
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

void create_mutant_directory(int set_index, int mutant_index, input_params& ip){
	char* dir_name = new char [50];
	
	cout << term->blue << "Creating set_" << set_index << "/" << mutant_dir_name[mutant_index] << " directory if necessary..." << term->reset << endl;
	if (ip.has_out_dir){
		sprintf(dir_name, "%s/set_%d/%s", ip.out_dir, set_index, mutant_dir_name[mutant_index]);	
	}
	else{
		sprintf(dir_name, "set_%d/%s", set_index, mutant_dir_name[mutant_index]);
	}
	if (mkdir(dir_name, 0775) != 0 && errno != EEXIST){
		cout << term->red << "Couldn't create '" << dir_name << "' directory!" << term->reset << endl;
		delete [] dir_name;
		exit(EXIT_FILE_WRITE_ERROR);
	}
	delete [] dir_name;
}

void print_debug_info(input_params& ip, embryo& em, parameters& pr){
	
}

void print_concentrations(input_params& ip, int set_index, int mutant_index, embryo& em){
	// allocate memory for directory and file names
	int num_files = NUM_KEEP_STATES + ip.num_print_states;
	char** filename = new char* [num_files];
	for (int i = 0; i < num_files; i++){
		filename[i] = new char [100];
	}
	char* dir_name = new char [50];
	
	// create directory name
	if (ip.has_out_dir){
		sprintf(dir_name, "%s/set_%d/%s", ip.out_dir, set_index, mutant_dir_name[mutant_index]);
	}
	else{
		sprintf(dir_name, "set_%d/%s", set_index, mutant_dir_name[mutant_index]);
	}
	
	// create file names
	create_concentrations_file_name(num_files, ip, dir_name , filename);
	
	// print out onto files
	for (int i = 0; i < num_files; i ++){
		print_one_state_concentrations(em, i, filename);
	}
	
	// delete directory and file names
	delete [] dir_name;
	for (int i = 0; i < num_files; i ++){
		delete [] filename[i];
	}
	delete [] filename;
}

void create_concentrations_file_name(int num_files, input_params& ip, char* dir_name, char** filename){
	sprintf(filename[KEEPMH1], "%s/%s.txt", dir_name, state_file_name[MH1]);
	sprintf(filename[KEEPMH7], "%s/%s.txt", dir_name, state_file_name[MH7]);
	int state_index;
	for (int i = NUM_KEEP_STATES; i < num_files; i ++){
		state_index = ip.print_states[i - NUM_KEEP_STATES];
		sprintf(filename[i], "%s/%s.txt", dir_name, state_file_name[state_index]);
	}
}

void print_one_state_concentrations(embryo& em, int state_index, char** filename){
	ofstream file_cons;
	open_file(&file_cons, filename[state_index], true);
	vector<int> * cons;
	int num_steps = (em.time_record)->size();
	for (int i = 0; i < num_steps; i++){
		file_cons << (em.time_record)->at(i) << ","; 
	}
	file_cons << "\n";
	for (int i = 0; i < em.num_cells; i++){
		cons = (em.cell_list[i]->cons_record)[state_index];
		for (int j = 0; j < num_steps; j++){
			file_cons << cons->at(j) << ",";
		}
		file_cons << "\n";
	}
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

void print_smooth_data (input_params& ip, embryo& em, peak_trough& pt, int cons_index, int set_index){
	// create smooth_file name
	char* filename = new char [100];
	if (ip.has_out_dir){
		sprintf(filename, "%s/set_%d/%s/%s_smoothed.txt", ip.out_dir, set_index, mutant_dir_name[WT],state_file_name[cons_index]);
	}
	else{
		sprintf(filename, "set_%d/%s/%s_smoothed.txt", set_index, mutant_dir_name[WT], state_file_name[cons_index]);
	}
	ofstream file_smooth;
	open_file(&file_smooth, filename, true);
	int num_steps = em.time_record->size() - (WINDOW_SIZE * 2) - 1;
	// write the smoothed data
	// write time first
	for (int i = 0; i < num_steps; i++){
		file_smooth << (em.time_record)->at(i + WINDOW_SIZE) << ",";
	}
	file_smooth << (em.time_record)->at(num_steps + WINDOW_SIZE) << "\n";
	// write smoothed concentrations of cells
	for (int i = 0; i < em.num_cells; i ++){
		double * smooth_record = pt.smooth_cons[i];
		for (int j = 0; j < num_steps; j++){
			file_smooth << smooth_record[j] << ",";
		}
		file_smooth << smooth_record[num_steps] << "\n";
	}
	
	delete [] filename;
	
	// create peak file 
	char* peaksfile = new char [100];
	if (ip.has_out_dir){
		sprintf(peaksfile, "%s/set_%d/%s/%s_peaks.txt", ip.out_dir, set_index, mutant_dir_name[WT],state_file_name[cons_index]);
	}
	else{
		sprintf(peaksfile, "set_%d/%s/%s_peaks.txt", set_index, mutant_dir_name[WT], state_file_name[cons_index]);
	}
	ofstream peaks;
	open_file(&peaks, peaksfile, true);
	for (int i = 0; i < pt.num_cells; i ++){
		vector<int> * cell_peaks = pt.peaks[i];
		int num_peaks = cell_peaks->size() - 1;
		for (int j = 0; j < (num_peaks); j++){
			peaks << cell_peaks->at(j) << ",";
		}
		peaks << cell_peaks->at(num_peaks) << "\n";
	}
	delete [] peaksfile;
	
	// create trough file
	char* troughsfile = new char [100];
	if (ip.has_out_dir){
		sprintf(peaksfile, "%s/set_%d/%s/%s_troughs.txt", ip.out_dir, set_index, mutant_dir_name[WT],state_file_name[cons_index]);
	}
	else{
		sprintf(peaksfile, "set_%d/%s/%s_troughs.txt", set_index, mutant_dir_name[WT], state_file_name[cons_index]);
	}
	ofstream troughs;
	open_file(&troughs, troughsfile, true);
	for (int i = 0; i < pt.num_cells; i ++){
		vector<int> * cell_troughs = pt.troughs[i];
		int num_troughs = cell_troughs->size() - 1;
		for (int j = 0; j < (num_troughs); j++){
			troughs << cell_troughs->at(j) << ",";
		}
		troughs << cell_troughs->at(num_troughs) << "\n";
	}
	delete [] troughsfile;
}

void print_features_data (input_params& ip, embryo& em, features& fts, int set_index, int mutant_index){
	char* filename = new char [100];
	if (ip.has_out_dir){
		sprintf(filename, "%s/set_%d/%s/features.txt", ip.out_dir, set_index, mutant_dir_name[mutant_index]);
	}
	else{
		sprintf(filename, "set_%d/%s/features.txt", set_index, mutant_dir_name[mutant_index]);
	}
	ofstream features_file;
	open_file(&features_file, filename, true);
	features_file.precision(32);
	// num_cells, num_bin
	features_file << fts.num_cells << "," << fts.num_bin << "\n";
	// print amplitude
	for (int i = 0; i < NUM_KEEP_STATES; i ++){
		for (int j = 0; j < (fts.num_cells - 1); j++){
			features_file << fts.avg_amplitude[i][j] << ",";
		}
		features_file << fts.avg_amplitude[i][fts.num_cells - 1] << "\n";
	}
	// print period
	for (int i = 0; i < NUM_KEEP_STATES; i++){
		for (int j = 0; j < (fts.num_cells - 1); j++){
			features_file << fts.avg_period[i][j] << ",";
		}
		features_file << fts.avg_period[i][fts.num_cells - 1] << "\n";
	}
	// Binned concentrations
	for (int i = 0 ; i < (fts.num_bin - 1); i++){
		features_file << fts.avg_cons[i] << ",";
	}
	features_file << fts.avg_cons[fts.num_bin - 1] << "\n";
	
	// binned intrinsic noise
	for (int i = 0; i < (fts.num_bin - 1); i++){
		features_file << fts.intrinsic[i] << ",";
	}
	features_file << fts.intrinsic[fts.num_bin - 1] << "\n";
	
	// binned extrinsic noise
	for (int i = 0; i < (fts.num_bin - 1); i++){
		features_file << fts.extrinsic[i] << ",";
	}
	features_file << fts.extrinsic[fts.num_bin - 1] << "\n";
	delete [] filename;
}

void print_perturb_rates(input_params& ip, rates& rs, parameters& pr, int set_index){
	char* filename = new char [100];
	if (ip.has_out_dir){
		sprintf(filename, "%s/set_%d/rates.txt", ip.out_dir, set_index);
	}
	else{
		sprintf(filename, "set_%d/rates.txt", set_index);
	}
	ofstream file_rates;
	open_file(&file_rates, filename, true);
	// first print num cells, num rates
	file_rates << rs.num_cells << "," << NUM_NETWORK_RATES << endl;
	// second print the base rates
	file_rates << pr.data[set_index][0];
	for (int i = 1; i < NUM_RATES; i ++){
		file_rates << "," << pr.data[set_index][i];
	}
	file_rates << endl;
	// third print the perturbation rates, each cell has its own line
	for (int i = 0; i < rs.num_cells; i ++){
		file_rates << rs.perturb_rates[i][0];
		for (int j = 1; j < NUM_NETWORK_RATES; j ++){
			file_rates << "," << rs.perturb_rates[i][j];
		}
		file_rates << endl;
	}
	// fourth print the actual rates, each cell has its own line
	for (int i = 0; i < rs.num_cells; i ++){
		file_rates << rs.data[i][0];
		for (int j = 1; j < NUM_NETWORK_RATES; j ++){
			file_rates << "," << rs.data[i][j];
		}
		file_rates << endl;
	}
}
