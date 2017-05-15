Zebrafish somitogenesis stochastic simulation
==================================

Table of contents
-----------------

| 0: Compatibility and system requirements  
|__ 0.0: Tested systems and portability  
|__ 0.1: MPI dependency  
|__ 0.2: SCons dependency  
|__ 0.3: Python and SciPy dependencies  
|__ 0.4: FFmpeg dependency  
| 1: Compilation  
|__ 1.0: Compiling with and without SCons  
|__ 1.1: Compilation options  
|__ 1.2: Compiling for MPI  
|__ 1.3: Random number generation discrepancies  
| 2: Running simulations  
|__ 2.0: Biological and computational description of the simulation  
|__ 2.1: Setting up a simulation  
|__ 2.2: Input and output  
|____ 2.2.0: Methods of input  
|____ 2.2.1: Parameter sets  
|____ 2.2.2: Perturbations and gradients  
|____ 2.2.3: Command-line arguments  
|____ 2.2.4: Input file formats  
|______ 2.2.4.0: Parameter sets format  
|______ 2.2.4.1: Ranges format  
|____ 2.2.5: Output file formats  
|______ 2.2.5.1: Concentrations text format  
|______ 2.2.5.2: Oscillation features format  
|____ 2.2.6: Piping in parameter sets from other applications  
|____ 2.2.7: Generating random parameter sets  
| 3: Finding parameter sets with SRES  
|__ 3.0: How SRES finds parameters  
|____ 3.0.0: How evolutionary algorithms work  
|____ 3.0.1: How SRES works  
|____ 3.0.2: How SRES-gradients works  
|__ 3.1: Setting up SRES  
|____ 3.1.0: Searching for parameters  
|____ 3.1.1: Searching for gradients  
|____ 3.1.2: Using MPI  
|__ 3.2: Input and output  
|____ 3.2.0: Methods of input  
|____ 3.2.1: Ranges  
|____ 3.2.2: Simulation calls  
|____ 3.2.3: Command-line arguments  
|____ 3.2.4: Input file formats  
|____ 3.2.5: Output file formats  
|__ 3.3: Modifying the code  
|____ 3.3.0: Adding command-line arguments  
|____ 3.3.1: Changing which arguments SRES sends to simulations  
|____ 3.3.2: Adding input and output files  
|____ 3.3.3: Modifying libSRES  
|____ 3.3.4: Other modifications
|_ 4: Refine ranges 
|__ 4.0: Overview
|__ 4.1: How to run the program
|___ 4.1.0: Command-line arguments
|___ 4.1.1: Example program calls
|_ 5: Authorship and licensing  

0: Compatibility and system requirements
----------------------------------------

**0.0: Tested systems and portability**

This package has been designed for Unix based machines. It has been tested on Fedora 18 with GCC version 4.7.2 and OSX 10.7 (Lion) with GCC 4.2.1 using 64-bit processors. Because this package's only external dependencies are the C++98 standard libraries, Python, and SciPy, most of its functionality is platform-independent. However, various file I/O functions rely on Unix system commands making full compatibility with Windows and other non-Unix operating systems unlikely.

***********************
**0.1: MPI dependency**

Due to the computationally intensive nature of the simulations, we recommend either a very powerful, multicore processor or a cluster of many processors so that multiple simulations can be run in parallel. To achieve paralellization, most applications in this package that run multiple simulations are parallelized using MPI (Message Passing Interface), which can be downloaded at http://www.open-mpi.org/software/ompi/v1.6/. While all applications have serial versions for testing purposes, replicating our results serially would take several weeks of computation time, making the installation of MPI essentially a requirement.

*************************
**0.2: SCons dependency**

The software constructor tool SCons 2.3.0 is used to simplify compilation. SCons can be downloaded at http://www.scons.org/download.php but this tool is not necessary and instructions for direct compilation are provided in Section 1.

**************************************
**0.3: Python and SciPy dependencies**

Scripts for creating figures and parsing data files require Python 2.7.3 and SciPy, which can be downloaded at http://www.python.org/download/ and http://www.scipy.org/install.html, respectively. These scripts are independent from the simulation code but are useful for parsing and visualizing results.

**************************
**0.4: FFmpeg dependency**

Movies can be made from tissue snapshots using FFmpeg 1.0.7, which can be downloaded at http://www.ffmpeg.org/download.html. If you have no need to create movies, FFmpeg's intallation can be skipped.

1: Compilation
--------------

**1.0: Compiling with and without SCons**

To compile an application in its default configuration, open a terminal window and navigate to the package's root directory. If SCons is installed on the machine, simply enter 'scons' to compile the source. If SCons cannot be installed on the machine, each application can be compiled manually by entering its associated g++ compilation statement:
* simulation: 'g++ -O2 -Wall -o simulation main.cpp init.cpp sim.cpp feats.cpp io.cpp memory.cpp debug.cpp'
* sres: 'g++ -O2 -Wall -o sres main.cpp init.cpp sres.cpp io.cpp memory.cpp'
* sensitivity: 'g++ -O2 -Wall -o sensitivity source/analysis.cpp source/init.cpp source/io.cpp source/memory.cpp finite-difference/finite-difference.cpp'

If g++ is not install on the machine, you need to install it or an equivalent compiler.

****************************
**1.1: Compilation options**

All applications come with at least three compilation options, 'profile', 'debug', and 'memtrack'. By entering 'scons profile=1', 'scons debug=1', or 'scons memtrack=1', the application is compiled with compile and link flags designed for profiling, debugging, and memory tracking, respectively. Profiling adds the '-pg' compile and link flags, which adds extra code that enables gprof profiling analysis. Debugging adds the '-g' compile flag, which adds extra code that enables GDB debugging. Memory tracking adds the '-D MEMTRACK' compile flag, which adds a custom macro indicating the program should track its heap memory allocation. For more information on these options, see Section 6.

**************************
**1.2: Compiling for MPI**

When MPI-compatible applications are compiled with 'scons mpi=1' then the mpic++ compiler is used instead of g++ and a custom macro indicating the program should parallelize with MPI is added via '-D MPI'.

***********************************************
**1.3: Random number generation discrepancies**

As of this publication, Linux and Mac use different versions of GCC whose standard random number generators happen to produce different random numbers even when given the same seed. Because the modeled system is robust, parameter sets should receive similar scores regardless of the initial seed but it is worth noting that an identically configured simulation may produce different results on different operating systems when random number generation is incorporated via perturbations.

2: Running simulations
----------------------

**2.0: Biological and computational description of the simulation**

*******************************
**2.0.0: The biological model**

This simulation is a one or two dimensional, gene-and-cell stochastic simulation of the presomitic mesoderm (PSM) in zebrafish embryos. It includes the genes Her1, Her7, and Delta. Delta performs all of the functionality of the DeltaC-Notch intercellular signaling pathway. Her1 and Her7 homodimerize and heterodimerize, with Her1Her1 binding to DNA regulatory region to repress _her1_, _her7_ and _delta_ mRNA transcription. The simulation uses discrete reactions to model mRNA and protein synthesis and degradation as well as dimer association, dissociation, and degradation. 
********************************
**2.0.1: Possible tissue sizes**

The cell tissue can be either two cells, a one dimensional chain of three or more cells, or a two dimensional hexagonal grid of 16 cells, pr a squre of 4 cells

**************************************
**2.0.2: Rates**

The simulation is controlled by a set of parameters, each representing a biological rate such as _her1_ mRNA synthesis. 

******************
**2.0.3: Mutants and scores**

Right now, parameter sets are tested with wild type conditions. The maximum score that a parameter set can get is 6 scores. The score returned to SRES is MAX_SCORE (specified in ./stochastic/source/macros.hpp) - model score (calculated by simulate_one_param_set in sim.cpp).

********************************
**2.1: Setting up a simulation**

After compiling (see Section 1), there should be an executable named "stochasti" in the stochastic directory. The simulation has many configurations but most have default values and do not necessarily have to be specified. To run a simulation, navigate to the stochastic directory in the terminal and enter "./stochastic ..." where "..." contains all of the command-line options. Once the simulation starts, it will not ask for further input and will print out various milestones when it is finished. 

*************************
**2.2: Input and output**

***************************
**2.2.0: Methods of input**

The only interaction with a simulation is through input you give it before it starts and output it gives you after it ends. There are four aspects of input:
* Parameter sets
* Command-line arguments

***********************************
**2.2.1: Parameter sets**

Parameter sets are inputted in one of three ways: an input file containing each set, the file descriptor of a pipe that will send them in, or an input file containing valid ranges for each parameter and the number of random sets to generate based off of these ranges. The file formats for parameter sets and ranges files are covered in Section 2.2.4 and piping is covered in Section 2.2.6.

*********************************
**2.2.2: Command-line arguments**

Each command-line argument can be entered in a short or long form. Short forms begin with a single dash (-) followed by a single letter. Long forms begin with a double dash (--) followed by a word or phrase. Lower-case letters are considered distinct from their upper-case equivalents. The short and long forms are equivalent - short forms are for convenience and long forms are for clarity. The following is a comprehensive list of all command-line options and their uses:

-i,  --input-params-file  	[filename]   : the relative filename of the parameter sets input file, default=none
-r,  --ranges-file        	[filename]   : the relative filename of the parameter ranges input file, default=none
-ps, --print-states	   		[string]	 : the space separated list of state indices that users want to keep track of over time. MH1, MH7, MD are default. Ex: -ps "4 5 6"
-pc, --print-cons         	[N/A]        : print concentration values to the specified output directory, default=false
-pf, --print-features		[N/A]		 : print noise, period and amplitude features of the simulation, default=false
-pd, --print-debug			[filename]	 : print all simulation information to debug
-od, --output-directory   	[directory]  : the relative directory where output of the program is saved, default=none
-ns, --num-parameter-sets 	[int]        : the number of parameters for which to simulate the model, min=1, default=1
-st, --total-time         	[int]        : the number of minutes to simulate before ending, min=1, default=600
-ss, --sim-seed           	[int]        : the seed to generate random numbers for simulation, min=1, default=generated from the time and process ID
-m,  --mutants            	[int]        : the space-separated list of mutants indices to run for each parameter set, default= list of all possible mutants
-pi, --pipe-in            	[file desc.] : the file descriptor to pipe data from (usually passed by the sampler), default=none
-po, --pipe-out           	[file desc.] : the file descriptor to pipe data into (usually passed by the sampler), default=none
-v,  --verbose            	[N/A]        : print detailed messages about the program and simulation state, default=unused
-q,  --quiet              	[N/A]        : hide the terminal output, default=unused
-nc, --num-cells          	[N/A]        : number of cells in the sytem to simulate. Default = 16
-cdg,--check-done-granularity[int] 	 : number of reactions fired in between two times that the program check whether it is done or not. Default = 120
-rg, --record-granularity	[int]		 : number of reactions fired in between two times that we record the concentrations of states we care about. Default = 30
We require that check-done-granularity and record-granularity are both at least 1, and that check-done-granularity is invisible by record-granularity
-l,  --licensing          [N/A]        : view licensing information (no simulations will be run)
-h,  --help               [N/A]        : view usage information (i.e. this)

Example: ./simulation -i parameters.csv --parameters 10 -m 2000 --no-color

*****************************
**2.2.4: Input file formats**

**********************************
**2.2.4.0: Parameter sets format**

A parameter sets file consists of a list of parameter sets. Each parameter set is placed on its own line, with each parameter separated by a comma. There is not a comma after the final parameter set. Each parameter is included and has a floating point or integral value. Blank lines and lines beginning with "#" are ignored. There must be at least one parameter set per file. There is no limit to the number of parameter sets allowed in a single file, although because every parameter set must be loaded into memory at once (to calculate the maximum delay size used to determine the size of various concentration levels structs ahead of time to avoid repeated reallocation), placing millions of parameter sets in a single file may exhaust your computer's RAM.

The following three lines represent an example file:
```
# This is a comment
63.049748,34.955167,0,42.993889,0,44.600336,0.429514,0.453927,0,0.279120,0,0.335253,45.227915,30.430273,0,46.391646,0,22.217098,0.157011,0.338070,0,0.307667,0,0.274870,0.019982,0.000370,0,0.026413,0,0.008877,0,0.027406,0,0,0,0,0.005284,0,0,0.219352,0.032981,0,0.051730,0,0.175033,0,0.116019,0,0,0,0,0.006626,0,0,0.278811,0.269840,0,0.279624,0,0.285013,0,0.238326,0,0,0,0,0.242239,0,0,11.017495,9.139970,0,0.000000,0,7.329024,0.426847,1.210969,0,0.596138,0,12.461520,361.441083,485.384224,0,0,0,0,499.577384
63.000647,33.757737,0,43.849951,0,47.332097,0.434420,0.455262,0,0.274844,0,0.346678,43.338772,30.019011,0,54.822609,0,25.281511,0.141475,0.315663,0,0.345098,0,0.269280,0.018546,0.003612,0,0.028153,0,0.008334,0,0.025200,0,0,0,0,0.011394,0,0,0.170959,0.041615,0,0.044836,0,0.237797,0,0.248760,0,0,0,0,0.017808,0,0,0.280718,0.310334,0,0.343655,0,0.210100,0,0.233876,0,0,0,0,0.214772,0,0,10.916983,9.443056,0,0.000000,0,7.742257,0.445980,1.035695,0,0.578762,0,12.446215,231.836670,477.034572,0,0,0,0,540.815524
```

**************************
**2.2.4.1: Ranges format**

A ranges files consists of a list of ranges for every parameter in the simulation. Each line contains a string identifying the parameter (used by humans only) and, after whitespace, a (lower bound, upper bound) pair that defines the range. Each pair is enclosed in brackets, [ and ], and the values in each pair are separated by a comma. Both values are floating point or integral and also nonnegative since no parameter represents a value that should be negative. The ranges are read sequentially, with the first range read applied to the first parameter and so on. The number of ranges matches the number of parameters exactly. Blank lines and lines beginning with "#" are ignored. Text after the closing bracket is ignored and can thus be used for further comments.

The following lines represent an excerpt from an example file, where "..." represents missing lines:
```
msh1 [30.0,65.0]
msh7 [30.0,65.0]
msh13 [30.0,65.0]

...

critph1h1 [150.0,900.0]
critph7h13 [150.0,900.0]
critpdelta [200.0,800.0]
```

******************************
**2.2.5: Output file formats**

*******************************
**2.2.5.0: Passed sets format**

A passed sets file consists of a list of parameter sets. Its format is identical to a parameter sets input file (see Section 2.2.4.0) although it never contains blank lines or comments (lines beginning with "#"). Only parameter sets that receive a perfect score are added to the passed sets file.

****************************************
**2.2.5.1: Concentrations text format**

Concentrations are stored as text unless. Text files are easily readable by humans. Each time '-pc' is specified in the command line before running stochastic, the following files are generated in ./<User_specified output directory if necessary>/set_<set_index>/<mutantName>/:
- mHer1_smoothed.txt
- mHer1_peaks.txt
- mHer1_troughs.txt
- mHer7_smoothed.txt
- mHer7_peaks.txt
- mHer7_troughs.txt
- features.txt
- mHer1.txt
- mHer7.txt

Users can specify where they want to store the model's ouput using flags: "-od" or "--output-directory"

*****************************************
**2.2.5.3: Oscillation features format**

An oscillation features file consists of a list of numbers representing various oscillatory features of each parameter set simulated. This data is stored in ./<User_specified output directory if necessary>/set_<set_index>/<mutantName>/features.txt. The following lines represent the format of such file.
num_cells, num_bins
amplitude_her1_mRNA cell 0, amplitude_her1_mRNA cell 1 ,....
amplitude_her7_mRNA cell 0, amplitude_her7_mRNA cell 1, ....
period_her1_mRNA cell0, period_her1_mRNA cell 1, ...
period_her1_mRNA cell 0 , period_her1_mRNA cell 1, ...
avg_her_concentrations bin 0, avg_her_concentrations bin1, avg_her_concentrations bin2, ..
avg_intrinsic_noise_concentrations bin 0, avg_intrinsic_noise_concentrations bin1, avg_intrinsic_noise_concentrations bin2, ...
avg_extrinsic_noise_concentrations bin 0, avg_extrinsic_noise_concentrations bin1, avg_extrinsic_noise_concentrations bin2, ...


***********************************************************
**2.2.6: Piping in parameter sets from other applications**

SRES parameter search program communicates with stochastic simulation through piping, the program has the ability to read parameter sets from one end of a pipe ("-pi") and write its resulting score to the other end ("-po"). This is accomplished by passing to the simulation the file descriptors identifying the two pipe ends via the command-line with -I or --pipe-in for the end to read from and -O or --pipe-out for the end to write to.

To send parameters to the program, write to the input pipe the number of parameters per set (as an int), the number of sets (as an int), and then each parameter of each set as a sequence of doubles. Sending the number of parameters per set is important because the program needs to know how many doubles to read in a row before moving on to the next set. After reading in the number of sets specified, the program stops reading from the pipe and begins simulating. After it finishes every given parameter set, it writes to the output pipe the maximum possible score a set can achieve (as an int) and then the score of each set as a sequence of ints. After writing all of the scores, the program closes the pipe.

*******************************************
**2.2.7: Generating random parameter sets**

Parameter sets can be generated by the program using a ranges file (via the command-line with -r or --ranges-file) and run on the fly. For details on how to construct a ranges file, see Section 2.2.4.3. To record any generated parameter sets, enter a filename to store them via the command-line with -P or --print-sets.

***************************

3: Finding parameter sets with SRES
-----------------------------------

**3.0: How Evolutionary Strategy finds parameters**

Evolution Strategy is used to find biologically realistic parameters by repeatedly running simulations with different parameter sets or configurations. While the theory behind evolutionary algorithms and the details of the library used are beyond the scope of this README, a brief explanation of them follows with the hope that it is enough to reproduce our results.

*******************************************
**3.0.0: How evolutionary algorithms work**

Evolutionary algorithms "evolve" parameters to perform a particular function better. A population is created, each member of which gets its own random set of parameters and performs a particular function that returns a fitness score. Members with the best fitness are allowed to continue to the next generation, at which point they are mutated (a parameter value is randomly modified by a small amount) and crossed over (two members create a new member with a combination of both parents' parameters). The new members perform the same function with their new parameters and new fitness scores are returned. This continues for a specified number of generations, hopefully evolving better sets of parameters.

*************************
**3.0.1: How SRES works**

SRES (Stochastically Ranked Evolution Strategy) uses the libSRES library, which can be downloaded at http://rich.yunda.org/uga/science/libSRES/. libSRES uses a stochastically ranked evolution strategy that, unlike the basic evolutionary strategy, stochastically ranks members when deciding which ones to propogate to the next generation. The full explanation can be found in the paper at https://notendur.hi.is/~tpr/software/sres/Tec311r.pdf, but in short, this algorithm allows there to be a chance that a worse performing member continues to the next generation instead of a better performing one. This prevents the algorithm from getting stuck in local maxima. If only the best members are used then once they find a range of parameter sets that perform better than all sets around them in the problem space, the algorithm has no incentive to branch out to find better maxima because every direction looks worse than its current location.

In our case, members of the population are simulations and the parameter sets picked, mutated, etc. are the usual parameter sets used in the simulation. SRES pipes in a parameter set to the simulation and uses the score it receives as its fitness score, replicating this process dozens of times for hundreds of generations. By evolving the parameter sets, simulations receive increasingly better scores on condition tests until they pass every test, at which point there can obviously be no more improvement.

***********************************
**3.0.2: How SRES-gradients works**

SRES-gradients is a variant on SRES that evolves gradients instead of parameters. When given start point, end point, and gradient factor ranges, the program searches within the ranges using a specified parameter set and tries to find the best gradient configuration. Multiple gradients can be given concurrently, although they all receive the same start point, end point, and factor. The program is intended to be used once a suitable parameter set has been found and can aid in finding the gradients required to meet further conditions that gradient-less simulations cannot achieve.

************************
**3.1: Setting up SRES**

After compiling (see Section 1), there should be an executable named "sres" in the sres directory. SRES has many configurations but most have default values and do not necessarily have to be specified. The only parameter that must be specified is the file to read ranges from (-r or --ranges-file). To begin evolving parameter sets, navigate to the sres directory in the terminal and enter "./sres ..." where "..." contains all of the command-line options. Once the process starts, it will not ask for further input and will print out various milestones and when it is finished. Each generation, the best fitness score and its associated parameter set are printed to the terminal. If the program cannot interpret a part of your input, it will let you know with as much specific information as possible; it will not silently fail.

***********************************
**3.1.0: Searching for parameters**

To use SRES to search for parameters, create or use an existing ranges file to select from. This ranges file determines how low and high each parameter may be set when searching for the best values. Call sres with this ranges file (via the command-line with -r or --ranges-file) and any other configurations desired. To send arguments to the simulation, enter -a or --arguments and then the arguments to pass. For example, "./sres -r input.ranges -a -x 2 -w 2 -y 1". Note that SRES automatically adds the required piping arguments, which take precedence over any given input file. We recommend running simulations in quiet mode (-q or --quiet) to avoid flooding the terminal with unnecessary output. The default arguments for SRES should be fine for most purposes but reducing the population size or number of generations is recommended if you're just testing out the functionality of the program since each generation takes a long time with a sufficiently sized population.

After each generation finishes, the best fitness score and its associated parameter set are printed to the terminal. The program continues until the given number of generations has been completed, regardless of the score each generation receives. You can specify to print sets that receive a good enough to an output file. Set the filename with -o or --good-sets-file and set the threshold with -G or --good-set-threshold where 0.0 is a perfect score and 1.0 is a complete failure.

********************
**3.1.1: Using MPI**

Due to how computationally expensive SRES is, running hundreds of simulations for hundreds of generations on a single processor is impractical. To reduce the runtime, SRES can be distributed among multiple processors via MPI (Message Passing Interface). A full explanation of MPI can be found at http://www.open-mpi.org/. Basically, MPI allows multiple processes to communicate via messages, allowing multiple processes to work on different parts of the same task and report their results back to a process that organizes the results and doles out more tasks.

SRES uses MPI to distribute each simulation (i.e. population member) to a different process, each of which is ideally run on a different processor. This enables each generation to complete much faster than it would serially and allows productive evolutions to take place via large populations for many generations. Because libSRES, the library SRES uses, has separate files for MPI and non-MPI versions, SRES must be compiled with or without MPI and cannot run both versions with one compilation. For details on how to compile with MPI, see Section 1.

There are two possible commands to run SRES with MPI: mpiexec and mpirun. mpiexec is the standardized command and is usually the better choice. mpirun is implementation specific and may not work as expected. Since the correct command is environment specific, we recommend asking your server / cluster / network administrator. Using the command of your choice, set the number of processes to use via the command-line with "-np X", where X is the number of processes. This should be passed to mpiexec/mpirun and is automatically communicated to SRES. Here is an example call:
```
mpiexec -np 24 ./sres -d 45 -r range2014.txt -g 1000 -p 200 -P 30 -G good_set.txt -a -v
```

*************************
**3.2: Input and output**

***************************
**3.2.0: Methods of input**

The only interaction with a run is through input you give it before it starts and output it gives you after it ends. There are three aspects of input:
* Ranges
* Simulation calls
* Command-line arguments

*****************
**3.2.1: Ranges**

A ranges file must be passed via the command-line with -r or --ranges-file each time SRES is called. This ranges file defines the parameter space SRES looks throughout. If SRES fails to find parameter sets that receive high enough scores, consider generating a more restricted ranges file based on the best results so far. This may enable SRES to focus on more relevant regions in the space. Be aware, however, that the parameter space may contain peaks in unexpected locations; an alternative to reducing the space is to expand it in the hopes that the highest peaks are outside your current ranges. There is no known formal methodology for picking ranges for such an interconnected system beyond adopting experimental results found _in vivo_.

***************************
**3.2.2: Simulation calls**

Each SRES run needs to know where the simulation executable is located. The location can be specified via the command-line with -f or --simulation but defaults are set inside ./sres/source/structs/inside input_params.simulation feature if no argument is given. Ensure that any executable you specify matches the ranges file and number of dimensions given and that you have permission to execute the file.

*********************************
**3.2.3: Command-line arguments**

Each command-line argument can be entered in a short or long form. Short forms begin with a single dash (-) followed by a single letter. Long forms begin with a double dash (--) followed by a word or phrase. Lower-case letters are considered distinct from their upper-case equivalents. The short and long forms are equivalent - short forms are for convenience and long forms are for clarity. The following is a comprehensive list of all command-line options and their uses:

-r, --ranges-file        [filename]   : the relative filename of the ranges input file, default=none
-f, --simulation         [filename]   : the relative filename of the simulation executable, default=../simulation/simulation
-d, --dimensions         [int]        : the number of dimensions (i.e. rate parameters) to explore, min=1, default=45
-P, --parent-population  [int]        : the population of parent simulations to use each generation, min=1, default=3
-p, --total-population   [int]        : the population of total simulations to use each generation, min=1, default=20
-g, --generations        [int]        : the number of generations to run before returning results, min=1, default=1750
-G, --good-set-threshold [float]      : the threshold of score for good set. If a set scores less than or equal to the threshold (score 0 : good, score 1: bad), the set gets printed to file. Default: 0 (only print perfect sets)
-t, --print-good-set     [filename]   : the relative filename of the file to print good sets out
-s, --seed               [int]        : the seed used in the evolutionary strategy (not simulations), min=1, default=time
-e, --printing-precision [int]        : how many digits of precision parameters should be printed with, min=1, default=6
-a, --arguments          [N/A]        : every argument following this will be sent to the simulation
-c, --no-color           [N/A]        : disable coloring the terminal output, default=unused
-v, --verbose            [N/A]        : print detailed messages about the program state
-q, --quiet              [N/A]        : hide the terminal output, default=unused
-l, --licensing          [N/A]        : view licensing information (no simulations will be run)
-h, --help               [N/A]        : view usage information (i.e. this)

*****************************
**3.2.4: Input file formats**

**************************
**3.2.4.0: Ranges format**

For details on the format for a ranges file, see Section 2.2.4.3.

******************************
**3.2.5: Output file formats**

***********************************
**3.2.5.0: Terminal output format**

SRES prints its latest results after each generation finishes. Each result consists of two lines. The first line contains two values: current generation and best fitness score. They are printed in the format "current generation: ?, best fitness: ?", where each question mark contains the values. Generation values are integral and start at 1. Fitness values are floating point where 0.0 is a perfect score and 1.0 is a complete failure. The second line contains the parameter set with the best fitness score so far. It is printed in the format "best individual: ?,?,?", with a value for each dimension / parameter (here represented as question marks).

Note that the precision of the printed parameters is not as high as the precision of the parameters stored in memory that are sent to the simulation. Therefore, if the system is not robust, rerunning a parameter set taken from the terminal output with the same configuration may return different results because of slightly different values. Because most biological systems are quite robust, if 10^-8 percent changes in parameters produces different results then there may be a flaw in the simulation used. The simulation presented in this package represents a robust system and does not produce different results because of reduced precision except in the corner cases in which the parameter is on the edge of acceptability and the change in precision pushes it over. Regardless, the output precision can be changed via the command-line with -e or --printing-precision (this affects terminal output and good sets files).

*****************************
**3.2.5.1: Good sets format**

A good sets file consists of a list of parameter sets. Each parameter set having score lower than good_score_threshold (specified using "-G") is print out into file specified after flag "-t". Each line in this file represents a good parameter sets. The fist number in each line is the score that it gets, the following comma-separated numbers are values of parameters in the set.


*******************************************
**4: Refining ranges (refine-ranges.py)**

*******************
**4.0: Overview**

refine-ranges.py refines existing ranges into more precise ones based on parameter sets. Possible biologically realistic ranges are quite wide. After running enough parameter sets that produce valid results, more precise ranges may be inferred by creating ranges based on the average values for each parameter in the sets. More precisely, by giving this script a list of good parameter sets, your current ranges, and a filename to place the new, more refined ranges, this script finds the mean and standard deviation of every parameter and creates a new range by taking the values X standard deviations from the mean, where X is a positive integer you specify. However, if the new range is not as precise as the one found in the given current ranges file, the current range is used instead.

*********************************
**4.1: How to run the program:

**4.1.0: Command-line arguments**

* -s, --sets           [filename] : the relative filename of the parameter sets to run the simulations with, required
* -c, --current-ranges [filename] : the relative filename of the ranges file to base new ranges on, required
* -n, --new-ranges     [filename] : the relative filename of the ranges file to stores the new ranges in, required
* -d, --standard-dev   [int]      : how many standard deviations to stray from each mean range, default=2
* -r, --round-to       [int]      : how many digits to round each new range to, default=5
* -h, --help           [N/A]      : view usage information (i.e. this)

********************************
**4.1.1: Example program calls**

Assuming there is a parameter sets file named good.params and a ranges file named current.ranges, the following command refines the ranges using the parameter sets, using two standard deviations from the mean and (by default) rounding to 5 digits:

```
python refine-ranges.py -s good.params -c current.ranges -n new.ranges -d 2
```

5: Authorship and licensing
---------------------------

This program is created by Ha Vu and Ahmet Ay, 2017.
