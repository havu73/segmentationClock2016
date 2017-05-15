#ifndef MACROS_HPP
#define MACROS_HPP

#include <math.h>

// STATES
#define MH1		0 	//mRNA Her1
#define MH7		1 	//mRNA Her7
#define MD		2 	//mRNA delta
#define PH1		3 	//protein her1
#define	PH7		4 	//protein her7
#define	PD		5 	//protein delta 
#define	PH11	6 	//dimer Her1- Her1
#define	PH17	7 	//dimer Her1- Her7
#define	PH77	8 	//dimer Her7- Her7


#define	NUM_STATES			9	//Total number of states in the system

// RECORD STATES
#define KEEPMH1 	0
#define KEEPMH7		1
#define NUM_RECORD_STATES 	2	// We only want to keep record of concentrations of MH1 and MH7 
// RATES
#define PSH1		0		//protein Her1 synthesis rate
#define PSH7		1		//protein Her7 synthesis rate
#define PSD			2		//protein delta systhesis rate

#define PDH1		3		//protein Her1 degradation rate
#define PDH7		4		//protein Her7 degradation rate 
#define PDD			5		//protein delta degradation rate

#define MSH1		6		//mRNA Her1 transcription rate
#define MSH7		7		//mRNA Her7 transcription rate 
#define MSD			8		//mRNA delta transcription rate

#define MDH1		9		//mRNA Her1 degradation rate
#define MDH7		10		//mRNA Her7 degradation rate 
#define MDD			11		//mRNA delta degradation rate

#define DAH11		12		//dimer H1- H1 association rate
#define DAH17		13		//dimer H1-H7 association rate
#define DAH77		14		//dimer H7- H7 association rate

#define DDH11		15		//dimer H1- H1 dissociation rate
#define DDH17		16		//dimer H1-H7 dissocation rate
#define DDH77		17		//dimer H7-H7 dissocation rate

#define PDH11		18		//dimer degradation rate H1-H1
#define PDH17		19		//dimer degradation rate H1-H7
#define PDH77		20		//dimer degradation rate H7-H7

#define NUM_CONSTANT_RATES		21 // The number of rates that are kept between rates and parameters structures

#define NMH1		21		//transcription time delay mRNA H1
#define NMH7		22		//transcription time delay mRNA H7
#define NMD			23		//transcription time delay mRNA delta

#define NPH1		24		//translation time delay H1
#define NPH7		25		//translation time delay H7
#define NPD			26		//translation time delay delta

#define CRITPDH1	27		// critical PD
#define CRITPH11H1	28		// critical PH11
#define CRITPDH7	29
#define CRITPH11H7	30
#define CRITPH11D	31

#define NUM_RATES	32		// Total number of rates in the systems

#define KAG1PN		27		//association constant for gene1 and protein NICD
#define	KAG1PH11	28		//association constant for gene1 and protein dimer H1 H1
#define KAG7PN		29		//association constant for gene7 and protein NICD
#define KAG7PH11	30		//association constant for gene7 and protein dimer H1 H1
#define KAGDPH11	31		//association constant for gene Delta and protein dimer H1 H1

#define KDG1PN		32		//dissociation constant for gene1 and protein NICD
#define KDG1PH11	33		//dissociation constant for gene1 and protein dimer H1 H1
#define KDG7PN		34		//dissociation constant for gene7 and protein NICD
#define KDG7PH11	35		//dissociation constant for gene7 and protein dimer H1 H1
#define KDGDPH11	36		//dissociation constant for gene Delta and protein dimer H1 H1

#define NUM_PARAMS	37		//total number of parameters in the system

// Testing parameters
#define LOW_BOUND_NUM_PEAKS	3
#define LOW_BOUND_NUM_AMP	3
#define LOW_BOUND_P2T		1.5
#define UP_BOUND_M2L		2
#define LOW_BOUND_PERIOD	27
#define UP_BOUND_PERIOD		33

// scoring system
#define MAX_SCORE		5
// NUM MUTANTS
#define NUM_MUTANTS 	1


// Exit statuses
#define EXIT_SUCCESS			0
#define EXIT_MEMORY_ERROR		1
#define EXIT_FILE_READ_ERROR	2
#define EXIT_FILE_WRITE_ERROR	3
#define EXIT_PIPE_READ_ERROR	4
#define EXIT_PIPE_WRITE_ERROR	5
#define EXIT_INPUT_ERROR		6

#endif

