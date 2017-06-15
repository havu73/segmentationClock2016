#ifndef MACROS_HPP
#define MACROS_HPP

// NEIGHBOR NETWORK PARAMETERS
#define SIXTEEN_NEIGHBORS	6	// each cell has maximum 6 neighbors
#define FOUR_NEIGHBORS 		3
#define TWO_NEIGHBORS		1

// STATES
#define MH1		0 	//mRNA Her1
#define MH7		1 	//mRNA Her7
#define MD		2 	//mRNA delta
#define PH1		3 	//protein her1
#define	PD		4 	//protein delta 
#define	PH11	5 	//dimer Her1- Her1
#define	G1		6 	//gene her1
#define	G1N		7 	//gene her1 bound by protein nicd
#define	G1PH11	8 	//gene her1 bound by dimer Her1-Her1
#define	G7		9 	//gene her7
#define	G7N		10 	//gene her7 bound by protein nicd
#define	G7PH11	11 	//gene her7 bound by dimer Her1-Her1
#define	GD		12 	//gene delta
#define	GDPH11	13 	//gene delta bound by dimer Her1-Her1

#define	NUM_STATES			14	//Total number of states in the system

// RATES
#define PSH1		0		//protein Her1 synthesis rate
#define PSD			1		//protein delta systhesis rate

#define PDH1		2		//protein Her1 degradation rate
#define PDD			3		//protein delta degradation rate

#define MSH1		4		//mRNA Her1 transcription rate
#define MSH7		5		//mRNA Her7 transcription rate 
#define MSD			6		//mRNA delta transcription rate

#define MDH1		7		//mRNA Her1 degradation rate
#define MDH7		8		//mRNA Her7 degradation rate 
#define MDD			9		//mRNA delta degradation rate

#define DAH11		10		//dimer H1- H1 association rate

#define DDH11		11		//dimer H1- H1 dissociation rate

#define PDH11		12		//dimer degradation rate H1-H1

#define NMH1		13		//transcription time delay mRNA H1
#define NMH7		14		//transcription time delay mRNA H7
#define NMD			15		//transcription time delay mRNA delta

#define NPH1		16		//translation time delay H1
#define NPD			17		//translation time delay delta


//rates that we will calculate based on what sres gave us
/* KonPH11 = Koff / (critPH11 ^ 2)
 * KonPD = Koff / (critPD)
 */ 
#define KAG1PN		18		//association constant for gene1 and protein NICD
#define	KAG1PH11	19		//association constant for gene1 and protein dimer H1 H1
#define KAG7PN		20		//association constant for gene7 and protein NICD
#define KAG7PH11	21		//association constant for gene7 and protein dimer H1 H1
#define KAGDPH11	22		//association constant for gene Delta and protein dimer H1 H1

#define KDG1PN		23		//dissociation constant for gene1 and protein NICD
#define KDG1PH11	24		//dissociation constant for gene1 and protein dimer H1 H1
#define KDG7PN		25		//dissociation constant for gene7 and protein NICD
#define KDG7PH11	26		//dissociation constant for gene7 and protein dimer H1 H1
#define KDGDPH11	27		//dissociation constant for gene Delta and protein dimer H1 H1

#define GRADIENT	28		// How much the end of the anterior's rates of NPH1 and NPD should be
#define NUM_NETWORK_RATES	29

#define PERTURB		29		// maximum variation each parameter in each cell should be 

#define NUM_RATES	30		//total number of rates in the system

// REACTIONS
#define RPSH1		0		//mh1 --> ph1
#define RPSD		1		//md -> pd 

#define RPDH1		2		//ph1 --> null
#define RPDD		3		//pd -> null
#define RPDH11		4		//ph11 --> null

#define RDAH11		5		//ph1 + ph1 --> ph11

#define RDDH11		6		//ph11 --> ph1 + ph1

#define RMDH1		7		//mh1 --> null
#define RMDH7		8		//mh7 -->null
#define RMDD		9		//md -> null

#define RMSH1		10		//g1 -> mh1
#define RMSH1N		11		//g1n -> mh1

#define RAG1PH11	12		//g1 + ph11 --> g1ph11
#define RDG1PH11	13		//g1ph11 --> g1 + ph11

#define RAG1N		14		//g1 + nicd -> g1n
#define RDG1N		15		//g1n -> g1 + nicd

#define RMSH7		16		//g7 -> mh7
#define RMSH7N		17		//g7n -> mh7

#define RAG7PH11	18		//g7 + ph11 --> g7ph11
#define RDG7PH11	19		//g7ph11 --> g7 + ph11

#define RAG7N		20		//g7 + nicd -> g7n
#define RDG7N		21		//g7n -> g7 + nicd

#define RMSD		22		//gd --> md

#define RAGDPH11	23		//gd + ph11 --> gdph11
#define RDGDPH11	24		//gdph11 --> gd + ph11

#define NUM_REACTIONS		25		//Total number of reactions

//WHAT STATES TO KEEP RECORD OF OVER TIME
#define KEEPMH1		0		// Index of MH1
#define KEEPMH7		1		// Index of MH7

#define NUM_KEEP_STATES		2

// MUTANTS and MUTANT SCORE
#define WT					0
#define WT_SCORE			7
#define DELTA_MUTANT		1
#define DELTA_SCORE			1
#define DAPT_MUTANT			2
#define DAPT_SCORE			0
#define NUM_MUTANTS			1	// should be 2

//MUTANTS SCORES
#define TOTAL_SC			7	// should be different

// TIME RELATED TO GRADIENTS
#define ABSOLUTE_RATE_TIME	60	// from 0 mins to ABSOLUTE_RATE_TIME, the protein synthesis delay rates should be kept as exactely as the input
								// But after that, after every 5 mins, these rates will be changed by a certain constants, to amke sure that
								// gradient from posterior to anterior are increasing/decreasing
								
// PARAMETERS RELATED TO THE TIMING OF MUTANTS
#define INDUCTION_DAPT_TIME	100
#define WITHDRAW_DAPT_TIME	190

// PARAMETERS FOR SMOOTHING DATA AND CHECKING PEAKS AND TROUGHS
#define WINDOW_SIZE			30 // each smoothed data point is equal to the average of data point and 200 points to the left and 200 points ot the right
#define WING_CHECK_SIZE_LEVEL_1		10 // each potential peak/trough needs to satisfy that its concentrations is higher / lower than the WING_CHECK_SIZE_LEVEL_1 data points to its left and right
#define WING_CHECK_SIZE_LEVEL_2		10	// after they pass the first test, the potential peaks/troughs will also need to pass the second test of the next WING_CHECK_SIZE_LEVEL_2 points
#define WING_CHECK_SIZE_LEVEL_3		30	// after they pass the first test, the potential peaks/troughs will also need to pass the second test of the next WING_CHECK_SIZE_LEVEL_3 points
#define NUM_MID_PT			3
#define NUM_LAST_PT			3
inline int low_num_pt(){  // in order to be considered to have sustained oscillation, need to have at least LOW_NUM_PT peaks and at lease LOW_NUM_PT troughs
	return NUM_MID_PT + NUM_LAST_PT + 1;
}

// Default user parameter values
#define	DEFAULT_NUM_BIN 	3
#define DEFAULT_NUM_SET		1
#define DEFAULT_TOTAL_TIME	310
#define DEFAULT_NUM_CELLS	4
#define	DEFAULT_RECORD_GRAN	0.1

// number of bins for calculating coefficent of variation squared
#define CVS_BIN		5
#define CVS_BIN_TIME	60

// PARAMETER FOR SLICING
#define MINUTE_PER_SLICE	5

// CORRELATION WT THRESHOLD
#define WT_CORR_SCORE		0.75

// CONDITIONS CHECKING PARAMETERS
#define LOWER_BOUND_PTT		1.5
#define UPPER_BOUND_MTL		1.5
#define WT_H1_AMP_LOW		10//35
#define WT_H1_AMP_HIGH		70//45
#define WT_H7_AMP_LOW		10//43
#define WT_H7_AMP_HIGH		70//51
#define WT_H1_ALT_AMP_LOW	35
#define WT_H1_ALT_AMP_HIGH	45
#define WT_H7_ALT_AMP_LOW	43
#define WT_H7_ALT_AMP_HIGH	51
#define WT_LOW_PER			22//27
#define WT_HIGH_PER			38//33
// Average Her mRNA levels in bins
#define	LOW_WT_HER_BIN1			25.06664788
#define UP_WT_HER_BIN1			37.59997181
#define	LOW_WT_HER_BIN2			61.93060553
#define UP_WT_HER_BIN2			92.89590829
#define	LOW_WT_HER_BIN3			98.56669283
#define	UP_WT_HER_BIN3			147.8500392
// Average intrinsic noise in bins
#define LOW_WT_IN_NOISE_BIN1	0.160721585
#define UP_WT_IN_NOISE_BIN1		0.241082377
#define	LOW_WT_IN_NOISE_BIN2	0.061040281
#define UP_WT_IN_NOISE_BIN2		0.091560421
#define	LOW_WT_IN_NOISE_BIN3	0.029502671
#define UP_WT_IN_NOISE_BIN3		0.044254006
// Average extrinsic noise in bins
#define	LOW_WT_EX_NOISE_BIN1	0.213986597
#define	UP_WT_EX_NOISE_BIN1		0.320979895
#define	LOW_WT_EX_NOISE_BIN2	0.167874357
#define	UP_WT_EX_NOISE_BIN2		0.251811535
#define LOW_WT_EX_NOISE_BIN3	0.136449288
#define UP_WT_EX_NOISE_BIN3		0.204673932

// PARAMETERS FOR TESTS IN DELTA MUTANT
// Ratio between delta noise / wildtype noise
#define DELTA_NOISE_WT_RATIO	1.4
// ranges for H1 and H7 spatial amplitude
#define DELTA_H1_AMP_LOW		12
#define DELTA_H1_AMP_HIGH		30
#define	DELTA_H7_AMP_LOW		18
#define DELTA_H7_AMP_HIGH		42

// Average Her mRNA levels in bins
#define LOW_DELTA_HER_BIN1		2.803207317
#define	UP_DELTA_HER_BIN1		11.16748999
#define	LOW_DELTA_HER_BIN2		7.950440396
#define	UP_DELTA_HER_BIN2		33.33296424
#define LOW_DELTA_HER_BIN3		13.18983534
#define	UP_DELTA_HER_BIN3		55.90542736

// Average intrinsic noise in bins
#define	LOW_DELTA_IN_NOISE_BIN1	0.130888597
#define	UP_DELTA_IN_NOISE_BIN1	0.397140778
#define	LOW_DELTA_IN_NOISE_BIN2	0.059928008
#define UP_DELTA_IN_NOISE_BIN2	0.203655053
#define LOW_DELTA_IN_NOISE_BIN3	0.031897972
#define UP_DELTA_IN_NOISE_BIN3	0.143711875

// Average extrinsic noise in bins
#define LOW_DELTA_EX_NOISE_BIN1	0.160419697
#define UP_DELTA_EX_NOISE_BIN1	0.512638225
#define LOW_DELTA_EX_NOISE_BIN2	0.198647086
#define UP_DELTA_EX_NOISE_BIN2	0.645999064
#define LOW_DELTA_EX_NOISE_BIN3	0.162110177
#define UP_DELTA_EX_NOISE_BIN3	0.603766761

// LIMITS OF MRNA and PROTEINs
#define MRNA_LIM				200
#define PROTEIN_LIM				3000


// EXIT STATUS
#define EXIT_SUCCESS			0
#define EXIT_MEMORY_ERROR		1
#define EXIT_FILE_READ_ERROR	2
#define EXIT_FILE_WRITE_ERROR	3
#define EXIT_PIPE_READ_ERROR	4
#define EXIT_PIPE_WRITE_ERROR	5
#define EXIT_INPUT_ERROR		6

#define MAX(x, y) ((x) < (y) ? (y) : (x))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
