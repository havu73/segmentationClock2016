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
#define	PH7		4 	//protein her7
#define	PD		5 	//protein delta 
#define	PH11	6 	//dimer Her1- Her1
#define	PH17	7 	//dimer Her1- Her7
#define	PH77	8 	//dimer Her7- Her7
#define	G1		9 	//gene her1
#define	G1N		10 	//gene her1 bound by protein nicd
#define	G1PH11	11 	//gene her1 bound by dimer Her1-Her1
#define	G7		12 	//gene her7
#define	G7N		13 	//gene her7 bound by protein nicd
#define	G7PH11	14 	//gene her7 bound by dimer Her1-Her1
#define	GD		15 	//gene delta
#define	GDPH11	16 	//gene delta bound by dimer Her1-Her1

#define	NUM_STATES			17	//Total number of states in the system

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

#define NMH1		21		//transcription time delay mRNA H1
#define NMH7		22		//transcription time delay mRNA H7
#define NMD			23		//transcription time delay mRNA delta

#define NPH1		24		//translation time delay H1
#define NPH7		25		//translation time delay H7
#define NPD			26		//translation time delay delta


// rates that we will get from sres
#define CRITPDH1	27		// critical PD-H1
#define CRITPH11H1	28		// critical PH11-H1
#define CRITPDH7	29		// critical PD-H7
#define CRITPH11H7	30		// critical PH11-H7
#define CRITPH11D	31		// critical PH11-D

//rates that we will calculate based on what sres gave us
/* KonPH11 = Koff / (critPH11 ^ 2)
 * KonPD = Koff / (critPD)
 */ 
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

#define NUM_RATES	37		//total number of rates in the system

// REACTIONS
#define RPSH1		0		//mh1 --> ph1
#define RPSH7		1		//mh7 --> ph7 
#define RPSD		2		//md -> pd 

#define RPDH1		3		//ph1 --> null
#define RPDH7		4		//ph7 --> null
#define RPDD		5		//pd -> null
#define RPDH11		6		//ph11 --> null
#define RPDH17		7		//ph17 --> null 
#define RPDH77		8		//ph77 --> null 

#define RDAH11		9		//ph1 + ph1 --> ph11
#define RDAH17		10		//ph1 + ph7 --> ph17
#define RDAH77		11		//ph7 + ph7 --> ph77

#define RDDH11		12		//ph11 --> ph1 + ph1
#define RDDH17		13		//ph17 --> ph1 + ph7 
#define RDDH77		14		//ph77 --> ph7 + ph7 

#define RMDH1		15		//mh1 --> null
#define RMDH7		16		//mh7 -->null
#define RMDD		17		//md -> null

#define RMSH1		18		//g1 -> mh1
#define RMSH1N		19		//g1n -> mh1

#define RAG1PH11	20		//g1 + ph11 --> g1ph11
#define RDG1PH11	21		//g1ph11 --> g1 + ph11

#define RAG1N		22		//g1 + nicd -> g1n
#define RDG1N		23		//g1n -> g1 + nicd

#define RMSH7		24		//g7 -> mh7
#define RMSH7N		25		//g7n -> mh7

#define RAG7PH11	26		//g7 + ph11 --> g7ph11
#define RDG7PH11	27		//g7ph11 --> g7 + ph11

#define RAG7N		28		//g7 + nicd -> g7n
#define RDG7N		29		//g7n -> g7 + nicd

#define RMSD		30		//gd --> md

#define RAGDPH11	31		//gd + ph11 --> gdph11
#define RDGDPH11	32		//gdph11 --> gd + ph11

#define NUM_REACTIONS		33		//Total number of reactions

//WHAT STATES TO KEEP RECORD OF OVER TIME
#define KEEPMH1		0		// Index of MH1
#define KEEPMH7		1		// Index of MH7

#define NUM_KEEP_STATES		2

// MUTANTS
#define WT					0
#define DELTA				1
#define NUM_MUTANTS			1	// should be 2

//MUTANTS SCORES
#define WT_SC				0	// should be different

#define TOTAL_SC			0	// should be different

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

// PARAMETERS FOR CHECKING CORRELATIONS
#define LOW_BOUND_WT_CORR 	0.7

// PARAMETERS FOR BINNING DATA
#define	DEFAULT_NUM_BIN 	3

// PARAMETER FOR SLICING
#define MINUTE_PER_SLICE	5

// CONDITIONS CHECKING PARAMETERS
#define UPPER_BOUND_HER		250
#define LOWER_BOUND_PTT		1.5
#define UPPER_BOUND_MTL		1.5
#define WT_H1_AMP_LOW		10//35
#define WT_H1_AMP_HIGH		70//45
#define WT_H7_AMP_LOW		10//43
#define WT_H7_AMP_HIGH		70//51
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

// Score
#define	MAX_SCORE 				6

// EXIT STATUS
#define EXIT_SUCCESS			0
#define EXIT_MEMORY_ERROR		1
#define EXIT_FILE_READ_ERROR	2
#define EXIT_FILE_WRITE_ERROR	3
#define EXIT_PIPE_READ_ERROR	4
#define EXIT_PIPE_WRITE_ERROR	5
#define EXIT_INPUT_ERROR		6

#endif
