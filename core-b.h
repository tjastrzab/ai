#ifndef _CORE_B_H_
#define _CORE_B_H_

#include "mpi.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>

#define MAX_A  26				/* max. number of symbols in alphabet */
#define MAX_WL 500				/* max. length of word */

#define TRUE   1
#define FALSE  0
#define UNDEFINED_STATE -55555
#define REMOVED_0_STATE -100000
#define REMOVED_0_WORD -500000
#define REMOVED_L_STATE 500000

#define REALLOC_SIZE 10
#define RW_REALLOC_SIZE 100

#define BUFFER_SIZE_0   7500000
#define BUFFER_SIZE_1   1000
#define BUFFER_SIZE_2   100

// End of array is marked with NULL
typedef char* word;
typedef int* number_set;
typedef char** word_set;

typedef struct word_set_with_size {
	word_set elem;
	int size;
} WSS;

struct st {
	number_set states;					/* states array */
	int size;						/* size of array of states */
	int word_no;					/* number of word in language and also reference to y */
} *T;							/* T - array of significant states for given word */
int sT;							/* size of T = size of L (at the beginning) */
word_set* y;					/* suffixes for states in T */

struct q {
	int exists;					/* exists = k means that state exists k times (k_max = sT) */
	int cur_size;
	int** rem;					/* removed states; size max: sT x MADFA_q */
	/* indexed with 'exists' */
	int* size;					/* rem lists sizes; size max: MADFA_q (start from index 1) */
} *S;							/* initial decomposition set */
int sS;							/* size of S; sS = MADFA_q */
int active_sS;					/* count of active (exists > 0) states */

// MADFA = Minimal Acyclic Deterministic Finite Automaton
char* MADFA_Sigma;				/* alphabet */
int Sigma_map[MAX_A];

int MADFA_q;					/* number of states */
word_set* MADFA_Q;				/* states */

int MADFA_f;					/* number of final states */
number_set MADFA_F;				/* final states (1st always = 0) */

number_set* MADFA_delta;		/* transition function; delta[state, number[letter]] = new_state */

WSS Lng;						/* input language */
WSS *LeftLng;					/* left Language for given state */
WSS *RightLng;					/* right language for given state */

int sD;							/* size of D */
int sP;							/* size of P */
number_set P;					/* decomposition states */
number_set D;					/* additional states */

number_set significant;			/* TRUE/FALSE array denoting state as significant */

int increment;
int* aux_array;					/* auxiliary array of size MADFA_q */
int* state_list;				/* useful while removing states for TPrime */
int sState_list;				/* size of state_list */

// Variables used in 'searching'
int* pos;						/* positions of combinations in D */
int* Sx;						/* states in the form of array */
int sSx;						/* size of Sx */

// Buffers for repeatability verification
int** buffer_1;					/* buffer of checked decompositions */
int** buffer_2;					/* buffer of accepted S & D pairs */
int buffer_idx_1;			/* buffer index */
int buffer_idx_2;			/* buffer index */

// Parameter count variables
long pairs_acc;				/* accepted S & D pairs */
long pairs_prop;			/* proposed S & D pairs */
long dec_checked;			/* decompositions checked */
long dec_to_check;			/* decompositions to be checked */
long dec_to_check_all;
long dec_to_check_rall;
long dec_to_check_rank;

// Timing variables
double aux_time;				/* auxiliary variable */
double sum_time;			/* total S & D processing time */
double start_time;
double meas_time_1, meas_time_2;/* time measurement variable */
double phase1_time, phase2_time;
double phase1_aux_time, phase2_aux_time;

int max_level;
int found_cnt;
int total_wrd_length;


FILE *test;						/* for testing */

int K;							/* cut-off threshold */
int level;					/* level of recursion */
char in_file[150];
char out_file[150];

int rank;
int comm_size;

number_set *all_Ps;
number_set final_st;

void init();
int check_P();
void create_D();
void remove_T();
void save_time();
void decompose();
void remove_MADFA();
void createLRLang();
void remove_buffers();
int check_if_checked();
int easy(word_set lang);
void create_state_list();
void remove_y(word_set* y);
int lambda_set(word_set Z);
int power(word_set wrd_st);
int minimalADFA(word_set X);
int catenation(WSS A, WSS B);
void remove_all_Ps_L1s_L2s();
void remove_state(int state);
int searching(int* D, int sD);
void restore_states(int state);
void save_max_level(char *type);
void define_significant_states();
void gather_all_decompositions();
void remove_S(struct q* S, int sS);
void remove_word_set(word_set wrd_st);
void verify_and_check_P(int *cacheNo);
int ends_with(char *wrd, char *ending);
int equal_sets(word_set A, word_set B);
void reduce_duplicate_decompositions();
int from_initial_to(int n, word_set Z);
WSS read_words_from_file(char* filename);
void assertion(int condition, char* msg);
int compare_strings(char **s1, char **s2);
word_set lq(char prefix, word_set wrd_st);
WSS copy_word_set(word_set wrd_st, int count);
int compare_ints(const void* a, const void* b);
WSS sum(WSS* words, int set_size, int sum_size);
void create_MADFA(word alphabet, word_set lang);
void find_alphabet(word_set lang, word alphabet);
WSS product(WSS* words, int set_size, int max_prod_size, int row);
void move_from(int from, int to, word lop, int lop_idx, word_set Z, int* z_idx);

#endif // !_CORE_H_
