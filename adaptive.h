#ifndef _ADAPTIVE_
#define _ADAPTIVE_

#include "core-a.h"
#include "adapt.h"

typedef struct state_words {
	int word_no;
	struct state_words *next;
} sw;
sw **states_to_words;

typedef struct state_with_lls {
	int lls;
	int state_no;
} LLS;

typedef struct rls_tracker {
	int level;
	number_set counts;
	struct rls_tracker *next;
	struct rls_tracker *prev;
} rls;
rls *rlsTail;

typedef struct min_tracker {
	int min;
	int row;
	struct min_tracker *next;
	struct min_tracker *prev;
} mt;
mt *mtTail;

counter *suffCounts;
int suffCountsMap[MAX_A];
counter *suffCountsPtr[MAX_A];						
int s_sc;
							
number_set actualRLS;			/* stores actual right lang. sizes */

int changed = TRUE;	// signals local significance change

long changedCount = 0;
long statesRemoved = 0;

double mem_time = 0.0;
double remove1_time = 0.0;
double remove2_time = 0.0;
double restore_time = 0.0;
double dec_sum_time = 0.0;
double dec_prd_time = 0.0;
double dec_cat_time = 0.0;

void create_ss();
void create_counts();
void create_T_and_y();
void create_T_and_y_dec(int level);
void restore_T_and_y_dec(int level);
int compare_lls(LLS *lls1, LLS *lls2);
int compare_ss(counter *ss1, counter *ss2);
int compare_counters(counter *f1, counter *f2);

#endif
