/*
* Basic scenario - no adaptive K
*/

#include "basic.h"

/* Creates initial state array T and suffix array y */
void create_Ty()
{
	word s;
	int i, j, k, m, current_st;

	sT = Lng.size;
	T = calloc(sT, sizeof(struct st));
	assertion(T != NULL, "Error 30");
	y = calloc(sT, sizeof(word_set));
	assertion(y != NULL, "Error 30a");
	for (i = 0; i < sT; i++)
	{
		/* computing 'size' of row T[i].states */
		j = 0;
		current_st = 0;
		for (s = Lng.elem[i]; *s != '\0'; s++)
		{
			if (significant[current_st])
			{
				j++;
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			j++;
		}
		T[i].states = calloc(j, sizeof(int));
		assertion(T[i].states != NULL, "Error 32");
		T[i].size = j;
	}
	for (i = 0; i < Lng.size; i++)
	{
		T[i].word_no = i; /* word number (also reference in y) */
		y[i] = calloc(T[i].size, sizeof(word));
		assertion(y[i] != NULL, "Error 32a");
		current_st = 0;
		/* k is the position in word */
		for (s = Lng.elem[i], j = 0, k = -1; *s != '\0'; s++, k++)
		{
			if (significant[current_st])
			{
				T[i].states[j] = current_st;
				y[i][j++] = Lng.elem[i] + k + 1;
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			T[i].states[j] = current_st;
			y[i][j++] = Lng.elem[i] + k + 1;
		}
	}
}

int check_P()
{
	WSS* words;
	WSS L1, L2;
	int i, prod_size_row;
	number_set *tmp_all_Ps;
	int sum_size = 0, max_prod_size = 0, ret_value = 0;

	for (i = 0; i < sP; i++)
	{
		sum_size += LeftLng[P[i]].size;
	}
	for (max_prod_size = Lng.size + 1, prod_size_row = -1, i = 0; i < sP; i++)
	{
		if (RightLng[P[i]].size < max_prod_size)
		{
			max_prod_size = RightLng[P[i]].size;
			prod_size_row = i;
		}
	}
	if (!check_if_checked())
	{
		return FALSE;
	}
	words = malloc(sP * sizeof(WSS));
	for (i = 0; i < sP; i++)
	{
		words[i] = RightLng[P[i]];
	}
	L2 = product(words, sP, max_prod_size, prod_size_row);
	for (i = 0; i < sP; i++)
	{
		words[i] = LeftLng[P[i]];
	}
	L1 = sum(words, sP, sum_size);
	free(words);
	dec_checked++;
	ret_value = !(lambda_set(L1.elem) || lambda_set(L2.elem));
	if (ret_value)
	{
		ret_value = catenation(L1, L2);
		if (ret_value)
		{
			found_cnt++;
			tmp_all_Ps = realloc(all_Ps, found_cnt * sizeof(number_set));
			assertion(tmp_all_Ps != NULL, "Error Ps");
			all_Ps = tmp_all_Ps;
			all_Ps[found_cnt - 1] = calloc(MADFA_q, sizeof(int));
			printf("Found %d: %d", rank, found_cnt);
			for (i = 0; i < sP; i++)
			{
				all_Ps[found_cnt - 1][P[i]] = 1;
				printf(" %d", P[i]);
			}
			printf("\n");
			fflush(stdout);
		}
	}
	free(L1.elem);
	free(L2.elem);
	return ret_value;
}

void decompose()
{
	int* state;
	int sState;
	word suffix;
	int **rem_ptr;
	int *size_ptr;
	int new_size = 0;
	word_set suffixes;
	int marked_for_skip;
	int i, j, k, l, min, list_no, r, row, tmp;

	level++;
	/* step 1. - decomposition below cut-off level */
	if (sState_list - active_sS <= K)
	{
		phase1_aux_time = MPI_Wtime() - phase1_aux_time;
		phase1_time += phase1_aux_time;
		create_D();
		phase2_aux_time = MPI_Wtime();
		searching(D, sD);
		phase2_aux_time = MPI_Wtime() - phase2_aux_time;
		phase2_time += phase2_aux_time;
		printf("Rt1t2: %d, %.5f, %.5f\n", rank, phase1_aux_time, phase2_aux_time);
		level--;
		phase1_aux_time = MPI_Wtime();
		return;
	}
	/* step 2. - searching row of T with min. number of states */
	assertion(sT > 0, "Error 37");
	for (row = -1, min = MADFA_q, i = 0; i < sT; i++)
	{
		if (T[i].word_no >= 0)
		{
			for (tmp = 0, j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] >= 0)
				{
					tmp++;
				}
			}
			if (tmp < min)
			{
				min = tmp;
				row = i;
			}
		}
	}
	if (row == -1)
	{
		level--;
		return;
	}
	assertion(row >= 0, "Error 38");
	if (min == 0) /* T[i] = Fi */
	{
		level--;
		return;
	}
	/* step 3. - searching recursively */
	/* 3a - restoring states from 'row' */
	state = malloc(MADFA_q * sizeof(int));
	assertion(state != NULL, "Error 39");
	for (k = 0, j = 0; j < T[row].size; j++)
	{
		if (T[row].states[j] >= 0)
		{
			state[k++] = T[row].states[j];
		}
	}
	sState = k;
	assertion(sState > 0, "Error 42");
	suffixes = y[T[row].word_no];
	T[row].word_no = ((T[row].word_no == 0) ? REMOVED_0_WORD : -T[row].word_no);
	/* Inserting states into S */
	for (r = 0; r < sState; r++)  /* r is the index in array 'state' */
	{
		assertion(state[r] >= 0, "Error 60");
		if (r != 0)
		{
			/* restoring removed states basing on TPrime, S, and state[r - 1] */
			tmp = state[r - 1];
			restore_states(tmp);
			assertion(S[tmp].exists > 0, "Error 40");
			S[tmp].size[S[tmp].exists--] = 0;			/* list out-of-date, retracting previous state */
			if (S[tmp].exists == 0)
			{
				active_sS--;
			}
			assertion(S[tmp].exists >= 0, "Error 40a");
		}
		list_no = ++S[state[r]].exists;				/* inserting current state */
		if (list_no == 1)
		{
			active_sS++;
		}
		assertion(list_no > 0, "Error 60a");
		assertion(list_no <= Lng.size, "Error 60b");
		/* finding position of state[r] in T[row].states[.] */
		for (i = 0; i < T[row].size; i++)
		{
			if (T[row].states[i] == state[r])
			{
				break;
			}
		}
		assertion(T[row].states[i] == state[r], "Error 44");
		suffix = suffixes[i];
		/* actually used states in TPrime */
		create_state_list();
		l = 0; /* index in removed states list */
		if (S[state[r]].cur_size <= list_no)
		{
			new_size = S[state[r]].cur_size + increment;
			rem_ptr = realloc(S[state[r]].rem, new_size * sizeof(number_set));
			assertion(rem_ptr != NULL, "Error 56");
			S[state[r]].rem = rem_ptr;
			for (j = S[state[r]].cur_size; j < new_size; j++)
			{
				S[state[r]].rem[j] = calloc(MADFA_q, sizeof(int));
				assertion(S[state[r]].rem[j] != NULL, "Error 56a");
			}
			free(S[state[r]].rem[0]);
			S[state[r]].rem[0] = NULL; /* 0th list always empty */
			new_size = S[state[r]].cur_size + increment;
			size_ptr = realloc(S[state[r]].size, new_size * sizeof(int));
			assertion(size_ptr != NULL, "Error 57");
			S[state[r]].size = size_ptr;
			S[state[r]].cur_size = new_size;
		}
		marked_for_skip = FALSE;
		for (i = 0; i < sState_list; i++) /* state removal */
		{
			assertion(significant[state_list[i]], "Error 48");
			assertion(state_list[i] >= 0, "Error 48a");
			if (bsearch(&suffix, RightLng[state_list[i]].elem, RightLng[state_list[i]].size, sizeof(word), compare_strings) == NULL)
			{
				remove_state(state_list[i]);
				S[state[r]].rem[list_no][l++] = state_list[i];
				if (S[state_list[i]].exists > 0)
				{
					marked_for_skip = TRUE;
					break;
				}
			}
		}
		S[state[r]].size[list_no] = l;
		sState_list -= l;
		if (marked_for_skip)
		{
			continue;
		}
		decompose();
	}
	assertion(sState - 1 >= 0, "Error 61");
	tmp = state[sState - 1];
	restore_states(tmp);
	assertion(S[tmp].exists > 0, "Error 61a");
	S[tmp].size[S[tmp].exists--] = 0;
	if (S[tmp].exists == 0)
	{
		active_sS--;
	}
	assertion(S[tmp].exists >= 0, "Error 61b");
	T[row].word_no = ((T[row].word_no == REMOVED_0_WORD) ? 0 : -T[row].word_no);
	free(state);
	level--;
}

/*
* Removes given state from  the list of states
* @param state			- state to be removed
*/
void remove_state(int state)
{
	int i, j, val;

	assertion(state >= 0, "Error 58");
	val = ((state > 0) ? -state : REMOVED_0_STATE);
	for (i = 0; i < sT; i++)
	{
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] == state)
				{
					T[i].states[j] = val;
					break;
				}
			}
		}
	}
}

/*
* Searches for decompositions set
* @param D			- set of additional states
* param sD			- size of D
*/
int searching(int* D, int sD)
{
	int cacheNo = 0;
	int sum_lls, min_rls, expSum;
	int i, j, k, m, r, cnt, is_in_buffer;
	
	for (j = 0, i = 0; i < sS; i++)
	{
		if (S[i].exists > 0)
		{
			Sx[j++] = i;
		}
	}
	sSx = j;
	pairs_prop++;
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		is_in_buffer = TRUE;
		for (j = 0; j < sSx; j++)
		{
			if (buffer_2[i][Sx[j]] != 1)
			{
				is_in_buffer = FALSE;
				break;
			}
		}
		if (is_in_buffer == TRUE)
		{
			for (j = 0; j < sD; j++)
			{
				if (buffer_2[i][D[j]] != 2)
				{
					is_in_buffer = FALSE;
					break;
				}
			}
		}
		if (is_in_buffer == TRUE)
		{
			for (j = 0, cnt = 0; j < MADFA_q; j++)
			{
				if (buffer_2[i][j] == 1)
				{
					cnt++;
				}
				if (cnt == sSx)
				{
					return FALSE;
				}
			}
		}
	}
	pairs_acc++;
	memset(buffer_2[buffer_idx_2], 0, MADFA_q * sizeof(int));
	for (i = 0; i < sSx; i++)
	{
		buffer_2[buffer_idx_2][Sx[i]] = 1;
	}
	for (i = 0; i < sD; i++)
	{
		buffer_2[buffer_idx_2][D[i]] = 2;
	}
	buffer_idx_2 = (++buffer_idx_2 % BUFFER_SIZE_2);
	aux_time = MPI_Wtime();
	memcpy(P, Sx, sSx * sizeof(int));
	dec_to_check_all += pow(2.0, sD);
	dec_to_check_rall += pow(2.0, K);
	sP = sSx;
	if (sSx != 0)
	{
		verify_and_check_P(&cacheNo);
	}
	for (k = 1; k <= sD; k++)
	{
		m = sSx;
		for (i = 0; i < k; i++)
		{
			pos[i] = i;
			P[m++] = D[i];
		}
		sP = m;
		verify_and_check_P(&cacheNo);
		do {
			j = k - 1;
			while (j >= 0)
			{
				if (pos[j] < sD - k + j)
				{
					pos[j]++;
					for (r = j + 1; r < k; r++)
					{
						pos[r] = pos[r - 1] + 1;
					}
					m = sSx;
					for (i = 0; i < k; i++)
					{
						P[m++] = D[pos[i]];
					}
					sP = m;
					verify_and_check_P(&cacheNo);
					break;
				}
				else
				{
					j--;
				}
			}
		} while (j >= 0);
	}
	save_time();
	return FALSE;
}

/*
* Restores removed states from T basing on S and 'state'
* @param state			- state to be removed
*/
void restore_states(int state)
{
	int list_no;		/* list of states to be restored */
	int i, j, k, s;		/* state to be restored */
	int val, new_val;

	assertion(state >= 0, "Error 59");
	list_no = S[state].exists;
	assertion(list_no > 0, "Error 59a");
	assertion(S[state].size[list_no] >= 0, "Error 59b");
	for (k = 0; k < S[state].size[list_no]; k++)
	{
		s = S[state].rem[list_no][k];	/* state to recover */
		assertion(s >= 0, "Error 59c");
		if (s > 0)
		{
			val = -s;
			new_val = s;
		}
		else
		{
			val = REMOVED_0_STATE;
			new_val = 0;
		}
		for (i = 0; i < sT; i++)
		{
			if (T[i].word_no >= 0)
			{
				for (j = 0; j < T[i].size; j++)
				{
					if (T[i].states[j] == val)
					{
						T[i].states[j] = new_val;
						break;
					}
				}
			}
		}
	}
}

int main(int argc, char **argv)
{
	int i;
	FILE* o_file;
	double total_time;
	char alphabet[MAX_A];

	init();

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	total_time = MPI_Wtime();
	start_time = MPI_Wtime();

	if (argc == 4)
	{
		K = atoi(argv[1]);
		strcpy(in_file, argv[2]);
		sprintf(out_file, "%s%d", argv[3], rank);
	}
	else
	{
		K = 10;
		strcpy(in_file, "we");
		sprintf(out_file, "%s%d", "wy", rank);
	}

	if (!(test = fopen("test", "a")))
	{
		printf("Error opening file test\n");
		MPI_Finalize();
		return 1;
	}
	Lng = read_words_from_file(in_file);
	if (Lng.size == 0)
	{
		remove_word_set(Lng.elem);
		MPI_Finalize();
		return 0;
	}
	if (easy(Lng.elem))
	{
		printf("Composite language, easily decomposable.");
		remove_word_set(Lng.elem);
		MPI_Finalize();
		return 0;
	}
	find_alphabet(Lng.elem, alphabet);
	create_MADFA(alphabet, Lng.elem);
	define_significant_states();
	createLRLang();
	create_Ty();
	increment = Lng.size / 5;
	S = calloc(MADFA_q, sizeof(struct q));
	assertion(S != NULL, "Error 55");
	sS = MADFA_q;
	for (i = 0; i < sS; i++)
	{
		S[i].exists = 0;
		S[i].rem = NULL;
		S[i].size = NULL;
		S[i].cur_size = 0;
	}
	buffer_1 = malloc(BUFFER_SIZE_1 * sizeof(number_set));
	assertion(buffer_1 != NULL, "Error 63");
	for (i = 0; i < BUFFER_SIZE_1; i++)
	{
		buffer_1[i] = calloc(MADFA_q, sizeof(int));
		assertion(buffer_1[i] != NULL, "Error 64");
	}
	buffer_2 = malloc(BUFFER_SIZE_2 * sizeof(number_set));
	assertion(buffer_2 != NULL, "Error 63a");
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		buffer_2[i] = calloc(MADFA_q, sizeof(int));
		assertion(buffer_2[i] != NULL, "Error 64a");
	}
	D = calloc(MADFA_q, sizeof(int));
	assertion(D != NULL, "Error 28");
	P = calloc(MADFA_q, sizeof(int));
	assertion(P != NULL, "Error 27");
	aux_array = malloc(MADFA_q * sizeof(int));
	assertion(aux_array != NULL, "Error 33");
	state_list = malloc(MADFA_q * sizeof(int));
	assertion(state_list != NULL, "Error 47");
	Sx = malloc(MADFA_q * sizeof(int));
	assertion(Sx != NULL, "Error 34");
	pos = malloc(MADFA_q * sizeof(int));
	assertion(pos != NULL, "Error 35");
	meas_time_1 = MPI_Wtime();
	create_state_list();
	phase1_aux_time = MPI_Wtime();
	decompose();
	phase1_aux_time = MPI_Wtime() - phase1_aux_time;
	phase1_time += phase1_aux_time;
	printf("Rt1: %d, %.5f\n", rank, phase1_aux_time);
	meas_time_1 = MPI_Wtime() - meas_time_1;
	meas_time_2 = MPI_Wtime();
	reduce_duplicate_decompositions();
	gather_all_decompositions();
	if (o_file = fopen(out_file, "a"))
	{
		fprintf(o_file, "No. of decompositions in %d = %d\n", rank, found_cnt);
	}
	else
	{
		printf("No. of decompositions in %d = %d\n", rank, found_cnt);
	}
	meas_time_2 = MPI_Wtime() - meas_time_2;
	remove_all_Ps_L1s_L2s();
	free(P);
	free(D);
	free(Sx);
	free(pos);
	free(aux_array);
	free(state_list);
	for (i = 0; i < MADFA_q; i++)
	{
		if (significant[i])
		{
			remove_word_set(LeftLng[i].elem);
			remove_word_set(RightLng[i].elem);
		}
	}
	free(LeftLng);
	free(RightLng);
	remove_MADFA();
	remove_y(y);
	remove_S(S, sS);
	remove_T();
	remove_buffers();
	free(significant);
	fclose(test);
	total_time = MPI_Wtime() - total_time;
	if (o_file)
	{
		fprintf(o_file, "Rank: %d; execution times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2, phase1_time, phase2_time);
		fprintf(o_file, "Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, dec_to_check_all, dec_to_check_rall);
		fclose(o_file);
	}
	else
	{
		printf("Rank: %d; execution times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2, phase1_time, phase2_time);
		printf("Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, dec_to_check_all, dec_to_check_rall);
	}
	MPI_Finalize();
	return 0;
}
