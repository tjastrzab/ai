#include "core-b.h"

void init()
{
	sum_time = 0.0;
	start_time = 0.0;
	phase1_time = 0.0;
	phase2_time = 0.0;
	phase1_aux_time = 0.0;
	phase2_aux_time = 0.0;
	
	level = 0;
	max_level = 0;
	comm_size = 1;
	found_cnt = 0;
	pairs_acc = 0;
	pairs_prop = 0;
	dec_checked = 0;
	buffer_idx_1 = 0;
	buffer_idx_2 = 0;
	dec_to_check = 0;
	total_wrd_length = 0;
	dec_to_check_all = 0;
	dec_to_check_rall = 0;
	dec_to_check_rank = 0;
}

/*
* Builds array of additional states D for function searching(...)
*/
void create_D()
{
	int i, j;

	memset(aux_array, 0, MADFA_q * sizeof(int));
	for (i = 0; i < sT; i++)
	{
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] >= 0)
				{
					aux_array[T[i].states[j]] = 1;
				}
			}
		}
	}
	for (i = 0, j = 0; i < MADFA_q; i++)
	{
		if (aux_array[i] && !S[i].exists) {
			D[j++] = i;
		}
	}
	sD = j;
}

/*
* Deallocates structure T
*/
void remove_T()
{
	int i;

	assertion(sT >= 0, "Error 49");
	if (T != NULL)
	{
		for (i = 0; i < sT; i++)
		{
			assertion(T[i].size >= 0, "Error 50");
			free(T[i].states);
		}
		free(T);
	}
}

/* Saves execution time */
void save_time()
{
	aux_time = MPI_Wtime() - aux_time;
	sum_time = sum_time + aux_time;
}

/* Frees the memory for MADFA */
void remove_MADFA()
{
	int i;

	free(MADFA_Q);
	free(final_st);
	for (i = 0; i < MADFA_q; i++)
	{
		free(MADFA_delta[i]);
	}
	free(MADFA_delta);
	free(MADFA_Sigma);
}

/* Creates left and right languages */
void createLRLang()
{
	int i;
	word_set Z;

	Z = calloc(Lng.size, sizeof(word));
	assertion(Z != NULL, "Error 13");
	for (i = 0; i < Lng.size; i++)
	{
		Z[i] = calloc(MADFA_q, sizeof(char));
		assertion(Z[i] != NULL, "Error 14");
	}
	LeftLng = calloc(MADFA_q, sizeof(WSS));
	assertion(LeftLng != NULL, "Error 15");
	RightLng = calloc(MADFA_q, sizeof(WSS));
	assertion(RightLng != NULL, "Error 16");
	for (i = 0; i < MADFA_q; i++)
	{
		if (significant[i])
		{
			LeftLng[i] = copy_word_set(Z, from_initial_to(i, Z));
			RightLng[i].elem = MADFA_Q[i];
			RightLng[i].size = power(RightLng[i].elem);
		}
        }

	for (i = 0; i < Lng.size; i++)
	{
		free(Z[i]);
	}
	free(Z);
}

/* Clears memory for buffers */
void remove_buffers()
{
	int i;

	for (i = 0; i < BUFFER_SIZE_1; i++)
	{
		free(buffer_1[i]);
	}
	for (i = 0; i < BUFFER_SIZE_2; i++)
	{
		free(buffer_2[i]);
	}
	free(buffer_1);
	free(buffer_2);
}

int check_if_checked()
{
	int i, j, cnt, is_in_buffer;

	for (i = 0; i < BUFFER_SIZE_1; ++i)
	{
		is_in_buffer = TRUE;
		for (j = 0; j < sP; ++j)
		{
			if (!buffer_1[i][P[j]])
			{
				is_in_buffer = FALSE;
				break;
			}
		}
		if (is_in_buffer)
		{
			for (j = 0, cnt = 0; j < MADFA_q; ++j)
			{
				cnt += buffer_1[i][j];
			}
			if (cnt == sP)
			{
				return FALSE;
			}
		}
	}
	memset(buffer_1[buffer_idx_1], 0, MADFA_q * sizeof(int));
	for (i = 0; i < sP; ++i)
	{
		buffer_1[buffer_idx_1][P[i]] = 1;  /* decomposition set */
	}
	buffer_idx_1 = (++buffer_idx_1 % BUFFER_SIZE_1);
	return TRUE;
}

/*
* Checks if language has the form (Xa^{-1})a
* @param X				- language
* @return				- 0 if language is not of the form (Xa^{-1})a, 1 otherwise
*/
int easy(word_set lang)
{
	int i = 1;
	char last_letter;

	if (lang[0][0] == '\0')
	{
		return FALSE;
	}
	last_letter = lang[0][strlen(lang[0]) - 1];
	while (lang[i] != NULL)
	{
		if (lang[i][strlen(lang[i]) - 1] != last_letter)
		{
			return FALSE;
		}
		++i;
	}
	return TRUE;
}

/*
* Finds different states in T
*/
void create_state_list()
{
	int i, j;

	memset(aux_array, 0, MADFA_q * sizeof(int));
	for (i = 0; i < sT; i++)
	{
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] >= 0)
				{
					aux_array[T[i].states[j]] = 1;
				}
			}
		}
	}
	for (j = 0, i = 0; i < MADFA_q; i++)
	{
		if (aux_array[i] != 0)
		{
			state_list[j++] = i;
		}
	}
	sState_list = j;
}


/*
* Deallocates suffixes y
* @param y				- set of suffixes
*/
void remove_y(word_set* y)
{
	int i;

	assertion(sT > 0, "Error 49a");
	assertion(y != NULL, "Error 49b");
	for (i = 0; i < sT; i++)
	{
		free(y[i]);
	}
	free(y);
}

/*
* Checks whether the set is lambda set
* @param Z				- set to be verified
* @return				- 0 if Z is not lambda, 1 otherwise
*/
int lambda_set(word_set Z)
{
	return (Z[0] != NULL) && (Z[1] == NULL) && (strlen(Z[0]) == 0);
}

/*
* Finds the power of a word set
* @param wrd_st			- word set
* @return					- power (count of words) in wrd_st
*/
int power(word_set wrd_st)
{
	int i = 0;

	while (wrd_st[i++] != NULL);
	return i - 1;
}

/*
* Builds the transition map (delta)
* @param X			- word set being the basis for MADFA
* @return			- next step in transition map
*/
int minimalADFA(word_set X)
{
	word_set U;
	int i, p, i_letter;

	MADFA_Q[MADFA_q] = X;
	if (X[0][0] == '\0')
	{
		MADFA_F[MADFA_f++] = MADFA_q;
	}
	p = MADFA_q++;
	i_letter = 0;
	while (MADFA_Sigma[i_letter])
	{
		U = lq(MADFA_Sigma[i_letter], X);
		if (U[0] != NULL)
		{
			i = 0;
			while (MADFA_Q[i] != NULL)
			{
				if (equal_sets(MADFA_Q[i], U))
				{
					break;
				}
				i++;
			}
			if (MADFA_Q[i])
			{
				MADFA_delta[p][i_letter] = i;
				remove_word_set(U);
			}
			else
			{
				MADFA_delta[p][i_letter] = minimalADFA(U);
			}
		}
		else
		{
			remove_word_set(U);
		}
		i_letter++;
	}
	return p;
}

/*
* Performs word sets catenation
* @param A, B		- structures representing word sets of known size
* @return			- result of comparison between language resulting from catenation and L
*/
int catenation(WSS A, WSS B)
{
	char *new_item;
	word found_item;
	int *words_found;	// stores 0 if word generated, 1 otherwise
	int i, j, idx, cnt, len;

	new_item = calloc(MAX_WL, sizeof(char));
	words_found = calloc(Lng.size, sizeof(int));
	for (i = 0, cnt = 0; i < A.size; ++i)
	{
		len = strlen(A.elem[i]);
		strcpy(new_item, A.elem[i]);
		for (j = 0; j < B.size; ++j)
		{
			strcpy(new_item + len, B.elem[j]);
			found_item = bsearch(&new_item, Lng.elem, Lng.size, sizeof(word), compare_strings);
			idx = (found_item - (char*)Lng.elem) / sizeof(word);
			if (!words_found[idx])
			{
				++cnt;
				words_found[idx] = 1;
			}
		}
	}
	free(new_item);
	free(words_found);
	return (cnt == Lng.size);
}

void remove_all_Ps_L1s_L2s()
{
	int i;

	for (i = 0; i < found_cnt; i++)
	{
		free(all_Ps[i]);
	}
	free(all_Ps);
}

/* Finds significant states */
void define_significant_states()
{
	int i, j, arrows_out;
	int n = strlen(MADFA_Sigma);

	significant = calloc(MADFA_q, sizeof(int));
	assertion(significant != NULL, "Error 12");
	for (i = 0; i < MADFA_q; i++)
	{
		arrows_out = 0;
		for (j = 0; j < n; j++)
		{
			if (MADFA_delta[i][j] != UNDEFINED_STATE)
			{
				arrows_out++;
			}
		}
		significant[i] = ((final_st[i] && (arrows_out > 0)) || (!final_st[i]/* && (arrows_out > 1)*/));
	}
}

/*
* Deallocates set of states S
* @param S				- set of states
* @param sS			- size of S
*/
void remove_S(struct q* S, int sS)
{
	int i, j;

	assertion(sS >= 0, "Error 49");
	if (S != NULL)
	{
		for (i = 0; i < sS; i++)
		{
			for (j = 0; j < S[i].cur_size; j++)
			{
				free(S[i].rem[j]);
			}
			free(S[i].size);
			free(S[i].rem);
		}
		free(S);
	}
}

/*
* Deallocates word set
* @param wrd_st			- word set to deallocate
*/
void remove_word_set(word_set wrd_st)
{
	int i;

	for (i = 0; wrd_st[i] != NULL; i++)
	{
		free(wrd_st[i]);
	}
	free(wrd_st);
}

/*
* Checks whether given buffer should be accepted and
* verifies if it constitutes decomposition set
*/
void verify_and_check_P(int *cacheNo)
{
	dec_to_check++;
	if (((*cacheNo)++ % comm_size) == rank)
	{
		dec_to_check_rank++;
		check_P();
	}
	if (*cacheNo == 1000)
	{
		*cacheNo = 0;
	}
}

int ends_with(char *wrd, char *ending)
{
	int i;
	int lenW = strlen(wrd) - 1;
	int lenE = strlen(ending) - 1;
	
	for (i = lenE; i >= 0; --i) 
	{
		if (wrd[lenW - lenE + i] != ending[i])
		{
			return FALSE;
		}
	}
	return TRUE;
}

/*
* Compares words sets
* @param A, B		- word sets to compare
* @return			- 0 for not equal sets, 1 for equal ones
*/
int equal_sets(word_set A, word_set B)
{
	int i = 0;

	while ((A[i] != NULL) && (B[i] != NULL))
	{
		if (strcmp(A[i], B[i]) != 0)
		{
			return FALSE;
		}
		++i;
	}
	return A[i] == B[i];
}

void reduce_duplicate_decompositions()
{
	int i, j, k;
	int new_cnt = 0;
	number_set *new_all_Ps;
	int is_in_buffer1, is_in_buffer2;

	if (found_cnt > 0)
	{
		new_all_Ps = calloc(found_cnt, sizeof(number_set));
		for (i = 0; i < found_cnt; i++)
		{
			is_in_buffer2 = FALSE;
			for (j = i + 1; j < found_cnt; j++)
			{
				is_in_buffer1 = TRUE;
				for(k = 0; k < MADFA_q; k++)
				{
					if(all_Ps[i][k] != all_Ps[j][k])					
					{
						is_in_buffer1 = FALSE;
						break;
					}
				}
				if(is_in_buffer1)
				{
					is_in_buffer2 = TRUE;
					break;
				}
			}
			if (!is_in_buffer2)
			{
				new_all_Ps[new_cnt] = calloc(MADFA_q, sizeof(int));
				memcpy(new_all_Ps[new_cnt], all_Ps[i], MADFA_q * sizeof(int));
				new_cnt++;
			}
		}
		remove_all_Ps_L1s_L2s();
		all_Ps = realloc(new_all_Ps, new_cnt * sizeof(number_set));
		assertion(all_Ps != NULL, "Error reduce Ps");
		found_cnt = new_cnt;
	}
}

/*
* Traces the path from initial state to given one
* @param n				- destination state
* @param Z				- output set of words on the path
* @return				- length of Z
*/
int from_initial_to(int n, word_set Z)
{
	int index = 0;
	char lop[MAX_WL];

	move_from(0, n, lop, 0, Z, &index);
	qsort(Z, index, sizeof(word), compare_strings);
	return index;
}

/*
* Reads words from file. In the 1st row there is
* number of words, next lines contain words.
* Words are unique.
* @param filename		- name of the file to open
* @return				- language
*/
WSS read_words_from_file(char* filename)
{
	FILE *file;
	WSS out_wss;
	int i, len, position = 0;
	char line[MAX_WL], buffer[BUFFER_SIZE_0];

	if (rank == 0)
	{
		file = fopen(filename, "r");
		assertion(file != NULL, "Error opening file");
		fscanf(file, "%d\n", &out_wss.size);
		MPI_Pack(&out_wss.size, 1, MPI_INT, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
		out_wss.elem = calloc(out_wss.size + 1, sizeof(word));
		assertion(out_wss.elem != NULL, "Error 11");
		for (i = 0; i < out_wss.size; i++)
		{
			line[0] = '\0';
			fscanf(file, "%s\n", line);
			out_wss.elem[i] = strdup(line);
			len = strlen(line) + 1;
			total_wrd_length += len;
			MPI_Pack(&len, 1, MPI_INT, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
			MPI_Pack(line, len, MPI_CHAR, buffer, BUFFER_SIZE_0, &position, MPI_COMM_WORLD);
			assertion(position < BUFFER_SIZE_0, "Error 66");
		}
		fclose(file);
		out_wss.elem[out_wss.size] = NULL;
		qsort(out_wss.elem, out_wss.size, sizeof(word), compare_strings);
		MPI_Bcast(buffer, BUFFER_SIZE_0, MPI_PACKED, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Bcast(buffer, BUFFER_SIZE_0, MPI_PACKED, 0, MPI_COMM_WORLD);
		MPI_Unpack(buffer, BUFFER_SIZE_0, &position, &out_wss.size, 1, MPI_INT, MPI_COMM_WORLD);
		out_wss.elem = calloc(out_wss.size + 1, sizeof(word));
		assertion(out_wss.elem != NULL, "Error 11a");
		for (i = 0; i < out_wss.size; i++)
		{
			MPI_Unpack(buffer, BUFFER_SIZE_0, &position, &len, 1, MPI_INT, MPI_COMM_WORLD);
			total_wrd_length += len;
			MPI_Unpack(buffer, BUFFER_SIZE_0, &position, line, len, MPI_CHAR, MPI_COMM_WORLD);
			out_wss.elem[i] = strdup(line);
		}
		out_wss.elem[out_wss.size] = NULL;
		qsort(out_wss.elem, out_wss.size, sizeof(word), compare_strings);
	}
	return out_wss;
}

/*
* Performs assertion check
* @param condition			- condition to verify
* @param msg				- message to store in case of error
*/
void assertion(int condition, char* msg)
{
	if (!condition)
	{
		fprintf(test, "%s\n", msg);
		fclose(test);
		MPI_Finalize();
		exit(1);
	}
}

/*
* Compares strings
* @param s1, s2	- strings to compare
* @return			- -1 if s1 < s2, 0 if s1 = s2, 1 if s1 > s2
*/
int compare_strings(char **s1, char **s2)
{
	return strcoll(*s1, *s2);
}

/*
* Finds the left quotient of word set

* @param prefix			- word prefix used for lq finding
* @param wrd_st			- word set
* @return					- word set of lq suffixes
*/
word_set lq(char prefix, word_set wrd_st)
{
	word_set out_wrd_st;
	int i, j = 0, count_w = 0;

	for (i = 0; wrd_st[i] != NULL; i++)
	if (wrd_st[i][0] == prefix)
	{
		++count_w;
	}
	out_wrd_st = calloc(count_w + 1, sizeof(word));
	assertion(out_wrd_st != NULL, "Error 2");
	for (i = 0; j < count_w; i++)
	{
		if (wrd_st[i][0] == prefix)
		{
			out_wrd_st[j++] = strdup(wrd_st[i] + 1);
		}
	}
	out_wrd_st[j] = NULL;
	qsort(out_wrd_st, count_w, sizeof(word), compare_strings);
	return out_wrd_st;
}

/*
* Copies word set of 'count' elements
* @param wrd_st			- word set to copy
* @param count				- power of wrd_st
* @return					- a copy of wrd_st
*/
WSS copy_word_set(word_set wrd_st, int count)
{
	int i;
	WSS out_copy;

	out_copy.elem = calloc(count + 1, sizeof(word));
	assertion(out_copy.elem != NULL, "Error 1");
	for (i = 0; i < count; i++)
	{
		out_copy.elem[i] = strdup(wrd_st[i]);
	}
	out_copy.elem[i] = NULL;
	out_copy.size = count;
	return out_copy;
}

/*
* Compares integers
* @param a, b		- integers to compare
* @return			- -1 if a < b, 0 if a = b, 1 if a > b
*/
int compare_ints(const void* a, const void* b)
{
	int* arg1 = (int*)a;
	int* arg2 = (int*)b;
	int diff = *arg1 - *arg2;

	return (diff < 0 ? -1 : (diff == 0 ? 0 : 1));
}

/*
* Finds the sum of set_size word sets
* @param words			- sets whose sum is to be found
* @param set_size		- size of the array of word sets
* @param sum_size		- size of resulting set
* @return				- sum of all sets
*/
WSS sum(WSS* words, int set_size, int sum_size)
{
	int i, j, m;
	WSS out_sum;

	out_sum.size = sum_size;
	out_sum.elem = calloc(sum_size + 1, sizeof(word));
	assertion(out_sum.elem != NULL, "Error 20");
	for (i = 0, m = 0; i < set_size; i++)
	{
		for (j = 0; j < words[i].size; j++)
		{
			out_sum.elem[m++] = words[i].elem[j];
		}
	}
	out_sum.elem[m] = NULL;
	return out_sum;
}

/*
* Builds whole MADFA (delta, F, Q, T)
* @param alphabet		- alphabet
* @param lang			- language
*/
void create_MADFA(word alphabet, word_set lang)
{
	int** delta_ptr;
	word_set* Q_ptr;
	int i, j, max_state_count = 0, alph_len;  /* maximum state count */

	MADFA_q = 0;
	MADFA_f = 0;
	alph_len = strlen(alphabet);
	MADFA_Sigma = strdup(alphabet);
	for (i = 0; lang[i] != NULL; i++)
	{
		max_state_count += strlen(lang[i]);
	}
	max_state_count += 1;
	MADFA_Q = calloc(max_state_count, sizeof(word_set));
	assertion(MADFA_Q != NULL, "Error 3");
	MADFA_F = calloc(max_state_count, sizeof(int));
	assertion(MADFA_F != NULL, "Error 4");
	MADFA_delta = calloc(max_state_count, sizeof(number_set));
	assertion(MADFA_delta != NULL, "Error 5");
	for (i = 0; i < max_state_count; i++)
	{
		MADFA_F[i] = -1;
		MADFA_Q[i] = NULL;
		MADFA_delta[i] = calloc(alph_len, sizeof(int));
		assertion(MADFA_delta[i] != NULL, "Error 6");
		for (j = 0; j < alph_len; j++)
		{
			MADFA_delta[i][j] = UNDEFINED_STATE;
		}
	}
	minimalADFA(lang);
	Q_ptr = realloc(MADFA_Q, MADFA_q * sizeof(word_set));
	assertion(Q_ptr != NULL, "Error 7");
	MADFA_Q = Q_ptr;
	delta_ptr = realloc(MADFA_delta, MADFA_q * sizeof(number_set));
	assertion(delta_ptr != NULL, "Error 9");
	MADFA_delta = delta_ptr;
	final_st = calloc(MADFA_q, sizeof(int));
	for (i = 0; i < MADFA_f; i++)
	{
		final_st[MADFA_F[i]] = 1;
	}
	free(MADFA_F);
}

/*
* Determines the alphabet
* @param lang			- language
* @param alphabet		- output array with alphabet
*/
void find_alphabet(word_set lang, word alphabet)
{
	word s;
	int i, j;

	memset(alphabet, '\0', MAX_A);
	memset(Sigma_map, MAX_A, MAX_A * sizeof(int));
	for (i = 0; lang[i] != NULL; i++)
	{
		for (s = lang[i]; *s != '\0'; s++)
		{
			alphabet[(*s) - 'a'] = *s;
		}
	}
	for (i = 1; i < MAX_A; i++)
	{
		if (alphabet[i] != '\0')
		{
			j = i - 1;
			while ((j >= 0) && (alphabet[j] == '\0'))
			{
				alphabet[j] = alphabet[j + 1];
				alphabet[j + 1] = '\0';
				j--;
			}
		}
	}
	for (i = 0; i < strlen(alphabet); ++i)
	{
		Sigma_map[alphabet[i] - 'a'] = i;
	}
}

/*
* Finds the product of set_size word sets
* @param words			- sets whose product is to be found
* @param set_size		- size of the array of word sets
* @param prod_size		- maximum size of product = min RightLng size
* @param row			- row containing prod_size elements
* @return				- product of all sets
*/
WSS product(WSS* words, int set_size, int prod_size, int row)
{
	word wrd;
	WSS out_prod;
	int i, j, m;
	int common_cnt, not_found;

	out_prod.elem = calloc(prod_size + 1, sizeof(word)); /* min + 1 - max. product size */
	assertion(out_prod.elem != NULL, "Error 19");
	for (i = 0, m = 0; i < prod_size; i++)
	{
		wrd = words[row].elem[i];
		for (j = 0, common_cnt = 0; j < set_size; j++)
		{
			if (j != row)
			{
				not_found = (bsearch(&wrd, words[j].elem, words[j].size, sizeof(word), compare_strings) == NULL);
				if (not_found)
				{
					break;
				}
				else
				{
					common_cnt++;
				}

			}
		}
		if (common_cnt == set_size - 1)
		{
			out_prod.elem[m++] = words[row].elem[i];
		}
	}
	out_prod.elem[m] = NULL;
	out_prod.size = m;
	return out_prod;
}

/*
* Traces the path between two states
* @param from				- source state
* @param to				- destination state
* @param lop				- letters on path
* @param lop_idx			- position in lop
* @param Z					- output set with words on path
* @param z_idx				- position in Z
*/
void move_from(int from, int to, word lop, int lop_idx, word_set Z, int* z_idx)
{
	if (from == to)
	{
		lop[lop_idx] = '\0';
		strcpy(Z[(*z_idx)], lop);
		(*z_idx)++;
	}
	else
	{
		int j;
		char letter;

		for (j = 0, letter = MADFA_Sigma[0]; letter != '\0'; letter = MADFA_Sigma[++j])
		{
			if (MADFA_delta[from][j] != UNDEFINED_STATE)
			{
				lop[lop_idx] = letter;
				move_from(MADFA_delta[from][j], to, lop, lop_idx + 1, Z, z_idx);
			}
		}
	}
}

void gather_all_decompositions()
{
	int *int_line;
	MPI_Status status;
	int cur_val, max_val, sum;
	number_set *tmp_all_Ps;
	number_set sizes;
	number_set found_counts;
	number_set displacements;
	char line[MAX_WL], *buffer, *all_buffer;
	int i, j, k, len, tmp, size, position = 0;
	double s_time, e_time, t_time1 = 0.0, t_time2 = 0.0, t_time3 = 0.0, t_time4 = 0.0;

	s_time = MPI_Wtime();
	sizes = calloc(comm_size, sizeof(int));
	found_counts = calloc(comm_size, sizeof(int));
	displacements = calloc(comm_size, sizeof(int));
	int_line = calloc(MADFA_q, sizeof(int));
	cur_val = (found_cnt * (MADFA_q)) * sizeof(int); // for all_Ps
	e_time = MPI_Wtime();
	t_time1 += (e_time - s_time);
	s_time = MPI_Wtime();
	max_val = 0;
	MPI_Gather(&found_cnt, 1, MPI_INT, found_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		sum = 0;
		for(i = 0; i < comm_size; i++)
		{
			sum += found_counts[i];
			sizes[i] = found_counts[i] * MADFA_q * sizeof(int);
			max_val += sizes[i];
			if(i > 0)
			{
				displacements[i] = displacements[i - 1] + sizes[i - 1];
			}
		}
	}
	e_time = MPI_Wtime();
	t_time2 += (e_time - s_time);
	s_time = MPI_Wtime();
	buffer = malloc(cur_val * sizeof(char));
	all_buffer = malloc(max_val * sizeof(char));
	e_time = MPI_Wtime();
	t_time1 += (e_time - s_time);
	s_time = MPI_Wtime();
	for(i = 0; i < found_cnt; i++)
	{
		MPI_Pack(all_Ps[i], MADFA_q, MPI_INT, buffer, cur_val, &position, MPI_COMM_WORLD);
	}
	e_time = MPI_Wtime();
	t_time3 += (e_time - s_time);
	s_time = MPI_Wtime();
	MPI_Gatherv(buffer, cur_val, MPI_PACKED, all_buffer, sizes, displacements, MPI_PACKED, 0, MPI_COMM_WORLD);
	e_time = MPI_Wtime();
	t_time4 += (e_time - s_time);
	s_time = MPI_Wtime();
	if(rank == 0)
	{
		found_cnt = 0;
		tmp_all_Ps = realloc(all_Ps, sum * sizeof(number_set));
		assertion(tmp_all_Ps != NULL, "Error Gather 3");
		all_Ps = tmp_all_Ps;
		for(i = 0; i < comm_size; i++)
		{
			if(found_counts[i] > 0)
			{
				position = 0;
				for(j = 0; j < found_counts[i]; j++)
				{
					all_Ps[found_cnt + j] = calloc(MADFA_q, sizeof(int));
					MPI_Unpack(&all_buffer[displacements[i]], sizes[i], &position, int_line, MADFA_q, MPI_INT, MPI_COMM_WORLD);
					memcpy(all_Ps[found_cnt + j], int_line, MADFA_q * sizeof(int));
				}
				found_cnt += found_counts[i];
			}
		}
		if(found_cnt > 0)
		{
			reduce_duplicate_decompositions();
			printf("Decomposition set(s) found: %d\n", found_cnt);
		}
		else
		{
			printf("Decomposition set NOT found\n");
		}
	}
	e_time = MPI_Wtime();
	t_time3 += (e_time - s_time);
	printf("Timings %d: %5.2f, %5.2f, %5.2f, %5.2f\n", rank, t_time1, t_time2, t_time3, t_time4);
	fflush(stdout);
	free(sizes);
	free(buffer);
	free(int_line);	
	free(all_buffer);
	free(found_counts);
	free(displacements);
}
