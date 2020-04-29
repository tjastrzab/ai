/*
* Scenario with Reduce & Verify & locally significant states - adaptive K
*/

#include "adaptive.h"

void create_ss()
{
	create_T_and_y();
}

void create_counts()
{
	word s;
	int i, j, s_cnt = 0;
	counter* suffCountsTmp;
	counter *suffixes = malloc(total_wrd_length * sizeof(counter));
	
	for(i = 0; i < MADFA_q; i++)
	{
		if(significant[i])
		{
			for(j = 0; j < RightLng[i].size; j++)
			{
				suffixes[s_cnt].word = RightLng[i].elem[j];
				suffixes[s_cnt].count = LeftLng[i].size;
				s_cnt++;
			}	
		}
	}
	qsort(suffixes, s_cnt, sizeof(counter), compare_counters);
	s_sc = 0;
	suffCounts = malloc(s_cnt * sizeof(counter));
	s = suffixes[0].word;
	suffCounts[s_sc].word = s;
	suffCounts[s_sc].level = 0;
	suffCounts[s_sc].count = suffixes[0].count;
	suffCounts[s_sc].states = calloc(MADFA_q, sizeof(int));
	suffCounts[s_sc].statesTmp = calloc(MADFA_q, sizeof(int));
	if(s[0] != '\0')
	{
		suffCountsPtr[s[0] - 'a'] = &suffCounts[s_sc];
		suffCountsMap[s[0] - 'a'] = 1;
	}
	for(i = 1; i < s_cnt; i++)
	{		
		if (strcmp(s, suffixes[i].word))
		{
			s_sc++;
			s = suffixes[i].word;
			suffCounts[s_sc].word = s;
			suffCounts[s_sc].level = 0;
			suffCounts[s_sc].count = suffixes[i].count;
			suffCounts[s_sc].states = calloc(MADFA_q, sizeof(int));
			suffCounts[s_sc].statesTmp = calloc(MADFA_q, sizeof(int));
			suffCountsMap[s[0] - 'a']++;
			if(suffCountsPtr[s[0] - 'a'] == NULL || strcmp(s, suffCountsPtr[s[0] - 'a']->word) < 0)
			{
				suffCountsPtr[s[0] - 'a'] = &suffCounts[s_sc];
			}
		}
		else
		{
				suffCounts[s_sc].count += suffixes[i].count;
		}
	}
	s_sc++;
	suffCountsTmp = realloc(suffCounts, s_sc * sizeof(counter));
	assertion(suffCountsTmp != NULL, "Error ssf");
	suffCounts = suffCountsTmp;
	free(suffixes);
}

/* Counts locally significant states during decomposition */
void create_T_and_y_dec(int level)
{
	WSS **ptr;
	rls *newRLS;
	counter search, *result;
	int i, j, k, i_idx, idx;

	changed = FALSE;
	newRLS = malloc(sizeof(rls));
	newRLS->level = level;
	newRLS->prev = rlsTail;
	newRLS->next = NULL;
	newRLS->counts = calloc(MADFA_q, sizeof(int));

	// Clear the states in suffix counts
	for(i = 0; i < s_sc; i++)
	{
		memset(suffCounts[i].statesTmp, 0, MADFA_q * sizeof(int));
	}	
	for (i = 0; i < sT; i++)
	{
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] >= 0)
				{
					result = T[i].suffixes[j];
					idx = T[i].states[j];
					switch(result->statesTmp[idx])
					{
						case 0:
							if ((long) (result->count * (long) actualRLS[idx]) < Lng.size)
							{
								++statesRemoved;
								changed = TRUE;
								state_list[idx]--;
								if(!state_list[idx])
								{
									sState_list--;
								}
								newRLS->counts[idx]++;
								T[i].actualSize--;
								if(T[i].actualSize < mtTail->min)
								{
									mtTail->min = T[i].actualSize;
									mtTail->row = i;
								}
								T[i].states[j] = -idx - REMOVED_L_STATE;
								result->statesTmp[idx] = T[i].states[j];
							}
							else
							{
								result->statesTmp[idx] = 1;
							}		
							break;
						case 1:
							break;
						default:
							++statesRemoved;
							changed = TRUE;
							state_list[idx]--;
							if(!state_list[idx])
							{
								sState_list--;
							}
							T[i].actualSize--;
							if(T[i].actualSize < mtTail->min)
							{
								mtTail->min = T[i].actualSize;
								mtTail->row = i;
							}
							T[i].states[j] = result->statesTmp[idx];
							break;
					}
				}
			}
		}
	}
	if (changed)
	{
		++changedCount;
		rlsTail->next = newRLS;
		rlsTail = newRLS;
		for (i = 0; i < MADFA_q; i++)
		{
			actualRLS[i] -= newRLS->counts[i];
		}
	}
	else
	{
		free(newRLS->counts);
		free(newRLS);
	}
}

void restore_T_and_y_dec(int level)
{
	int i, j, k, idx, i_idx;
	counter search, *result;

	changed = FALSE;
	if(rlsTail->level == level)
	{
		changed = TRUE;
		for(i = 0; i < s_sc; i++)
		{
			memset(suffCounts[i].statesTmp, 0, MADFA_q * sizeof(int));
		}
		for (i = 0; i < MADFA_q; i++)
		{
			actualRLS[i] += rlsTail->counts[i];
		}
		rlsTail = rlsTail->prev;
		free(rlsTail->next->counts);
		free(rlsTail->next);
		rlsTail->next = NULL;
		
		for (i = 0; i < sT; i++)
		{
			if (T[i].word_no >= 0)
			{
				for (j = 0; j < T[i].size; j++)
				{
					if (T[i].states[j] != PERM_REMOVED_STATE && T[i].states[j] <= -REMOVED_L_STATE)
					{
						result = T[i].suffixes[j];				
						idx = -T[i].states[j] - REMOVED_L_STATE;
						if(result->statesTmp[idx] == 0)
						{
							if ((long) (result->count * (long) actualRLS[idx]) >= Lng.size)
							{
								T[i].states[j] = idx;
								state_list[idx]++;
								if(state_list[idx] == 1)
								{
									sState_list++;
								}
								T[i].actualSize++;
								result->statesTmp[idx] = idx;
							}
							else
							{
								result->statesTmp[idx] = T[i].states[j];
							}		
						}
						else if(result->statesTmp[idx] > 0)
						{
							T[i].states[j] = result->statesTmp[idx];
							state_list[idx]++;
							if(state_list[idx] == 1)
							{
								sState_list++;
							}
							T[i].actualSize++;
						}
					}
				}
			}
		}
	}
}

int compare_lls(LLS *lls1, LLS *lls2)
{
	int lls_lls1 = lls1->lls;
	int lls_lls2 = lls2->lls;
	int diff = lls_lls1 - lls_lls2;

	return (diff < 0 ? 1 : (diff == 0 ? 0 : -1));
}
int compare_ss(counter *ss1, counter *ss2)
{
	int result = compare_strings(&ss1->word, &ss2->word);

	if (!result)
	{
		result = compare_ints(&ss1->count, &ss2->count);
	}
	return result;
}

void create_T_and_y()
{
	word s;
	sw *swPtr;
	int i, j, k, current_st;
	counter countTmp, *foundCountTmp;
	
	sT = Lng.size;
	T = malloc(Lng.size * sizeof(struct st));
	assertion(T != NULL, "error T");
	for (i = 0; i < Lng.size; i++)
	{
		T[i].word_no = i;
		T[i].actualSize = 0;
		T[i].size = strlen(Lng.elem[i]) + 1;
		T[i].states = calloc(T[i].size, sizeof(int));
		T[i].suffixes = calloc(T[i].size, sizeof(counter*));

		current_st = 0;
		// For each suffix
		for (s = Lng.elem[i], k = 0; *s != '\0'; s++, k++)
		{
			// If it appears in significant state
			if (significant[current_st])
			{
				countTmp.word = s;
     			        assertion(countTmp.word != NULL, "err ctay 1");
				foundCountTmp = bsearch(&countTmp, suffCountsPtr[countTmp.word[0] - 'a'], 
					suffCountsMap[countTmp.word[0] - 'a'], sizeof(counter), compare_counters);
				switch(foundCountTmp->states[current_st])
				{
					case 0:
						// If cannot become a valid division point
						if ((current_st != 0 || foundCountTmp->count > 1) && 
							(long) (foundCountTmp->count * (long) RightLng[current_st].size) >= Lng.size)
						{
							actualRLS[current_st]++;
							T[i].states[k] = current_st;
							T[i].actualSize++;
							swPtr = malloc(sizeof(sw));
							swPtr->word_no = i;
							swPtr->next = states_to_words[current_st];
							states_to_words[current_st] = swPtr;
							foundCountTmp->states[current_st] = 1;
							T[i].suffixes[k] = foundCountTmp;
						}
						else
						{
							T[i].states[k] = PERM_REMOVED_STATE;
							foundCountTmp->states[current_st] = PERM_REMOVED_STATE;
						}		
						break;
					case PERM_REMOVED_STATE:
						T[i].states[k] = PERM_REMOVED_STATE;
						break;
					default:
						swPtr = malloc(sizeof(sw));
						swPtr->word_no = i;
						swPtr->next = states_to_words[current_st];
						states_to_words[current_st] = swPtr;
						T[i].states[k] = current_st;
						T[i].actualSize++;
						T[i].suffixes[k] = foundCountTmp;
						break;
				}
			}
			else
			{
				T[i].states[k] = PERM_REMOVED_STATE;
			}
			current_st = MADFA_delta[current_st][Sigma_map[(*s) - 'a']];
		}
		if (significant[current_st])
		{
			countTmp.word = s;
	 		assertion(countTmp.word != NULL, "err ctay 1");
			foundCountTmp = &suffCounts[0];
			switch(foundCountTmp->states[current_st])
			{
				case 0:
					// If cannot become a valid division point
					if ((current_st != 0 || foundCountTmp->count > 1) && 
						(long) (foundCountTmp->count * (long) RightLng[current_st].size) >= Lng.size)
					{
						actualRLS[current_st]++;
						T[i].actualSize++;
						swPtr = malloc(sizeof(sw));
						swPtr->word_no = i;
						swPtr->next = states_to_words[current_st];
						states_to_words[current_st] = swPtr;
						T[i].states[k] = current_st;
						foundCountTmp->states[current_st] = 1;
						T[i].suffixes[k] = foundCountTmp;
					}
					else
					{
						T[i].states[k] = PERM_REMOVED_STATE;
						foundCountTmp->states[current_st] = PERM_REMOVED_STATE;
					}		
					break;
				case PERM_REMOVED_STATE:
					T[i].states[k] = PERM_REMOVED_STATE;
					break;
				default:
					swPtr = malloc(sizeof(sw));
					swPtr->word_no = i;
					swPtr->next = states_to_words[current_st];
					states_to_words[current_st] = swPtr;
					T[i].actualSize++;
					T[i].states[k] = current_st;
					T[i].suffixes[k] = foundCountTmp;
					break;
			}
		}
		else
		{
			T[i].states[k] = PERM_REMOVED_STATE;
		}
		if(T[i].actualSize < mtTail->min)
		{
			mtTail->min = T[i].actualSize;
			mtTail->row = i;
		}
	}
}

int compare_counters(counter *f1, counter *f2)
{
	return compare_strings(&f1->word, &f2->word);
}

int check_P()
{
	WSS* words;
	WSS L1, L2;
	int i, prod_size_row;
	number_set *tmp_all_Ps;
	int sum_size = 0, max_prod_size = 0, ret_value = 0;
	double s_time, e_time;

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
	if ((long) (sum_size * (long) max_prod_size) < Lng.size)
	{
		return FALSE;
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
	s_time = MPI_Wtime();
	L2 = product(words, sP, max_prod_size, prod_size_row);
	e_time = MPI_Wtime();
	dec_prd_time += e_time - s_time;
	if ((long) (sum_size * (long) L2.size) < Lng.size)
	{
		free(L2.elem);
		free(words);
		return FALSE;
	}
	for (i = 0; i < sP; i++)
	{
		words[i] = LeftLng[P[i]];
	}
	s_time = MPI_Wtime();
	L1 = sum(words, sP, sum_size);
	e_time = MPI_Wtime();
	dec_sum_time += e_time - s_time;
	free(words);
	dec_checked++;
	ret_value = !(lambda_set(L1.elem) || lambda_set(L2.elem));
	if (ret_value)
	{
		s_time = MPI_Wtime();
		ret_value = catenation(L1, L2);
		e_time = MPI_Wtime();
		dec_cat_time += e_time - s_time;
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

/*
* Decompose - creates decompositions set
*/
void decompose()
{
	int sState;
	word suffix;
	word realSuffix;
	int **rem_ptr;
	int *size_ptr;
	mt *newMtTail;
	int new_size = 0;
	int marked_for_skip, removed;
	int i, j, k, l, min, list_no, r, row, tmp;
	double s_time, e_time;	

	level++;
	K = adaptK(sState_list - active_sS, K);
	/* step 1. - decomposition below cut-off level */
	if (sState_list - active_sS <= K)
	{
		phase1_aux_time = MPI_Wtime() - phase1_aux_time;
		phase1_time += phase1_aux_time;
		create_D();
		phase2_aux_time = MPI_Wtime();
		searching(D, sD);
		phase2_aux_time = MPI_Wtime() - phase2_aux_time;
		printf("Rt1t2: %d, %.5f, %.5f\n", rank, phase1_aux_time, phase2_aux_time);
		phase2_time += phase2_aux_time;
		level--;
		phase1_aux_time = MPI_Wtime();
		return;
	}
	/* step 2. - searching row of T with min. number of states */
	assertion(sT > 0, "Error 37");
	row = mtTail->row;
	min = mtTail->min;
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
	suffix = Lng.elem[T[row].word_no];
	T[row].word_no = ((T[row].word_no == 0) ? REMOVED_0_WORD : -T[row].word_no);
	/* Inserting states into S */
	for (r = 0; r < T[row].size; r++)  /* r is the index in array 'state' */
	{
		if(T[row].states[r] >= 0)
		{
			assertion(T[row].states[r] >= 0, "Error 60");
			list_no = ++S[T[row].states[r]].exists;				/* inserting current state */
			if (list_no == 1)
			{
				active_sS++;
			}
			assertion(list_no > 0, "Error 60a");
			assertion(list_no <= Lng.size, "Error 60b");
			/* actually used states in TPrime */
			s_time = MPI_Wtime();
			l = 0; /* index in removed states list */
			if (S[T[row].states[r]].cur_size <= list_no)
			{
				new_size = S[T[row].states[r]].cur_size + increment;
				rem_ptr = realloc(S[T[row].states[r]].rem, new_size * sizeof(number_set));
				assertion(rem_ptr != NULL, "Error 56");
				S[T[row].states[r]].rem = rem_ptr;
				for (j = S[T[row].states[r]].cur_size; j < new_size; j++)
				{
					S[T[row].states[r]].rem[j] = calloc(MADFA_q, sizeof(int));
					assertion(S[T[row].states[r]].rem[j] != NULL, "Error 56a");
				}
				free(S[T[row].states[r]].rem[0]);
				S[T[row].states[r]].rem[0] = NULL; /* 0th list always empty */
				new_size = S[T[row].states[r]].cur_size + increment;
				size_ptr = realloc(S[T[row].states[r]].size, new_size * sizeof(int));
				assertion(size_ptr != NULL, "Error 57");
				S[T[row].states[r]].size = size_ptr;
				S[T[row].states[r]].cur_size = new_size;
			}
			e_time = MPI_Wtime();
			mem_time += (e_time - s_time);
			removed = FALSE;
			marked_for_skip = FALSE;
			realSuffix = suffix + r;
			s_time = MPI_Wtime();
			newMtTail = malloc(sizeof(mt));
			newMtTail->min = min;
			newMtTail->row = -1;
			newMtTail->next = NULL;
			newMtTail->prev = mtTail;
			mtTail->next = newMtTail;			
			mtTail = newMtTail;
			if(!T[row].suffixes[r]->level)
			{
				T[row].suffixes[r]->level = level;	
				for (i = 0; i < MADFA_q; i++) /* state removal */
				{
					if(state_list[i] > 0 && !T[row].suffixes[r]->states[i])
					{
						removed = TRUE;
						remove_state(i);
						S[T[row].states[r]].rem[list_no][l++] = i;
						if (S[i].exists > 0)
						{
							marked_for_skip = TRUE;
							break;
						}
					}
				}
			}
			e_time = MPI_Wtime();
			remove1_time += (e_time - s_time);
			S[T[row].states[r]].size[list_no] = l;
			sState_list -= l;
			if (marked_for_skip)
			{
				s_time = MPI_Wtime();
				/* restoring removed states basing on TPrime, S, and state[r - 1] */
				tmp = T[row].states[r];
				restore_states(tmp);
				assertion(S[tmp].exists > 0, "Error 40");
				S[tmp].size[S[tmp].exists--] = 0;			/* list out-of-date, retracting previous state */
				if (S[tmp].exists == 0)
				{
					active_sS--;
				}
				assertion(S[tmp].exists >= 0, "Error 40a");
				restore_T_and_y_dec(level);
				mtTail = mtTail->prev;
				free(mtTail->next);
				if(T[row].suffixes[r]->level == level)
				{
					T[row].suffixes[r]->level = 0;
				}
				e_time = MPI_Wtime();
				restore_time += (e_time - s_time);
				continue;
			}
			// Something removed now or previously
			if (removed || changed)
			{	
				s_time = MPI_Wtime();
				create_T_and_y_dec(level);
				e_time = MPI_Wtime();
				remove2_time += (e_time - s_time);
			}
			if(mtTail->row == -1)
			{
				for(i = 0; i < sT; i++)
				{
					if(T[i].word_no >= 0 && T[i].actualSize == mtTail->min)
					{
						mtTail->row = i;
						break;
					}
				}				
			}
			if(mtTail->row == -1)
			{
				mtTail->min = MADFA_q;
				for(i = 0; i < sT; i++)
				{
					if(T[i].word_no >= 0 && T[i].actualSize < mtTail->min)
					{
						mtTail->row = i;
						mtTail->min = T[i].actualSize;
					}
				}				
			}
			decompose();
			s_time = MPI_Wtime();
			/* restoring removed states basing on TPrime, S, and state[r - 1] */
			tmp = T[row].states[r];
			restore_states(tmp);
			assertion(S[tmp].exists > 0, "Error 40");
			S[tmp].size[S[tmp].exists--] = 0;			/* list out-of-date, retracting previous state */
			if (S[tmp].exists == 0)
			{
				active_sS--;
			}
			assertion(S[tmp].exists >= 0, "Error 40a");
			restore_T_and_y_dec(level);
			mtTail = mtTail->prev;
			free(mtTail->next);
			if(T[row].suffixes[r]->level == level)
			{
				T[row].suffixes[r]->level = 0;
			}
			e_time = MPI_Wtime();
			restore_time += (e_time - s_time);
		}
	}
	T[row].word_no = ((T[row].word_no == REMOVED_0_WORD) ? 0 : -T[row].word_no);
	level--;
}

/*
* Removes given state from  the list of states
* @param state			- state to be removed
*/
void remove_state(int state)
{
	sw *swPtr;
	int i, j, val;
	counter search, *result;

	assertion(state >= 0, "Error 58");
	val = ((state > 0) ? -state : REMOVED_0_STATE);
	state_list[state] = -state_list[state];
	swPtr = states_to_words[state];
	while(swPtr != NULL)
	{
		i = swPtr->word_no;
		if (T[i].word_no >= 0)
		{
			for (j = 0; j < T[i].size; j++)
			{
				if (T[i].states[j] == state)
				{
					result = T[i].suffixes[j];
					// Decrease by 1
					result->count--;
					T[i].states[j] = val;
					T[i].actualSize--;
					if(T[i].actualSize < mtTail->min)
					{
						mtTail->min = T[i].actualSize;
						mtTail->row = i;
					}
					break;
				}
			}
		}
		swPtr = swPtr->next;
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
	LLS *sizes_lls = NULL, *sums_lls = NULL;
	int int_st, int_en, step, cur_pos, st_pos = -1;
	double s_time, e_time, s_time1, e_time1, t_time1 = 0.0;

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
	sP = sSx;
	printf("RS: %d; %d, %d, %lu\n", rank, sD, K, dec_checked);
	dec_to_check_all += pow(2.0, sD);
	dec_to_check_rall += pow(2.0, K);
	if (sSx != 0)
	{
		s_time = MPI_Wtime();
		verify_and_check_P(&cacheNo);
		e_time = MPI_Wtime();
		t_time1 += e_time - s_time;
	}
	for (i = 0, sum_lls = 0, min_rls = Lng.size; i < sP; i++)
	{
		sum_lls += LeftLng[P[i]].size;
		if (RightLng[P[i]].size < min_rls)
		{
			min_rls = RightLng[P[i]].size;
		}
	}
	if (min_rls == Lng.size) // means that only q0 is in sP or sP is empty
	{
		for (i = 0, min_rls = 0; i < sD; i++)
		{
			if (RightLng[D[i]].size > min_rls)
			{
				min_rls = RightLng[D[i]].size;
			}
		}
		if (min_rls == 0)
		{
			min_rls = Lng.size;
		}
	}
	expSum = ceil((double)Lng.size / min_rls) - sum_lls;
	if (sD != 0)
	{
		sums_lls = malloc(sD * sizeof(LLS));
		sizes_lls = malloc(sD * sizeof(LLS));
		for (i = 0; i < sD; i++)
		{
			sizes_lls[i].state_no = i;
			sizes_lls[i].lls = LeftLng[D[i]].size;
		}
		qsort(sizes_lls, sD, sizeof(LLS), compare_lls);
		sums_lls[0].lls = sizes_lls[0].lls;
		sums_lls[0].state_no = sizes_lls[0].state_no;
		for (i = 1; i < sD; i++)
		{
			sums_lls[i].state_no = sizes_lls[i].state_no;
			sums_lls[i].lls = sums_lls[i - 1].lls + sizes_lls[i].lls;
		}
		int_st = 0;
		int_en = sD;
		step = (int_en - int_st) >> 1;
		cur_pos = step;
		while (step != 0)
		{
			if (sums_lls[cur_pos].lls < expSum)
			{
				int_st = cur_pos;
				step = (int_en - int_st) >> 1;
				cur_pos += step;
			}
			else if (sums_lls[cur_pos].lls > expSum)
			{
				int_en = cur_pos;
				step = (int_en - int_st) >> 1;
				cur_pos -= step;
			}
			else
			{
				st_pos = cur_pos + 1;
				break;
			}
		}
		if (st_pos == -1)
		{
			if (int_st == 0)
			{
				for (i = 0; i < cur_pos; i++)
				{
					if (sums_lls[i].lls >= expSum)
					{
						st_pos = i + 1;
						break;
					}
				}
				if (st_pos == -1)
				{
					st_pos = cur_pos + 1;
				}
			}
			else if (int_en == sD)
			{
				for (i = cur_pos + 1; i < sD; i++)
				{
					if (sums_lls[i].lls >= expSum)
					{
						st_pos = i + 1;
						break;
					}
				}
				if (st_pos == -1)
				{
					st_pos = sD + 1;
				}
			}
			else
			{
				st_pos = cur_pos + 1;
			}
		}
	}
	else
	{
		st_pos = 1;
	}
	free(sums_lls);
	free(sizes_lls);
	s_time = MPI_Wtime();
	for (k = st_pos; k <= sD; k++)
	{
		m = sSx;
		for (i = 0; i < k; i++)
		{
			pos[i] = i;
			P[m++] = D[i];
		}
		sP = m;
		s_time1 = MPI_Wtime();
		verify_and_check_P(&cacheNo);
		e_time1 = MPI_Wtime();
		t_time1 += e_time1 - s_time1;
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
					s_time1 = MPI_Wtime();
					verify_and_check_P(&cacheNo);
					e_time1 = MPI_Wtime();
					t_time1 += e_time1 - s_time1;
					break;
				}
				else
				{
					j--;
				}
			}
		} while (j >= 0);
	}
	e_time = MPI_Wtime();
	save_time();
	printf("RE: %d; %d, %d, %lu\n", rank, sD, K, dec_checked);
	printf("RTime: %d; %5.2f, %5.2f\n", rank, t_time1, e_time - s_time);
	return FALSE;
}

/*
* Restores removed states from T basing on S and 'state'
* @param state			- state to be removed
*/
void restore_states(int state)
{
	sw *swPtr;
	int list_no;		/* list of states to be restored */
	int i, j, k, s;		/* state to be restored */
	int val, new_val;
	counter search, *result;

	assertion(state >= 0, "Error 59");
	list_no = S[state].exists;
	assertion(list_no > 0, "Error 59a");
	assertion(S[state].size[list_no] >= 0, "Error 59b");
	sState_list += S[state].size[list_no];
	for (k = 0; k < S[state].size[list_no]; k++)
	{
		s = S[state].rem[list_no][k];	/* state to recover */
		assertion(s >= 0, "Error 59c");
		state_list[s] = -state_list[s];
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
		swPtr = states_to_words[s];
		while(swPtr != NULL)
		{
			i = swPtr->word_no;
			if (T[i].word_no >= 0)
			{
				for (j = 0; j < T[i].size; j++)
				{
					if (T[i].states[j] == val)
					{
						T[i].states[j] = new_val;
						result = T[i].suffixes[j];
						// Increase by 1
						result->count++;
						T[i].actualSize++;
						break;
					}
				}
			}
			swPtr = swPtr->next;
		}
	}
}

int main(int argc, char **argv)
{
	int i;
	sw *swPtr;
	FILE *o_file, *proc_file;
	int proc_name_len;
	double total_time;
	char *line;
	char alphabet[MAX_A];
	char proc_name[MPI_MAX_PROCESSOR_NAME], buf[8192];
	double s_time, e_time;
	double s_time1, s_time2, e_time1, e_time2;

	init();
	initAdapt();

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Get_processor_name(proc_name, &proc_name_len);
	printf("Processor %d, %s\n", rank, proc_name);
	proc_file = fopen("/proc/self/stat", "r");
	fread(buf, sizeof(char), 8192, proc_file);
	fclose(proc_file);

	line = strtok(buf, " ");
	for (i = 1; i < 38; i++)
	{
		line = strtok(NULL, " ");
	}
	line = strtok(NULL, " ");
	printf("R: %d, P: %s, PID: %d\n", rank, proc_name, atoi(line));

	total_time = MPI_Wtime();
	start_time = MPI_Wtime();
	s_time = MPI_Wtime();
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

	KPlus1 = 2 * K;//K + 5;
	KPlus2 = 3 * K;//K + 15;
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
	increment = Lng.size / 5;
	rlsTail = malloc(sizeof(rls));
	rlsTail->prev = NULL;
	rlsTail->next = NULL;
	rlsTail->level = -1;
	rlsTail->counts = NULL;
	mtTail = malloc(sizeof(mt));
	mtTail->min = MADFA_q;
	mtTail->row = -1;
	mtTail->next = NULL;
	mtTail->prev = NULL;
	states_to_words = calloc(MADFA_q, sizeof(sw*));
	actualRLS = calloc(MADFA_q, sizeof(int));
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
	create_counts();
	create_ss();
	e_time = MPI_Wtime();
	meas_time_1 = MPI_Wtime();
	sState_list = 0;
	for(i = 0; i < MADFA_q; i++)
	{
		state_list[i] = LeftLng[i].size * actualRLS[i];
		if(state_list[i])
		{
			sState_list++;
		}
	}
	phase1_aux_time = MPI_Wtime();
	decompose();
	phase1_aux_time = MPI_Wtime() - phase1_aux_time;
	phase1_time += phase1_aux_time;
	printf("Rt1: %d, %.5f\n", rank, phase1_aux_time);
	meas_time_1 = MPI_Wtime() - meas_time_1;
	meas_time_2 = MPI_Wtime();
	s_time1 = MPI_Wtime();
	reduce_duplicate_decompositions();
	e_time1 = MPI_Wtime();
	s_time2 = MPI_Wtime();
	gather_all_decompositions();
	e_time2 = MPI_Wtime();
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
		swPtr = states_to_words[i];
		while(swPtr != NULL)
		{
			states_to_words[i] = swPtr->next;
			free(swPtr);
			swPtr = states_to_words[i];
		}
	}
	free(states_to_words);
	free(LeftLng);
	free(RightLng);
	remove_MADFA();
	remove_S(S, sS);
	remove_T();
	free(mtTail);
	free(rlsTail);
	remove_buffers();
	free(actualRLS);
	free(suffCounts);
	free(significant);
	fclose(test);
	total_time = MPI_Wtime() - total_time;
	if (o_file)
	{
		fprintf(o_file, "Rank: %d; execution times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2, phase1_time, phase2_time);
		fprintf(o_file, "Rank: %d; measurement times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, remove1_time, remove2_time, restore_time, mem_time, e_time - s_time, e_time1 - s_time1, e_time2 - s_time2, dec_sum_time, dec_prd_time, dec_cat_time);
		fprintf(o_file, "Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu, %d, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, statesRemoved, changedCount, maxK, dec_to_check_all, dec_to_check_rall);
		fclose(o_file);
	}
	else
	{
		printf("Rank: %d; execution times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, total_time, meas_time_1, meas_time_2, phase1_time, phase2_time);
		printf("Rank: %d; measurement times: %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f\n", rank, remove1_time, remove2_time, restore_time, mem_time, e_time - s_time, e_time1 - s_time1, e_time2 - s_time2, dec_sum_time, dec_prd_time, dec_cat_time);
		printf("Rank: %d; counts: %lu, %lu, %lu, %lu, %lu, %lu, %lu, %d, %lu, %lu\n", rank, pairs_prop, pairs_acc, dec_to_check, dec_to_check_rank, dec_checked, statesRemoved, changedCount, maxK, dec_to_check_all, dec_to_check_rall);
	}
	MPI_Finalize();
	return 0;
}
