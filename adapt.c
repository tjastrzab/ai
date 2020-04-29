#include "adapt.h"

void initAdapt()
{
	maxK = 0;
	prevCnt = -1; 
	prevConstCnt = 0;
	recCalls = REC_CALLS;
}

int adaptK(int diff, int K)
{
	if (prevCnt == -1)
	{
		prevCnt = diff;
	}
	else
	{
		if (prevCnt == diff)
		{
			prevConstCnt++;
			if (prevConstCnt == recCalls)
			{
				K++;
				prevConstCnt = 0;
			}
		}
		else
		{
			prevConstCnt = 0;
			prevCnt = diff;
		}
	}
	if (K == KPlus1 && recCalls == REC_CALLS)
	{
		recCalls += REC_CALLS;
	}
	else if (K == KPlus2 && recCalls == 2 * REC_CALLS)
	{
		recCalls += 2 * REC_CALLS;
	}
	if (K > maxK)
	{
		maxK = K;
	}
	return K;
}
