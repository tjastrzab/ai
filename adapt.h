#ifndef _ADAPT_H_
#define _ADAPT_H_

#define REC_CALLS 150

int KPlus1, KPlus2, maxK;
int prevCnt, prevConstCnt, recCalls;

void initAdapt();
int adaptK(int diff, int K);

#endif
