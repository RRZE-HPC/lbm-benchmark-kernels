#ifndef __PADDING_H__
#define __PADDING_H__

typedef struct PadInfo_
{
	int nEntries;
	int * Modulus;
	int * Offset;
} PadInfo;

PadInfo * PadInfoFromStr(const char * padStr);

int PadCells(int nCells, int cellSizeBytes, PadInfo ** padInfo);
int PadCellsAndReport(int nCells, int cellSizeBytes, PadInfo ** padInfoPtr);

void PadInfoPrint(PadInfo * padInfo, FILE * f, const char * prefix);
void PadInfoFree(PadInfo * padInfo);

#endif // __PADDING_H__
