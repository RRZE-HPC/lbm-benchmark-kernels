#include "Base.h"
#include "Kernel.h"
#include "Padding.h"

#include <errno.h>
#include <limits.h>

// Generates a PadInfo struct from padStr. If padStr is NULL then a NULL
// pointer is returned.
//
// Padding string is either "auto" or in the format of
// <modulus>+<remainder>(,<modulus>+<remainder>)*.
//
// If modulus is 0 then the remainder will just be added and no modulo will be
// applied.
//
// Returned PadInfo structs must be freed with PadInfoFree.

PadInfo * PadInfoFromStr(const char * padStr)
{
	if (padStr == NULL) {
		return NULL;
	}
	else if (!strncmp("auto", padStr, 4)) {

		PadInfo * padInfo = (PadInfo *)malloc(sizeof(PadInfo) + 4 * sizeof(int));
		Assert(padInfo != NULL);

		padInfo->nEntries = 2;
		padInfo->Modulus = (int *)(padInfo + 1);
		padInfo->Offset = ((int *)(padInfo + 1)) + 2;

		// Intel TLB 2 MiB pages
		padInfo->Modulus[0] = 2 * 1024 * 1024 * 8;
		padInfo->Offset[0]  = 2 * 1024 * 1024 + 512;

		// Intel L2, 256 sets
		padInfo->Modulus[1] = 256 * 64;
		padInfo->Offset[1]  = 128;

		return padInfo;
	}
	else if (!strncmp("no", padStr, 2) || !strncmp("off", padStr, 3)) {
		PadInfo * padInfo = (PadInfo *)malloc(sizeof(PadInfo));
		Assert(padInfo != NULL);

		padInfo->nEntries = 0;
		padInfo->Modulus  = NULL;
		padInfo->Offset   = NULL;

		return padInfo;
	}

	// Count number of commas. Number of commas + 1 gives us the number of
	// padding entries.
	const char * tmp = padStr;
	int nCommas = 0;

	while (*tmp != 0x00) {
		if (*tmp == ',') ++nCommas;
		++tmp;
	}

	// Number of padding entries = nCommas + 1.
	size_t padInfoSize = sizeof(PadInfo) + (nCommas + 1) * 2 * sizeof(int);
	PadInfo * padInfo = (PadInfo *)malloc(padInfoSize);
	Assert(padInfo != NULL);

	memset(padInfo, 0x00, padInfoSize);

	padInfo->Modulus = (int *)(padInfo + 1);
	padInfo->Offset = ((int *)(padInfo + 1)) + nCommas + 1;

	tmp = padStr;

	int padEntryIndex = 0;
	int entriesSeen   = 0;

	int modulusSeen   = 0;
	int offsetSeen    = 0;

	// TODO: parsing is currently a mess. We assume a fixed format, which
	// should be easier to parse now.
	while (*tmp != 0x00) {

		if (*tmp == ',') {
			// Check if modulus was seen in this section, then move to the next.

			if (!modulusSeen) {
				fprintf(stderr, "ERROR: modulus missing before next section in padding string.\n");
				exit(1);
			}

			modulusSeen = 0;
			offsetSeen  = 0;
			++padEntryIndex;
			++tmp;
		}
		else if (*tmp == '+') {
			// Check if modulus was seen in this section and no offset yet.

			if (!modulusSeen) {
				fprintf(stderr, "ERROR: modulus missing before offset in padding string.\n");
				exit(1);
			}

			if (offsetSeen) {
				fprintf(stderr, "ERROR: offset is only allowed to be specified once per section in padding string.\n");
				exit(1);
			}

			++tmp;
		}
		else if (*tmp >= '0' && *tmp <= '9') {
			// Parse number and check for all the errors.

			char * endPtr;

			errno = 0;

			int value = strtol(tmp, &endPtr, 10);

			if ((value == LONG_MIN || value == LONG_MAX) && errno != 0) {
				fprintf(stderr, "ERROR: over- or underflow in modulus/offset of padding string.\n");
				exit(1);
			}
			else if (value == 0 && errno != 0) {
				fprintf(stderr, "ERROR: error parsing modulus/offset of padding string: %d - %s\n", errno, strerror(errno));
				exit(1);
			}
			else if (tmp == endPtr) {
				// No digits found: empty string or invalid pad string.
				fprintf(stderr, "ERROR: modulus or offset missing in padding string.\n");
				exit(1);
			}

			if (value < 0) {
				fprintf(stderr, "ERROR: modulus and offset must be >= 0 in padding string.\n");
				exit(1);
			}

			if (!modulusSeen) {
				if (offsetSeen) {
					fprintf(stderr, "ERROR: offset already seen, but not modulus in padding string.\n");
					exit(1);
				}

				Verify(padEntryIndex < (nCommas + 1));

				padInfo->Modulus[padEntryIndex] = value;

				modulusSeen = 1;
				++entriesSeen;
			}
			else {
				if (offsetSeen) {
					fprintf(stderr, "ERROR: offset is only allowed to be specified once per section in padding string.\n");
					exit(1);
				}

				Verify(padEntryIndex < (nCommas + 1));

				padInfo->Offset[padEntryIndex] = value;

				offsetSeen = 1;
			}
			// No increment of tmp needed, as endPtr points already to the first
			// character which is not part of the number.
			tmp = endPtr;
		}
		else {
			fprintf(stderr, "ERROR: padding string contains invalid characters.\n");
			exit(1);
		}

	}

	if (entriesSeen != nCommas + 1) {
		fprintf(stderr, "ERROR: did not find all padding entries.\n");
		exit(1);
	}

	for (int i = 0; i < nCommas + 1; ++i) {
		if (padInfo->Offset[i] >= padInfo->Modulus[i]) {
			fprintf(stderr, "ERROR: offset in padding entry %d is equal or larger than modulus.\n", i);
			exit(1);
		}
	}

	padInfo->nEntries = entriesSeen;

	return padInfo;
}

int PadCells(int nCells, int cellSizeBytes, PadInfo ** padInfoPtr)
{
	Assert(padInfoPtr != NULL);

	// If the padInfo is NULL determine if we use the "auto" configuration.
	const int defaultAutoPadding = 1;


	if (*padInfoPtr == NULL) {
		if (!defaultAutoPadding) {
			return nCells;
		}

		*padInfoPtr = PadInfoFromStr("auto");
	}

	PadInfo * padInfo = *padInfoPtr;

	Assert(padInfo->nEntries >= 0);

	for (int i = 0; i < padInfo->nEntries; ++i) {
		Assert(padInfo->Modulus[i] >=  0);
		Assert(padInfo->Offset[i]  >= 0);

		int nPadCells = 0;

		if (padInfo->Modulus[i] == 0) {
			// When the modulus is just zero then we only add the offset.
			nPadCells = padInfo->Offset[i] * cellSizeBytes;
		}
		else {
			int nModCells    = padInfo->Modulus[i] / cellSizeBytes;

			if (nModCells == 0) {
				fprintf(stderr, "ERROR: modulus of %d byte in padding entry %d becomes zero for PDF size %d byte.\n",
						padInfo->Modulus[i], i, cellSizeBytes);
				exit(1);
			}

			int nOffsetCells = padInfo->Offset[i]  / cellSizeBytes;
			int nRemainder   = nCells % nModCells;

			nPadCells = (nOffsetCells + nModCells - nRemainder) % nModCells;
		}

		nCells += nPadCells;
	}

	return nCells;
}

void PadInfoPrint(PadInfo * padInfo, FILE * f, const char * prefix)
{
	Assert(padInfo != NULL);
	Assert(padInfo->nEntries >= 0);
	Assert(f != NULL);
	Assert(prefix != NULL);

	for (int i = 0; i < padInfo->nEntries; ++i) {
		fprintf(f, "%sm: %10d b  o: %10d b   m: %f KiB  o: %f KiB   m: %f MiB  o: %f MiB\n",
			prefix,
			padInfo->Modulus[i],                   padInfo->Offset[i],
			padInfo->Modulus[i] / 1024.0,          padInfo->Offset[i] / 1024.0,
			padInfo->Modulus[i] / 1024.0 / 1024.0, padInfo->Offset[i] / 1024.0 / 1024.0);
	}

	return;
}

void PadInfoFree(PadInfo * padInfo)
{
	if (padInfo != NULL) {
		free(padInfo);
	}

	return;
}

// Returns the new padded cell size and reports if padding was performed and how.

int PadCellsAndReport(int nCells, int cellSizeBytes, PadInfo ** padInfoPtr)
{
	// Apply padding.
	int nNewCells = PadCells(nCells, sizeof(PdfT), padInfoPtr);

	if (nCells != nNewCells) {
		printf("# padding info:\n");
		PadInfoPrint(*padInfoPtr, stdout, "#                        ");
		// int nPaddedCells = nNewCells - nCells;
		// printf("# padding per dir.:      %10d nodes (%f MiB)\n", nPaddedCells, nPaddedCells / 1024.0 / 1024.0 * sizeof(PdfT));
		// printf("# padding total:         %10d nodes (%f MiB)\n", 19 * nPaddedCells, 19 * nPaddedCells / 1024.0 / 1024.0 * sizeof(PdfT));

		int nPadCells = nNewCells - nCells;

		printf("#\n# padding %d nodes with %d nodes, %f MiB per direction, %f MiB in total\n",
			nCells,
			nPadCells,
			nPadCells * sizeof(PdfT) / 1024.0 / 1024.0,
			nPadCells * sizeof(PdfT) / 1024.0 / 1024.0 * 19);
	}
	else {
		printf("# padding info:          no padding used\n");
	}

	return nNewCells;
}
