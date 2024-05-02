/**********************************************************
ChromSizesExt  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 05/02/2024
-------------------------
***********************************************************/

#pragma once

#include "Features.h"

// 'ChromSizesExt' provides functionality extention for the ChromSizes
class ChromSizesExt : public ChromSizes
{
	mutable chrid	_treatedCnt;	// number of treated chromosomes

	// Gets chrom's effective (treated) real length: a double length for autosomes, a single somatic
	//	@param it: ChromSizes iterator
	chrlen SetEffLength(cIter it) const { return Data(it).SetEffDefined(Chrom::IsAutosome(CID(it))); }

public:
	// Creates and initializes an instance.
	//	@param gName: reference genome directory
	//	@param printMsg: true if print message about chrom.sizes generation (in case of reference genome)
	ChromSizesExt(const char* gName, bool printMsg, const char* sPath)
		: _treatedCnt(0), ChromSizes(gName, printMsg, sPath, true) {}

	// Gets chrom's defined effective (treated) length
	//	@param it: ChromSizes iterator
	chrlen DefEffLength(cIter it) const;

	// Gets count of treated chromosomes.
	chrid TreatedCount() const { return _treatedCnt; }

	// Sets actually treated chromosomes according template and custom chrom
	//	@param all: if true then all chroms are pre-approved
	//	@param templ: template bed or NULL
	//	@returns: number of treated chromosomes
	chrid	SetTreatedChroms(bool all, const Features* const templ);

	// Prints threated chroms short names
	void	PrintTreatedChroms() const;
};
