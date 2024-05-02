#include "ChromSizesExt.h"
#include "ChromSeq.h"

chrlen ChromSizesExt::DefEffLength(cIter it) const
{
	if (Data(it).Defined)	return Data(it).Defined;	// def.eff. length is initialized
	if (ChromSeq::LetGaps)	return SetEffLength(it);	// initialize def.eff. length by real size
	// initialize def.eff. length by chrN.region file
	ChromDefRegions rgns(RefName(CID(it)));
	if (rgns.Empty())		return SetEffLength(it);
	return Data(it).Defined = rgns.DefLength() << int(Chrom::IsAutosome(CID(it)));
}

chrid ChromSizesExt::SetTreatedChroms(bool all, const Features* const templ)
{
	// _treatedCnt is 0 after constructor

	if (Chrom::IsSetByUser())
		if (templ) {
			chrid userChr = Chrom::UserCID();
			if (templ->FindChrom(userChr))
				GetIter(userChr)->second.Treated = ++_treatedCnt;
		}
		else
			GetIter(Chrom::UserCID())->second.Treated = ++_treatedCnt;
	else
		for (Iter it = Begin(); it != End(); it++)
			_treatedCnt += 
				it->second.Treated = all || (templ && templ->FindChrom(CID(it)));

	return _treatedCnt;
}

void ChromSizesExt::PrintTreatedChroms() const
{
	if (TreatedCount() == ChromCount()) {
		cout << " all";
		return;
	}
	/*
	* sequential IDs printed as range: <first-inrange>'-'<last in range>
	* detached IDs or ranges are separated by comma
	*/
	chrid cID = 0, cIDlast = 0;		// current cid, last printed cid
	chrid unprintedCnt = 0;
	bool prFirst = true;	// true if first chrom in range is printed
	cIter itLast;
	auto getSep = [&unprintedCnt](chrid cnt) { return unprintedCnt > cnt ? '-' : COMMA; };
	auto printChrom = [](char sep, chrid cID) { dout << sep << Chrom::Mark(cID); };
	auto printNextRange = [&](chrid lim, chrid nextcID) {
		if (cID != cIDlast) printChrom(getSep(lim), cID);
		printChrom(COMMA, nextcID);
	};

	//== define last treated it
	for (cIter it = cBegin(); it != cEnd(); it++)
		if (IsTreated(it))	itLast = it;

	//== print treated chrom
	for (cIter it = cBegin(); it != itLast; it++) {
		if (!IsTreated(it))		continue;
		if (prFirst)
			printChrom(SPACE, cIDlast = CID(it)),
			prFirst = false;
		else
			if (CID(it) - cID > 1) {
				printNextRange(1, CID(it));
				cIDlast = CID(it);
				unprintedCnt = 0;
			}
			else
				unprintedCnt++;
		cID = CID(it);
	}

	// print last treated chrom
	cIDlast = CID(itLast);
	if (!prFirst && cIDlast - cID > 1)
		printNextRange(0, cIDlast);
	else
		printChrom(prFirst ? SPACE : getSep(0), cIDlast);
}
