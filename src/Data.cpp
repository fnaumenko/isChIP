#include "Data.h"
#ifdef _ISCHIP
	#include "isChIP.h"
#elif defined _BIOCC
	#include "Calc.h"
#endif	//_BIOCC


/************************ class Bed ************************/

static const char* Per = " per ";

Bed::Ambig::Msg Bed::Ambig::_Msgs [] = {
	{ NULL, "duplicated",	"duplicated", },
	{ NULL, "crossed",		"is intersected with previous" },
	{ NULL, "adjacent",		"is adjacent with previous" },
	{ NULL, "covered",		"is fully covered by previous" },
	{ NULL, "too short",	"too short" },
	{ NULL, "different size of", "different size of read" },
	{ NULL, "filtered by low score", "filtered by score" },
	{ NULL, "negligible",	"negligible chromosome" },
};

const char*	Bed::Ambig::_ActionMsgs[] = {
	"accepted",
	"joined",
	"skipped",
	"execution aborted"
};

const Bed::Ambig::ReportCase Bed::Ambig::_Actions[] = {
	&Bed::Ambig::Accept,
	&Bed::Ambig::Handle,
	&Bed::Ambig::Omit,
	&Bed::Ambig::OmitQuiet,
	&Bed::Ambig::Abort
};

// Prints entity name(s)
//	@cnt: number of entity
void Bed::Ambig::PrintEntityName(chrlen cnt) const
{
	dout << _entityName;
	if(cnt>1)	dout <<  's';
}

// Prints entities count
//	@cID: chromosome's ID or Chrom::UnID if all
//	@cnt: number of entity
void Bed::Ambig::PrintEntityCount(chrid cID, chrlen cnt) const
{
	dout << cnt << BLANK;
	PrintEntityName(cnt);
	if( cID != Chrom::UnID )	
		dout << Per << Chrom::TitleName(cID);
	//_infoPrinted = true;
}

// Creates an instance with omitted COVER, SHORT, SCORE and NEGL cases,
// and BedF cases by default:
// omitted DUPL cases and handled CROSS and ADJAC 
Bed::Ambig::Ambig (bool alarm, bool printFileName, const string& entityName,
	eAction dupl, eAction crossANDadjac, eAction diffsz) :
	_alarm(alarm),
	_printFileName(printFileName),
	_unsortedItems(false),
	_firstAlarm(true),
	_infoPrinted(false),
	_count(CHRLEN_UNDEF),
	_entityName(entityName),
	_file(NULL)
{
	memset(_cases, 0, _CasesCnt*sizeof(Case));	// initialize by 0
	_cases[DUPL].Type = dupl;
	_cases[CROSS].Type = _cases[ADJAC].Type = crossANDadjac;
	_cases[COVER].Type = _cases[SHORT].Type = _cases[SCORE].Type = OMIT;
	_cases[DIFFSZ].Type = diffsz;
}

// Print given ambiguity as alarm
//	@ambig: given ambiguity
void Bed::Ambig::PrintLineAlarm(eCase ambig) const
{
	if( _alarm ) {
		_infoPrinted = true;
		if( _firstAlarm )	{ dout << EOL;	_firstAlarm = false; }
		Err(_Msgs[ambig].LineAlarm + BLANK + _entityName + MSGSEP_BLANK, FileRecordNumb()).
			Warning(Message(ambig));
	}
}

// Outputs case statistics
//	@ambig: ambiguity's case
//	@allCnt: total count of ambiguities
//	@total: if true then prints total warning case
void Bed::Ambig::PrintCaseStat(eCase ambig, chrlen allCnt, bool total) const
{
	//chrlen cnt = _stats[ambig];	// count of case's ambiguities
	chrlen cnt = _cases[ambig].Count;	// count of case's ambiguities
	if( !cnt )	return;
	const char* totalAlarm = _Msgs[ambig].TotalAlarm;
	if(total)	dout << WARNING;
	else		dout << TAB;
	dout<< cnt
		<< sPercent(ULLONG(cnt), ULLONG(allCnt), 4) << BLANK
		<< _Msgs[ambig].StatInfo << BLANK;
	PrintEntityName(cnt);
	if(UnsortedItems())		dout << " arisen after sorting";
	if(totalAlarm)	dout << BLANK << totalAlarm;
	dout << "; " << Message(ambig);
	if(totalAlarm)	dout << '!';
	dout << EOL;
}

// Outputs statistics.
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@title: string at the beginning
//	@printItemCnt: if true then print count of items
//	@itemCnt: count of all items (features/reads) from file
//	@savedItemCnt: count of saved items (features/reads)
void Bed::Ambig::PrintStatistics(
	chrid cID, const char* title, bool printItemCnt, chrlen itemCnt, chrlen savedItemCnt) 
{
	if( title )	{
		dout << title << MSGSEP_BLANK;
		_infoPrinted = true;
	}
	if( printItemCnt )
		PrintEntityCount(cID, itemCnt);
	if( AmbigCount() ) {
		if( printItemCnt )	dout << ", from which\n";
		for(BYTE i=0; i<_CasesCnt; i++)
			PrintCaseStat(static_cast<eCase>(i), itemCnt);

		// print ambiguities of negligible chroms if all chroms are readed
		if( cID == Chrom::UnID ) {
			// calculate ambigs for negligible chroms:
			// rest of difference between in & out features minus count of ambigs
			_cases[NEGL].Count = itemCnt - savedItemCnt - AmbigCount();
			// correct (add) accepted features
			for(BYTE i=0; i<_CasesCnt-1; i++)	// loop excepting NEGL chroms case
				if( _cases[i].Type == Ambig::ACCEPT )	_cases[NEGL].Count += _cases[i].Count;
			
			PrintCaseStat(NEGL, itemCnt);
		}
		// print total remained entities
		dout<< TAB << Total << " remained: " << savedItemCnt
			<< sPercent(ULLONG(savedItemCnt), ULLONG(itemCnt), 4) << BLANK;
		PrintEntityName(savedItemCnt);
		_infoPrinted = true;
	}
	if( _infoPrinted )	dout << EOL;
	fflush(stdout);
}

// Prints total warning if some ambiguous are occurs; needs for silent and laconic modes
//	@printFull: true if full info should be printed; otherwise count of items only
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@printItemCnt: if true then print count of items
//	@itemCnt: count of all items (features/reads) from file
void Bed::Ambig::PrintInfo(bool printFull, chrid cID, bool printItemCnt, chrlen itemCnt) const
{
	if(printItemCnt)		PrintEntityCount(cID, itemCnt);
	if(!printFull || !AmbigCount())		return;
	if(printItemCnt)		dout << EOL;
	for(BYTE i=0; i<_CasesCnt; i++)
		PrintCaseStat(static_cast<eCase>(i), itemCnt, true);
	_infoPrinted = true;
	fflush(stdout);
}

// Adds statistics and print given ambiguity as alarm (if permitted)
//	@ambig: given ambiguity
//	return: treatment code: 1 - accept, 0 - handle, -1 - omit
int Bed::Ambig::TreatCase(eCase ambig)
{
	_cases[ambig].Count++;
	return (this->*_Actions[_cases[ambig].Type])(ambig);
}

// Gets count of ambiguities
chrlen Bed::Ambig::AmbigCount() const {
	if( _count == CHRLEN_UNDEF ) {
		_count = 0;
		for(BYTE i=0; i<_CasesCnt; i++)	_count += _cases[i].Count;
	}
	return _count;
}

// Initializes Region with checking
static void InitRegn(Region& rgn, long start, long end, long cLen, const Bed::Ambig& ambig)
{
	if(start < 0 || end < 0)	ambig.ThrowExcept(Err::BP_NEGPOS);
	if(start >= end)			ambig.ThrowExcept(Err::BP_BADEND);
	if(cLen > 0 && end > cLen)	ambig.ThrowExcept(Err::BP_EXCEED);
	rgn.Init(start, end);
}

// Initializes new instance by bed-file name.
//	@fName: name of bed-file
//	@ambig: ambiguous filters
//	@cSizes: chrom sizes to control the chrom length exceedeng; if NULL, no control
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
//	@info: type of feature ambiguties that should be printed
//	@printItemCnt: if true then print count of items
//	@return: count of features/reads in file
ULONG Bed::Init	(const string& fName, Ambig& ambig, const ChromSizes* cSizes,
	bool abortInvalid, eInfo info, bool printItemCnt)
{
	ULONG	cntLines = 0,		// count of lines beginning with 'chr'
			currInd = 0,		// current index in Feature's/Read's container.
								// Needed to avoid excess virtual method given current container size
			initFeatrsSize;		// initial size of _items
	bool	cAll = Chrom::StatedAll();	// true if all chroms are stated by user
	TabFile file(fName, _FieldsCnt, abortInvalid, TxtFile::READ, Chrom::Abbr, TabFile::Comment, false);
	Timer	timer;

	if( ambig.PrintFileName() )	{ dout << fName << MSGSEP_BLANK; fflush(stdout); }
	try {
		if( file.IsBad() )
			Err(file.ErrCode(), ambig.PrintFileName() ? sBACK : fName).Throw();
		ambig.SetFile(file);
		if( !file.GetFirstLine(&initFeatrsSize) )
			if( file.ErrCode() == Err::NONE )	// file is right, but no any features
				Err(Err::TF_EMPTY, ambig.PrintFileName() ? sBACK : fName, ItemTitle(true)).Throw();
			else								// wrong first line
				ambig.ThrowExcept(file.ErrCode());
		
		ULONG 	firstInd = 0;	// first index in feature's container for current chromosome
		Region	rgn;			// current feature positions
		chrlen	prevStart=0,	// start previous feature positions
				cLen = 0;		// the length of chromosome
		chrid	cCurrID,		// current chromosome's ID
				cNextID;		// next chromosome's ID
		bool	needSortChrom = false;
		char cName[Chrom::MaxShortNameLength + 1];	// chrom number (only) buffer

		Reserve(cAll ? Chrom::Count : 1);
		ReserveItemContainer(initFeatrsSize);
		cName[0] = 0;	// empty initial chrom name
		cCurrID = Chrom::ID(file.ChromName());		// first chrom ID
		if(cSizes)	cLen = cSizes->Size(cCurrID);

		do {
			if( strcmp(cName, file.ChromName()) ) {	// next chromosome?
				// chrom's name may be long such as 'chrY_random'
				cNextID = Chrom::ID(file.ChromName());
				if( cNextID == Chrom::UnID ) {		// negligible next chromosome?
					cCurrID = cNextID;
 					continue;
				}
				// now we can save this chrom name
				strcpy(cName, file.ChromName());
				if( cAll ) {			// are all chroms as parameter defined?
					if( currInd != firstInd	)		// have been features for this chrom saved?
						// save current chrom which features have been saved already
						AddVal(cCurrID, ChromItemsInd(firstInd, currInd));
					if( cNextID < cCurrID && cNextID != Chrom::M )	// unsorted chrom?
						needSortChrom = true;
				}
				else {		// single chromosome as parameter is defined
					// features in a single defined chrom were saved;
					// the chrom proper will be saved after loop
					if(rgn.End)			
						break;
					if(needSortChrom)
						Err("is unsorted. Option --chr " + Chrom::Name(Chrom::StatedID())
						+ " is forbidden", fName).
							Throw();
					if(cNextID != Chrom::StatedID())	continue;
				}
				cCurrID = cNextID;
				firstInd = currInd;
				if(cSizes)	cLen = cSizes->Size(cCurrID);
				InitRegn(rgn, file.LongField(1), file.LongField(2), cLen, ambig);
				cntLines++;
			}
			else { 
				if( cCurrID == Chrom::UnID )				// negligible chrom?
					continue;
				if( !cAll && cCurrID != Chrom::StatedID() )	// undefined chrom?
					continue;
				InitRegn(rgn, file.LongField(1), file.LongField(2), cLen, ambig);
				cntLines++;
				if( rgn.Start < prevStart )			// unsorted feature?
					ambig.ItemsAreUnsorted();
				// handle positions for the same chromosome only
				if( !CheckLastPos(rgn, ambig) )
					continue;
			}
			if( AddPos(rgn,	file.StrField(3), file.FloatField(4), file.StrField(5)) )
				currInd++;
			else
				ambig.TreatCase(ambig.SCORE);
			prevStart = rgn.Start;
		}
		while( file.GetLine() );

		if( file.IsBad() )	// line may be NULL because of wrong
			ambig.ThrowExcept(file.ErrCode());

		if( rgn.End && currInd ) {		// some features for as minimum one valid chrom were saved
			if( cCurrID != Chrom::UnID )	// valid last chrom?
				// save last chrom. Its features have been saved already.
				AddVal(cCurrID, ChromItemsInd(firstInd, currInd));
			if( initFeatrsSize/currInd > 2 )	
				ShrinkItemContainer();
			if( needSortChrom )
				Sort();
			if( ambig.UnsortedItems() )
				SortItems(ambig);
		}
	}
	catch(Err &err) { ThrowError(err, abortInvalid); }
	SetAllItemsCount();
	if( !_isBad ) {
		if(info == iSTAT)		// print statistics
			ambig.PrintStatistics(Chrom::StatedID(), 
				ambig.PrintFileName() ? NULL : fName.c_str(),	// fname is printed yet or not
				printItemCnt, cntLines, AllItemsCount());
		else
			ambig.PrintInfo(info == iNOTE, Chrom::StatedID(), printItemCnt, cntLines);
		_EOLPrinted = ambig.IsAddInfoPrinted();
	}
	if(Timer::Enabled && !_EOLPrinted)	dout << BLANK;
	timer.Stop(true, _EOLPrinted);
	if(!_EOLPrinted && (printItemCnt || ambig.PrintFileName())) {
		dout << EOL; 
		_EOLPrinted = true;
		//fflush(stdout);	
	}
	if( !_isBad && !currInd ) {
		string specify = ItemTitle(true);
		if(!cAll)	specify += Per + Chrom::TitleName(Chrom::StatedID());
		Err err(Err::TF_EMPTY, fName, specify);
		ThrowError(err, abortInvalid);
	}
	return cntLines;
}

//#ifdef _DENPRO
//// Shifts item's positions to collaps the 'holes' between regions.
////	@cID: chromosome's ID
////	@regns: valid (defined) regions
//void Bed::ShrinkByID(chrid cID, const Regions &regns)
//{
//	chrlen	rgCnt = regns.Count(),
//			iCnt = At(cID).ItemsCount(),	// count of reads/features
//			shift = 0,	// shift start position to left (accumulative regions start positions)
//			rgEnd,		// current region's stop position
//			k=0;		// index of reads/features
//
//	for(chrlen i=0; i<rgCnt; i++) {
//		shift += regns[i].Start;
//		rgEnd = regns[i].End;
//		for(; k<iCnt; k++)
//			if( !DecreasePos(cID, k, shift, rgEnd) )
//				break;
//	}
//	if(k > iCnt)
//		Err("item outside last region", "Bed::Shrink()").Throw();
//}
//#endif	// _DENPRO

#ifdef DEBUG
void Bed::PrintChrom() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::AbbrName(CID(it)) << TAB
		<< it->second.FirstInd << TAB << it->second.LastInd << MSGSEP_TAB
		<< it->second.ItemsCount() << TAB << ItemTitle() << 's' << EOL;
}

void Bed::Print(chrlen itemCnt) const
{
	chrlen i, iCnt;
	cout << "Bed's ";
	if( itemCnt )	cout << "first " << itemCnt << BLANK;
	cout << ItemTitle() << "s:\n";
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		iCnt = itemCnt ?
			(itemCnt > ItemsCount(CID(it)) ? ItemsCount(CID(it)) : it->second.FirstInd + itemCnt) - 1:
			it->second.LastInd;
		for(i=it->second.FirstInd; i<=iCnt; i++) {
			cout << Chrom::AbbrName(CID(it)) << TAB;
			PrintItem(i);
		}
	}
}
#endif

/************************ end of class Bed ************************/

#if !defined _ISCHIP && !defined _WIGREG

/************************ class BedR ************************/
const string BedR::_ItemTitle = "read";
const string BedR::_ItemTitles = BedR::_ItemTitle + 's';

// Checks the element for the new potential start/end positions for all possible ambiguous.
//	@it: iterator reffering to the compared element
//	@rgn: checked start/stop positions
//	@ambig: possible ambiguities: different Read length, duplicated Reads
//  return: true if item should be accepted; otherwise false
bool BedR::CheckItemPos(ItemsIter it, const Region& rgn, Ambig& ambig)
{
	if( !_readLen )		_readLen = rgn.Length() - 1;	// initialize Read length once
	else if( _readLen != rgn.Length() - 1				// different Read length?
	&& ambig.TreatCase(ambig.DIFFSZ) < 0)	return false;

	if( rgn.Start == it->Pos							// duplicated Read?
	&& ambig.TreatCase(ambig.DUPL) < 0 )	return false;
	// adjacent & crossed Reads are not checked: it's a normal cases
	// covered Reads are not checked: it's impossible case
	return true;
}

const string NotStated = " is not stated";

// Adds Read to the container.
//	@rgn: start/end Read positions
//	@rName: Read name field
//	@score: Read score
//	@strand: Read strand
//	return: true if Read was added successfully
bool BedR::AddPos(const Region& rgn, const char* rName, float score, const char* strand) {
	if( score <= _minScore )	return false;				// pass Read with under-threshhold score
#ifdef _VALIGN
	if( _rNameType == Read::nmUndef ) {	// first Read: define type of Read's name
		const char* sNumVal = strchr(rName, *Read::NmNumbDelimiter);
		if( sNumVal ) {
			if( strchr(rName, Read::NmSuffMate1[0]) )	_paired = true;
			_rNameType = *(sNumVal+1) == *(Read::NmNumbDelimiter+1) ?
				Read::nmNumb : Read::nmPos;
		}
		else	Err(Err::BR_RNAME, StrEmpty, "delimiter COLON" + NotStated).Throw();
	}
	if( !(rName = Chrom::FindNumb(rName)) )
		Err(Err::BR_RNAME, StrEmpty, Chrom::Title + NotStated).Throw();
	chrid cID = Chrom::ID(rName);
	// pass chrom number with unknown capacity and possible Number Delimiter
	rName = strchr(rName+1, *Read::NmNumbDelimiter) + 1;	
	if( _paired && _rNameType == Read::nmPos && *strand == Read::Strand[1]) {
		const char* posDelim = strchr(rName, Read::NmPosDelimiter);
		if( posDelim )		rName = posDelim + 1;
		else
			Err(Err::BR_RNAME, StrEmpty, "paired-end delimiter '-'" + NotStated).Throw();
	}
	_items.push_back(Read(rgn.Start, cID, size_t(atol(rName)), readscr(score)));
#else
	_items.push_back(Read(rgn.Start));
#endif	// _VALIGN
	if(score > _maxScore)	_maxScore = score;
	return true;
}

// Decreases Read's start position without checkup indexes.
//	@cID: chromosome's ID
//	@rInd: index of read
//	@shift: decrease read's start position on this value
//	@rgEnd: region's end position to control Read location
//	@return: true if Read is insinde this region
//bool BedR::DecreasePos(chrid cID, chrlen rInd, chrlen shift, chrlen rgEnd)
//{
//	chrlen start = Item(cID, rInd);	// start position of Read
//	//if(start < shift) {
//	//	string msgFileLine = "read " + IntToStr(start) + " outside region";
//	//	Err(msgFileLine.c_str(), "Region end="+IntToStr(rgEnd)).Throw();
//	//}
//	if(start + _readLen > rgEnd)
//		return false;
//	_items[At(cID).FirstInd + rInd] -= shift;
//	return true;
//}

/************************ end of class BedR ************************/

#endif	// !_ISCHIP && !_WIGREG

#ifndef _VALIGN
/************************ class BedF ************************/
const string BedF::_ItemTitle = "feature";
const string BedF::_ItemTitles = BedF::_ItemTitle + 's';

// Sets new end position on the feature if necessary.
//	@it: iterator reffering to the feature which end position may be corrected
//	@end: potential new end position
//	@treatCaseRes: result of treatment this ambiguity
//	return: true if ambiguity is permitted (feature is valid)
bool BedF::CorrectItemsEnd(ItemsIter it, chrlen end, int treatCaseRes) {
	if(treatCaseRes > 0)	return true;
	if(!treatCaseRes)				// treatment: merge or join features
		it->End = end;
	return false;
}

// Checks the element for the new potential start/end positions for all possible ambiguous.
// * duplicated features
// * short feature
// * adjacent features
// * crossed features
//	@it: iterator reffering to the compared element
//	@rgn: checked start/stop positions
//	@ambig: possible ambiguities
//  return: true if item should be accepted; otherwise false
bool BedF::CheckItemPos(ItemsIter it, const Region& rgn, Ambig& ambig)
{
	// features should be sorted!
#ifdef _BIOCC
	if( _uniformFtrLen )
		if(!_ftrLen)	_ftrLen = rgn.Length();	// initialize feature's length once
		else 			_uniformFtrLen = (_ftrLen == rgn.Length());
#endif
	Region currRgn = *it;
	if(rgn == currRgn)					// duplicated feature?
		return ambig.TreatCase(ambig.DUPL) >= 0;
#ifdef _ISCHIP
	if(rgn.Length() < _minFtrLen)		// short feature?
		return ambig.TreatCase(ambig.SHORT) >= 0;
#endif
	if(currRgn.Adjoin(rgn))				// adjacent feature?
		return CorrectItemsEnd(it, rgn.End, ambig.TreatCase(ambig.ADJAC));
	if(currRgn.Cover(rgn))				// covering feature?
		return ambig.TreatCase(ambig.COVER) >= 0;
	if(currRgn.Cross(rgn))				// crossed feature?
		return CorrectItemsEnd(it, rgn.End, ambig.TreatCase(ambig.CROSS));
	return true;
}

// Adds feature to the container.
//	return: true if Read was added successfully
bool BedF::AddPos(const Region& rgn, const char* name, float score, const char* strand)
{
	_items.push_back(Featr(rgn
#ifdef _ISCHIP
		, score));
	if(score > _maxScore)	_maxScore = score;
#else
		));
#endif
	return true;
	//if(end - start > _maxFtrLen) {	_maxFtrLen = end - start; _maxFtrStart = start; }
}

// Decreases Feature's positions without checkup indexes.
//	@cID: chromosome's ID
//	@fInd: index of Feature
//	@shift: decrease Feature's positions on this value
//	@rgEnd: region's end position to control feature location
//	@return: true if Feature is insinde this region
//bool BedF::DecreasePos(chrid cID, chrlen fInd, chrlen shift, chrlen rgEnd) {
//	//if(Item(cInd, fInd).Start < shift)
//	//	Err("feature outside region", "Region end="+IntToStr(rgEnd)).Throw();
//	if(Item(cID, fInd).End > rgEnd) 
//		return false;
//	chrlen ind = At(cID).FirstInd + fInd;
//	_items[ind].Start -= shift;
//	_items[ind].End -= shift;
//	return true;
//}

// Gets chromosome's treated length:
// a double length for numeric chromosomes or a single for named.
//	@it: chromosome's iterator
//	@multiplier: 1 for numerics, 0 for letters
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
ULONG BedF::FeaturesTreatLength(cIter it, BYTE multiplier, float fLen) const
{
	ChromItemsInd cII = it->second;
	ULONG	res = 0;
	for(chrlen i=cII.FirstInd; i<=cII.LastInd; i++)
		res += _items[i].Length() + int(2*fLen);
	return res << multiplier;
}

#ifdef _ISCHIP
// Scales defined score through all features to the part of 1.
void BedF::ScaleScores ()
{
	ItemsIter fit;
	for(cIter cit=Begin(); cit!=End(); cit++)
		for(fit=ItemsBegin(cit); fit!=ItemsEnd(cit); fit++)
			fit->Score /= _maxScore;	// if score is undef then it become 1
}
#else	// NO _ISCHIP
// Copies feature coordinates to external Regions.
void BedF::FillRegions(chrid cID, Regions& regn) const {
	const ChromItemsInd& cii = At(cID);
	regn.Reserve(cii.LastInd - cii.FirstInd + 1);
	//vector<Featr>::const_iterator itEnd = _items.end() + cii.LastInd + 1;
	//for(vector<Featr>::const_iterator it=_items.begin() + cii.FirstInd; it!=itEnd; it++)
	//	regn.AddRegion(it->Start, it->End);
	regn.Copy(_items, cii.FirstInd, cii.LastInd);
}
#endif	// _ISCHIP

// Extends or shrinks all chromosomes features positions on the fixed length in both directions.
//	@extLen: distance on which start should be decreased, end should be increased
//	or inside out if it os negative
//	@printStats: true if statistics should be printed
void BedF::Extend(int extLen, bool printStats)
{
	if( !extLen )	return;
	ULONG	rmvCnt = 0;		// counter of removed items in current chrom
	ItemsIter fit, fitLast;
	Ambig ambig(false, true, _ItemTitle);

	for(Iter cit=Begin(); cit!=End(); cit++) {
		fit		= _items.begin() + (cit->second.FirstInd -= rmvCnt);
		fitLast	= _items.begin() + (cit->second.LastInd  -= rmvCnt);
		fit->Start -= extLen;
		fit->End += extLen;
		for(rmvCnt=0, fit++; fit<=fitLast; fit++) {
			fit->Start -= extLen;
			fit->End += extLen;
			if( ( fit->Start < (fit-1)->End	&& ambig.TreatCase(ambig.CROSS) >= 0 )
			||( fit->Start == (fit-1)->End	&& ambig.TreatCase(ambig.ADJAC) >= 0 ) )
			{	// merge crossing/adjacent features
				(fit-1)->End = fit->End;
				_items.erase(fit--);
				cit->second.LastInd--;
				fitLast--;
				rmvCnt++;
			}
		}
	}
	if( printStats && ambig.AmbigCount() )
		ambig.PrintStatistics(
			ChromsCount()==1 ? CID(Begin()) : Chrom::UnID,
			"during expanding", false, _itemsCnt - rmvCnt, _itemsCnt
		);
	_itemsCnt -= rmvCnt;
}

// Checks whether all features length exceed gien length, throws exception otherwise.
//	@len: given control length
//	@lenDefinition: control length definition to print in exception message
//	@sender: exception sender to print in exception message
void BedF::CheckFeaturesLength(chrlen len, const string lenDefinition, const char* sender)
{
	ItemsIter fit;
	for(cIter cit=Begin(); cit!=End(); cit++)
		for(fit=ItemsBegin(cit); fit!=ItemsEnd(cit); fit++)
			if(fit->Length()<len)
				Err("Feature size " + NSTR(fit->Length()) +
					" is less than stated " + lenDefinition + sBLANK + NSTR(len),sender).Throw();
}

/************************ end of class BedF ************************/
#endif	// _VALIGN

/************************ class Nts ************************/
#define CNT_DEF_NT_REGIONS	10

// Copy current readed line to the nucleotides buffer.
//	@line: current readed line
//	@lineLen: current readed line length
void Nts::CopyLine(const char* line, chrlen lineLen)
{
	//ifIsBadReadPtr(line, lineLen))			cout << "BAD src PTR: lineLen = " << lineLen << endl;
	//if(IsBadReadPtr(_nts + _len, lineLen))	cout << "BAD dst PTR: _len = " << _len << endl;
	memcpy(_nts + _len, line, lineLen);
	_len += lineLen;
}

// Creates a newNts instance
//	@fName: file name
//	@minGapLen: minimal length which defines gap as a real gap
//	@fillNts: if true fill nucleotides and def regions, otherwise def regions only
//  @letN: if true then include 'N' on the beginning and on the end 
void Nts::Init(const string& fName, short minGapLen, bool fillNts, bool letN) 
{
	_nts = NULL;
	_defRgns.Reserve(CNT_DEF_NT_REGIONS);
	FaFile::Pocket pocket(_defRgns, minGapLen);
	FaFile file(fName, pocket);

	_len = pocket.ChromLength();
	if( fillNts ) {
		try { _nts = new char[_len]; }
		catch(const bad_alloc&) { Err(Err::F_MEM, fName).Throw(); }
		_len = 0;	// is accumulated while reading. Should be restore at the end
	}
	if( fillNts && !minGapLen && letN )	// fill nts without defRegions. First line is readed yet
		for(const char* line = file.Line(); line; line = file.GetLine() )
			CopyLine(line, file.LineLength());
	else if( minGapLen || !letN ) {		// fill nts and defRegions. First line is readed yet
		for(const char* line = file.Line(); line; line = file.GetLine(pocket) )
			if( fillNts )
				CopyLine(line, file.LineLength());
		pocket.CloseAddN();		// close def regions
		_cntN = pocket.CountN();
	}

	if( file.IsBad() ) {
		if(_nts)	delete [] _nts;
		_nts = NULL;
		Err(file.ErrCode(), fName).Throw();
	}
	// set _commonDefRgn
	if( !letN && _defRgns.Count() > 0 ) {
		_commonDefRgn.Start = _defRgns.FirstStart();
		_commonDefRgn.End = _defRgns.LastEnd();
	}
	else {
		_commonDefRgn.Start = 0;
		_commonDefRgn.End = _len-1;
	}
	//cout << "defRgns.Count(): " << _defRgns.Count() << EOL;
	//cout << "Start: " << _commonDefRgn.Start << "\tEnd: " << _commonDefRgn.End << EOL;
}

#if defined(_FILE_WRITE) && defined(DEBUG)
#define FA_LINE_LEN	50	// length of wrtied lines

void Nts::Write(const string & fName, const char *chrName) const
{
	FaFile file(fName, chrName);
	chrlen i, cnt = _len / FA_LINE_LEN;
	for(i=0; i<cnt; i++)
		file.AddLine(_nts + i * FA_LINE_LEN, FA_LINE_LEN);
	file.AddLine(_nts + i * FA_LINE_LEN, _len % FA_LINE_LEN);
	file.Write();
}
#endif	// DEBUG

/************************ end of class Nts ************************/

/************************  class ChromFiles ************************/
const string GenomeFileMsg(chrid cID) {
	return " genome file" + 
		( cID == Chrom::UnID ? "s" : " for given " + Chrom::TitleName(cID));
}

// Fills external vector by chrom IDs relevant to file's names found in given directory.
//	@files: empty external vector of file's names
//	@gName: name of .fa files directory or single .fa file
//	@cID: chromosomes ID that sould be treated, or Chrom::UnID if all
//	return: count of filled chrom IDs
//	Method first searches chroms among .fa files.
//	If there are not .fa files or there are not .fa file for given cID,
//	then searches among .fa.gz files
BYTE ChromFiles::GetChromIDs(vector<string>& files, const string& gName, chrid cID)
{
	if( !FS::GetFiles(files, gName, _ext) )		return 0;

	string name;
	chrid	cid,				// chrom ID relevant to current file in files
			i,					// index of current file in files
			wrongNamesCnt = 0;	// counter of additional chromosomes ("1_random" etc)
	int		prefixLen;			// length of prefix of chrom file name
	BYTE	extLen = _ext.length();
	BYTE	cnt = BYTE(files.size());

	// remove additional names and sort listFiles
	for(i=0; i<cnt; i++) {
		prefixLen = CommonPrefixLength(files[i], extLen);
		if( prefixLen >= 0 ) {				// right chrom file name 
			// set to name a chromosome name
			if( !_prefixName.size() )		// not initialized yet
				_prefixName = files[i].substr(0, prefixLen);
			name = files[i].substr(prefixLen, files[i].length() - prefixLen - extLen);
			// filter additional names
			cid = Chrom::ID(name);
			if( cid != Chrom::UnID ) {		// "pure" chrom's name
				if( cID == Chrom::UnID ) {
					// add '0' to a single numeric for correct sorting
					if( isdigit(name[0]) && (name.length() == 1 || !isdigit(name[1])) )
						name.insert(0, 1, '0');
				}
				else 
					if( cID != cid ) { name = StrEmpty; wrongNamesCnt++; }
			}
			else { name = StrEmpty; wrongNamesCnt++; }	// additional chrom's name
		}
		else { name = StrEmpty; wrongNamesCnt++; }		// some other .fa[.gz] file, not chrom
		files[i] = name;
	}
	sort(files.begin(), files.end());
	if( wrongNamesCnt ) {	// are there any additional chromosome's names?
		// remove empty names - they are at the beginning due to sorting
		files.erase(files.begin(), files.begin() + wrongNamesCnt);
		cnt = files.size();
	}
	// remove added '0' from a single numeric
	for(i=0; i<cnt; i++)
		if( files[i][0] == '0' )
			files[i].erase(0,1);
	return cnt;
}

// Creates and initializes an instance.
//	@gName: name of .fa files directory or single .fa file. If single file, then set Chrom::StatedID()
//	@getAll: true if all chromosomes should be extracted
ChromFiles::ChromFiles(const string& gName, bool extractAll)
	: _ext(FaFile::Ext), _extractAll(extractAll)
{
	vector<string> listFiles;
	BYTE cnt = 1;	// number of chroms readed from pointed location 

	if( FS::IsDirExist(gName.c_str()) ) {	// gName is a directory
		if( !GetChromIDs(listFiles, gName, Chrom::StatedID()) ) {	// fill IDs from .fa files
			_ext += ZipFileExt;
			if( !GetChromIDs(listFiles, gName, Chrom::StatedID()) // fill IDs from .fa.gz files
			&& Chrom::StatedAll() )					
				Err( Err::MsgNoFiles("*", FaFile::Ext), gName ).Throw();
		}
		_path = FS::MakePath(gName);
		cnt = listFiles.size();
		if( !cnt )	Err("no" + GenomeFileMsg(Chrom::StatedID()), gName).Throw();
	}
	else {		// gName is a file. It exists since it had been checked in main()
		BYTE extLen = _ext.length();
		if( FS::HasGzipExt(gName) ) {
			_ext += ZipFileExt;
			extLen += ZipFileExt.length();
		}
		string fName = FS::ShortFileName(gName);
		int prefixLen = CommonPrefixLength(fName, extLen);
		chrid  cid = Chrom::ID(fName.c_str(), prefixLen);
		Chrom::SetStatedID(cid);	// do not check Chrom::StatedID for UnID
									// since it had been checked in main()
		if( Chrom::StatedID() != Chrom::UnID && Chrom::StatedID() != cid )
			Err("wrong" + GenomeFileMsg(Chrom::StatedID()), gName).Throw();
		listFiles.push_back( Chrom::Name(cid) );
		_prefixName = fName.substr(0, prefixLen);
		_path = FS::DirName(gName, true);
	}
	// fill attributes
	Reserve(cnt);
	for(chrid i=0; i<cnt; i++)
		AddVal(Chrom::ID(listFiles[i]), listFiles[i]);
}

// Returns length of common prefix before abbr chrom name of all file names
//	@fName: full file name
//	@extLen: length of file name's extention
//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
inline int	ChromFiles::CommonPrefixLength(const string & fName, BYTE extLen)
{
	// a short file name without extention
	return Chrom::PrefixLength(	fName.substr(0, fName.length() - extLen).c_str());
}

// Returns full file name or first full file name by default
//	@cID: chromosome's ID
const string ChromFiles::FileName(chrid cID) const {
	if( !cID )	cID = FirstChromID();
	return FullCommonName() + Chrom::Name(cID) + _ext;
}

#ifdef _ISCHIP

// Sets actually treated chromosomes indexes and sizes according bed.
//	@bed: template bed. If NULL, set all chromosomes
//	return: count of treated chromosomes
chrid ChromFiles::SetTreated(const Bed* const bed)
{
	chrid	cnt = 0;
	LLONG	sz;
	string	fname;
	Iter it = Begin();
	bool isZipped = FS::HasGzipExt(FileName(CID(it)));

	for(; it!=End(); it++)
		if( _extractAll || !bed || bed->FindChrom(CID(it)) ) {
			cnt++;
			fname = FileName(CID(it));
			sz = isZipped ?
				FS::UncomressSize(fname.c_str()) :
				FS::Size(fname.c_str());
			if( sz < 0 )	Err(Err::F_OPEN, fname).Throw();
			it->second.FileLen = chrlen(sz); 
		}
	return cnt;
}

// Gets count of treated chromosomes.
chrid ChromFiles::TreatedCount() const
{
	if( _extractAll )	return ChromsCount();
	chrid res = 0;
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if( TREATED(it)() )
			res++;
	return res;
}

#endif	// _ISCHIP

// Gets the total length of files (dupl==false) or nucleotides (dupl==true) with EOLs
//  dupl: true for duplicated numeric's chromosomes
//ULLONG ChromFiles::TotalLength(const bool dupl) const
//{
//	ULLONG res = 0;
//	for(BYTE i=0; i<_cnt; i++)
//		res += (dupl ? _attrs[i].Multiplier : 1) * _attrs[i].FileSize;
//	return res;
//}

#ifdef DEBUG
void ChromFiles::Print() const
{
	cout << "ChromFiles: count of chroms: " << int(ChromsCount()) << endl;
	cout << "chrom\tNumeric\tFileLen\n";
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::AbbrName(CID(it)) << TAB
			 << int(it->second.Numeric) << TAB
			 << it->second.FileLen << EOL;
}
#endif	// DEBUG
/************************  end of class ChromFiles ************************/

/************************ class ChromSizes ************************/
const string ChromSizes::Ext = ".sizes";

// Initializes instance from file.
//	@fName: name of file.sizes
void ChromSizes::Init (const string& fName)
{
	TabFile file(fName, 2, true, TxtFile::READ, Chrom::Abbr);
	chrid cID;
	ULONG cntLines;
	const char* currLine = file.GetFirstLine(&cntLines);
	Reserve(chrid(cntLines));
	for(; currLine; currLine = file.GetLine())
		// skip random chromosomes
		if( (cID=Chrom::IDbyAbbrName(file.StrField(0))) != Chrom::UnID )
			AddVal(cID, file.LongField(1));
}

// Saves instance to file
//	@fName: full file name
void ChromSizes::Write(const string fName)
{
	LineFile file(fName, TAB);

	file.BeginWrite(Chrom::MaxNamedPosLength+1);

	// needed to sort because of unoredictable order of chr size adding
#ifdef _NO_HASHTABLE
	Sort();
	for(cIter it=cBegin(); it!=cEnd(); it++)
#else
	vector<cSize> cSizes(cBegin(), cEnd());
	sort(cSizes.begin(), cSizes.end(), SizeCompare);
	for(vector<cSize>::iterator it=cSizes.begin(); it!=cSizes.end(); it++)
#endif
		file.WriteLine(Chrom::AbbrName(CID(it)), it->second);

	file.Write();
}

// Creates a new instance by chrom files.
//	@cFiles: chrom files
//	@printReport: if true then print report about generation/addition size file to dout
// Reads an existing chrom sizes file if it exists, otherwise creates new instance.
// Cheks and adds chrom if it is absent.
// Saves instance to file if it is changed.
ChromSizes::ChromSizes (const ChromFiles& cFiles, bool printReport)
{
	const string fName = cFiles.Path() + FS::LastSubDirName(cFiles.FileName()) + ".chrom" + Ext;
	bool updated = !FS::IsFileExist(fName.c_str());	// false, if file exists
	bool dontCheck = updated;
	const char* report = updated ? "Generate " : "Redefine ";
	Timer tm;

	if(updated)	Reserve(cFiles.ChromsCount());	// file doesn't exist, create it
	else		Init(fName);

	for(ChromFiles::cIter it=cFiles.cBegin(); it!=cFiles.cEnd(); it++)
		if(dontCheck || !this->FindChrom(CID(it))) {
			if(printReport) {
				dout << report << Chrom::Title << " sizes file...";
				tm.Start();
				printReport = false;
				//fflush(dout);
			}
			AddValFromFile(CID(it), cFiles);
			updated = true;
		}
	if(updated) {
		Write(fName);
		dout << MsgDone;
		tm.Stop(false);
	}
}

// Gets total size of genome
//genlen ChromSizes::GenSize() const
//{
//	if( !_gsize )
//		for(cIter it=cBegin(); it!=cEnd(); it++)
//			_gsize += it->second;
//	return _gsize;
//}

//#ifdef _BIOCC
// Gets miminal size of chromosome
//chrlen ChromSizes::MinSize() const
//{
//	if( !_minsize ) {
//		cIter it=cBegin();
//		_minsize = it->second;
//		for(it++; it!=cEnd(); it++)
//			if( _minsize > it->second )
//				_minsize = it->second;
//	}
//	return _minsize;
//}
//#endif	// _BIOCC

#ifdef DEBUG
void ChromSizes::Print() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::TitleName(CID(it)) << TAB << it->second << EOL;
}
#endif	// DEBUG

/************************ end of class ChromSizes ************************/

#if defined _DENPRO || defined _BIOCC

/************************ class FileList ************************/
#ifdef OS_Windows
// Returns true if 'name' is file's pattern
bool	IsFilePattern	(const char* name)
{
	return strchr(name, '*') != NULL || strchr(name, '?') != NULL;
}

string GetPath	(const LPCTSTR name)
{
	const char* pch = strrchr(name, '/');
	return pch ? string(name, pch-name+1) : StrEmpty;
}

// Fills vector by filenames according name template.
// Works only if _UNICODE is not defined.
//	@files: vector of file names
//	@templ: name template, which can include '*' and '?' marks
void FillFilesByTemplate(vector<string>& files, const LPCTSTR templ)
{
	string path = GetPath(templ);
	WIN32_FIND_DATA ffd;
 
	// check directory and estimate listFileNames capacity
	HANDLE hFind = FindFirstFile( templ, &ffd );
	if( hFind == INVALID_HANDLE_VALUE )		
		Err("bad file or content", templ).Throw();
	if( files.capacity() == 0 ) {		// count files to reserve capacity
		short count = 1;
		for(; FindNextFile(hFind, &ffd); count++);
		files.reserve(count);
		hFind = FindFirstFile( templ, &ffd );
	}
	// fill the list
	do	files.push_back(path + string(ffd.cFileName));	//  works only if _UNICODE isn't defined
	while (FindNextFile(hFind, &ffd));
	FindClose(hFind);
}
#endif	// OS_Windows

FileList::FileList(char* files[], short cntFiles) : _files(NULL), _memRelease(true)
{
#ifdef OS_Windows
	short i;
	bool hasTemplate;
	for(i=0; i<cntFiles; i++)
		if( hasTemplate = IsFilePattern(files[i]) )
			break;
	if( hasTemplate ) {
		// First we fill vector of file names, second initialize _files by this vector.
		// It needs because files[] may contain file names and file template (pattern) as well,
		// so in case of Windows we don't know the total amount of files beforehand.
		vector<string> tmpFiles;
		if( cntFiles > 1 )
			tmpFiles.reserve(cntFiles);
		// else if it's a single template name, capacity will be reserved in FillFilesByTemplate(),
		// or if it's a single common name, capacity will not be reserved at all

		for(i=0; i<cntFiles; i++)
			if( IsFilePattern(files[i]) )
				FillFilesByTemplate(tmpFiles, files[i]);
			else		// real list of files
				tmpFiles.push_back(files[i]);

		_files = new char*[_count=tmpFiles.size()];
		for(i=0; i<_count; i++) {
			_files[i] = (char*)malloc(tmpFiles[i].length()+1);
			strcpy(_files[i], tmpFiles[i].c_str());
		}
	}
	else 
#endif	// OS_Windows
	{
		_files = files;
		_count = cntFiles;
		_memRelease = false;
	}
}

FileList::FileList(const char* fName) : _files(NULL), _memRelease(true)
{
	TabFile file(fName, 1, true, TxtFile::READ, NULL, TabFile::Comment);
	ULONG cntLines;
	char *dstStr;
	const char *srcStr;
	const char *currLine = file.GetFirstLine(&cntLines);
	vector<char*> tmpFiles;	// temporary vector because cntLines is not proof, but estimated capacity
	
	_count = 0;
	if(currLine == NULL)	return;		// no real lines
	tmpFiles.reserve(cntLines);
	while(currLine!=NULL) {
		srcStr = file.StrField(0);
		dstStr = (char*)malloc(strlen(srcStr)+1);
		strcpy(dstStr, srcStr);
		tmpFiles.push_back(dstStr);
		_count++;
		currLine=file.GetLine();
	}
	_files = new char*[_count];
	for(short i=0; i<_count; i++)
		_files[i] = tmpFiles[i];
}

FileList::~FileList()
{
	if( _files && _memRelease ) {
		for(short i=0; i<_count; i++)
			free(_files[i]);
		delete [] _files;
	}
}

#ifdef DEBUG
void FileList::Print() const
{
	if( _files )
		for(short i=0; i<_count; i++)
			cout << _files[i] << EOL;
	else
		cout << "Empty\n";
}
#endif	// DEBUG
/************************ end of class FileList ************************/

/************************ class ChromRegions ************************/
const string ChromRegions::_FileExt = ".region";

// Creates an instance from file 'chrN.regions', if it exists.
// Otherwise from .fa file then writes it to file 'chrN.regions'
// File 'chrN.regions' is placed in @gName directory
//	@commName: full common name of .fa files
//	@cID: chromosome's ID
//	@minGapLen: minimal length which defines gap as a real gap
ChromRegions::ChromRegions(const string& commName, chrlen cID, short minGapLen)
{
	// get the name of regions file. cID is checked already in GenomeRegions()
	string fName = commName + Chrom::Name(cID);
	string regionFileName = fName + DOT + NSTR(minGapLen) + _FileExt;

	// get data
	if( FS::IsFileExist(regionFileName.c_str())		// chrN.region already exist?
	&& Read(regionFileName) == minGapLen			// read it: has the same minGapLen?
	&&	Count() )									// are there regions in chrN.region?
		return;
	// read chrN.fa to define regions
	string faFileName = fName + FaFile::Ext;
	if( !FS::IsFileExist(faFileName.c_str()) ) {
		faFileName += ZipFileExt;
		if( !FS::IsFileExist(faFileName.c_str()) )
			Err(Err::MsgNoFiles(FS::ShortFileName(fName), FaFile::Ext),
				FS::DirName(fName, false)).Throw();
	}
	Copy(Nts(faFileName, minGapLen, true).DefRegions());
	Write(regionFileName, minGapLen);
}

/************************ end of class ChromRegions ************************/

/************************ class GenomeRegions ************************/

//void GenomeRegions::Init(const ChromSizes* cSizes)
//{
//	if( Chrom::StatedAll() ) {
//		Reserve(cSizes->ChromsCount());
//		for(ChromSizes::cIter it=cSizes->cBegin(); it != cSizes->cEnd(); it++)
//			AddClass(CID(it), Regions(0, it->second));
//	}
//	else
//		AddClass(Chrom::StatedID(), Regions(0, cSizes->Size(Chrom::StatedID())));
//
//}

GenomeRegions::GenomeRegions(const char* gName, const ChromSizes* cSizes, short minGapLen)
	: _minGapLen(minGapLen), _singleRgn(FS::HasExt(gName, ChromSizes::Ext))
{
	if(_singleRgn) {
		cSizes = new ChromSizes(gName);
		// initialize instance from chrom sizes
		if( Chrom::StatedAll() ) {
			Reserve(cSizes->ChromsCount());
			for(ChromSizes::cIter it=cSizes->cBegin(); it != cSizes->cEnd(); it++)
				AddClass(CID(it), Regions(0, it->second));
		}
		else
			AddClass(Chrom::StatedID(), Regions(0, cSizes->Size(Chrom::StatedID())));
	}
	else {
		const ChromFiles cFiles(gName);
		cSizes = new ChromSizes(cFiles);
		// initialize instance from chrom files
		_commonName = cFiles.FullCommonName();
	}
}

#ifdef _BIOCC
// Gets total genome's size: for represented chromosomes only
genlen GenomeRegions::GenSize() const
{
	genlen gsize = 0; 
	for(cIter it=cBegin(); it!=cEnd(); it++)
		gsize += Size(it);
	return gsize;
}

// Gets miminal size of chromosome: for represented chromosomes only
chrlen GenomeRegions::MinSize() const
{
	cIter it=cBegin();
	chrlen	minsize = Size(it);
	for(it++; it!=cEnd(); it++)
		if( minsize > Size(it) )
			minsize = Size(it);
	return	minsize;
}
#endif	// _BIOCC

#ifdef DEBUG
void GenomeRegions::Print() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout<< Chrom::TitleName(CID(it))
			<< TAB << it->second.FirstStart() 
			<< TAB << Size(it) << EOL;
}
#endif	// DEBUG
/************************ end of class GenomeRegions ************************/
#endif	// _DENPRO || _BIOCC
