/**********************************************************
Data.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 10.04.2019
-------------------------
Provides common data functionality
***********************************************************/

#include "Data.h"
#ifdef _ISCHIP
	#include "isChIP.h"
#elif defined _BIOCC
	#include "Calc.h"
#endif	//_BIOCC
#include <fstream>	// to write simple files without _FILE_WRITE


static const char* Per = " per ";

/************************ class Obj ************************/

Obj::Spotter::Msg Obj::Spotter::_Msgs [] = {
	{ NULL, "duplicated",	"duplicated", },
	{ NULL, "crossed",		"is intersected with previous" },
	{ NULL, "adjacent",		"is adjacent with previous" },
	{ NULL, "covered",		"is fully covered by previous" },
	{ NULL, "too short",	"too short" },
	{ NULL, "different size of", "different size of read" },
	{ NULL, "filtered by low score", "filtered by score" },
	{ NULL, "chrom exceeding",	"position exceeds chromosome length" },
	{ NULL, "negligible",	"negligible chromosome" },
};
//const BYTE Obj::Spotter::_CasesCnt = sizeof(Obj::Spotter::_Msgs)/sizeof(Obj::Spotter::Msg);

const char*	Obj::Spotter::_ActionMsgs[] = {
	"accepted",
	"joined",
	"omitted",
	"omitted",
	"execution aborted"
};

const Obj::Spotter::ReportCase Obj::Spotter::_Actions[] = {
	&Obj::Spotter::Accept,
	&Obj::Spotter::Handle,
	&Obj::Spotter::Omit,
	&Obj::Spotter::OmitQuiet,
	&Obj::Spotter::Abort
};

// Gets count of ambiguities
chrlen Obj::Spotter::Count() const
{
	if( _count == CHRLEN_UNDEF ) {
		_count = 0;
		for(BYTE i=0; i<_CasesCnt; i++)	_count += _cases[i].Count;
	}
	return _count;
}

// Print given spotter as line alarm
//	@spotter: given spotter
void Obj::Spotter::PrintLineAlarm(eCase spotter) const
{
	if( _alarm ) {
		if( !_alarmPrinted )	{ dout << EOL;	_alarmPrinted = true; }
		_file->ThrowLineWarning(
			_Msgs[spotter].LineAlarm + BLANK + ItemTitle() + SepCl,
			Message(spotter));
	}
}

// Outputs case statistics
//	@spotter: spotter's case
//	@allCnt: total count of ambiguities
//	@total: if true then prints total warning case
void Obj::Spotter::PrintCaseStat(eCase spotter, chrlen allCnt, bool total) const
{
	chrlen cnt = _cases[spotter].Count;	// count of case's ambiguities
	if( !cnt )	return;
	const char* totalAlarm = _Msgs[spotter].TotalAlarm;
	if(total)	dout << Notice;
	else		dout << TAB;
	dout<< cnt
		<< sPercent(ULLONG(cnt), ULLONG(allCnt), 4) << BLANK
		<< _Msgs[spotter].StatInfo << BLANK
		<< ItemTitle(cnt);
	//if(unsortedItems)		dout << " arisen after sorting";
	if(totalAlarm)	dout << BLANK << totalAlarm;
	dout << SepSCl << Message(spotter);
	if(totalAlarm)	dout << '!';
	dout << EOL;
}

const char* ACCEPTED = " accepted";

// Prints accepted items with specifying chrom
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@prAcceptItems: if true then prints number of accepted items
//	@itemCnt: count of accepted items after treatment
void Obj::Spotter::PrintItems(chrid cID, bool prAcceptItems, long itemCnt) const
{
	if(prAcceptItems) {
		dout << itemCnt;
		dout << ACCEPTED;
	}
	if(_info > Obj::iNM) {
		dout << BLANK << ItemTitle(itemCnt);
		if( cID != Chrom::UnID )	dout << Per << Chrom::ShortName(cID);
	}
}

// Creates an instance with omitted COVER, SHORT, SCORE and NEGL cases,
// and BedF cases by default:
// omitted DUPL cases and handled CROSS and ADJAC 
Obj::Spotter::Spotter (eInfo info, bool alarm, FT::eTypes format,
	eAction dupl, eAction crossANDadjac, eAction diffsz) :
	_info(info),
	_alarm(alarm),
	_fType(format),
	_alarmPrinted(info <= Obj::iLAC),	// if LAC, do not print EOL at the first time
	_count(CHRLEN_UNDEF),
	_file(NULL),
#ifndef _ISCHIP
	_treatcID(vUNDEF),
#endif
	unsortedItems(false),
	wasPrinted(false),
	chrLen(CHRLEN_UNDEF)
{
	memset(_cases, 0, _CasesCnt*sizeof(Case));	// initialize by 0
	_cases[DUPL].Action = dupl;
	_cases[CROSS].Action = _cases[ADJAC].Action = crossANDadjac;
	_cases[COVER].Action = _cases[SHORT].Action = _cases[SCORE].Action = OMIT;
	_cases[DIFFSZ].Action = diffsz;
	_cases[EXCEED].Action = OMIT;
}

// Initializes given Region by second and third current reading line positions, with validating
//	@rgn: Region that should be initialized
//	@prevStart: previous start position
//	return: true if Region was initialized successfully
bool Obj::Spotter::InitRegn(Region& rgn, chrlen prevStart)
{
	// dismiss check for negative
	//long start = _file->LongField(1);
	//if(start < 0)		ThrowExcept(Err::B_NEGPOS);
	//rgn.Set(start, _file->LongField(2));

	rgn.Set( _file->LongField(1), _file->LongField(2));
	if(rgn.Invalid())	ThrowExcept(Err::B_BADEND);
	if(rgn.End>chrLen && TreatCase(EXCEED)<0)	return false;	// exceeding chrom
	if(rgn.End < prevStart)		unsortedItems = true;
	return true;
}

// Prints statistics.
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@title: string at the beginning; if NULL then this instance is used while initialization and don't feeds line
//	@totalItemCnt: count of all items
//	@acceptItemCnt: count of accepted items after treatment
//	The last line never ends with EOL 
void Obj::Spotter::Print(chrid cID, const char* title, ULONG totalItemCnt, ULONG acceptItemCnt)
{
	if(_info <= Obj::iLAC || !totalItemCnt)		return;
	bool noSpotters = totalItemCnt == acceptItemCnt;

	if(wasPrinted)	dout << EOL;	// "sorting" was printed
	if( title )	{		// additional mode: after extension
		if(_info < Obj::iEXT || noSpotters)		return;	// no ambigs
		dout << "    " << title << SepCl;
		if(_info==Obj::iEXT)
			PrintItems(cID, true, acceptItemCnt);
		dout << EOL;
	}
	else {				// main mode: addition to file name
		bool printAccept = _info==Obj::iEXT && !noSpotters;		// print accepted items

		if(_info > Obj::iNM) {
			dout << SepCl << totalItemCnt;
			if(!noSpotters)	dout << BLANK << Total;
		}
		if(printAccept)		dout << SepCm;
		PrintItems(cID, printAccept, acceptItemCnt);
		wasPrinted = true;
	}
	if(Count() && _info == Obj::iSTAT) {
		if(!title)		dout << SepCm << "from which" << EOL;
		for(BYTE i=0; i<_CasesCnt; i++)
			PrintCaseStat(static_cast<eCase>(i), totalItemCnt);

		// print ambiguities of negligible chroms if all chroms are readed
		if( cID == Chrom::UnID ) {
			// calculate ambigs for negligible chroms:
			// rest of difference between in & out features minus count of ambigs
			_cases[NEGL].Count = totalItemCnt - acceptItemCnt - Count();
			// correct (add) accepted features
			for(BYTE i=0; i<_CasesCnt-1; i++)	// loop excepting NEGL chroms case
				if( _cases[i].Action == Spotter::ACCEPT )	_cases[NEGL].Count += _cases[i].Count;
			
			PrintCaseStat(NEGL, totalItemCnt);
		}
		// print total remained entities
		dout<< TAB << Total << ACCEPTED << SepCl << acceptItemCnt
			<< sPercent(ULLONG(acceptItemCnt), ULLONG(totalItemCnt), 4) << BLANK
			<< ItemTitle(acceptItemCnt);
		wasPrinted = true;
	}
	fflush(stdout);
}

#ifndef _ISCHIP
// Remember treated chrom.
void Obj::Spotter::SetTreatedChrom(chrid cid)
{
	if(_treatcID == -1)					_treatcID = cid;
	else if(_treatcID != Chrom::UnID)	_treatcID = Chrom::UnID;
}
#endif

// EOL rules.
// appearance of info on the screen while Init() and others | ended by EOL
// -----------------------------------------------------------------------
//	exception throwing from file.LongField() etc.	-
//	exception throwing from ChildInit()				-
//	line alarm warning								+
//	file title (silent mode)						-
//	file title (nonesilent mode)					-
//	file title + item info (nonesilent mode)		-
//	statistics										-
//	statistics after extention						-
//	


// Throws exception or warning message
//	@err: input error
//	@abortInvalid: if true, throws exception, otherwise throws warning
void Obj::ThrowError(Err &err, bool abortInvalid)
{
	_isBad = true;
	_EOLneeded = true;
	err.Throw(abortInvalid, false);
}

// Initializes new instance by by tab file name.
//	@title: title printed before file name
//	@fName: name of file
//	@spotter: ambiguities
//	@addObj: auxiliary object using while initializing
//	@isInfo: true if file info should be printed
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
void Obj::Init	(const char* title, const string& fName, Spotter& spotter, void* addObj,
	bool isInfo, bool abortInvalid)
{
	// in case of bioCC the list is checked already,
	// so better is to add a parameter to Init(),
	// but it is not worth it
	FT::CheckType(fName.c_str(), spotter.FileType(), true, abortInvalid);
	dchrlen items;
	Timer	timer(isInfo);

	if( isInfo ) {
		if(title)	dout << title << BLANK;
		dout << fName;
		fflush(stdout); 
		_EOLneeded = true;
	}
	try {
		TabFile file(fName, FT::FileParams(spotter.FileType()), abortInvalid, !isInfo);
		spotter.SetFile(file);
		items = InitChild(spotter, addObj);
	}
	catch(Err& err) {	// intercept an exception to manage _isBad and aborting if invalid
		ThrowError(err, abortInvalid);
	}
#ifndef _ISCHIP
	if(abortInvalid)
		spotter.KeepTreatedChrom();	// save treated chrom for primary
#endif
	if( !_isBad ) {
		if( !items.second ) {		// no items for given chrom
			string sender = isInfo ? strEmpty : fName;
			string specify = ItemTitle(true);
			if(!Chrom::NoCustom())	specify += Per + Chrom::ShortName(Chrom::CustomID());
			Err err(Err::TF_EMPTY, sender.c_str(), specify);
			ThrowError(err, abortInvalid);
		}
		spotter.Print(Chrom::CustomID(), NULL, items.first, items.second);
	}
	if(timer.IsEnabled())	dout << BLANK;
	timer.Stop(true, false);
	PrintEOL(spotter.wasPrinted);	// || spotter.IsAlarmPrinted());
}

// Prints EOL if needs.
//	@printEOL: true if EOL should be printed
void Obj::PrintEOL(bool printEOL)
{
	if(printEOL || _EOLneeded)	dout << EOL;
	_EOLneeded = false;
}

/******************** end of class Obj *********************/

/************************ class Bed ************************/

// Checks if chromosome is uniq and adds it to the container
//	@cID: chroms id
//	@firstInd: first item index
//	@lastInd: last item index
//	@fname: file name for exception message
void Bed::AddChrom(chrid cID, chrlen firstInd, chrlen lastInd, const string& fname)
{
	if(FindChrom(cID))
		Err(ItemTitle(true) + " are not consolidated on chromosomes", fname).Throw(true);
	AddVal(cID, ChromItemsInd(firstInd, lastInd));
}

// Initializes instance from tab file
//	@spotter: ambiguities
//	@pcSizes: chrom sizes
//	return: numbers of all and initialied items for given chrom
dchrlen Bed::InitChild	(Spotter& spotter, void* pcSizes)
{
	ULONG	initSize;		// initial size of _items
	TabFile& file = spotter.File();
	ChromSizes* cSizes = (ChromSizes*)pcSizes;

	if( !file.GetFirstLine(&initSize) )	return make_pair(0, 0);
		
	chrlen 	firstInd = 0,	// first index in feature's container for current chromosome
			cntLines = 0,	// count of lines beginning with 'chr'
			currInd	= 0,	// current index in Feature's/Read's container.
							// Needed to avoid excess virtual method given current container size
			prevStart = 0;	// start previous feature positions
	Region	rgn;			// current feature positions
	chrid	cCurrID,		// current chromosome's ID
			cNextID;		// next chromosome's ID
	bool	unsortChroms = false;
	bool	skipChrom = false;
	const bool	cAll = Chrom::NoCustom();		// true if all chroms are permitted
	char	cMark[Chrom::MaxMarkLength + 1];	// chrom mark (only) buffer

	Reserve(cAll ? Chrom::Count : 1);
	ReserveItemContainer(initSize);
	cMark[0] = 0;				// empty initial chrom name
	cCurrID = file.ChromID();	// first chrom ID

	do {
		if( strcmp(cMark, file.ChromMark()) ) {	// next chromosome?
			// chrom's name may be long such as 'chrY_random'
			cNextID = file.ChromID();
			if( cNextID == Chrom::UnID ) {		// negligible next chromosome?
				cCurrID = cNextID;
				skipChrom = true;
 				continue;
			}
			// now we can save this chrom name
			strcpy(cMark, file.ChromMark());
			skipChrom = false;
			if(cAll) {			// are all chroms specified?
				if(currInd != firstInd)		// have been features for this chrom saved?
					AddChrom(cCurrID, firstInd, currInd, file.CondFileName());
				if(cNextID < cCurrID)		// unsorted chrom?
					unsortChroms = true;
			}
			else {		// single chromosome is specified
				// features in a single defined chrom have already been saved;
				// the chrom proper will be saved after loop
				if(rgn.End)			break;
				if(unsortChroms)
					Err("is unsorted. Single " + Chrom::ShortName(Chrom::CustomID()) + 
						" extraction is allowed only for sorted file",
						file.CondFileName()).Throw(true);
				if(cNextID != Chrom::CustomID()) { skipChrom = true; continue; }
			}
			cCurrID = cNextID;
#ifndef _ISCHIP
			spotter.SetTreatedChrom(cCurrID);
#endif
			firstInd = currInd;
			if(cSizes)	spotter.chrLen = cSizes->At(cCurrID).Size();
			cntLines++;
			if(!spotter.InitRegn(rgn, 0))			continue;
		}
		else {		// the same chrom
			if(skipChrom)							continue;
			cntLines++;
			if(!spotter.InitRegn(rgn, prevStart))	continue;
			if(!CheckLastPos(rgn, spotter))			continue;	// check positions for the same chrom only
		}
		if(AddPos(rgn, file))	currInd++;
		else					spotter.TreatCase(spotter.SCORE);
		prevStart = rgn.Start;
	}
	while( file.GetLine() );

	if( rgn.End && currInd ) {			// some features for at least one valid chrom were saved
		if( cCurrID != Chrom::UnID )	// is last chrom valid?
			AddChrom(cCurrID, firstInd, currInd, file.CondFileName());
		if( initSize/currInd > 2 )	ShrinkItemContainer();
		if( unsortChroms )		Sort();
		if( spotter.unsortedItems ) {
			bool prInfo = spotter.Info() > Obj::iNONE;
			if(prInfo) {
				Err(ItemTitle(false) + " sorting...", file.CondFileName()).Throw(false, false);
				spotter.wasPrinted = true;
			}
			SortItems(spotter);
			if(prInfo)	dout << Done;// << EOL;
		}
	}
	SetAllItemsCount();
	
	return make_pair(cntLines, AllItemsCount());
}

// Prints items name and count, adding chrom name if the instance holds only one chrom
void Bed::PrintItemCount() const
{
	size_t iCnt = AllItemsCount();
	dout << iCnt << BLANK << ItemTitle(iCnt>1);
	if(ChromCount()==1)		dout << Per << Chrom::TitleName(CID(cBegin()));
	dout << EOL;
}

//#ifdef _DENPRO
//// Shifts item's positions to collaps the 'holes' between regions.
////	@cID: chromosome's ID
////	@regns: valid (defined) regions
//void Bed::ShrinkByID(chrid cID, const DefRegions &regns)
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
		<< it->second.FirstInd << TAB << it->second.LastInd << SepClTab
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

// Checks the element for the new potential start/end positions for all possible ambiguous.
//	@it: iterator reffering to the compared element
//	@rgn: checked start/stop positions
//	@spotter: possible ambiguities: different Read length, duplicated Reads
//  return: true if item should be accepted; otherwise false
bool BedR::CheckItemPos(ItemsIter it, const Region& rgn, Spotter& spotter)
{
	if(_readLen != rgn.Length())					// different Read length?
		if(!_readLen)	_readLen = rgn.Length();	// initialize Read length once
		else if(spotter.TreatCase(spotter.DIFFSZ) < 0)	return false;

	if( rgn.Start == it->Pos()							// duplicated Read?
	&& spotter.TreatCase(spotter.DUPL) < 0 )	return false;
	// adjacent & crossed Reads are not checked: there are common cases
	// covered Reads are not checked: it's impossible case
	return true;
}

const string NotStated = " is not stated";

#ifdef _VALIGN

const string isChIP = "isChIP";

// Sets chrom ID and returns initial position from Read name
//	@rName: Read's name
//	@cID: pointer to the chrom ID to set value
//	@strand: Read's strand
//	return: Read's initial position
chrlen	BedR::GetInitPos(const char* rName, chrid* cID, bool strand)
{
	const char* name = Chrom::FindMark(rName);
	if(name) {
		*cID = Chrom::ID(name);
		if(_paired && !strand)
			return chrlen(atol(strchr(name + 1, Read::NmPos2Delimiter) + 1) - _readLen);
		return chrlen(atol(strchr(name + 1, Read::NmPos1Delimiter) + 1));
	}
	return SetAlien(cID);
}

// Sets chrom ID and returns initial number from Read name
//	@rName: Read's name
//	@cID: pointer to the chrom ID to set value
//	return: Read's initial number
chrlen	BedR::GetInitNumb(const char* rName, chrid* cID)
{
	const char* name = Chrom::FindMark(rName);
	if(name) {
		*cID = Chrom::ID(name);
		return chrlen(atol(strchr(name + 1, Read::NmPos1Delimiter) + 1));
	}
	return SetAlien(cID);
}

// Sets undefined chrom ID and returns undefined position/number; for bed-file isn't generated by isChIP
//	@cID: pointer to the chrom ID to set undefined value
//	return: undefined value
chrlen	BedR::SetAlien	(chrid* cID)
{
	*cID = Chrom::UnID;
	return CHRLEN_UNDEF;
}

#endif

// Adds Read to the container.
//	@rgn: Region with mandatory fields
//	@file: file to access to additionally fields
//	return: true if Read was added successfully
bool BedR::AddPos(const Region& rgn, TabFile& file)
{
#ifdef _BEDR_EXT
	readscr score = file.FloatFieldValid(4);
	if(score <= _minScore)	return false;				// pass Read with under-threshhold score
	if(score > _maxScore)	_maxScore = score;
#endif
#ifdef _VALIGN
	chrid cID;
	chrlen numb;
	const char* rName = file.StrFieldValid(3);
	bool strand = *file.StrFieldValid(5) == Read::Strands[0];

	if(_rNameType == Read::nmPos)			numb = GetInitPos(rName, &cID, strand);
	//else if(_rNameType == Read::nmNumb)		numb = GetInitNumb(rName, &cID);
	else if(_rNameType == Read::nmAlien)	numb = SetAlien(&cID);
	else	// Read::nmUndef
		if( strncmp(rName, isChIP.c_str(), isChIP.length()) ) {	// bed-file isn't generated by isChIP
			_rNameType = Read::nmAlien;
			numb = SetAlien(&cID);
		}
		else {
			// initialize Read length once; 
			// it's doing in BedR::CheckItemPos(), but after the first call of BedR::AddPos()
			if(!_readLen)	_readLen = rgn.Length();
			rName += isChIP.length() + 1;
			if( Chrom::FindMark(rName) ) {	// valid Read name
				if(strchr(rName, Read::NmPos2Delimiter))	_paired = true;
				_rNameType = Read::nmPos;
				numb = GetInitPos(rName, &cID, strand);
			}
			else {							// wrong Read name: no 'chr' substring
				_rNameType = Read::nmAlien;
				numb = SetAlien(&cID);
			}
		}
	_items.push_back(Read(rgn.Start, cID, numb, strand, score));
#elif defined _FRAGPRO
	const char* rName = strchr(file.StrFieldValid(3), DOT);
	if(!rName || !strchr(rName, '/'))
		file.ThrowLineExcept("Inappropriate read's name: should be '<name>.<uniq_number>/<index>'");
	_items.push_back(Read(rgn.Start, atoi(++rName), *file.StrFieldValid(5) == PLUS));
#else
	_items.push_back(Read(rgn.Start));
#endif
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

#ifdef _BEDF
/************************ class BedF ************************/

// Sets new end position on the feature if necessary.
//	@it: iterator reffering to the feature which end position may be corrected
//	@end: potential new end position
//	@treatCaseRes: result of treatment this spotter
//	return: true if spotter is permitted (feature is valid)
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
//	@spotter: possible ambiguities
//  return: true if item should be accepted; otherwise false
bool BedF::CheckItemPos(ItemsIter it, const Region& rgn, Spotter& spotter)
{
#ifdef _BIOCC
	if( _unifLen )	// check if features have equel length
		if(_fLen)	_unifLen = abs(_fLen-long(rgn.Length())) <= 10;	// consider the difference 10 as a threshold
		else		_fLen = rgn.Length();	// initialize feature's length once
#endif
	Region& currRgn = *it;
	if(rgn == currRgn)					// duplicated feature?
		return spotter.TreatCase(spotter.DUPL) >= 0;
#ifdef _ISCHIP
	if(rgn.Length() < _minFtrLen)		// short feature?
		return spotter.TreatCase(spotter.SHORT) >= 0;
#endif
	if(currRgn.Adjoin(rgn))				// adjacent feature?
		return CorrectItemsEnd(it, rgn.End, spotter.TreatCase(spotter.ADJAC));
	if(currRgn.Cover(rgn))				// covering feature?
		return spotter.TreatCase(spotter.COVER) >= 0;
	if(currRgn.Cross(rgn))				// crossed feature?
		return CorrectItemsEnd(it, rgn.End, spotter.TreatCase(spotter.CROSS));
	return true;
}

// Adds feature to the container
//	@rgn: Region with mandatory fields
//	@file: file to access to additionally fields
//	return: true if Read was added successfully
bool BedF::AddPos(const Region& rgn, TabFile& file)
{
#ifdef _ISCHIP
	readscr score = file.FloatFieldValid(4);
	_items.push_back(Featr(rgn, score));
	if(score > _maxScore)	_maxScore = score;
#else
	_items.push_back(Featr(rgn));
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

// Gets count of features for chromosome or all by default.
//	@cID: chromosome's ID
//	return: count of features for existed chrom or 0, otherwise count of all features
chrlen BedF::FeatureCount(chrid cID) const
{
	if(cID==Chrom::UnID)	return AllItemsCount();
	cIter it = GetIter(cID);
	return it != cEnd() ? FeatureCount(it) : 0;
}

// Gets chromosome's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@it: chromosome's iterator
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
chrlen BedF::EnrRegLength(cIter it, BYTE multiplier, float fLen) const
{
	ChromItemsInd cII = it->second;
	ULONG	res = 0;
	for(chrlen i=cII.FirstInd; i<=cII.LastInd; i++)
		res += _items[i].Length() + int(2*fLen) - 2;	// -2 to consider BS engagement
	return res << multiplier;
}

// Gets chrom's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@cID: chromosome's ID
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
//	return: chrom's total enriched regions length, or 0 if chrom is absent
chrlen BedF::EnrRegLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegLength(it, multiplier, fLen) : 0;
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

// Gets count of treated chroms
chrid BedF::TreatedChromsCount() const
{
	chrid res = 0;
	for(cIter cit=cBegin(); cit!=cEnd(); cit++)	// loop through chroms
		if(TREATED(cit))		res++;
	return res;
}

#else	// NO _ISCHIP


// Copies feature coordinates to external DefRegions.
void BedF::FillRegions(chrid cID, Regions& regn) const
{
	const ChromItemsInd& cii = At(cID);
	regn.Reserve(cii.LastInd - cii.FirstInd + 1);
	//vector<Featr>::const_iterator itEnd = _items.end() + cii.LastInd + 1;
	//for(vector<Featr>::const_iterator it=_items.begin() + cii.FirstInd; it!=itEnd; it++)
	//	regn.AddRegion(it->Start, it->End);
	regn.Copy(_items, cii.FirstInd, cii.LastInd);
}
#endif	// _ISCHIP

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

// Extends all features positions on the fixed length in both directions.
// If extended feature starts from negative, or ends after chrom length, it is fitted.
//	@extLen: distance on which start should be decreased, end should be increased
//	or inside out if it os negative
//	@cSizes: chrom sizes
//	@info: displayed info
//	return: true if instance have been changed
bool BedF::Extend(int extLen, const ChromSizes* cSizes, eInfo info)
{
	if( !extLen )	return false;
	chrlen	rmvCnt;		// counter of removed items in current chrom
	chrlen	allrmvCnt = 0;
	Iter cit;
	ItemsIter fit, fitLast;
	chrlen	cLen = 0;			// chrom length
	Spotter spotter(info, false, FT::BED);

	for(cit=Begin(); cit!=End(); cit++) {	// loop through chroms
		if(cSizes)	cLen = cSizes->At(CID(cit)).Size();
		fit		= ItemsBegin(cit);
		fitLast	= ItemsEnd(cit);
		fit->Extend(extLen, cLen);			// first item
		rmvCnt = 0;
		for(fit++; fit!=fitLast; fit++) {
			fit->Extend(extLen, cLen);		// next item: compare to previous
			if( ( fit->Start < (fit-1)->End	&& spotter.TreatCase(spotter.CROSS) >= 0 )
			|| ( fit->Start == (fit-1)->End && spotter.TreatCase(spotter.ADJAC) >= 0 ) )
			{	// merge crossing/adjacent features
				rmvCnt++;
				(fit-rmvCnt)->End = fit->End;
				fit->Start = UNDEFINED;		// mark item as removed
			}
			else {	allrmvCnt += rmvCnt; rmvCnt=0;	}
		}
	}
	if(rmvCnt) {		// get rid of items merked as removed 
		vector<Featr> newItems;

		newItems.reserve(_itemsCnt - allrmvCnt);
		allrmvCnt = 0;
		// filling newItems by valid values and correct their indexes in the chrom
		for(cit=Begin(); cit!=End(); cit++) {	// loop through chroms
			fit		= ItemsBegin(cit);
			fitLast	= ItemsEnd(cit);
			for(rmvCnt=0; fit!=fitLast; fit++)
				if(fit->Start == UNDEFINED)		rmvCnt++;	// skip removed item
				else			newItems.push_back(*fit);
			cit->second.FirstInd -= allrmvCnt;				// correct indexes
			cit->second.LastInd  -= (allrmvCnt += rmvCnt);
		}
		// replace instance items with new ones
		_items.clear(); 
		_items = newItems;
	}
	spotter.Print(ChromCount()==1 ? CID(Begin()) : Chrom::UnID,
		"after extension", _itemsCnt, _itemsCnt - allrmvCnt);
	PrintEOL(spotter.wasPrinted);
	_itemsCnt -= allrmvCnt;
	return true;
}

// Checks whether all features length exceed gien length, throws exception otherwise.
//	@len: given control length
//	@lenDefinition: control length definition to print in exception message
//	@sender: exception sender to print in exception message
void BedF::CheckFeaturesLength(chrlen len, const string lenDefinition, const char* sender)
{
	ItemsIter fit, fitLast;
	for(cIter cit=Begin(); cit!=End(); cit++) {
		fitLast = ItemsEnd(cit);
		for(fit=ItemsBegin(cit); fit!=fitLast; fit++)
			if(fit->Length()<len)
				Err("Feature size " + NSTR(fit->Length()) +
					" is less than stated " + lenDefinition + sBLANK + NSTR(len),sender).Throw();
	}
}

/************************ end of class BedF ************************/
#endif	// _BEDF

/************************  class ChromSizes ************************/

string ChromSizes::ext;		// FA files real extention

// Initializes chrom sizes from file
void ChromSizes::InitChrSizes(const string& fName)
{
	TabFile file(fName, TxtFile::READ, 2, 2, cNULL, Chrom::Abbr);
	chrid cID;
	ULONG lineCnt;
	// check already done
	if(file.GetFirstLine(&lineCnt)) {
		Reserve(chrid(lineCnt));
		do	// fill by skipping 'random' chromosomes
			if( (cID=Chrom::ValidatedIDbyAbbrName(file.StrField(0))) != Chrom::UnID )
				AddVal(cID, file.LongField(1));
		while(file.GetLine());
	}
#ifdef _ISCHIP
	_treatedCnt = Count();
#endif
}

// Returns length of common prefix before abbr chrom name of all file names
//	@fName: full file name
//	@extLen: length of file name's extention
//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
inline int	ChromSizes::CommonPrefixLength(const string & fName, BYTE extLen)
{
	// a short file name without extention
	return Chrom::PrefixLength(	fName.substr(0, fName.length() - extLen).c_str());
}

// Saves chrom sizes to file
//	@fName: full file name
void ChromSizes::Write(const string& fName) const
{
	ofstream file;

	file.open (fName.c_str(), ios_base::out);
	for(cIter it=cBegin(); it!=cEnd(); it++)
		file << Chrom::AbbrName(CID(it)) << TAB << it->second.Size() << EOL;
	file.close();
}

// Fills external vector by chrom IDs relevant to file's names found in given directory.
//	@cIDs: filling vector of chrom's IDs
//	@gName: path to reference genome
//	@cID: chromosome that sould be treated, or Chrom::UnID if not stated
//	return: count of filled chrom's IDs
BYTE ChromSizes::GetChromIDs(vector<chrid>& cIDs, const string& gName, chrid cID)
{
	vector<string> files;
	if( !FS::GetFiles(files, gName, ext) )		return 0;

	chrid	cid;				// chrom ID relevant to current file in files
	int		prefixLen;			// length of prefix of chrom file name
	BYTE	extLen = ext.length();
	chrid	cnt = chrid(files.size());

	cIDs.reserve(cnt);
	// remove additional names and sort listFiles
	for(chrid i=0; i<cnt; i++) {
		if( (prefixLen = CommonPrefixLength(files[i], extLen)) < 0 )		// right chrom file name
			continue;
		// set to name a chromosome name
		//if( !_prefixName.size() )		// not initialized yet
		//	_prefixName = files[i].substr(0, prefixLen);
		// filter additional names
		cid = Chrom::ValidatedID(files[i].substr(prefixLen, files[i].length() - prefixLen - extLen));
		if( cid != Chrom::UnID 		// "pure" chrom's name
		&& cID == Chrom::UnID )
			cIDs.push_back(cid);
	}
	sort(cIDs.begin(), cIDs.end());
	return cIDs.size();
}

// Creates and initializes an instance
//	@gName: reference genome directory or chrom.sizes file
//	@printMsg: true if print message about chrom.sizes generation (in case of reference genome)
ChromSizes::ChromSizes(const string& gName, bool printMsg) //: _ext(FT::RealExt(FT::FA))
{
	ext = FT::RealExt(FT::FA);
	if( FS::IsDirExist(gName.c_str()) ) {	// gName is a directory
		_path = FS::MakePath(gName);
		const string fname = _path + FS::LastSubDirName(_path) + ".chrom.sizes";
		bool csExist = FS::IsFileExist(fname.c_str());
		vector<chrid> cIDs;		// chrom's ID list
		BYTE cnt;				// number of chroms readed from pointed location 
		// try to fill IDs from .fa files in any case, just to set _ext
		if( !(cnt = GetChromIDs(cIDs, gName, Chrom::CustomID())) ) {
			ext += ZipFileExt;
			// try to fill IDs from .fa.gz files
			if( !csExist && !(cnt = GetChromIDs(cIDs, gName, Chrom::CustomID()) )
			&& Chrom::NoCustom() )					
				Err( Err::MsgNoFiles("*", FT::RealExt(FT::FA)), gName ).Throw();
		}
		if(csExist) { InitChrSizes(fname); return; }
		// chrom.sizes file doesn't exist; generate it
#ifdef _ISCHIP
		_treatedCnt = cnt;
#endif
		Reserve(cnt);
		for(vector<chrid>::const_iterator it = cIDs.begin(); it != cIDs.end(); it++)
			AddVal(*it, FaFile(RefName(*it) + ext).ChromLength());
		Write(fname);
		if(printMsg)
			dout << FS::ShortFileName(fname) << " created\n", fflush(stdout);	// std::endl is unacceptable
	}
	else {		// gName is a chrom.sizes file. Its existence has been verified in main()
		_path = strEmpty;
		InitChrSizes(gName);
	}
}

#ifdef _ISCHIP

//// Gets chrom's defined effective (treated) length
////	@it: ChromSizes iterator
//chrlen ChromSizes::DefEffLength(cIter it) const
//{
//	if(RefSeq::LetGaps)		return it->second.Size() << Autosome(CID(it));
//	ChromDefRegions rgns(RefName(CID(it)));
//	return (rgns.Empty() ? it->second.Size() : (rgns.LastEnd() - rgns.FirstStart()))
//		<< Autosome(CID(it));
//}
//
//// Sets actually treated chromosomes according template and custom chrom
////	@templ: template bed or NULL
////	return: number of treated chromosomes
//chrid ChromSizes::SetTreated(bool statedAll, const Bed* const templ)
//{
//	_treatedCnt = 0;
//
//	for(Iter it = Begin(); it!=End(); it++)
//		if( it->second._treated = Chrom::IsCustom(CID(it)) 
//		&& (statedAll || (!templ || templ->FindChrom(CID(it)))) )
//			_treatedCnt++;
//	return _treatedCnt;
//}
//
//// Prints threated chroms short names
//void ChromSizes::PrintTreatedChroms() const
//{
//	bool next = false;
//	for(cIter it=cBegin(); it!=cEnd(); it++)
//		if(IsTreated(it)) {
//			if(next)	dout << SepCm;
//			dout << Chrom::Mark(CID(it));
//			next = true;
//		}
//}

#elif defined _BIOCC

// Gets total size of genome
genlen ChromSizes::GenSize() const
{
	if( !_gsize )
		for(cIter it=cBegin(); it!=cEnd(); it++)
			_gsize += it->second.Size();
	return _gsize;
}

#endif	// _ISCHIP

#ifdef DEBUG
void ChromSizes::Print() const
{
	cout << "ChromSizes: count: " << int(ChromCount()) << endl;
	cout << "ID\tchrom\t";
//#ifdef _ISCHIP
//		cout << "autosome\t";
//#endif
		cout << "size\n";
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		cout << int(CID(it)) << TAB << Chrom::AbbrName(CID(it)) << TAB;
//#ifdef _ISCHIP
//		cout << int(Autosome(CID(it))) << TAB;
//#endif
		cout << it->second.Size() << EOL;
	}
}
#endif	// DEBUG

/************************  end of ChromSizes ************************/

#ifdef _ISCHIP
/************************  ChromSizesExt ************************/
// Gets chrom's defined effective (treated) length
//	@it: ChromSizes iterator
chrlen ChromSizesExt::DefEffLength(cIter it) const
{
	if(RefSeq::LetGaps)		return it->second.Size() << Autosome(CID(it));
	ChromDefRegions rgns(RefName(CID(it)));
	return (rgns.Empty() ? it->second.Size() : (rgns.LastEnd() - rgns.FirstStart()))
		<< Autosome(CID(it));
}

// Sets actually treated chromosomes according template and custom chrom
//	@templ: template bed or NULL
//	return: number of treated chromosomes
chrid ChromSizesExt::SetTreated(bool statedAll, const Bed* const templ)
{
	_treatedCnt = 0;

	for(Iter it = Begin(); it!=End(); it++)
		if( it->second._treated = Chrom::IsCustom(CID(it)) 
		&& (statedAll || (!templ || templ->FindChrom(CID(it)))) )
			_treatedCnt++;
	return _treatedCnt;
}

// Prints threated chroms short names
void ChromSizesExt::PrintTreatedChroms() const
{
	bool next = false;
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if(IsTreated(it)) {
			if(next)	dout << SepCm;
			dout << Chrom::Mark(CID(it));
			next = true;
		}
}

/************************  end of ChromSizesExt ************************/
#endif	// _ISCHIP

/************************ class RefSeq ************************/

bool RefSeq::LetGaps = true;	// if true then include gaps at the edges of the ref chrom while reading
bool RefSeq::StatGaps = false;	// if true sum gaps for statistic output

// Initializes instance and/or chrom's defined regions
//	@fName: file name
//	@rgns: chrom's defined regions: ripe or new
//	@fill: if true fill sequence and def regions, otherwise def regions only
//	return: true if chrom def regions are stated
bool RefSeq::Init(const string& fName, ChromDefRegions& rgns, bool fill) 
{
	_seq = NULL;
	bool getN = StatGaps || LetGaps || rgns.Empty();	// if true then chrom def regions should be recorded
	FaFile file(fName, rgns.Empty() ? &rgns : NULL);

	_len = file.ChromLength();
	if(fill) {
		try { _seq = new char[_len]; }
		catch(const bad_alloc&) { Err(Err::F_MEM, fName.c_str()).Throw(); }
		const char* line = file.Line();		// First line is readed by FaFile()
		chrlen linelen;
		_len = 0;
		do	memcpy(_seq + _len, line, linelen = file.LineLength()),
			_len += linelen;
		while(line = file.GetLine());
	}
	else if (getN)	while(file.GetLine());	// just to fill chrom def regions
	file.CLoseReading();	// only makes sense if chrom def regions were filled
	_len -= Read::Len;
	return getN;
}

#if defined _ISCHIP || defined _VALIGN

// Creates an stub instance (for sampling cutting)
RefSeq::RefSeq(const ChromSizes& cSizes) : _seq(NULL)
{
	_effDefRgn.Set(0, _len=cSizes[0].Size()-Read::Len);
}

RefSeq::RefSeq(const string& refName)
{
	ChromDefRegions rgns(refName);	// initizlized by file or new (empty)

	if( Init(refName + ChromSizes::RefExt(), rgns, true) && !rgns.Empty() )
		_effDefRgn.Set(rgns.FirstStart(), rgns.LastEnd());
	else
		_effDefRgn.Set(0, Length());
	_gapLen = rgns.GapLen();
}

#elif defined _DENPRO || defined _BIOCC

// Creates an empty instance and fills chrom's defined regions
//	@fName: FA file name
//	@rgns: new chrom's defined regions
//	@minGapLen: minimal length which defines gap as a real gap
RefSeq::RefSeq(const string& fName, ChromDefRegions& rgns, short minGapLen)
{
	Init(fName, rgns, false);
	rgns.Combine(minGapLen);
}

#endif
//#if defined _FILE_WRITE && defined DEBUG
//#define FA_LINE_LEN	50	// length of wrtied lines
//
//void RefSeq::Write(const string & fName, const char *chrName) const
//{
//	FaFile file(fName, chrName);
//	chrlen i, cnt = _len / FA_LINE_LEN;
//	for(i=0; i<cnt; i++)
//		file.AddLine(_seq + i * FA_LINE_LEN, FA_LINE_LEN);
//	file.AddLine(_seq + i * FA_LINE_LEN, _len % FA_LINE_LEN);
//	file.Write();
//}
//#endif	// DEBUG

/************************ end of class RefSeq ************************/

#if defined _DENPRO || defined _BIOCC

/************************ DefRegions ************************/

DefRegions::DefRegions(const ChromSizes& cSizes, chrlen minGapLen)
	: _cSizes(cSizes), _minGapLen(minGapLen)
#ifdef _BIOCC
	, _singleRgn(true)
#endif
{
	if(_cSizes.IsExplicit()) {
		// initialize instance from chrom sizes
		if( Chrom::NoCustom() ) {
			Reserve(cSizes.ChromCount());
			for(ChromSizes::cIter it=cSizes.cBegin(); it != cSizes.cEnd(); it++)
				AddElem(CID(it), Regions(0, it->second.Size()));
		}
		else
			AddElem(Chrom::CustomID(), Regions(0, cSizes[Chrom::CustomID()].Size()));
	}
}

// Gets chrom regions by chrom ID; lazy for real chrom regions
const Regions& DefRegions::operator[] (chrid cID)
{
	if(FindChrom(cID))	return At(cID);
	string name = _cSizes.RefName(cID);
	ChromDefRegions rgns(name, _minGapLen);
	if(rgns.Empty())	// file with def regions doesn't exist?
		RefSeq(name + _cSizes.RefExt(), rgns, _minGapLen);
	return AddElem(cID, rgns);
}

#ifdef _BIOCC
// Gets total genome's size: for represented chromosomes only
genlen DefRegions::GenSize() const
{
	genlen gsize = 0; 
	for(cIter it=cBegin(); it!=cEnd(); it++)
		gsize += Size(it);
	return gsize;
}

// Gets miminal size of chromosome: for represented chromosomes only
chrlen DefRegions::MinSize() const
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
void DefRegions::Print() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout<< Chrom::TitleName(CID(it))
			<< TAB << it->second.FirstStart() 
			<< TAB << Size(it) << EOL;
}
#endif	// DEBUG
/************************ DefRegions: end ************************/
#endif	// _DENPRO || _BIOCC

#if defined _ISCHIP || defined _FRAGPRO
/************************ FragFreq ************************/

// Print statistics
//	@s: print stream
void FragFreq::Print(dostream& s) const
{
	const char* sFragLen = "FragLen";
	float avr = 0;
	ULLONG cnt = 0;
	map<chrlen,ULONG>::const_iterator it;

	//for(it=_freq.begin(); it!=_freq.end(); it++) {
	for(it=begin(); it!=end(); it++) {
		cnt += it->second;
		avr += it->first * it->second;
	}
	s << sFragLen << " average: " << avr/cnt;

	s << EOL << sFragLen << "\tfrequency\n";
	for(it=begin(); it!=end(); it++)
		s << it->first << TAB << it->second << EOL;
}

/************************ FragFreq: end ************************/
#endif	// _ISCHIP