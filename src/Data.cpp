#include "Data.h"
#ifdef _ISCHIP
	#include "isChIP.h"
#elif defined _BIOCC
	#include "Calc.h"
#endif	//_BIOCC

static const char* Per = " per ";

/************************ class Obj ************************/

Obj::Ambig::Msg Obj::Ambig::_Msgs [] = {
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
//const BYTE Obj::Ambig::_CasesCnt = sizeof(Obj::Ambig::_Msgs)/sizeof(Obj::Ambig::Msg);

const char*	Obj::Ambig::_ActionMsgs[] = {
	"accepted",
	"joined",
	"omitted",
	"omitted",
	"execution aborted"
};

const Obj::Ambig::ReportCase Obj::Ambig::_Actions[] = {
	&Obj::Ambig::Accept,
	&Obj::Ambig::Handle,
	&Obj::Ambig::Omit,
	&Obj::Ambig::OmitQuiet,
	&Obj::Ambig::Abort
};

// Gets count of ambiguities
chrlen Obj::Ambig::Count() const
{
	if( _count == CHRLEN_UNDEF ) {
		_count = 0;
		for(BYTE i=0; i<_CasesCnt; i++)	_count += _cases[i].Count;
	}
	return _count;
}

// Print given ambiguity as line alarm
//	@ambig: given ambiguity
void Obj::Ambig::PrintLineAlarm(eCase ambig) const
{
	if( _alarm ) {
		if( !_alarmPrinted )	{ dout << EOL;	_alarmPrinted = true; }
		_file->ThrowLineWarning(
			_Msgs[ambig].LineAlarm + BLANK + ItemTitle() + SepCl,
			Message(ambig));
	}
}

// Outputs case statistics
//	@ambig: ambiguity's case
//	@allCnt: total count of ambiguities
//	@total: if true then prints total warning case
void Obj::Ambig::PrintCaseStat(eCase ambig, chrlen allCnt, bool total) const
{
	chrlen cnt = _cases[ambig].Count;	// count of case's ambiguities
	if( !cnt )	return;
	const char* totalAlarm = _Msgs[ambig].TotalAlarm;
	if(total)	dout << Notice;
	else		dout << TAB;
	dout<< cnt
		<< sPercent(ULLONG(cnt), ULLONG(allCnt), 4) << BLANK
		<< _Msgs[ambig].StatInfo << BLANK
		<< ItemTitle(cnt);
	//if(unsortedItems)		dout << " arisen after sorting";
	if(totalAlarm)	dout << BLANK << totalAlarm;
	dout << SepSCl << Message(ambig);
	if(totalAlarm)	dout << '!';
	dout << EOL;
}

const char* ACCEPTED = " accepted";

// Prints accepted items with specifying chrom
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@prAcceptItems: if true then prints number of accepted items
//	@itemCnt: count of accepted items after treatment
void Obj::Ambig::PrintItems(chrid cID, bool prAcceptItems, long itemCnt) const
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
Obj::Ambig::Ambig (eInfo info, bool alarm, FT::eTypes format,
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
	chrLen(0)
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
bool Obj::Ambig::InitRegn(Region& rgn, chrlen prevStart)
{
	long start = _file->LongField(1);
	long end = _file->LongField(2);
	if(start < 0)		ThrowExcept(Err::BP_NEGPOS);
	if(start >= end)	ThrowExcept(Err::BP_BADEND);
	if(chrLen && chrlen(end)>chrLen && TreatCase(EXCEED)<0)	return false;	// exceeding chrom
	
	rgn.Init(start, end);
	if(start < prevStart)			// unsorted feature?
		unsortedItems = true;
	return true;
}

// Prints statistics.
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@title: string at the beginning; if NULL then this instance is used while initialization and don't feeds line
//	@totalItemCnt: count of all items
//	@acceptItemCnt: count of accepted items after treatment
//	The last line never ends with EOL 
void Obj::Ambig::Print(chrid cID, const char* title, ULONG totalItemCnt, ULONG acceptItemCnt)
{
	if(_info <= Obj::iLAC || !totalItemCnt)		return;
	bool noAmbigs = totalItemCnt == acceptItemCnt;

	if(wasPrinted)	dout << EOL;	// "sorting" was printed
	if( title )	{		// additional mode: after extension
		if(_info < Obj::iEXT || noAmbigs)		return;	// no ambigs
		dout << "    " << title << SepCl;
		if(_info==Obj::iEXT)
			PrintItems(cID, true, acceptItemCnt);
		dout << EOL;
	}
	else {				// main mode: addition to file name
		bool printAccept = _info==Obj::iEXT && !noAmbigs;		// print accepted items

		if(_info > Obj::iNM) {
			dout << SepCl << totalItemCnt;
			if(!noAmbigs)	dout << BLANK << Total;
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
				if( _cases[i].Action == Ambig::ACCEPT )	_cases[NEGL].Count += _cases[i].Count;
			
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
void Obj::Ambig::SetTreatedChrom(chrid cid)
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
//	@ambig: ambiguities
//	@addObj: auxiliary object using while initializing
//	@isInfo: true if file info should be printed
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
void Obj::Init	(const char* title, const string& fName, Ambig& ambig, void* addObj,
	bool isInfo, bool abortInvalid)
{
	// in case of bioCC the list is checked already,
	// so better is to add a parameter to Init(),
	// but it is not worth it
	FT::CheckType(fName.c_str(), ambig.FileType(), true, abortInvalid);
	dchrlen items;
	Timer	timer(isInfo);

	if( isInfo ) {
		if(title)	dout << title << BLANK;
		dout << fName;
		fflush(stdout); 
		_EOLneeded = true;
	}
	try {
		TabFile file(fName, FT::FileParams(ambig.FileType()), abortInvalid, !isInfo);	//, false);
		ambig.SetFile(file);
		items = InitChild(ambig, addObj);
	}
	catch(Err& err) {	// intercept an exception to manage _isBad and aborting if invalid
		ThrowError(err, abortInvalid);
	}
#ifndef _ISCHIP
	if(abortInvalid)
		ambig.KeepTreatedChrom();	// save treated chrom for primary
#endif
	if( !_isBad ) {
		if( !items.second ) {		// no items for given chrom
			string sender = isInfo ? strEmpty : fName;
			string specify = ItemTitle(true);
			if(!Chrom::StatedAll())	specify += Per + Chrom::ShortName(Chrom::StatedID());
			Err err(Err::TF_EMPTY, sender.c_str(), specify);
			ThrowError(err, abortInvalid);
		}
		ambig.Print(Chrom::StatedID(), NULL, items.first, items.second);
	}
	if(timer.IsEnabled())	dout << BLANK;
	timer.Stop(true, false);
	PrintEOL(ambig.wasPrinted);	// || ambig.IsAlarmPrinted());
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
//	@file: file to output message
void Bed::AddChrom(chrid cID, chrlen firstInd, chrlen lastInd, const TabFile& file)
{
	if( FindChrom(cID) )
		file.ThrowExcept(ItemTitle(true) + " are not consolidated on chromosomes");
	AddVal(cID, ChromItemsInd(firstInd, lastInd));
	//cout << "add " << Chrom::AbbrName(cID) << TAB << int(cID) << EOL;
}

// Initializes instance from tab file
//	@ambig: ambiguities
//	@pcSizes: chrom sizes
//	return: numbers of all and initialied items for given chrom
dchrlen Bed::InitChild	(Ambig& ambig, void* pcSizes)
{
	ULONG	initSize;		// initial size of _items
	TabFile& file = ambig.File();
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
	bool	unsortedChroms = false;
	bool	cAll = Chrom::StatedAll();				// true if all chroms are stated by user
	bool	skipChrom = false;
	char	cName[Chrom::MaxMarkLength + 1];	// chrom number (only) buffer

	Reserve(cAll ? Chrom::Count : 1);
	ReserveItemContainer(initSize);
	cName[0] = 0;	// empty initial chrom name
	cCurrID = Chrom::ID(file.ChromName());		// first chrom ID

	do {
		if( strcmp(cName, file.ChromName()) ) {	// next chromosome?
			// chrom's name may be long such as 'chrY_random'
			cNextID = Chrom::ID(file.ChromName());
			if( cNextID == Chrom::UnID ) {		// negligible next chromosome?
				cCurrID = cNextID;
				skipChrom = true;
 				continue;
			}
			// now we can save this chrom name
			strcpy(cName, file.ChromName());
			skipChrom = false;
			if(cAll) {			// are all chroms specified?
				if( currInd != firstInd	)		// have been features for this chrom saved?
					AddChrom(cCurrID, firstInd, currInd, file);	// chrom's features have already been saved
				if( cNextID < cCurrID && cNextID != Chrom::M )	// unsorted chrom?
					unsortedChroms = true;
			}
			else {		// single chromosome is specified
				// features in a single defined chrom have already been saved;
				// the chrom proper will be saved after loop
				if(rgn.End)			break;
				if(unsortedChroms)
					file.ThrowExcept(
						"is unsorted. Option --chr " + Chrom::Mark(Chrom::StatedID()) + " is forbidden");
				if(cNextID != Chrom::StatedID()) { skipChrom = true; continue; }
			}
			cCurrID = cNextID;
#ifndef _ISCHIP
			ambig.SetTreatedChrom(cCurrID);
#endif
			firstInd = currInd;
			if(cSizes)	ambig.chrLen = cSizes->Size(cCurrID);
			cntLines++;
			if(!ambig.InitRegn(rgn, 0))				continue;
		}
		else {		// the same chrom
			if(skipChrom)	continue;
			cntLines++;
			if(!ambig.InitRegn(rgn, prevStart))		continue;
			if(!CheckLastPos(rgn, ambig))			continue;	// check positions for the same chrom only
		}
		if(AddPos(rgn, file))	currInd++;
		else					ambig.TreatCase(ambig.SCORE);
		prevStart = rgn.Start;
	}
	while( file.GetLine() );

	if( rgn.End && currInd ) {			// some features for at least one valid chrom were saved
		if( cCurrID != Chrom::UnID )	// is last chrom valid?
			AddChrom(cCurrID, firstInd, currInd, file);	// Last chrom's features have already been saved
		if( initSize/currInd > 2 )	ShrinkItemContainer();
		if( unsortedChroms )		Sort();
		if( ambig.unsortedItems ) {
			bool prInfo = ambig.Info() > Obj::iNONE;
			if(prInfo) {
				file.PrintMsg(ItemTitle(false) + " sorting...", false);
				ambig.wasPrinted = true;
			}
			SortItems(ambig);
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
//	@ambig: possible ambiguities: different Read length, duplicated Reads
//  return: true if item should be accepted; otherwise false
bool BedR::CheckItemPos(ItemsIter it, const Region& rgn, Ambig& ambig)
{
	if(_readLen != rgn.Length())					// different Read length?
		if(!_readLen)	_readLen = rgn.Length();	// initialize Read length once
		else if(ambig.TreatCase(ambig.DIFFSZ) < 0)	return false;

	if( rgn.Start == it->Pos()							// duplicated Read?
	&& ambig.TreatCase(ambig.DUPL) < 0 )	return false;
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
	*cID = CHRID_UNDEF;
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
#endif	// _VALIGN
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
#ifdef _BIOCC
	if( _unifLen )	// check if features have equel length
		if(_fLen)	_unifLen = abs(_fLen-long(rgn.Length())) <= 10;	// consider the difference 10 as a threshold
		else		_fLen = rgn.Length();	// initialize feature's length once
#endif
	Region& currRgn = *it;
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


// Copies feature coordinates to external Regions.
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
	Ambig ambig(info, false, FT::BED);

	for(cit=Begin(); cit!=End(); cit++) {	// loop through chroms
		if(cSizes)	cLen = cSizes->Size(CID(cit));
		fit		= ItemsBegin(cit);
		fitLast	= ItemsEnd(cit);
		fit->Extend(extLen, cLen);			// first item
		rmvCnt = 0;
		for(fit++; fit!=fitLast; fit++) {
			fit->Extend(extLen, cLen);		// next item: compare to previous
			if( ( fit->Start < (fit-1)->End	&& ambig.TreatCase(ambig.CROSS) >= 0 )
			|| ( fit->Start == (fit-1)->End && ambig.TreatCase(ambig.ADJAC) >= 0 ) )
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
	ambig.Print(ChromCount()==1 ? CID(Begin()) : Chrom::UnID,
		"after extension", _itemsCnt, _itemsCnt - allrmvCnt);
	PrintEOL(ambig.wasPrinted);
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

/************************ class Nts ************************/
#define CNT_DEF_NT_REGIONS	10

bool Nts::LetN = true;		// if true then include 'N' at the edges of the ref chrom while reading
bool Nts::StatN = false;	// if true count 'N' for statistic output

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
void Nts::Init(const string& fName, short minGapLen, bool fillNts) 
{
	_nts = NULL;
	const char* line;
	bool allowN = StatN || !LetN;	// if true then 'N' should be recorded
	_defRgns.Reserve(CNT_DEF_NT_REGIONS);
	FaFile::Pocket pocket(_defRgns, minGapLen);
	FaFile file(fName, pocket);

	_cntN = 0;
	_len = pocket.ChromLength();
	if(fillNts) {
		try { _nts = new char[_len]; }
		catch(const bad_alloc&) { Err(Err::F_MEM, fName.c_str()).Throw(); }
		_len = 0;	// is accumulated while reading. Should be restore at the end
	}
	// First line is readed yet
	if(fillNts && !allowN)			// fill nts with N and without defRegions.
		for(line = file.Line(); line; line = file.GetLine() )
			CopyLine(line, file.LineLength());
	else if(minGapLen || allowN) {	// fill nts and defRegions.
		for(line = file.Line(); line; line = file.GetLine(pocket) )
			if(fillNts)		CopyLine(line, file.LineLength());
		pocket.CloseAddN();			// close def regions
		_cntN = pocket.CountN();
	}
	// set _commonDefRgn
	if( allowN && _defRgns.Count() > 0 ) {
		_commonDefRgn.Start = _defRgns.FirstStart();
		_commonDefRgn.End = _defRgns.LastEnd();
	}
	else {
		_commonDefRgn.Start = 0;
		_commonDefRgn.End = _len;
	}
	//cout << "defRgns.Count(): " << _defRgns.Count() << EOL;
	//cout << "Start: " << _commonDefRgn.Start << "\tEnd: " << _commonDefRgn.End << EOL;
}

//#if defined _FILE_WRITE && defined DEBUG
//#define FA_LINE_LEN	50	// length of wrtied lines
//
//void Nts::Write(const string & fName, const char *chrName) const
//{
//	FaFile file(fName, chrName);
//	chrlen i, cnt = _len / FA_LINE_LEN;
//	for(i=0; i<cnt; i++)
//		file.AddLine(_nts + i * FA_LINE_LEN, FA_LINE_LEN);
//	file.AddLine(_nts + i * FA_LINE_LEN, _len % FA_LINE_LEN);
//	file.Write();
//}
//#endif	// DEBUG

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
					if( cID != cid ) { name = strEmpty; wrongNamesCnt++; }
			}
			else { name = strEmpty; wrongNamesCnt++; }	// additional chrom's name
		}
		else { name = strEmpty; wrongNamesCnt++; }		// some other .fa[.gz] file, not chrom
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
#ifdef _ISCHIP
	, _treatedCnt(0)
#endif
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
		listFiles.push_back( Chrom::Mark(cid) );
		_prefixName = fName.substr(0, prefixLen);
		_path = FS::DirName(gName, true);
	}
	// fill attributes
	Reserve(cnt);
	for(vector<string>::const_iterator it = listFiles.begin(); it != listFiles.end(); it++)
		AddChrom(*it);
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
	return FullCommonName() + Chrom::Mark(cID) + _ext;
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
			if( sz < 0 )	Err(Err::F_OPEN, fname.c_str()).Throw();
			it->second._fileLen = chrlen(sz); 
		}
	return cnt;
}

// Gets count of treated chromosomes.
chrid ChromFiles::TreatedCount() const
{
	if(!_treatedCnt) {
		if( _extractAll )
			_treatedCnt = ChromCount();
		else
			for(cIter it=cBegin(); it!=cEnd(); it++)
				if( TREATED(it)() )
					_treatedCnt++;
		}
	return _treatedCnt;
}

// Prints threated chroms short names
void ChromFiles::PrintTreatedNames() const
{
	bool next = false;
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if( TREATED(it)() ) {
			if(next)	dout << SepCm;
			dout << Chrom::Mark(CID(it));
			next = true;
		}
}

// Gets first treated chrom ID
chrid ChromFiles::FirstTreatedChromID() const
{
	if(_extractAll)	return FirstChromID();
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if( TREATED(it)() )	return CID(it);
	return FirstChromID();
}

// Gets first treated chrom ID
chrlen ChromFiles::FirstTreatedFileLength() const
{
	if(_extractAll)	return cBegin()->second._fileLen;
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if( TREATED(it)() )	return it->second._fileLen;
	return 0;
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
	cout << "ChromFiles: count of chroms: " << int(ChromCount()) << endl;
	cout << "chrom\tNumeric\tFileLen\n";
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::AbbrName(CID(it)) << TAB
#ifdef _ISCHIP
		<< int(it->second.Numeric()) << TAB
#endif
		<< it->second._fileLen << EOL;
}
#endif	// DEBUG
/************************  end of class ChromFiles ************************/

/************************ class ChromSizes ************************/
const string ChromSizes::Ext = ".sizes";

// Initializes instance from file.
//	@fName: name of file.sizes
void ChromSizes::Init (const string& fName)
{
	TabFile file(fName, TxtFile::READ, 2, 2, '\0', Chrom::Abbr);
	chrid cID;
	ULONG cntLines;
	// no needs to check since aborting invalid file is set
	const char* currLine = file.GetFirstLine(&cntLines);
	Reserve(chrid(cntLines));
	for(; currLine; currLine = file.GetLine())
		// skip random chromosomes
		if( (cID=Chrom::IDbyAbbrName(file.StrField(0))) != Chrom::UnID )
			AddVal(cID, file.LongField(1));
}

#ifdef _FILE_WRITE
// Saves instance to file
//	@fName: full file name
void ChromSizes::Write(const string fName)
{
	TxtOutFile file(fName);

	file.SetLineBuff(Chrom::MaxNamedPosLength+1);

	// needed to sort because of unoredictable order of chr size adding
#ifdef _NO_MAP
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
#endif

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
	bool print = true;

	if(updated)	Reserve(cFiles.ChromCount());	// file doesn't exist, create it
	else		Init(fName);

	for(ChromFiles::cIter it=cFiles.cBegin(); it!=cFiles.cEnd(); it++)
		if(dontCheck || !this->FindChrom(CID(it))) {
			if(printReport && print) {
				dout << report << Chrom::Title << " sizes file...";
				fflush(stdout);
				print = false;		// don't print title next time
				tm.Start();
			}
			AddValFromFile(CID(it), cFiles);
			updated = true;
		}
#ifdef _FILE_WRITE
	if(updated) {
		Write(fName);
		if(printReport) {
			dout << Done << TAB;
			tm.Stop(true, false);
			dout << EOL;
			fflush(stdout);
		}
	}
#endif
}

// Gets total size of genome
genlen ChromSizes::GenSize() const
{
	if( !_gsize )
		for(cIter it=cBegin(); it!=cEnd(); it++)
			_gsize += it->second;
	return _gsize;
}

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

/************************ DefRegions ************************/

const string DefRegions::DefRegionsFromFile::_FileExt = ".region";

// Creates an instance from file 'chrN.regions', if it exists.
// Otherwise from .fa file then writes it to file 'chrN.regions'
// File 'chrN.regions' is placed in @gName directory
//	@commName: full common name of .fa files
//	@cID: chromosome's ID
//	@minGapLen: minimal length which defines gap as a real gap
DefRegions::DefRegionsFromFile::DefRegionsFromFile(const string& commName, chrlen cID, short minGapLen)
{
	// get the name of regions file. cID is checked already in DefRegions()
	string fName = commName + Chrom::Mark(cID);
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
	Copy(Nts(faFileName, minGapLen).DefRegionsFromFile());
	Write(regionFileName, minGapLen);
}

//void DefRegions::Init(const ChromSizes* cSizes)
//{
//	if( Chrom::StatedAll() ) {
//		Reserve(cSizes->ChromCount());
//		for(ChromSizes::cIter it=cSizes->cBegin(); it != cSizes->cEnd(); it++)
//			AddElem(CID(it), Regions(0, it->second));
//	}
//	else
//		AddElem(Chrom::StatedID(), Regions(0, cSizes->Size(Chrom::StatedID())));
//
//}

DefRegions::DefRegions(const char* gName, ChromSizes** cSizes, short minGapLen)
	: _minGapLen(minGapLen), _singleRgn(FS::HasExt(gName, ChromSizes::Ext))
{
	if(_singleRgn) {
		ChromSizes* cSzs = new ChromSizes(gName);
		// initialize instance from chrom sizes
		if( Chrom::StatedAll() ) {
			Reserve(cSzs->ChromCount());
			for(ChromSizes::cIter it=cSzs->cBegin(); it != cSzs->cEnd(); it++)
				AddElem(CID(it), Regions(0, it->second));
		}
		else
			AddElem(Chrom::StatedID(), Regions(0, cSzs->Size(Chrom::StatedID())));
		*cSizes = cSzs;
	}
	else {
		const ChromFiles cFiles(gName);
		*cSizes = new ChromSizes(cFiles);
		// initialize instance from chrom files
		_commonName = cFiles.FullCommonName();
	}
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