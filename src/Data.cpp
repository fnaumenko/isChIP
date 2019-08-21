/**********************************************************
Data.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 22.06.2019
-------------------------
Provides common data functionality
***********************************************************/

#include "Data.h"
#ifdef _ISCHIP
	#include "isChIP.h"
#endif	//_BIOCC
#include <fstream>	// to write simple files without _FILE_WRITE


static const char* Per = " per ";

/************************ class Obj ************************/

Obj::Spotter::Msg Obj::Spotter::_Msgs [] = {
	{ NULL, "duplicated",			"duplicated", },
	{ NULL, "crossed",				"is intersected with previous" },
	{ NULL, "adjacent",				"is adjacent with previous" },
	{ NULL, "covered",				"is fully covered by previous" },
	{ NULL, "too short",			"too short" },
	{ NULL, "different size of",	"different size of read" },
	{ NULL, "filtered by low score","filtered by score" },
	{ NULL, "chrom exceeding",		"position exceeds chromosome length" },
	{ NULL, "negligible",			"negligible chromosome" },
};
//const BYTE Obj::Spotter::_CasesCnt = sizeof(Obj::Spotter::_Msgs)/sizeof(Obj::Spotter::Msg);

const char*	Obj::Spotter::_ActionMsgs[] = {
	"accepted",
	"joined",
	"rejected",
	"rejected",
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
		//Err(_Msgs[spotter].LineAlarm + BLANK + ItemTitle() + SepCl, _file->LineNumbToStr().c_str()).
		Err(_Msgs[spotter].LineAlarm + BLANK + ItemTitle() + SepCl, _file->LineNumbToStr().c_str()).
			Warning(Message(spotter));
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
		<< sPercent(ULLONG(cnt), ULLONG(allCnt), 4, 0, true) << BLANK
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
// and Features cases by default:
// omitted DUPL cases and handled CROSS and ADJAC 
Obj::Spotter::Spotter (eInfo info, bool alarm, FT::fType format,
	eAction dupl, eAction crossANDadjac, eAction diffsz) :
	_info(info),
	_alarm(alarm),
	_fType(format),
	_alarmPrinted(info <= Obj::iLAC),	// if LAC, do not print EOL at the first time
	_count(CHRLEN_UNDEF),
	_file(NULL),
#ifdef _BIOCC
	_treatcID(vUNDEF),
#endif
	unsortedItems(false),
	wasPrinted(false),
	chrLen(CHRLEN_UNDEF)
{
	noCheck = info != Obj::iSTAT && dupl == ACCEPT && diffsz == ACCEPT && !alarm;
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

	rgn.Set( _file->ItemStart(), _file->ItemEnd());
	if(rgn.Invalid())	ThrowExcept(Err::B_BADEND);
	if(rgn.End>chrLen && TreatCase(EXCEED)<0)	return false;	// exceeding chrom
	if(rgn.End < prevStart)		unsortedItems = true;
	return true;
}

// Prints statistics.
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@title: string at the beginning; if NULL then this instance is used while initialization and don't feeds line
//	@itemCnts: pair of count of all items AND count of accepted items 
//	The last line never ends with EOL 
void Obj::Spotter::Print(chrid cID, const char* title, const p_ulong& itemCnts)
{
	if(_info <= Obj::iLAC || !itemCnts.first)		return;
	bool noSpotters = itemCnts.first == itemCnts.second;

	if(wasPrinted)	dout << EOL;	// "sorting" was printed
	if( title )	{		// additional mode: after extension
		if(_info < Obj::iEXT || noSpotters)		return;	// no ambigs
		dout << "    " << title << SepCl;
		if(_info==Obj::iEXT)
			PrintItems(cID, true, itemCnts.second);
		dout << EOL;
	}
	else {				// main mode: addition to file name
		bool printAccept = _info==Obj::iEXT && !noSpotters;		// print accepted items

		if(_info > Obj::iNM) {
			dout << SepCl << itemCnts.first;
			if(!noSpotters)	dout << BLANK << Total;
		}
		if(printAccept)		dout << SepCm;
		PrintItems(cID, printAccept, itemCnts.second);
		wasPrinted = true;
	}
	if(Count() && _info == Obj::iSTAT) {
		if(!title)		dout << SepCm << "from which" << EOL;
		for(BYTE i=0; i<_CasesCnt; i++)
			PrintCaseStat(static_cast<eCase>(i), itemCnts.first);

		// print ambiguities of negligible chroms if all chroms are readed
		if( cID == Chrom::UnID ) {
			// calculate ambigs for negligible chroms:
			// rest of difference between in & out features minus count of ambigs
			_cases[NEGL_CHR].Count = itemCnts.first - itemCnts.second - Count();
			// correct (add) accepted features
			for(BYTE i=0; i<_CasesCnt-1; i++)	// loop excepting NEGL chroms case
				if( _cases[i].Action == Spotter::ACCEPT )	_cases[NEGL_CHR].Count += _cases[i].Count;
			
			PrintCaseStat(NEGL_CHR, itemCnts.first);
		}
		// print total remained entities
		dout<< TAB << Total << ACCEPTED << SepCl << itemCnts.second
			<< sPercent(ULLONG(itemCnts.second), ULLONG(itemCnts.first), 4, 0, true) << BLANK
			<< ItemTitle(itemCnts.second);
		wasPrinted = true;
	}
	fflush(stdout);
}

#ifdef _BIOCC
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


// Initializes new instance by by tab file name.
//	@title: title printed before file name
//	@fName: name of file
//	@spotter: ambiguities
//	@cSizes: chrom sizes to control chrom length exceedeing
//	@isInfo: true if file info should be printed
//	@abortInval: true if invalid instance shold be completed by throwing exception
//	@scoreInd: index of 'score' filed; is set for FBED only
void Obj::Init	(const char* title, const string& fName, Spotter& spotter,
	ChromSizes& cSizes, bool isInfo, bool abortInval, BYTE scoreInd)
{
	p_ulong items;
	Timer	timer(isInfo);

	if(isInfo) {
		if(title)	dout << title << BLANK;
		dout << fName, fflush(stdout), _EOLneeded = true;
	}
	try {
#ifdef _BAM
		if( spotter.FileType() == FT::BAM) {
			BamInFile file(fName, !isInfo);
			spotter.SetFile(file);
			Chrom::SetRelativeID();		// irrespective of initial setting
#ifdef _FRAGDIST
			Chrom::Validate(file.GetHeaderText());	// validate all chroms ID
#else
			if(!cSizes.IsFilled())	
				cSizes.Init(file.GetHeaderText(), file.ChromCount());

#endif
			items = InitDerived(spotter, cSizes);
		}
		else 
#endif
#ifdef _BIOCC
		if( spotter.FileType() == FT::WIG) {
			TabFile file(fName, spotter.FileType(), abortInval, !isInfo);
			spotter.SetFile((DataFile&)file);
			items = InitDerived(spotter, cSizes);
		}
		else 
#endif
		{
			BedInFile file(fName, spotter.FileType(), scoreInd, abortInval, !isInfo);
			spotter.SetFile(file);
			items = InitDerived(spotter, cSizes);
		}
	}
	catch(Err& err) {	// intercept an exception to manage _isBad and aborting if invalid
		err.Throw(abortInval, !(_isBad = _EOLneeded = true));
	}
#ifdef _BIOCC
	if(abortInval)
		spotter.KeepTreatedChrom();	// save treated chrom for primary
#endif
	if( !_isBad ) {
		if( !items.second ) {		// no items for given chrom
			string specify = ItemTitle(_isBad = _EOLneeded = true);
			if(!Chrom::NoCustom())	specify += Per + Chrom::ShortName(Chrom::CustomID());
			Err(Err::TF_EMPTY, (isInfo ? strEmpty : fName).c_str(), specify).Throw(abortInval, false);
		}
		spotter.Print(Chrom::CustomID(), NULL, items);
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

/************************ class BaseItems ************************/

// Checks if chromosome is uniq and adds it to the container
//	@cID: chroms id
//	@firstInd: first item index
//	@lastInd: last item index
//	@fname: file name for exception message
void BaseItems::AddChrom(chrid cID, chrlen firstInd, chrlen lastInd, const string& fname)
{
	if(FindChrom(cID))
		Err(ItemTitle(true) + " are not consolidated on chromosomes", fname).Throw(true);
	AddVal(cID, ItemIndexes(firstInd, lastInd));
}

// Initializes instance from tab file
//	@spotter: spotter to control ambiguities
//	@cSizes: chrom sizes to control chrom length exceedeing
//	return: numbers of all and number of initialized items
p_ulong BaseItems::InitDerived (Spotter& spotter, const ChromSizes& cSizes)
{
	DataFile& file = spotter.File();
	ULONG	itemCnt = file.ItemCount();

	if( !itemCnt )		return make_pair(0, 0);
		
	chrlen	cntLines  = 0,	// count of lines beginning with 'chr'
			prevStart = 0,	// start previous feature positions
			firstInd  = 0,	// first index in feature's container for current chromosome
			currInd	  = 0;	// current index in Feature's/Read's container.
							// Needed to avoid excess virtual method given current container size
	Region	rgn;			// current feature positions
	chrid	currCID = Chrom::UnID,	// current chromosome's ID
			nextCID;				// next chromosome's ID
	bool	unsorted = false;		// true if chroms are unsorted
	bool	skipChrom = false;

	Reserve(Chrom::NoCustom() ? Chrom::Count : 1);
	ReserveItemContainer(itemCnt);
	
	while(file.GetNextItem()) {
		if( file.IsNextChrom() ) {
			if((nextCID = file.SetNextChrom()) == Chrom::UnID)	continue;	// negligible next chrom
			if(Chrom::NoCustom()) {			// are all chroms specified?
				if(currInd != firstInd)		// have been features for this chrom saved?
					AddChrom(currCID, firstInd, currInd, file.CondFileName());
				if(nextCID < currCID)	unsorted = true;		// unsorted chrom
			}
			else {		// single chrom is specified; items have already been saved,
						// the chrom itself will be saved after loop
				if(rgn.End)			break;
				if(unsorted)
					Err("is unsorted. Single " + Chrom::ShortName(Chrom::CustomID()) + 
						" extraction is allowed only for sorted file", file.CondFileName()).Throw(true);
				if(skipChrom = nextCID != Chrom::CustomID())	continue;
			}
			currCID = nextCID;
			firstInd = currInd;
#ifdef _BIOCC
			spotter.SetTreatedChrom(currCID);
#endif
			if(cSizes.IsFilled())	spotter.chrLen = cSizes[currCID];
			cntLines++;
			if(!spotter.InitRegn(rgn, 0))			continue;
		}
		else {		// the same chrom
			if(skipChrom)							continue;
			cntLines++;
			if(!spotter.InitRegn(rgn, prevStart))	continue;
			if( !CheckLastPos(rgn, spotter))			continue;	// check positions for the same chrom only
		}
#ifdef _VALIGN
		if(AddItem(rgn, file))	currInd++;		// check score
		else	spotter.TreatCase(spotter.SCORE);
#else
		AddItem(rgn, file);		currInd++;		// don't check score
#endif
		prevStart = rgn.Start;
	};

	if( rgn.End && currInd ) {			// some features for at least one valid chrom were saved
		if( currCID != Chrom::UnID )	// is last chrom valid?
			AddChrom(currCID, firstInd, currInd, file.CondFileName());
		if( itemCnt/currInd > 2 )	ShrinkItemContainer();
		if( unsorted )		Sort();		// sort chroms
		if( spotter.unsortedItems ) {
			const bool prInfo = spotter.Info() > Obj::iNONE;
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
void BaseItems::PrintItemCount() const
{
	size_t iCnt = AllItemsCount();
	dout << iCnt << BLANK << ItemTitle(iCnt>1);
	if(ChromCount()==1)		dout << Per << Chrom::TitleName(CID(cBegin()));
	dout << EOL;
}

//#ifdef _READDENS
//// Shifts item's positions to collaps the 'holes' between regions.
////	@cID: chromosome's ID
////	@regns: valid (defined) regions
//void BaseItems::ShrinkByID(chrid cID, const DefRegions &regns)
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
//		Err("item outside last region", "BaseItems::Shrink()").Throw();
//}
//#endif	// _READDENS

#ifdef DEBUG
void BaseItems::PrintChrom() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::AbbrName(CID(it)) << TAB
		<< Data(it).FirstInd << TAB << Data(it).LastInd << SepClTab
		<< Data(it).ItemsCount() << TAB << ItemTitle() << 's' << EOL;
}

void BaseItems::Print(chrlen itemCnt) const
{
	chrlen i, iCnt;
	cout << "BaseItems's ";
	if( itemCnt )	cout << "first " << itemCnt << BLANK;
	cout << ItemTitle() << "s:\n";
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		iCnt = itemCnt ?
			(itemCnt > ItemsCount(CID(it)) ? ItemsCount(CID(it)) : Data(it).FirstInd + itemCnt) - 1:
			Data(it).LastInd;
		for(i=Data(it).FirstInd; i<=iCnt; i++) {
			cout << Chrom::AbbrName(CID(it)) << TAB;
			PrintItem(i);
		}
	}
}
#endif

/************************ end of class BaseItems ************************/

#if !defined _ISCHIP && !defined _WIGREG

/************************ class Reads ************************/

// Checks the element for the new potential start/end positions for all possible ambiguous.
// * different read length
// * duplicated read
//	@rgn: checked start/stop positions
//	@it: iterator reffering to the compared element
//	@spotter: possible ambiguities: different Read length, duplicated Reads
//  return: true if item should be accepted; otherwise false
bool Reads::CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter)
{
	if(_readLen != rgn.Length())					// different Read length?
		if(!_readLen)	_readLen = rgn.Length();	// initialize Read length once
		//else if(spotter.TreatCase(spotter.DIFFSZ) < 0)	return false;
		else spotter.TreatCase(spotter.DIFFSZ);

	if(!spotter.noCheck && !spotter.unsortedItems
	//if( rgn.Start == it->Pos
	&& rgn.Start == it->Pos
	&& _strand == spotter.File().ItemStrand()
	&& spotter.TreatCase(spotter.DUPL) < 0 )	return false;	// duplicated Read?
	// adjacent & crossed Reads are not checked: there are common cases
	// covered Reads are not checked: it's impossible case
	return true;
}

const string NotStated = " is not stated";

// Adds Read to the container.
//	@rgn: Region with mandatory fields
//	@file: file to access to additionally fields
//	return: true if Read was added successfully
bool Reads::AddItem(const Region& rgn, DataFile& file)
{
#if defined _VALIGN || defined _FRAGDIST
#ifdef _FRAGDIST
	if(_isPEonly && !file.ItemIsPaired())
		Err("only paired-end reads are acceptable", file.CondFileName()).Throw();
#endif
	bool strand = file.ItemStrand();
	const char* name = file.ItemName();
	const char* pos = strrchr(name, DOT);	// last DOT
	if(!pos || !isdigit(*(++pos)))
		Err("Cannot find number in the read's name. It should be '*.<number>'",
			file.LineNumbToStr().c_str()).Throw();
#ifdef _FRAGDIST
	_items.push_back( Read(rgn.Start, atol(pos), file.ItemStrand()) );
#else
	const char* cid = strstr(name, Chrom::Abbr);
	if(!cid)
		Err("Cannot find chrom in the read's name. It should be '<chr>*.<number>'",
			file.LineNumbToStr().c_str()).Throw();
	float score = file.ItemScore();
	if(score <= _minScore)	return false;				// pass Read with under-threshhold score
	if(score > _maxScore)	_maxScore = score;

	_items.push_back( Read(rgn.Start, Chrom::IDbyAbbrName(cid), atol(pos), file.ItemStrand(), score) );
#endif	// _FRAGDIST #else
#else	
	_items.push_back( Read(rgn.Start) );
#endif	//  _VALIGN || _FRAGDIST
	_strand = file.ItemStrand();
	return true;
}


// Decreases Read's start position without checkup indexes.
//	@cID: chromosome's ID
//	@rInd: index of read
//	@shift: decrease read's start position on this value
//	@rgEnd: region's end position to control Read location
//	@return: true if Read is insinde this region
//bool Reads::DecreasePos(chrid cID, chrlen rInd, chrlen shift, chrlen rgEnd)
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

/************************ end of class Reads ************************/

#endif	// !_ISCHIP && !_WIGREG

#ifdef _FEATURES
/************************ class Features ************************/

// Sets new end position on the feature if necessary.
//	@rgn: current feature
//	@end: potential new end position
//	@treatCaseRes: result of treatment this spotter
//	return: true if spotter is permitted (feature is valid)
bool Features::CorrectItemsEnd(Region& rgn, chrlen end, int treatCaseRes) {
	if(treatCaseRes > 0)	return true;
	if(!treatCaseRes)		rgn.End = end;		// treatment: merge or join features
		
	return false;
}

// Checks the element for the new potential start/end positions for all possible ambiguous.
// * duplicated features
// * short feature
// * adjacent features
// * crossed features
//	@rgn: checked start/stop positions
//	@it: iterator reffering to the compared element
//	@spotter: possible ambiguities
//  return: true if item should be accepted; otherwise false
bool Features::CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter)
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
		return CorrectItemsEnd(currRgn, rgn.End, spotter.TreatCase(spotter.ADJAC));
	if(currRgn.Cover(rgn))				// covering feature?
		return spotter.TreatCase(spotter.COVER) >= 0;
	if(currRgn.Cross(rgn))				// crossed feature?
		return CorrectItemsEnd(currRgn, rgn.End, spotter.TreatCase(spotter.CROSS));
	return true;
}

// Adds feature to the container
//	@rgn: Region with mandatory fields
//	@file: file to access to additionally fields
//	return: true if Read was added successfully
bool Features::AddItem(const Region& rgn, DataFile& file)
{
#ifdef _ISCHIP
	float score = file.ItemScore();
	_items.push_back(Featr(rgn, score));
	if(score > _maxScore)	_maxScore = score;
#else
	_items.push_back(rgn);
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
//bool Features::DecreasePos(chrid cID, chrlen fInd, chrlen shift, chrlen rgEnd) {
//	//if(Item(cInd, fInd).Start < shift)
//	//	Err("feature outside region", "Region end="+IntToStr(rgEnd)).Throw();
//	if(Item(cID, fInd).End > rgEnd) 
//		return false;
//	chrlen ind = At(cID).FirstInd + fInd;
//	_items[ind].Start -= shift;
//	_items[ind].End -= shift;
//	return true;
//}

// Gets chromosome's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@it: chromosome's iterator
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
chrlen Features::EnrRegLength(cIter it, BYTE multiplier, float fLen) const
{
	ItemIndexes cII = it->second.Data;
	ULONG	res = 0;
	for(chrlen i=cII.FirstInd; i<=cII.LastInd; i++)
		res += _items[i].Length() + int(2*fLen);
	return res << multiplier;
}

// Gets chrom's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@cID: chromosome's ID
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
//	return: chrom's total enriched regions length, or 0 if chrom is absent
chrlen Features::EnrRegLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegLength(it, multiplier, fLen) : 0;
}

#ifdef _ISCHIP
// Scales defined score through all features to the part of 1.
void Features::ScaleScores ()
{
	ItemsIter fit;
	for(cIter cit=Begin(); cit!=End(); cit++)
		for(fit=ItemsBegin(cit); fit!=ItemsEnd(cit); fit++)
			fit->Score /= _maxScore;	// if score is undef then it become 1
}

#else	// NO _ISCHIP


// Copies feature coordinates to external DefRegions.
void Features::FillRegions(chrid cID, Regions& regn) const
{
	const ItemIndexes& cii = At(cID).Data;
	regn.Reserve(cii.LastInd - cii.FirstInd + 1);
	//vector<Featr>::const_iterator itEnd = _items.end() + cii.LastInd + 1;
	//for(vector<Featr>::const_iterator it=_items.begin() + cii.FirstInd; it!=itEnd; it++)
	//	regn.AddRegion(it->Start, it->End);
	regn.Copy(_items, cii.FirstInd, cii.LastInd);
}
#endif	// _ISCHIP

// Return min feature length
chrlen Features::GetMinFeatureLength() const
{
	cItemsIter fit;
	chrlen len, minLen = CHRLEN_MAX;

	for(cIter cit=cBegin(); cit!=cEnd(); cit++) {	// loop through chroms
		const cItemsIter fitLast = cItemsEnd(cit);
		for(fit=cItemsBegin(cit); fit!=fitLast; fit++)
			if((len = fit->Length()) < minLen)	minLen = len;
	}
	return minLen;
}

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

// Return min distance between features boundaries
chrlen Features::GetMinDistance() const
{
	cItemsIter fit;
	chrlen minDist = CHRLEN_MAX, dist, end;			// current feature's end position

	for(cIter cit=cBegin(); cit!=cEnd(); cit++) {	// loop through chroms
		const cItemsIter fitLast = cItemsEnd(cit);
		fit	= cItemsBegin(cit);
		end = fit->End;
		for(fit++; fit!=fitLast; fit++) {
			dist = fit->Start - end;
			if(dist < minDist)	minDist = dist;
			end = fit->End;
		}
	}
	return minDist;
}

// Extends all features positions on the fixed length in both directions.
// If extended feature starts from negative, or ends after chrom length, it is fitted.
//	@extLen: distance on which start should be decreased, end should be increased
//	or inside out if it os negative
//	@cSizes: chrom sizes
//	@info: displayed info
//	return: true if instance have been changed
bool Features::Extend(chrlen extLen, const ChromSizes& cSizes, eInfo info)
{
	if( !extLen )	return false;
	chrlen	rmvCnt;		// counter of removed items in current chrom
	chrlen	allrmvCnt = 0;
	Iter cit;
	ItemsIter fit;
	chrlen	cLen = 0;			// chrom length
	Spotter spotter(info, false, FT::BED);

	for(cit=Begin(); cit!=End(); cit++) {	// loop through chroms
		if(cSizes.IsFilled())	cLen = cSizes[CID(cit)];
		fit	= ItemsBegin(cit);
		const ItemsIter fitLast	= ItemsEnd(cit);
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
			const ItemsIter fitLast	= ItemsEnd(cit);
			for(rmvCnt=0; fit!=fitLast; fit++)
				if(fit->Start == UNDEFINED)		rmvCnt++;	// skip removed item
				else			newItems.push_back(*fit);
			Data(cit).FirstInd -= allrmvCnt;				// correct indexes
			Data(cit).LastInd  -= (allrmvCnt += rmvCnt);
		}
		// replace instance items with new ones
		_items.clear(); 
		_items = newItems;
	}
	spotter.Print(ChromCount()==1 ? CID(Begin()) : Chrom::UnID,
		"after extension", make_pair(_itemsCnt, _itemsCnt - allrmvCnt));
	PrintEOL(spotter.wasPrinted);
	_itemsCnt -= allrmvCnt;
	return true;
}

// Checks whether all features length exceed given length, throws exception otherwise.
//	@len: given control length
//	@lenDefinition: control length definition to print in exception message
//	@sender: exception sender to print in exception message
void Features::CheckFeaturesLength(chrlen len, const string& lenDefinition, const char* sender)
{
	ItemsIter fit, fitLast;
	for(cIter cit=Begin(); cit!=End(); cit++) {
		fitLast = ItemsEnd(cit);
		for(fit=ItemsBegin(cit); fit!=fitLast; fit++)
			if(fit->Length()<len) {
				ostringstream oss;
				oss << "Feature size " << fit->Length() << " is less than stated "
					<< lenDefinition << sBLANK << len;
				Err(oss.str(), sender).Throw();
			}
	}
}

#ifdef _BIOCC
// Gets count of features for chromosome or all by default.
//	@cID: chromosome's ID
//	return: count of features for existed chrom or 0, otherwise count of all features
chrlen Features::Count(chrid cID) const
{
	if(cID==Chrom::UnID)	return AllItemsCount();
	cIter it = GetIter(cID);
	return it != cEnd() ? Count(it) : 0;
}
#endif	// _BIOCC
/************************ end of class Features ************************/
#endif	// _FEATURES

/************************  class ChromSizes ************************/

//string ChromSizes::ext;		// FA files real extention

// Returns length of common prefix before abbr chrom name of all file names
//	@fName: full file name
//	@extLen: length of file name's extention
//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
inline int	ChromSizes::CommonPrefixLength(const string & fName, BYTE extLen)
{
	// a short file name without extention
	return Chrom::PrefixLength(	fName.substr(0, fName.length() - extLen).c_str());
}

// Initializes chrom sizes from file
void ChromSizes::Read(const string& fName)
{
	TabFile file(fName, TxtFile::READ, 2, 2, cNULL, Chrom::Abbr);
	chrid cID;
	ULONG lineCnt;
	// check already done
	if(file.GetFirstLine(&lineCnt)) {
		Reserve(chrid(lineCnt));
		do	// fill by skipping 'random' chromosomes
			if( (cID=Chrom::ValidateIDbyAbbrName(file.StrField(0))) != Chrom::UnID )
				AddValue(cID, ChromSize(file.LongField(1)));
		while(file.GetLine());
	}
}

// Saves chrom sizes to file
//	@fName: full file name
void ChromSizes::Write(const string& fName) const
{
	ofstream file;

	file.open (fName.c_str(), ios_base::out);
	for(cIter it=cBegin(); it!=cEnd(); it++)
		file << Chrom::AbbrName(CID(it)) << TAB << Length(it) << EOL;
	file.close();
}

// Fills external vector by chrom IDs relevant to file's names found in given directory.
//	@cIDs: filling vector of chrom's IDs
//	@gName: path to reference genome
//	return: count of filled chrom's IDs
BYTE ChromSizes::GetChromIDs(vector<chrid>& cIDs, const string& gName)
{
	vector<string> files;
	if( !FS::GetFiles(files, gName, _ext) )		return 0;

	chrid	cid;				// chrom ID relevant to current file in files
	int		prefixLen;			// length of prefix of chrom file name
	BYTE	extLen = _ext.length();
	chrid	cnt = chrid(files.size());
	
	cIDs.reserve(cnt);
	sort(files.begin(), files.end());
	// remove additional names and sort listFiles
	for(chrid i=0; i<cnt; i++) {
		if( (prefixLen = CommonPrefixLength(files[i], extLen)) < 0 )		// right chrom file name
			continue;
		// filter additional names
		cid = Chrom::ValidateID(files[i].substr(prefixLen, files[i].length() - prefixLen - extLen));
		if(cid != Chrom::UnID) 		// "pure" chrom's name
			cIDs.push_back(cid);
	}
	sort(cIDs.begin(), cIDs.end());
	return cIDs.size();
}

// Initializes the paths
//	@gPath: reference genome directory
//	@sPath: service directory
//	@prMsg: true if print message about service fodler and chrom.sizes generation
void ChromSizes::SetPath(const string& gPath, const char* sPath, bool prMsg)
{
	_gPath = FS::MakePath(gPath);
	if(sPath && !FS::CheckDirExist(sPath, false))
		_sPath = FS::MakePath(sPath);
	else
		if(FS::IsDirWritable(_gPath.c_str()))	_sPath = _gPath;
		else {
			_sPath = strEmpty;
			if(prMsg)
				Err("reference folder closed for writing and service folder is not pointed.\n").
					Warning("Service files will not be saved!");
		}
}

// Creates and initializes an instance
//	@gName: reference genome directory or chrom.sizes file
//	@sPath: service directory
//	@prMsg: true if print message about service fodler and chrom.sizes generation
ChromSizes::ChromSizes(const char* gName, const char* sPath, bool prMsg)
{
	_ext = _gPath = _sPath = strEmpty;
	
	if(gName)
		if( FS::IsDirExist(FS::CheckedFileDirName(gName)) ) {	// gName is a directory
			_ext = FT::Ext(FT::FA);
			SetPath(gName, sPath, prMsg);
			const string cName = _sPath + FS::LastDirName(gName) + ".chrom.sizes";
			const bool csExist = FS::IsFileExist(cName.c_str());
			chrid cnt;				// number of chroms readed from pointed location 
			vector<chrid> cIDs;		// chrom's ID fill list

			// fill list with inizialised chrom ID and set _ext
			if( !(cnt = GetChromIDs(cIDs, gName)) ) {				// fill list from *.fa
				_ext += ZipFileExt;			// if chrom.sizes exists, get out - we don't need a list
				if( !csExist && !(cnt = GetChromIDs(cIDs, gName)) )	// fill list from *.fa.gz
					Err( Err::MsgNoFiles("*", FT::Ext(FT::FA)), gName ).Throw();
			}

			if(csExist)		Read(cName);
			else {							// generate chrom.sizes
				Reserve(cnt);
				for(vector<chrid>::const_iterator it = cIDs.begin(); it != cIDs.end(); it++)
					AddValue(*it, ChromSize(FaFile(RefName(*it) + _ext).ChromLength()));
				if(IsServAvail())	Write(cName);
				if(prMsg)
					dout << FS::ShortFileName(cName) << BLANK 
					<< (IsServAvail() ? "created" : "generated") << EOL,
					fflush(stdout);			// std::endl is unacceptable
			}
		}
		else	Read(gName);		// gName is a chrom.sizes file
	else if(sPath)
		_gPath = _sPath = FS::MakePath(sPath);	// initialized be service dir; _ext is empty!
	// else instance remains empty

	Chrom::SetCustomID();
}

#if defined _READDENS || defined _BIOCC || defined _VALIGN
// Initializes empty instance by SAM header data
//	@cCnt: chroms count
void ChromSizes::Init(const string& samHeader, chrid cCnt)
{
	chrid cID;

	Reserve(cCnt);
	for(const char* header = samHeader.c_str();
		header = strstr(header, Chrom::Abbr);
		header = strchr(header, EOL) + 7)
	{
		cID = Chrom::ValidateIDbyAbbrName(header);
		header = strchr(header, TAB) + 4;
		AddValue(cID, atol(header));
	}
	Chrom::SetCustomID();
}
#endif

#ifdef _BIOCC

// Gets total size of genome
genlen ChromSizes::GenSize() const
{
	if( !_gsize )
		for(cIter it=cBegin(); it!=cEnd(); _gsize += Data(it++).Real);
	return _gsize;
}

#endif	// _BIOCC

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
//		cout << int(IsAutosome(CID(it))) << TAB;
//#endif
		cout << Length(it) << EOL;
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
	if(Data(it).Defined)	return Data(it).Defined;	// def.eff. length is initialized
	if(RefSeq::LetGaps)		return SetEffLength(it);	// initialize def.eff. length by real size
	// initialize def.eff. length by chrN.region file
	ChromDefRegions rgns(RefName(CID(it)));
	if(rgns.Empty())		return SetEffLength(it);
	return Data(it).Defined = rgns.DefLength() << int(IsAutosome(CID(it)));
}

// Sets actually treated chromosomes according template and custom chrom
//	@templ: template bed or NULL
//	return: number of treated chromosomes
chrid ChromSizesExt::SetTreated(bool statedAll, const BaseItems* const templ)
{
	_treatedCnt = 0;

	for(Iter it = Begin(); it!=End(); it++)
		_treatedCnt += 
			(it->second.Treated = Chrom::IsCustom(CID(it)) 
			&& (statedAll || (!templ || templ->FindChrom(CID(it)))));
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

// Creates and fills new instance
RefSeq::RefSeq(chrid cID, const ChromSizes& cSizes)
{
	ChromDefRegions rgns(cSizes.ServName(cID));	// read from file or new (empty)

	if( Init(cSizes.RefName(cID) + cSizes.RefExt(), rgns, true)	&& !rgns.Empty() )
		_effDefRgn.Set(rgns.FirstStart(), rgns.LastEnd());
	else
		_effDefRgn.Set(0, Length());
	_gapLen = rgns.GapLen();
}

#elif defined _READDENS || defined _BIOCC

// Creates an empty instance and fills chrom's defined regions
//	@fName: FA file name with extension
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

#if defined _READDENS || defined _BIOCC

/************************ DefRegions ************************/

DefRegions::DefRegions(ChromSizes& cSizes, chrlen minGapLen)
	: _cSizes(cSizes), _minGapLen(minGapLen)
#ifdef _BIOCC
	, _singleRgn(true)
#endif
{
	if(cSizes.IsExplicit()) {
		// initialize instance from chrom sizes
		if( Chrom::NoCustom() ) {
			Reserve(cSizes.ChromCount());
			for(ChromSizes::cIter it=cSizes.cBegin(); it != cSizes.cEnd(); it++)
				AddElem(CID(it), Regions(0, cSizes[CID(it)]));
		}
		else
			AddElem(Chrom::CustomID(), Regions(0, cSizes[Chrom::CustomID()]));
	}
}

// Gets chrom regions by chrom ID; lazy for real chrom regions
const Regions& DefRegions::operator[] (chrid cID)
{
	if(FindChrom(cID))	return At(cID).Data;
	ChromDefRegions rgns(_cSizes.ServName(cID), _minGapLen);
	if(rgns.Empty())		// file with def regions doesn't exist?
	{
		_cSizes.IsFilled();
		const string ext = _cSizes.RefExt();
		if(!ext.length())	// no .fa[.gz] file, empty service dir: _cSizes should be initialized by BAM
			return AddElem(cID, Regions(0, _cSizes[cID])).Data;
			//Err(Err::F_NONE, (_cSizes.ServName(cID) + ChromDefRegions::Ext).c_str()).Throw();
		RefSeq(_cSizes.RefName(cID) + ext, rgns, _minGapLen);
	}
	return AddElem(cID, rgns).Data;
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
			<< TAB << Data(it).FirstStart() 
			<< TAB << Size(it) << EOL;
}
#endif	// DEBUG
/************************ DefRegions: end ************************/
#endif	// _READDENS || _BIOCC

#if defined _ISCHIP || defined _FRAGDIST
/************************ FragFreq ************************/

// Set smm SMMbase (subset length) and clear it
void FragFreq::SMA::SetSize(BYTE halfBase)
{ 
	_size = 2 * halfBase + 1; 
	_sum = _count = 0;
	while(_q.size()) _q.pop();		// clear queque
}

//	Add element and return average
float FragFreq::SMA::Push(ULONG x)
{
	_count++;
	_sum += x;
	_q.push(x);
	if(_q.size() > _size) {	_sum -= _q.front();	_q.pop(); }
	return _count < _size ? 0 : (float)_sum / _size;
}

// Set smm SMMbase (subset length) and clear it
//	@end: 'end' iterator of external collection
void FragFreq::SMM::SetSize(BYTE halfBase, citer end)
{
	_end = end;
	_v.clear();
	_v.reserve(_size = 2 * (_middle = halfBase) + 1);
}

//	Add element and return median
ULONG FragFreq::SMM::Push(citer it)
{
	_v.clear();
	for(size_t i=0; i< _size; i++) {
		_v.push_back(it->second);
		if(++it == _end)
			Err("SMM splicer has reached the end of the collection").Throw();
	}
	sort(_v.begin(), _v.end());
	return _v[_middle];
}

// Calculate and print called lognormal distribution parameters
void FragFreq::CalcDistrParams(dostream& s) const
{
	const float K = 2;	//1.5;		// ratio of the summit height to height of the measuring point
	BYTE	base, SMAbase, SMMbase;	// spliners bases
	int		i = 0;					// multi-purpose counter
	chrlen	x, maxX = 0;			// current point, point with max height (Mode)
	ULONG	y, y0 = 0,				// current height, previous height
			maxY = 0, scatt = 0;	// max height, spread of y values (irregularity)
	float	fy, fy0, halfX,		 	// final point with half height
			fmaxY = 0;				// max height
	citer	it=begin();				// current iterator
	bool	cond = true;			// loop condition
	SMA		sma;
	SMM		smm;

	// STEP 1. Set bases and identify scanning limit (quarter)
	for(ULONG diff, quarter=0; cond; y0=y, it++)
		if((y = it->second) > maxY) {	// increasing part
			if(y < y0 && (diff = y0 - y) > scatt)	scatt = diff;
			quarter = (maxY = y) / 4;
			i++;						// count increasing steps
		}
		else {							// decreasing part
			if(y > y0 && (diff = y - y0) > scatt)	scatt = diff;
			cond = y > quarter;			// loop until 1/4 of summit height
		}
	// set bases
	SMAbase = 100 * float(scatt)/maxY + 2;
	if(i < SMAbase) {
		Err("distribution is cropped to the left; parameters may be inaccurate").Warning();
		SMAbase = i;
	}
	SMMbase = (bool(scatt) + 1) * SMAbase;	// SMAbase for smooth distribution, otherwise 2*SMAbase

	//s << "STEP1: " << maxY << TAB << scatt << TAB << Percent(scatt, maxY) << '%'
	//  << "\tSMAbase = " << int(SMAbase) << "  SMMbase = " << int(SMMbase) << EOL;

	// STEP 2. Splice summit and find Mode
	sma.SetSize(SMAbase);
	if(scatt)	smm.SetSize(base = SMMbase, end());
	else		base = 0;
	for(it=begin(), cond=true; cond; it++) {
		x = it->first - SMAbase + base;
		fy = sma.Push(scatt ? smm.Push(it) : it->second);	// spline by SMA after SMM or just SMA
		if(fy >= fmaxY)		maxX = x, fmaxY = fy;
		else				cond = fy > fmaxY * 4 / 5;
		//s << x << TAB << fy << EOL;
	}
	//s << "\nSUMMIT: " << maxX << TAB << fmaxY << EOL;
	
	// STEP 3. Splice right part and find Mean
	fmaxY /= 2;								// now fmaxY is equal to middle!
	sma.SetSize(SMMbase);
	for(i=0; i<SMMbase; i++)	it--;		// decrease iterator to start with valuable point
	for(; it!=end(); it++)
		if(fy = sma.Push(it->second)) {		// spline
			x = it->first - SMMbase;
			//s << x << TAB << fy << EOL;
			if(fy > fmaxY)	fy0 = fy, halfX = float(x);
			else {
				halfX = fy < fmaxY ? x - (fmaxY - fy)/(fy0 - fy) : x;	// proportional x or x
				break;
			}
		}
	//s << "\nMIDDLE: " << halfX << TAB << fmaxY << EOL;
	
	// STEP 4. Calculate mean and sigma
	fmaxY = log(K*maxX/halfX);		// reuse fmaxY: logarifm of 2*mean / point with half height
	fy = log((float)maxX);			// reuse fy: logarifm of Mode
	fy0 = log(halfX);				// reuse fy0: logarifm of point with half height
	float mean = (fy * fmaxY + (fy0*fy0 - fy*fy)/2) / (fmaxY + fy0 - fy);
	float sigma = sqrt(mean - fy);

	s << "\ncalled lognormal parameters:"
	  << "\nmean: " << mean << "\tsigma: " << sigma
	  << "\nMode: " << maxX << "\t Mean: " << (exp(mean + sigma*sigma/2)) << EOL;
	//s << "half: " << halfX << TAB << halfY << EOL;
}

FragFreq::FragFreq(const char* fName)
{
	TabFile file(fName, TxtFile::READ, 2, 2);
	for(int x; file.GetLine();)
		if(x = file.IntField(0))	// 0 if not integer
			(*this)[x] = file.IntField(1);
}

// Calculate and print dist params
//	@s: print stream
void FragFreq::Print(dostream& s, bool prDistr) const
{
	if(empty())		s << "empty distribution\n";
	else {
		CalcDistrParams(s);
		if(!prDistr)	return;

		const chrlen maxLen = INT_MAX/10;

		s << "\nFragLen\tfrequency\n";
		for(citer it=begin(); it!=end(); it++) {
			if(it->first > maxLen)	break;
			s << it->first << TAB << it->second << EOL;
		}
	}
}

/************************ FragFreq: end ************************/
#endif	// _ISCHIP