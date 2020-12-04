/**********************************************************
Data.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 3.12.2020
-------------------------
Provides common data functionality
***********************************************************/

#include "Data.h"
#ifdef _ISCHIP
	#include "isChIP.h"
#endif	//_BIOCC
#include <fstream>	// to write simple files without _FILE_WRITE
#include <memory>	// uniq_ptr

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
//	@ecase: spotter's case
void Obj::Spotter::PrintLineAlarm(eCase ecase) const
{
	if( _alarm ) {
		if( !_alarmPrinted )	{ dout << LF;	_alarmPrinted = true; }
		Err(_Msgs[ecase].LineAlarm + SPACE + ItemTitle() + SepCl, _file->LineNumbToStr().c_str()).
			Warning(Message(ecase));
	}
}

// Outputs case statistics
//	@ecase: spotter's case
//	@allCnt: total count of ambiguities
//	@total: if true then prints total warning case
void Obj::Spotter::PrintCaseStat(eCase ecase, chrlen allCnt, bool total) const
{
	chrlen cnt = _cases[ecase].Count;	// count of case's ambiguities
	if( !cnt )	return;
	const char* totalAlarm = _Msgs[ecase].TotalAlarm;
	if(total)	dout << "NOTICE: ";
	else		dout << TAB;
	dout<< cnt
		<< sPercent(ULLONG(cnt), ULLONG(allCnt), 4, 0, true) << SPACE
		<< _Msgs[ecase].StatInfo << SPACE
		<< ItemTitle(cnt);

	//if(unsortedItems)		dout << " arisen after sorting";
	if(totalAlarm)	dout << SPACE << totalAlarm;
	dout << SepSCl << Message(ecase);
	if(totalAlarm)	dout << '!';
	dout << LF;
}

const char* ACCEPTED = " accepted";

// Prints accepted items with specifying chrom
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@prAcceptItems: if true then prints number of accepted items
//	@estItemCnt: count of accepted items after treatment
void Obj::Spotter::PrintItems(chrid cID, bool prAcceptItems, long estItemCnt) const
{
	if (prAcceptItems || cID != Chrom::UnID)	dout << estItemCnt;
	if (prAcceptItems)							dout << ACCEPTED;
	if(_info > Obj::eInfo::NM) {
		dout << SPACE << ItemTitle(estItemCnt);
		if( cID != Chrom::UnID )	dout << Per << Chrom::ShortName(cID);
	}
}

// Creates an instance with omitted COVER, SHORT, SCORE and NEGL cases,
// and Features cases by default:
// omitted DUPL cases and handled CROSS and ADJAC 
Obj::Spotter::Spotter (FT::eType format, eInfo info, bool alarm,
	eAction dupl, eAction crossANDadjac, eAction diffsz) :
	_info(info),
	_alarm(alarm),
	_fType(format),
	_alarmPrinted(info <= Obj::eInfo::LAC),	// if LAC, do not print LF at the first time
	_count(CHRLEN_UNDEF),
	_file(NULL),
#ifdef _BIOCC
	_treatcID(vUNDEF),
#endif
	unsortedItems(false),
	hasPrinted(false),
	chrLen(CHRLEN_UNDEF)
{
	noCheck = info != Obj::eInfo::STAT && dupl == eAction::ACCEPT && diffsz == eAction::ACCEPT && !alarm;
	memset(_cases, 0, _CasesCnt*sizeof(Case));	// initialize by 0
	_cases[DUPL].Action = dupl;
	_cases[CROSS].Action = _cases[ADJAC].Action = crossANDadjac;
	_cases[COVER].Action = _cases[SHORT].Action = _cases[SCORE].Action = eAction::OMIT;
	_cases[DIFFSZ].Action = diffsz;
	_cases[EXCEED].Action = eAction::OMIT;
}

// Initializes given Region by second and third current reading line positions, with validating
//	@rgn: Region that should be initialized
//	return: true if Region was initialized successfully
bool Obj::Spotter::InitRegn(Region& rgn)
{
	// dismiss check for negative
	//long start = _file->LongField(1);
	//if(start < 0)		ThrowExcept(Err::B_NEGPOS);
	//rgn.Set(start, _file->LongField(2));

	chrlen end = _file->ItemEnd();
	if (end < rgn.Start)		
		unsortedItems = true;
	rgn.Set( _file->ItemStart(), end);
	if(rgn.Invalid())	ThrowExcept(Err::B_BADEND);
	return (end <= chrLen || TreatCase(EXCEED) >= 0);		// false if exceeding chrom

	//if(end>chrLen && TreatCase(EXCEED)<0)	return false;	// exceeding chrom
	//return true;
}

// Prints statistics.
//	@cID: readed chromosome's ID or Chrom::UnID if all
//	@title: string at the beginning; if NULL then this instance is used while initialization and don't feeds line
//	@itemCnts: pair of count of all items AND count of accepted items 
//	The last line never ends with LF 
void Obj::Spotter::Print(chrid cID, const char* title, const p_ulong& itemCnts)
{
	if(_info <= Obj::eInfo::LAC || !itemCnts.first)		return;
	bool noSpotters = itemCnts.first == itemCnts.second;

	//if(hasPrinted)	dout << LF;	// "sorting" was printed
	if( title )	{		// additional mode: after extension
		if(_info < Obj::eInfo::STD || noSpotters)		return;	// no ambigs
		dout << "    " << title << SepCl;
		if(_info==Obj::eInfo::STD)
			PrintItems(cID, true, itemCnts.second);
		dout << LF;
	}
	else {				// main mode: addition to file name
		//bool printAccept = itemCnts.first != itemCnts.second
		//|| (_info==Obj::eInfo::STD && !noSpotters);		// print accepted items
		bool printAccept = _info == Obj::eInfo::STD && !noSpotters;		// print accepted items

		if(_info > Obj::eInfo::NM ) {
			dout << SepCl;
			if (cID == Chrom::UnID) {
				dout << itemCnts.first;
				if (!noSpotters)	dout << SPACE << sTotal;
			}
			if (printAccept)		dout << SepCm;
			PrintItems(cID, printAccept, itemCnts.second);
			hasPrinted = true;
		}
	}
	if(Count() && _info == Obj::eInfo::STAT) {
		if(!title)		dout << SepCm << "from which" << LF;
		for(BYTE i=0; i<_CasesCnt; i++)
			PrintCaseStat(static_cast<eCase>(i), itemCnts.first);

		// print ambiguities of negligible chroms if all chroms are readed
		if( cID == Chrom::UnID ) {
			// calculate ambigs for negligible chroms:
			// rest of difference between in & out features minus count of ambigs
			_cases[NEGL_CHR].Count = itemCnts.first - itemCnts.second - Count();
			// correct (add) accepted features
			for(BYTE i=0; i<_CasesCnt-1; i++)	// loop excepting NEGL chroms case
				if( _cases[i].Action == Spotter::eAction::ACCEPT )	_cases[NEGL_CHR].Count += _cases[i].Count;
			
			PrintCaseStat(NEGL_CHR, itemCnts.first);
		}
		// print total remained entities
		dout<< TAB << sTotal << ACCEPTED << SepCl << itemCnts.second
			<< sPercent(ULLONG(itemCnts.second), ULLONG(itemCnts.first), 4, 0, true) << SPACE
			<< ItemTitle(itemCnts.second);
		hasPrinted = true;
	}
	//fflush(stdout);
}

#ifdef _BIOCC
// Remember treated chrom.
void Obj::Spotter::SetTreatedChrom(chrid cid)
{
	if(_treatcID == -1)					_treatcID = cid;
	else if(_treatcID != Chrom::UnID)	_treatcID = Chrom::UnID;
}
#endif

// LF rules.
// appearance of info on the screen while Init() and others | ended by LF
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
//	@scoreNumb: number of 'score' filed (for FBED and BedGraph)
void Obj::Init	(const char* title, const string& fName, Spotter& spotter,
	ChromSizes& cSizes, bool isInfo, bool abortInval, BYTE scoreNumb)
{
	p_ulong items;
	Timer	timer(isInfo);

	if(isInfo) {	// print title and file name
		if(title)	dout << title << SPACE;
		dout << fName;	fflush(stdout);	_EOLneeded = true;
	}
	try {
#ifdef _BAM
		if( spotter.FileType() == FT::eType::BAM) {
			BamInFile file(fName, !isInfo);
			spotter.SetFile(file);
			Chrom::SetRelativeID();		// irrespective of initial setting
#ifdef _FRAGDIST
			Chrom::Validate(file.GetHeaderText());	// validate all chroms ID
#else
			if(!cSizes.IsFilled())	
				cSizes.Init(file.GetHeaderText());
#endif
			items = InitDerived(spotter, cSizes);
		}
		else 
#endif	//_BAM
		{
			BedInFile file(fName, spotter.FileType(), scoreNumb, abortInval, !isInfo);
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
	if(timer.IsEnabled())	dout << SPACE;
	timer.Stop(true, false);
	PrintEOL(false);
	}

// Prints LF if needs and flash stdout
//	@printEOL: true if LF should be printed explicitly
void Obj::PrintEOL(bool printEOL)
{
	if(printEOL || _EOLneeded)	dout << LF;
	fflush(stdout);
	_EOLneeded = false;
}

/******************** end of class Obj *********************/

/************************ class BaseItems ************************/

// Checks if items are initialized, chromosome is uniq, and adds it to the instance
//	@cID: chroms id
//	@firstInd: first item index
//	@spotter: spotter to close items
//	@fname: file name for exception message
void BaseItems::AddChrom(chrid cID, chrlen firstInd, const Spotter& spotter)
{
	FillChromItems(spotter);
	chrlen lastInd = ItemsCount();
	if (firstInd == lastInd)	return;	// uninitialized items
	if (FindChrom(cID))					// this chrom is already saved in Chroms
		Err(ItemTitle(true) + " are not consolidated on chromosomes", spotter.File().CondFileName()).Throw(true);
	AddVal(cID, ItemIndexes(firstInd, lastInd + FinishItems(spotter)));
}

// Initializes generalized BED instance from tab file
//	@spotter: spotter to control ambiguities
//	@cSizes: chrom sizes to control chrom length exceedeing
//	return: numbers of all and number of initialized items
// It's separated because it's invoked not only in BaseItems::InitDerived, but also in Cover::InitDerived
p_ulong BaseItems::InitBed(Spotter& spotter, const ChromSizes& cSizes)
{
	DataInFile& file = spotter.File();
	ULONG	estItemCnt = file.EstItemCount();

	if( !estItemCnt )		return make_pair(0, 0);
		
	chrlen	cntLines = 0;		// count of lines beginning with 'chr'
	chrlen	cntAccLines = 0;	// count of accepted lines beginning with 'chr'
	chrlen	firstInd = 0;		// first index in feature's container for current chromosome
	Region	rgn;				// current item positions
	chrid	cID = Chrom::UnID;	// current chromosome's ID
	//bool	unsorted = false;	// true if chroms are unsorted
	bool	skipChrom = false;

	ReserveItems(estItemCnt);

	for (; file.GetNextItem(); ++cntLines) {
		if( file.GetNextChrom() ) {
			chrid nextCID = file.GetChrom();
			if (nextCID == Chrom::UnID)		continue;
			if(Chrom::NoCustom()) {					// are all chroms specified?
				// Last item index is equal of total number of recorded items.
				// In the first pass cID==UnID, but it doesn't added because firstInd==ItemsCount()==0
				AddChrom(cID, firstInd, spotter);
				//if(cID != Chrom::UnID && nextCID < cID)	unsorted = true;		// unsorted chrom
			}
			else {				// single chrom is specified
				if (rgn.End)	// rigion is initialized: items for the specified chrom are already saved
					break;		// the chrom itself will be saved after loop
				//if(unsorted)
				//	Err("is unsorted. Single " + Chrom::ShortName(Chrom::CustomID()) + 
				//		" extraction is allowed only for sorted file", file.CondFileName()).Throw(true);
				if(skipChrom = nextCID != Chrom::CustomID())	
					continue;
			}
			cID = nextCID;
			firstInd = ItemsCount();
			rgn.Start = 0;			// clear for the new chromosome to avoid false sorting
#ifdef _WIG
			spotter.lastEnd = 0;	// clear for the new chromosome (BEDGRAPH)
#endif
#ifdef _BIOCC
			spotter.SetTreatedChrom(cID);
#endif
			if(cSizes.IsFilled())	spotter.chrLen = cSizes[cID];
			if(!spotter.InitRegn(rgn))			continue;
		}
		else {		// the same chrom
			if(skipChrom)						continue;
			if(!spotter.InitRegn(rgn))			continue;
			if( !CheckLastPos(rgn, spotter))	continue;	// check positions for the same chrom only
		}
#ifdef _VALIGN
		if (!AddItem(rgn, spotter))
			spotter.TreatCase(spotter.SCORE);	// check score
		else	cntAccLines++;
#else
		cntAccLines += AddItem(rgn, spotter);				// don't check score
#endif
	}
	// save last chrom
	if( rgn.End && cntAccLines) {	// some features for at least one valid chrom were saved
		if (cID != Chrom::UnID)		// is last chrom valid?
			AddChrom(cID, firstInd, spotter);
		SortIfNecessary(spotter, estItemCnt);
	}
	//cout << " est/fact: " << float(estItemCnt) / ItemsCount() << SPACE;
	return make_pair(cntLines, cntAccLines);
}

// Prints items name and count, adding chrom name if the instance holds only one chrom
//	@prLF: if true then print line feed
void BaseItems::PrintItemCount(bool prLF) const
{
	size_t iCnt = ItemsCount();
	dout << iCnt << SPACE << ItemTitle(iCnt>1);
	if(ChromCount()==1)		dout << Per << Chrom::TitleName(CID(cBegin()));
	if(prLF)	dout << LF;
}

#ifdef _DEBUG
void BaseItems::PrintChrom() const
{
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout << Chrom::AbbrName(CID(it)) << TAB
		<< Data(it).FirstInd << TAB << Data(it).LastInd << SepClTab
		<< Data(it).ItemsCount() << TAB << ItemTitle() << 's' << LF;
}

void BaseItems::Print(const char* title, chrlen estItemCnt) const
{
	chrlen i, iCnt;
	cout << "BaseItems's ";
	if( estItemCnt )	cout << title << SPACE << estItemCnt << SPACE;
	cout << ItemTitle() << "s:\n";
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		iCnt = estItemCnt ?
			(estItemCnt > ItemsCount(CID(it)) ? ItemsCount(CID(it)) : Data(it).FirstInd + estItemCnt) - 1:
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
//	@spotter: temporary values & ambiguities
//	return: true if Read was added successfully
bool Reads::AddItem(const Region& rgn, Spotter& spotter)
{
	CheckStrand(spotter);
	if (_readLen != rgn.Length())					// different Read length?
		if (!_readLen)	_readLen = rgn.Length();	// initialize Read length once
		//else if(spotter.TreatCase(spotter.DIFFSZ) < 0)	return false;
		else spotter.TreatCase(spotter.DIFFSZ);

	_strand = spotter.File().ItemStrand();
#if defined _VALIGN || defined _FRAGDIST
	DataInFile& file = spotter.File();
#ifdef _FRAGDIST
	if (_isPEonly && !file.ItemIsPaired())
		Err("only paired-end reads are acceptable", file.CondFileName()).Throw();
#endif
	const char* name = file.ItemName();
	const char* numb = strrchr(name, DOT);	// "number" position, beginning with last DOT
	if(!numb || !isdigit(*(++numb)))
		Err("Cannot find number in the read's name. It should be '<name>.<number>'",
			file.LineNumbToStr().c_str()).Throw();
#ifdef _FRAGDIST
	_items.push_back( Read(rgn, atol(numb), _strand) );
#else	// _VALIGN
	//const char* cid = strstr(name, Chrom::Abbr);
	const char* cid = strchr(name, Read::NmDelimiter);
	if(cid)		cid = strstr(cid + 1, Chrom::Abbr);		// to be sure that cid pointed to 'chr'
	if(!cid)
		Err("Cannot find chrom in the read's name. It should be '<chr>*.<number>'",
			file.LineNumbToStr().c_str()).Throw();
	const char* pos = strchr(cid, Read::NmPos1Delimiter);	// "position" position, begining with last ':'
	if (!pos || !isdigit(*(++pos)))
		Err("Cannot find position in the read's name. It should be '<name>:<pos>.<number>'",
			file.LineNumbToStr().c_str()).Throw();
	float score = file.ItemValue();
	if(score <= _minScore)	return false;				// pass Read with under-threshhold score
	if(score > _maxScore)	_maxScore = score;

	_items.push_back( Read(rgn, atol(numb), _strand, atol(pos), score, Chrom::IDbyAbbrName(cid)) );
#endif	// _FRAGDIST #else
#else	
	_items.push_back( Read(rgn) );
#endif	//  _VALIGN || _FRAGDIST
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

	// Creates new instance from bed-file with output info.
	// Invalid instance will be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: name of bed-file
	//	@cSizes: chrom sizes to control the chrom length exceedeng; if NULL, no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@abortInval: true if invalid instance should abort excecution
	//	@alarm: true if warning messages should be printed 
	//	@PEonly: true if only paired-end reads are acceptable
	//	@acceptDupl: true if duplicates are acceptable 
	//	@minScore: score threshold (Reads with score <= minScore are skipping)
Reads::Reads(const char* title, const string& fName, ChromSizes& cSizes, eInfo info,
	bool printfName, bool abortInval, bool alarm, bool PEonly, bool acceptDupl, int minScore)
	: _readLen(0), _strand(true)
#if defined _VALIGN || defined _FRAGDIST
	, _isPEonly(PEonly), _paired(false)
#endif
#ifdef _VALIGN
	, _minScore(float(minScore)), _maxScore(0)
#endif
{
	FT::eType type = FT::GetType(fName.c_str());
	if (type == FT::eType::BED)	type = FT::eType::ABED;
	else if (type != FT::eType::BAM)
		Err("wrong extension", printfName ? fName.c_str() : NULL).Throw(abortInval);
	Spotter spotter(type, info, alarm,
		acceptDupl ? Spotter::eAction::ACCEPT : Spotter::eAction::OMIT_SILENT,	// duplicated reads
		Spotter::eAction::ACCEPT,		// crossed & adjacent reads: typical
		//ignoreDiffSize ? Spotter::OMIT : Spotter::ABORTING	// different Read size
		//Spotter::OMIT			// omit different Read size
		Spotter::eAction::ACCEPT			// accept different Read size
	);
	Init(title, fName, spotter, cSizes, info > eInfo::LAC || printfName, abortInval);
}

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
// Items<> abstract method implementation.
//	@rgn: checked start/stop positions
//	@it: iterator reffering to the compared element
//	@spotter: possible ambiguities
//  return: true if item should be accepted; otherwise false
bool Features::CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter)
{
	Region& currRgn = *it;
	if(rgn == currRgn)					// duplicated feature?
		return spotter.TreatCase(spotter.DUPL) >= 0;
#ifdef _ISCHIP
	if(rgn.Length() < _minFtrLen)		// short feature?
		return spotter.TreatCase(spotter.SHORT) >= 0;
#endif
	if(currRgn.Adjoin(rgn))				// adjacent feature?
		return CorrectItemsEnd(currRgn, rgn.End, spotter.TreatCase(spotter.ADJAC));
	if(currRgn.BaseCover(rgn))				// covering feature?
		return spotter.TreatCase(spotter.COVER) >= 0;
	if(currRgn.Cross(rgn))				// crossed feature?
		return CorrectItemsEnd(currRgn, rgn.End, spotter.TreatCase(spotter.CROSS));
	return true;
}

// Adds feature to the container
//	@rgn: Region with mandatory fields
//	@spotter: temporary values & ambiguities
//	return: true if Read was added successfully
bool Features::AddItem(const Region& rgn, Spotter& spotter)
{
#ifdef _BIOCC
	if (_isStrandPres < 0)
		_isStrandPres = spotter.File().IsItemHoldStrand();	// initialize once
#endif
#ifdef _ISCHIP
	float score = spotter.File().ItemValue();
	_items.push_back(Featr(rgn, score));
	if(score > _maxScore)	_maxScore = score;
#else
	_items.push_back(rgn);
#endif
	return true;
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
			fit->Value /= _maxScore;	// if score is undef then it become 1
}

#else	// NO _ISCHIP


// Copies feature coordinates to external DefRegions.
void Features::FillRegions(chrid cID, Regions& regn) const
{
	const ItemIndexes& cii = At(cID).Data;
	regn.Reserve(cii.ItemsCount());
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
		const cItemsIter fitLast = ItemsEnd(cit);
		for(fit=ItemsBegin(cit); fit!=fitLast; fit++)
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
		const cItemsIter fitLast = ItemsEnd(cit);
		fit	= ItemsBegin(cit);
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
	chrlen	rmvCnt = 0, allrmvCnt = 0;	// counters of removed items in current chrom and all removed items
	chrlen	cLen = 0;					// chrom length
	Iter cit;
	ItemsIter fit;
	Spotter spotter(FT::eType::BED, info, false);

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
	size_t itemsCnt = _items.size();
	if(rmvCnt) {		// get rid of items merked as removed 
		vector<Featr> newItems;
		newItems.reserve(itemsCnt - allrmvCnt);
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
		"after extension", make_pair(itemsCnt, itemsCnt - allrmvCnt));
	PrintEOL(spotter.hasPrinted);
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

/************************ end of class Features ************************/
#endif	// _FEATURES

/************************  class ChromSizes ************************/

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
	TabFile file(fName, FT::eType::CSIZE);
	chrid cID;
	// check already done
	while (file.GetNextLine())
		if ((cID = Chrom::ValidateIDbyAbbrName(file.StrField(0))) != Chrom::UnID)
			AddValue(cID, ChromSize(file.LongField(1)));
}

// Saves chrom sizes to file
//	@fName: full file name
void ChromSizes::Write(const string& fName) const
{
	ofstream file;

	file.open (fName.c_str(), ios_base::out);
	for(cIter it=cBegin(); it!=cEnd(); it++)
		file << Chrom::AbbrName(CID(it)) << TAB << Length(it) << LF;
	file.close();
}

// Fills external vector by chrom IDs relevant to file's names found in given directory.
//	@cIDs: filling vector of chrom's IDs
//	@gName: path to reference genome
//	return: count of filled chrom's IDs
chrid ChromSizes::GetChromIDs(vector<chrid>& cIDs, const string& gName)
{
	vector<string> files;
	if( !FS::GetFiles(files, gName, _ext) )		return 0;

	chrid	cid;				// chrom ID relevant to current file in files
	int		prefixLen;			// length of prefix of chrom file name
	chrid	extLen = BYTE(_ext.length());
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
	return chrid(cIDs.size());
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
		if(FS::IsDirWritable(_gPath.c_str()))	
			_sPath = _gPath;
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
//	checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
ChromSizes::ChromSizes(const char* gName, const char* sPath, bool prMsg, bool checkGRef)
{
	_ext = _gPath = _sPath = strEmpty;
	
	if (gName)
		if (FS::IsDirExist(FS::CheckedFileDirName(gName))) {	// gName is a directory
			_ext = FT::Ext(FT::eType::FA);
			SetPath(gName, sPath, prMsg);
			const string cName = _sPath + FS::LastDirName(gName) + FT::Ext(FT::eType::CSIZE);
			const bool csExist = FS::IsFileExist(cName.c_str());
			vector<chrid> cIDs;		// chrom's ID fill list

			// fill list with inizialised chrom ID and set _ext
			if (!GetChromIDs(cIDs, gName)) {				// fill list from *.fa
				_ext += ZipFileExt;			// if chrom.sizes exists, get out - we don't need a list
				if (!csExist && !GetChromIDs(cIDs, gName))	// fill list from *.fa.gz
					Err(Err::MsgNoFiles("*", FT::Ext(FT::eType::FA)), gName).Throw();
			}

			if (csExist)	Read(cName);
			else {							// generate chrom.sizes
				for (vector<chrid>::const_iterator it = cIDs.begin(); it != cIDs.end(); it++)
					AddValue(*it, ChromSize(FaFile(RefName(*it) + _ext).ChromLength()));
				if (IsServAvail())	Write(cName);
				if (prMsg)
					dout << FS::ShortFileName(cName) << SPACE
					<< (IsServAvail() ? "created" : "generated") << LF,
					fflush(stdout);			// std::endl is unacceptable
			}
		}
		else {
			if (checkGRef)	Err("is not a directory", gName).Throw();
			Read(gName);		// gName is a chrom.sizes file
			_sPath = FS::DirName(gName, true);
		}
	else if(sPath)
		_gPath = _sPath = FS::MakePath(sPath);	// initialized be service dir; _ext is empty!
	// else instance remains empty

	Chrom::SetCustomID();
}

#if defined _READDENS || defined _BIOCC || defined _VALIGN
// Initializes empty instance by SAM header data
void ChromSizes::Init(const string& samHeader)
{
	chrid cID;

	for(const char* header = samHeader.c_str();
		header = strstr(header, Chrom::Abbr);
		header = strchr(header, LF) + 7)
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
		cout << Length(it) << LF;
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
			&& (statedAll || !templ || templ->FindChrom(CID(it))));
	return _treatedCnt;
}

inline void PrintChromID(char sep, chrid cID) { dout << sep << Chrom::Mark(cID); }

// Prints threated chroms short names, starting with SPACE
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

	//== define last treated it
	for (cIter it = cBegin(); it != cEnd(); it++)
		if (IsTreated(it))	itLast = it;

	//== print treated chrom
	for(cIter it=cBegin(); it!=cEnd(); it++)
		if(IsTreated(it)) {
			if (it == itLast) {
				char sep;
				if (CID(it) - cID > 1) {
					if (cID != cIDlast)
						PrintChromID(unprintedCnt > 1 ? '-' : COMMA, cID);		// last chrom in the last range
					sep = COMMA;
				}
				else 
					sep = prFirst ? 
						SPACE :								// single
						(unprintedCnt >= 1 ? '-' : COMMA);	// last
				PrintChromID(sep, CID(it));
				break;
			}
			if (prFirst)
				PrintChromID(SPACE, cIDlast = CID(it));
			else 
				if (CID(it) - cID > 1) {
					if (cID != cIDlast)
						PrintChromID(unprintedCnt > 1 ? '-' : COMMA, cID);
					PrintChromID(COMMA, cIDlast = CID(it));
					unprintedCnt = 0;
				}
				else
					unprintedCnt++;
			cID = CID(it);
			prFirst = false;
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
		while(line = file.NextGetLine());
	}
	else if (getN)	while(file.NextGetLine());	// just to fill chrom def regions
	file.CLoseReading();	// only makes sense if chrom def regions were filled
	_len -= Read::FixedLen;
	return getN;
}

#if defined _ISCHIP || defined _VALIGN

// Creates and fills new instance
RefSeq::RefSeq(chrid cID, const ChromSizes& cSizes)
{
	_ID = cID;
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

void DefRegions::Init()
{
	if(IsEmpty())
		if (_cSizes.IsFilled()) {
			// initialize instance from chrom sizes
			if (Chrom::NoCustom())
				for (ChromSizes::cIter it = _cSizes.cBegin(); it != _cSizes.cEnd(); it++)
					AddElem(CID(it), Regions(0, _cSizes[CID(it)]));
			else
				AddElem(Chrom::CustomID(), Regions(0, _cSizes[Chrom::CustomID()]));
			//_isEmpty = false;
		}
}

// Gets chrom regions by chrom ID; lazy for real chrom regions
const Regions& DefRegions::operator[] (chrid cID)
{
	if(FindChrom(cID))	return At(cID).Data;
	ChromDefRegions rgns(_cSizes.ServName(cID), _minGapLen);
	if(rgns.Empty())		// file with def regions doesn't exist?
	{
		//_cSizes.IsFilled();
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

#ifdef _DEBUG
void DefRegions::Print() const
{
	cout << "DefRegions:\n";
	for(cIter it=cBegin(); it!=cEnd(); it++)
		cout<< Chrom::TitleName(CID(it))
			<< TAB << Data(it).FirstStart() 
			<< TAB << Size(it) << LF;
}
#endif	// _DEBUG
/************************ DefRegions: end ************************/
#endif	// _READDENS || _BIOCC


#if defined _ISCHIP || defined _FRAGDIST
/************************ LenFreq ************************/

static const string sDDistrib = "Distribution";
const char* LenFreq::sTitle[] = { "Normal ", "Lognormal " };
const string LenFreq::sSpec[] = {
	strEmpty,
	sDDistrib + " is smooth",
	sDDistrib + " is modulated",
	sDDistrib + " is even",
	sDDistrib + " is cropped to the left",
	sDDistrib + " is heavily cropped to the left",
	sDDistrib + " looks slightly defective on the left",
	sDDistrib + " looks defective on the left"
};
const string LenFreq::sParams = "parameters";
const string LenFreq::sInaccurateParams = sParams + " may be inaccurate";

//	Add element and return average
float LenFreq::SMA::Push(ULONG x)
{
	_sum += x;
	_q.push(x);
	const fraglen size = fraglen(_q.size());
	if (size < _size)	return 0;
	if (size > _size)	_sum -= _q.front(),	_q.pop();
	return (float)_sum / _size;
}

//	Add element and return median
ULONG LenFreq::SMM::Push(citer it)
{
	_v.clear();
	for (fraglen i = 0; i < _size; i++)	_v.push_back(it++->second);
	sort(_v.begin(), _v.end());
	return _v[_middle];
}

//#define _DEBUG

const float SDPI = float(sqrt(3.1415926 * 2));		// square of doubled Pi

// array of lambdas treated like a regular function and assigned to a function pointer
double (*Distrs[])(float, float, fraglen, float) = {
	[](float mean, float sigma, fraglen x, float dsigma2) ->
		double { return exp(-pow(((x - mean) / sigma), 2) / 2) / (sigma * SDPI); },		// normal
	[](float mean, float sigma, fraglen x, float dsigma2) ->
		double { return exp(-pow((log(x) - mean), 2) / dsigma2) / (sigma * SDPI * x); }	// lognormal
};

// Returns estimated slicing base
//	@pointCnt: returned count of actually scanning points
//	return: estimated half base, or 0 in case of degenerate distribution
fraglen LenFreq::GetBase(chrlen& pointCnt)
{
	const int CutoffFrac = 100;	// fraction of the maximum height below which scanning stops on the first pass
	pointCnt = 0;				// count of scanning points
	ULONG	cutoffY = 0;			// Y-value below which scanning stops on the first pass
	fraglen base = 1, halfX = 0;
	USHORT peakCnt = 0;
	bool up = false;
	auto it = begin();
	spoint p0(*it), p;			// previous, current point
	spoint pMin(0, 0), pMax(pMin), pMMax(pMin);	// current, previous, maximum point
	SMA	sma(base);
	vector<spoint> extr;		// local extremums

	//== define pMMax and halfX
	extr.reserve(20);
	for (it++; it != end(); p0 = p, pointCnt++, it++) {
		p.first = it->first;
		p.second = ULONG(sma.Push(it->second));
		//cout << p.first << TAB << p.second << LF;
		if (p.second > p0.second) {		// increasing part
			if (!up) {					// treat pit
				extr.push_back(pMax); pMin = p0; up = true;
			}
		}
		else {							// decreasing part
			if (up) {					// treat peak
				extr.push_back(pMin); pMax = p0; up = false;
				
				if (p0.second > pMMax.second) {
					pMMax = p0;
					cutoffY = pMMax.second / CutoffFrac;
				}
				peakCnt++;
			}
			if (peakCnt && p.second >= pMMax.second / 2)
				halfX = p.first;
			if (p.second < cutoffY) {
				extr.push_back(pMax);
				break;
			}
		}
	}
	if (!halfX || pMMax.second - pMin.second <= 4 ) {	// why 4? 5 maybe enough to identify a peak?
		_spec = eSpec::EVEN;
		return 0;
	}
#ifdef _DEBUG
	cout << "pMMax: " << pMMax.first << TAB << pMMax.second <<LF;
#endif
	//== define splined max point
	pMMax = make_pair(0, 0);
	for (spoint p : extr) {
		if (p.second > pMMax.second)	pMMax = p;
		//cout << p.first << TAB << p.second << LF;
	}

	//== define if sequence is modulated
	auto itv = extr.begin();	// always 0,0
	itv++;						// always 0,0 as well
	p0 = *(++itv);				// first point in sequence
	pMin = pMMax;
	pMax = make_pair(0, 0);
	int i = 1;
	bool isDip = false, isPeakAfterDip = false;

	// from now odd is always dip, even - peak, last - peak
	// looking critical dip in extr
	//s << p0.first << TAB << p0.second << LF;
	for (itv++; itv != extr.end(); i++, itv++) {
		p = *itv;	// we don't need point, using the Y-coordinate is enough. Çoint is used for debugging
		//cout << p.first << TAB << p.second << TAB;
		//if(i % 2)		cout << float(p0.second - p.second) / pMMax.second << "\tdip\n";
		//else			cout << float(p.second - p0.second) / pMMax.second << "\tpeak\n";
		if (i % 2)		// dip
			isDip = float(p0.second - p.second) / pMMax.second > 0.3;
		else 			// peak
			if (isPeakAfterDip = isDip && float(p.second - p0.second) / pMMax.second > 0.1)
				break;
		p0 = p;
	}
	// set  
	if(isPeakAfterDip)	_spec = eSpec::MODUL;

	if (!halfX)		return smoothBase;
	fraglen diffX = halfX - pMMax.first;
#ifdef _DEBUG
	base = float(diffX) * (isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F));
	cout << "isPeakAfterDip: " << isPeakAfterDip << "\thalfX: " << halfX << "\tdiffX: " << diffX << "\tbase: " << base << LF;
	return base;
#else
	return float(diffX) * (isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F));
#endif
}

// Defines key pointers
//	@base: splining base
//	@summit: returned X,Y coordinates of spliced (smoothed) summit
//	@fillSpline: true if fill splining curve (container) to fill
//	@spline: splining curve (container) to fill
//	return: key pointers: X-coord of highest point, X-coord of right middle hight point
LenFreq::fpair LenFreq::GetKeyPoints(fraglen base, point& summit, bool fillSpline, vector<point>& spline) const
{
	const bool isSmooth = base == smoothBase;	// true if input curve is smooth
	fraglen baseSMM = base;
	fraglen shift = isSmooth ? 0 : baseSMM;		// to correct SMM forward shift
	auto it = begin();
	auto End = end();
	point p0(*it), p(0, 0);						// current, previous, maximum point
	SMA	sma(base);
	SMM smm(baseSMM);							// have tried 2*baseSMM, no real difference

	summit.first = 0, summit.second = 0;
	for (short i = smm.Base(); i; i--)	End--;	// reduce End because of SMM forward reading
	//advance(it, size() - smm.Base());			// slower
	for (; it != End; it++) {
		p.first = it->first - base + shift;		// minus SMA base back shift plus SMM base forward shift
		p.second = sma.Push(isSmooth ? it->second : smm.Push(it));	// spline smooth input by SMA only
		//p.first = it->first - base;			// SMA base back shift
		//p.second = sma.Push(it->second);		// spline by SMA
		//p.first = it->first + shift;			// SMM base forward shift
		//p.second = smm.Push(it);				// spline by SMM
		if (fillSpline)		spline.push_back(p);
		if (p.second >= summit.second)	summit = p;
		else if (p.second < summit.second / hRatio)	// to debug print change hRatio to 10
			break;
		p0 = p;
	}
	return make_pair(
		float(summit.first),							// summit X
		p0.first + p0.second / (p.second + p0.second)	// final point with half height; proportional X
	);
}

// Calls distribution parameters
//	@type: NORM or LNORM
//	@keypts: key pointers: X-coord of highest point, X-coord of right middle hight point
//	return: mean & sigma
LenFreq::fpair LenFreq::GetParams(eType type, const fpair& keypts) const
{
	if (keypts.first) {
		if (type == eType::NORM)
			return make_pair(
				keypts.first,													// mean
				sqrt(pow(keypts.second - keypts.first, 2) / (2 * log(hRatio)))	// sigma
			);
		// lognorm
		const float A = log(keypts.first * hRatio / keypts.second);	// logarifm of 2*mean / middle height
		const float lgM = log(keypts.first);						// logarifm of Mode
		const float lgH = log(keypts.second);						// logarifm middle height
		const float mean = (lgM * A + (lgH * lgH - lgM * lgM) / 2) / (A + lgH - lgM);

		return make_pair(mean, sqrt(mean - lgM));
	}
	return make_pair(0, 0);
}

// Compares this sequence with calculated one by given mean & sigma, and returns PCC
//	@type: NORM or LNORM
//	@params: mean & sigma
//	@peakX: X-coordinate of summit
//	@full: if true then correlate from the beginning, otherwiase from summit
//	return: Pearson correlation coefficient,
//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
float LenFreq::CalcPCC(eType type, const fpair& params, fraglen peakX, bool full) const
{
	const float mean = params.first;
	const float sigma = params.second;
	const float dsigma2 = 2 * sigma * sigma;
	const float cutoffY = float(Distrs[int(type)](mean, sigma, peakX, dsigma2) / 1000);	// break when Y became less then 0.1% of max value
	unsigned int cnt = 0;			// count of points
	double	a, b;					// y coordinate (value) of the original and calculated sequence
	double	sumA = 0, sumA2 = 0;	// original and calculated sum of values
	double	sumB = 0, sumB2 = 0;	// original and calculated sum of squares of values
	double	sumAB = 0;				// sum of products of original and calculated values

	// one pass PCC calculation
	for (auto it = begin(); it != end(); it++) {
		if (!full && it->first < peakX)	continue;
		b = Distrs[int(type)](mean, sigma, it->first, dsigma2);
		if (it->first > peakX && b < cutoffY) break;
		sumA += a = it->second;
		sumB += b;
		sumA2 += a * a;
		sumB2 += b * b;
		sumAB += a * b;
		cnt++;
	}
	return float((sumAB * cnt - sumA * sumB) / 
		sqrt((sumA2 * cnt - sumA * sumA) * (sumB2 * cnt - sumB * sumB)));
}

// Returns X-coordinate of the middle height of the left branch
float LenFreq::GetLeftMiddleHalf(vector<point>& spline, const point& summit) const
{
	float maxY = summit.second / 2;

	point p(0, 0);
	float y0 = 0;
	for (auto it = spline.cbegin(); it != spline.cend(); it++)
		if ((p=*it).second > maxY)	break;
		else						y0 = p.second;
	return p.first + y0 / (p.second + y0);		// final point with half height; proportional x
}

// Prints PCC, mean and sigma
//	@s: print stream
//	@type: NORM or LNORM
//	@pcc: printed PCC
void LenFreq::PrintBaseParams(dostream& s, eType type, float pcc) const
{
	s << LF << sTitle[int(type)] << sDistrib << SepCl;
	s << "PCC = " << setprecision(4) << pcc << LF;
	s << sMean << SepCl << _params.first << TAB << sSigma << SepCl << _params.second << LF;
}

// Calculates and print called lognormal distribution parameters
//	@s: print stream
//	@eType: type of distribution
//	@callNorm: call normal parameters anyway
void LenFreq::CallParams(dostream& s, eType type, bool callNorm)
{
	chrlen pointCnt = 0;
	fraglen base = GetBase(pointCnt);
	if (_spec == eSpec::EVEN) {
		Err(Spec(_spec) + SepSCl + sParams + " are not called").Throw(false);
		return;
	}
	if (base == smoothBase)
		s << Spec(eSpec::SMOOTH) << LF;
	fpair keyPts;			// X-coord of highest point, X-coord of right middle hight point
	float Pcc = 0;
	fpair params;			// current mean, sigma
	float leftMH = 0;		// X-coordinate of the middle height of the left branch
	eType dt = type;
	point summit;
	vector<point> spline;	// saved spline to test for normal distribution
	bool fillSpline = false;

	if (type == eType::UNDEF || callNorm) {
		spline.reserve(pointCnt);
		dt = eType::LNORM;
		fillSpline = true;
	}
#ifdef _DEBUG
	fillSpline = true;
	BYTE i = 0;
#endif
	BYTE failCnt = 0;	// count of failure iterations (steps), when PCC is less then previous
	fraglen base0 = 0;
	const bool full = true;
	BYTE failCntLim = 2;

	for (; base; base--) {
		fpair keypts = GetKeyPoints(base, summit, fillSpline, spline);
		params = GetParams(dt, keypts);
		float pcc = CalcPCC(dt, params, summit.first, full);
#ifdef _DEBUG
		s << "base: " << setw(2) << base << "  summitX: " << keypts.first << "\tpcc: " << pcc;
		if (pcc > Pcc)	s << "\t>";
		s << LF;
		i++;
#endif
		if (pcc > Pcc) {
			Pcc = pcc;
			keyPts = keypts;
			_params = params;
			failCnt = 0;
			base0 = base;
		}
		else {
			if (pcc > 0)	failCnt++;		// negative PCC is possible in rare cases
			if (failCnt > 2)	break;
		}
	}
#ifdef _DEBUG
	//s << "spline:\nlength\tfrequency\tbase = " << base << LF;
	//for(spoint p: spline)
	//	if(p.second)	// zero at the begining
	//		s << p.first << TAB << p.second << LF;
	s << "\nbase: " << base0 << "\tsteps: " << int(i) << LF;
#endif
	if (type == eType::UNDEF)
		if (summit.first - begin()->first < SMA::Base(base0)
		|| begin()->second / summit.second > 0.95)
			Err(Spec(eSpec::HCROP) + SepSCl + sInaccurateParams).Warning();
		else if (_spec == eSpec::MODUL)
			s << Spec(_spec) << LF;
		else if(begin()->second / summit.second > 0.5)
			Err(Spec(eSpec::CROP)).Warning();
		else {
			float pcc = CalcPCC(dt, params, summit.first, !full);
#ifdef _DEBUG
			s << "summit: " << keyPts.first << "\tPCCsummit: " << pcc << "\tdiff PCC: " << pcc - Pcc << LF;
#endif
			float diffPCC = pcc - Pcc;
			if (diffPCC > 0.01)
				Err(Spec(eSpec::DEFECT) + SepSCl + sInaccurateParams).Warning();
			else if (diffPCC > 0.002)
				Err(Spec(eSpec::SDEFECT)).Warning();
			leftMH = GetLeftMiddleHalf(spline, summit);
		}

	if (type != eType::NORM) {
		PrintBaseParams(s, eType::LNORM, Pcc);
		// Mode: keypts.first is exactly equal to exp(mean - sigma * sigma)
		s << "Mode: " << keyPts.first << "\t Mean: " << (exp(_params.first + _params.second * _params.second / 2)) << LF;
	}
#ifdef _DEBUG
	if(leftMH)	s << ">> leftMH: " << leftMH << TAB << keyPts.first - leftMH << " / " << keyPts.second - keyPts.first << " = " << (keyPts.first - leftMH) / (keyPts.second - keyPts.first) << LF;
#endif
	if (type == eType::NORM
	|| callNorm
	|| type == eType::UNDEF && leftMH && (keyPts.first - leftMH) / (keyPts.second - keyPts.first) > 0.72) {
		_params = GetParams(eType::NORM, keyPts);
		Pcc = CalcPCC(eType::NORM, _params, summit.first, true);
		PrintBaseParams(s, eType::NORM, Pcc);
	}
}

// Constructor by file name; for test porpose
//	@fname: name of file, paired-end alignment
LenFreq::LenFreq(const char* fName)
{
	TabFile file(fName, FT::eType::DIST);
	for(int x; file.GetNextLine();)
		if(x = file.IntField(0))	// IntField(0) returns 0 if zero field is not an integer
			(*this)[x] = file.IntField(1);
}

// Calculate and print dist params
//	@s: print stream
//	@type: type of distribution
//	@callNorm: call normal parameters anyway
//	@prDistr: if true then print distribution additionally
void LenFreq::Print(dostream& s, eType type, bool callNorm, bool prDistr)
{
	if(empty())		s << "empty " << sDistrib << LF;
	else
		if (IsDegenerate())
			s << "Degenerate " << sDistrib << " (only " << size() << " points)\n";
		else
			CallParams(s, type, callNorm);
		if (prDistr) {
			const chrlen maxLen = INT_MAX / 10;

			s << "\nOriginal " << sDistrib << COLON << "\nlength\tfrequency\n";
			for (auto it = begin(); it != end(); it++) {
				if (it->first > maxLen)	break;
				s << it->first << TAB << it->second << LF;
			}
		}
		fflush(stdout);		// when called from a package
}

/************************ LenFreq: end ************************/
#endif	// _ISCHIP
