/**********************************************************
Data.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides common data functionality
***********************************************************/

#include "Data.h"
#include <fstream>	// to write simple files without _FILE_WRITE

#ifdef _FEATURES
/************************ class Features ************************/

// Adds chrom to the instance
//	@cID: chrom
//	@cnt: count of chrom items
void Features::AddChrom(chrid cID, size_t cnt)
{
	if (!cnt)	return;
	const chrlen lastInd = _items.size();
	AddVal(cID, ItemIndexes(lastInd - cnt, lastInd));
}

// treats current item
//	return: true if item is accepted
bool Features::operator()()
{
	if (_file->IsJoined()) {
		_items.back().End = _file->ItemEnd();
		return false;
	}
#ifdef _ISCHIP
	float score;
	if(_uniScore) _maxScore = score = 1;
	else {
		score = _file->ItemValue();
		if (score < 0)	_uniScore = _maxScore = score = 1;		// score is undefined in input data
		else if (score > _maxScore)	_maxScore = score;
	}
	_items.emplace_back(_file->ItemRegion(), score);
#else
	_items.emplace_back(_file->ItemRegion());
#endif
	return true;
}

// Gets the sum length of all chromosome's features
//	@cit: chromosome's iterator
chrlen Features::FeaturesLength(cIter cit) const
{
	chrlen res = 0;
	ForChrItems(cit, [&res](cItemsIter it) { res += it->Length(); } );
	return res;
}

// Gets chrom's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@cID: chromosome's ID
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
//	return: chrom's total enriched regions length, or 0 if chrom is absent
chrlen Features::EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegnLength(it, multiplier, fLen) : 0;
}

#ifdef _ISCHIP
// Scales defined score through all features to the part of 1.
void Features::ScaleScores()
{
	for (auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);

		for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
			it->Value /= _maxScore;		// if score is undef then it become 1
	}
}
#else	// NO _ISCHIP

// Copies feature coordinates to external DefRegions.
void Features::FillRegions(chrid cID, Regions& regn) const
{
	const auto& data = At(cID).Data;
	const auto itEnd = ItemsEnd(data);

	regn.Reserve(ItemsCount(data));
	for(auto it = ItemsBegin(data); it!=itEnd; it++)
		regn.Add(*it);
}
#endif	// _ISCHIP

// Return min feature length
chrlen Features::GetMinFeatureLength() const
{
	chrlen len, minLen = CHRLEN_MAX;

	ForAllItems( [&](cItemsIter it) { if ((len = it->Length()) < minLen) minLen = len; } );
	return minLen;
}

// Return min distance between features boundaries
chrlen Features::GetMinDistance() const
{
	chrlen dist, minDist = CHRLEN_MAX;

	for (const auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		for (chrlen end = it++->End; it != itEnd; end = it++->End)
			if ((dist = it->Start - end) < minDist)	minDist = dist;
	}
	return minDist;
}

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

// Extends all features positions on the fixed length in both directions.
// If extended feature starts from negative, or ends after chrom length, it is fitted.
//	@extLen: distance on which start should be decreased, end should be increased
//	or inside out if it os negative
//	@cSizes: chrom sizes
//	@action: action for overlapping features
//	return: true if instance have been changed
bool Features::Extend(chrlen extLen, const ChromSizes& cSizes, UniBedInFile::eAction action)
{
	if (!extLen)	return false;
	size_t	cRmvCnt = 0, tRmvCnt = 0;	// counters of removed items in current chrom and total removed items

	for (auto& c : Container()) {							// loop through chroms
		const chrlen cLen = cSizes.IsFilled() ? cSizes[c.first] : 0;	// chrom length
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		it->Extend(extLen, cLen);			// first item
		cRmvCnt = 0;
		for (it++; it != itEnd; it++) {
			it->Extend(extLen, cLen);						// next item: compare to previous
			if (it->Start <= prev(it)->End)				// overlapping feature
				if (action == UniBedInFile::eAction::JOIN) {
					cRmvCnt++;
					it->Start = UNDEFINED;					// mark item as removed
					(it - cRmvCnt)->End = it->End;
				}
				else if (action == UniBedInFile::eAction::ACCEPT)
					tRmvCnt += cRmvCnt,
					cRmvCnt = 0;
				else if (action == UniBedInFile::eAction::ABORT) {
					//Err("overlapping feature with an additional extension of " + to_string(extLen)).Throw(false, true);
					dout << "overlapping feature with an additional extension of " << extLen << LF;
					return false;
				}
				else if (prev(it)->Start != UNDEFINED)		// OMIT: unmarked item
					cRmvCnt++,
					it->Start = UNDEFINED;	// mark item as removed
		}
	}
	if (cRmvCnt) {		// get rid of items marked as removed 
		vector<Featr> newItems;
		newItems.reserve(_items.size() - tRmvCnt);
		tRmvCnt = 0;
		for (auto& c : Container()) {							// loop through chroms
			ItemIndexes& data = c.second.Data;
			const auto itEnd = ItemsEnd(data);
			cRmvCnt = 0;
			for (auto it = ItemsBegin(data); it != itEnd; it++)
				if (it->Start == UNDEFINED)		cRmvCnt++;	// skip removed item
				else			newItems.emplace_back(*it);
			data.FirstInd -= tRmvCnt;				// correct indexes
			data.LastInd -= (tRmvCnt += cRmvCnt);
		}
		_items.swap(newItems);
	}
	return true;
}

// Checks whether all features length exceed given length, throws exception otherwise.
//	@len: given control length
//	@lenDef: control length definition to print in exception message
//	@sender: exception sender to print in exception message
void Features::CheckFeaturesLength(chrlen len, const string& lenDef, const char* sender) const
{
	ForAllItems( [&](cItemsIter it) {
		if (it->Length() < len) {
			ostringstream oss("Feature size ");
			oss << it->Length() << " is less than stated " << lenDef << SPACE << len;
			Err(oss.str(), sender).Throw();
		}
		} );
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
	TabFile file(fName, FT::eType::CSIZE);	// file check already done

	while (file.GetNextLine()) {
		chrid cID = Chrom::ValidateIDbyAbbrName(file.StrField(0));
		if (cID != Chrom::UnID)
			AddValue(cID, ChromSize(file.LongField(1)));
	}
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
//	@customChrOpt: id of 'custom chrom' option
//	@prMsg: true if print message about service fodler and chrom.sizes generation
//	@sPath: service directory
//	checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
ChromSizes::ChromSizes(const char* gName, BYTE customChrOpt, bool prMsg, const char* sPath, bool checkGRef)
{
	_ext = _gPath = _sPath = strEmpty;
	
	Chrom::SetCustomOption(customChrOpt);
	if (gName) {
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
				for(chrid cid : cIDs)
					AddValue(cid, ChromSize(FaFile(RefName(cid) + _ext).ChromLength()));
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
		Chrom::SetCustomID();
	}
	else if(sPath)
		_gPath = _sPath = FS::MakePath(sPath);	// initialized be service dir; _ext is empty!
	// else instance remains empty
}

// Initializes ChromSizes by SAM header
void ChromSizes::Init(const string& headerSAM)
{
	if (!IsFilled())
		Chrom::ValidateIDs(headerSAM,
			[this](chrid cID, const char* header) { AddValue(cID, atol(header));  }
		);
}

//void ChromSizesInit(ChromSizes* cs, const string& headerSAM)
//{
//	if (!cs->IsFilled())
//		Chrom::ValidateIDs(headerSAM,
//			[cs](chrid cID, const char* header) { cs->AddValue(cID, atol(header));  });
//};

// Gets total size of genome
genlen ChromSizes::GenSize() const
{
	if( !_gsize )
		for(cIter it=cBegin(); it!=cEnd(); _gsize += Length(it++));
	return _gsize;
}

#ifdef _DEBUG
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
	return Data(it).Defined = rgns.DefLength() << int(Chrom::IsAutosome(CID(it)));
}

// Sets actually treated chromosomes according template and custom chrom
//	@templ: template bed or NULL
//	return: number of treated chromosomes
chrid ChromSizesExt::SetTreated(bool statedAll, const Features* const templ)
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
	_seq = nullptr;
	bool getN = StatGaps || LetGaps || rgns.Empty();	// if true then chrom def regions should be recorded
	FaFile file(fName, rgns.Empty() ? &rgns : nullptr);

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
	_len;
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
			if (Chrom::IsCustom())
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
		RefSeq rs(_cSizes.RefName(cID) + ext, rgns, _minGapLen);
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
	for(DefRegions::cIter it=cBegin(); it!=cEnd(); it++)
		cout<< Chrom::TitleName(CID(it))
			<< TAB << Data(it).FirstStart() 
			<< TAB << Size(it) << LF;
}
#endif	// _DEBUG
/************************ DefRegions: end ************************/
#endif	// _READDENS || _BIOCC


#if defined _ISCHIP || defined _CALLDIST
/************************ LenFreq ************************/

const float SDPI = float(sqrt(3.1415926 * 2));		// square of doubled Pi

const float LenFreq::DParams::UndefPCC = -1;
const float LenFreq::lghRatio = float(log(LenFreq::hRatio));	// log of ratio of the summit height to height of the measuring point
const string LenFreq::sParams = "parameters";
const string LenFreq::sInaccurate = " may be inaccurate";
const char* LenFreq::sTitle[] = { "Norm", "Lognorm", "Gamma" };
const string LenFreq::sSpec[] = {
	"is degenerate",
	"is smooth",
	"is modulated",
	"is even",
	"is cropped to the left",
	"is heavily cropped to the left",
	"looks slightly defective on the left",
	"looks defective on the left"
};

// Add last value and pop the first one (QUEUE functionality)
void LenFreq::MW::PushVal(ULONG x)
{
	move(begin() + 1, end(), begin());
	*(end() - 1) = x;
}

//	Add value and return average
float LenFreq::MA::Push(ULONG x)
{
	_sum += x - *begin();
	PushVal(x);
	return float(_sum) / size();
}

//	Add value and return median
ULONG LenFreq::MM::GetMedian(ULONG x)
{
	PushVal(x);
	copy(begin(), end(), _ss.begin());
	sort(_ss.begin(), _ss.end());
	return _ss[size() >> 1];		// mid-size
}

// Constructor
//	@base: moving window half-length;
//	if 0 then empty instance (initialized by empty Push() method)
LenFreq::MM::MM(fraglen base) : MW(base)
{
	if (base)	_ss.insert(_ss.begin(), size(), 0);
	else		_push = &MM::GetMedianStub;
}

// Returns two constant terms of the distrib equation of type, supplied as an index
//	@p: distrib params: mean/alpha and sigma/beta
//	@return: two initialized constant terms of the distrib equation
//	Implemeted as an array of lambdas treated like a regular function and assigned to a function pointer.
fpair (*InitEqTerms[])(const fpair& p) = {
	[](const fpair& p) -> fpair { return { p.second * SDPI, 0}; },							// normal
	[](const fpair& p) -> fpair { return { p.second * SDPI, 2 * p.second * p.second}; },	// lognormal
	[](const fpair& p) -> fpair { return { p.first - 1, float(pow(p.second, p.first)) }; }	// gamma
};

// Returns y-coordinate by x-coordinate of the distrib of type, supplied as an index
//	@p: distrib params: mean/alpha and sigma/beta
//	@x: x-coordinate
//	@eqTerms: two constant terms of the distrib equation
double (*Distrs[])(const fpair& p, fraglen x, const fpair& eqTerms) = {
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// normal
		exp(-pow(((x - p.first) / p.second), 2) / 2) / eqTerms.first; },
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// lognormal
		exp(-pow((log(x) - p.first), 2) / eqTerms.second) / (eqTerms.first * x); },
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// gamma
		pow(x, eqTerms.first) * exp(-(x/ p.second)) / eqTerms.second; }
};

// Returns distribution mode
//	@p: distrib params: mean/alpha and sigma/beta
float(*GetMode[])(const fpair& p) = {
	[](const fpair& p) -> float { return 0; },									// normal
	[](const fpair& p) -> float { return exp(p.first - p.second * p.second); },	// lognormal
	[](const fpair& p) -> float { return (p.first - 1) * p.second; }			// gamma 
};

// Returns distribution Mean (expected value)
//	@p: distrib params: mean/alpha and sigma/beta
float(*GetMean[])(const fpair& p) = {
	[](const fpair& p) -> float { return 0; },										// normal
	[](const fpair& p) -> float { return exp(p.first + p.second * p.second / 2); },	// lognormal
	[](const fpair& p) -> float { return p.first * p.second; }						// gamma 
};

// Calls distribution parameters by consecutive distribution type as index
//	@keypts: key points: X-coord of highest point, X-coord of right middle hight point
//	@p: returned params: mean(alpha) & sigma(beta)
//	Defined as a member of the class only to use the private short 'lghRatio'
void (*LenFreq::SetParams[])(const fpair& keypts, fpair& p) = {
	[](const fpair& keypts, fpair& p) {	//** normal
		p.first = keypts.first;														// mean
		p.second = float(sqrt(pow(keypts.second - p.first, 2) / lghRatio / 2));		// sigma
	},
	[](const fpair& keypts, fpair& p) {	//** lognormal
		const float lgM = log(keypts.first);	// logarifm of Mode
		const float lgH = log(keypts.second);	// logarifm of middle height
		
		p.first = (lgM * (lghRatio + lgM - lgH) + (lgH * lgH - lgM * lgM) / 2) / lghRatio;
		p.second = sqrt(p.first - lgM);
	},
	[](const fpair& keypts, fpair& p) {	//** gamma
		p.second = (keypts.second - keypts.first * (1 + log(keypts.second / keypts.first))) / lghRatio;
		p.first = (keypts.first / p.second) + 1;
	}
};

// Prints restored distr parameters
//	@s: print stream
//	@maxPCC: masimum PCC to print relative PCC percentage
void LenFreq::AllDParams::QualDParams::Print(dostream& s, float maxPCC) const
{
	if (IsSet()) {
		s << Title << TAB;
		if (dParams.IsUndefPcc())
			s << "parameters cannot be called";
		else {
			s << setprecision(5) << dParams.PCC << TAB;
			if (maxPCC) {		// print percent to max PCC
				if (maxPCC != dParams.PCC)	s << setprecision(3) << 100 * ((dParams.PCC - maxPCC) / maxPCC) << '%';
				s << TAB;
			}
			s << setprecision(4) << dParams.Params.first << TAB << dParams.Params.second << TAB;
			const dtype type = GetDType(Type);
			float Mode = GetMode[type](dParams.Params);
			if (Mode)		// for normal Mode is equal to 0
				s << Mode << TAB << GetMean[type](dParams.Params);
		}
		s << LF;
	}
}

// Returns true if distribution parameters set in sorted instance
bool LenFreq::AllDParams::IsSetInSorted(eCType ctype) const
{
	for (const QualDParams& dp : _allParams)
		if (dp.Type == ctype)
			return dp.IsSet();
	return false;
}

// Returns number of distribution parameters set in sorted instance
int LenFreq::AllDParams::SetCntInSorted() const
{
	int cnt = 0;
	for (const QualDParams& dp : _allParams)
		cnt += dp.IsSet();
	return cnt;
}

// Sorts in PCC descending order
void LenFreq::AllDParams::Sort()
{
	if (!_sorted) {
		sort(_allParams.begin(), _allParams.end(),
			[](const QualDParams& dp1, const QualDParams& dp2) -> bool
			{ return dp1.dParams > dp2.dParams; }
		);
		_sorted = true;
	}
}

// Default constructor
LenFreq::AllDParams::AllDParams()
{
	int i = 0;
	for (QualDParams& dp : _allParams)
		dp.SetTitle(i++);
	//for (int i = 0; i < DTCNT; i++)
	//	_allParams[i].SetTitle(i);
}

// Returns parameters of distribution with the highest (best) PCC
//	@dParams: returned PCC, mean(alpha) & sigma(beta)
//	return: consecutive distribution type with the highest (best) PCC
LenFreq::dtype LenFreq::AllDParams::GetBestParams(DParams& dParams)
{
	Sort();
	const QualDParams& QualDParams = _allParams[0];
	dParams = QualDParams.dParams;
	return GetDType(QualDParams.Type);
}

// Prints sorted distibutions params
//	@s: print stream
void LenFreq::AllDParams::Print(dostream& s)
{
	static const char* N[] = { sMean, sSigma };
	static const char* G[] = { "alpha", "beta" };
	static const char* P[] = { "p1", "p2" };
	static const char* a[] = { "*", "**" };
	const bool notSingle = SetCntInSorted() > 1;	// more then 1 output distr type
	float maxPCC = 0;
	bool note = false;

	Sort();			// should be sorted before by PrintTraits(), but just in case
	s << LF << "\t PCC\t";
	if (notSingle)
		s << "relPCC\t",
		maxPCC = _allParams[0].dParams.PCC;
	if (!IsSetInSorted(eCType::GAMMA))	s << N[0] << TAB << N[1];
	else if (note = notSingle)			s << P[0] << a[0] << TAB << P[1] << a[1];
	else								s << G[0] << TAB << G[1];
	if(notSingle || !IsSetInSorted(eCType::NORM))
		s << "\tmode\texp.val";
	s << LF;
	for (const QualDParams& params : _allParams)	params.Print(s, maxPCC);
	if (note) {
		s << LF;
		for (int i = 0; i < 2; i++)
			s << setw(3) << a[i] << P[i] << " - " << N[i] << ", or "
			<< G[i] << " for " << LenFreq::sTitle[GetDType(eCType::GAMMA)] << LF;
	}
}

// Returns estimated moving window half-length ("base")
//	return: estimated base, or 0 in case of degenerate distribution
fraglen LenFreq::GetBase()
{
	const int CutoffFrac = 100;	// fraction of the maximum height below which scanning stops on the first pass
	ULONG	cutoffY = 0;			// Y-value below which scanning stops on the first pass
	fraglen halfX = 0;
	USHORT peakCnt = 0;
	bool up = false;
	auto it = begin();
	spoint p0(*it), p;			// previous, current point
	spoint pMin(0, 0), pMax(pMin), pMMax(pMin);	// current, previous, maximum point
	MA	sma(1);
	vector<spoint> extr;		// local extremums

	//== define pMMax and halfX
	extr.reserve(20);
	for (it++; it != end(); p0 = p, it++) {
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
	if (!halfX || pMMax.second - pMin.second <= 4 ) {	// why 4? 5 maybe enough to identify a peak
		//_spec = eSpec::EVEN;
		Err(Spec(_spec) + SepSCl + sParams + " are not called").Throw(false);	// even distribution
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
	fraglen base = fraglen((isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F)) * diffX);
	cout << "isPeakAfterDip: " << isPeakAfterDip << "\thalfX: " << halfX << "\tdiffX: " << diffX << "\tbase: " << base << LF;
	return base;
#else
	return fraglen(float(diffX) * (isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F)));
#endif
}

// Defines key points
//	@base: moving window half-length
//	@summit: returned X,Y coordinates of spliced (smoothed) summit
//	return: key points: X-coord of highest point, X-coord of right middle hight point
fpair LenFreq::GetKeyPoints(fraglen base, point& summit) const
{
	const fraglen baseSMM = base <= smoothBase ? 0 : base;
	point p0(*begin()), p(0, 0);			// previous, current point
	MA ma(base);
	MM mm(baseSMM);

	summit.second = 0;
#ifdef _DEBUG
	_spline.clear();
	fpair keyPts(0, 0);
#endif
	for (const value_type& f : *this) {
		p.first = f.first - base - baseSMM;		// X: minus MA & MM base back shift
		p.second = ma.Push(mm.Push(f.second));	// Y: splined
#ifdef _DEBUG
		// *** to print splicing individually
		//p.first = f.first - base;		p.second = sma.Push(f.second);	// spline by MA
		//p.first = f.first - baseSMM;	p.second = smm.Push(f.second);	// spline by MM
		if (_fillSpline)	_spline.push_back(p);
#endif
		if (p.second >= summit.second)	summit = p;
		else {
			if (p.second < summit.second / hRatio) {
#ifdef _DEBUG
				if (!keyPts.first) {
					keyPts.first = float(summit.first);
					keyPts.second = p0.first + p0.second / (p.second + p0.second);
				}
			}
			if (p.second < summit.second / (hRatio * 5)) {
#endif
				break;
			}
			p0 = p;
		}
	}

#ifdef _DEBUG
	return keyPts;
#else
	return fpair(
		float(summit.first),							// summit X
		p0.first + p0.second / (p.second + p0.second)	// final point with half height; proportional X
	);
#endif
}

// Compares this sequence with calculated one by given mean & sigma, and returns PCC
//	@type: consecutive distribution type
//	@dParams: returned PCC, input mean(alpha) & sigma(beta)
//	@Mode: X-coordinate of summit
//	@full: if true then correlate from the beginning, otherwiase from summit
//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
void LenFreq::CalcPCC(dtype type, DParams& dParams, fraglen Mode, bool full) const
{
	const fpair eqTerms = InitEqTerms[type](dParams.Params);		// two constant terms of the distrib equation
	const double cutoffY = Distrs[type](dParams.Params, Mode, eqTerms) / 1000;	// break when Y became less then 0.1% of max value
	double	sumA = 0, sumA2 = 0;	// sum, sum of squares of original values
	double	sumB = 0, sumB2 = 0;	// sum, sum of squares of calculated values
	double	sumAB = 0;				// sum of products of original and calculated values
	UINT	cnt = 0;				// count of points

	// one pass PCC calculation
	dParams.SetUndefPcc();
	for (const value_type& f : *this) {
		if (!full && f.first < Mode)		continue;
		const double b = Distrs[type](dParams.Params, f.first, eqTerms);	// y-coordinate (value) of the calculated sequence
		if (isNaN(b))						return;
		if (f.first > Mode && b < cutoffY)	break;
		const double a = f.second;											// y-coordinate (value) of the original sequence
		sumA += a;
		sumB += b;
		sumA2 += a * a;
		sumB2 += b * b;
		sumAB += a * b;
		cnt++;
	}
	float pcc = float((sumAB * cnt - sumA * sumB) /
		sqrt((sumA2 * cnt - sumA * sumA) * (sumB2 * cnt - sumB * sumB)));
	if (!isNaN(pcc))	dParams.PCC = pcc;
}

// Calculates distribution parameters
//	@type: consecutive distribution type
//	@keyPts: key points: X-coord of highest point, X-coord of right middle hight point
//	@dParams: returned PCC, mean(alpha) & sigma(beta)
//	@Mode: X-coordinate of summit
void LenFreq::SetPCC(dtype type, const fpair& keypts, DParams& dParams, fraglen Mode) const
{
	SetParams[type](keypts, dParams.Params);
	CalcPCC(type, dParams, Mode);
}

// Calculates and print called distribution parameters
//	@type: consecutive distribution type
//	@base: moving window half-length
//	@summit: returned X,Y coordinates of spliced (smoothed) summit
void LenFreq::CallParams(dtype type, fraglen base, point& summit)
{
	const BYTE failCntLim = 2;	// max count of base's decreasing steps after which PCC is considered only decreasing
	BYTE failCnt = 0;			// counter of base's decreasing steps after which PCC is considered only decreasing
	point summit0;				// temporary summit
	DParams dParams0, dParams;	// temporary, final  PCC & mean(alpha) & sigma(beta)
#ifdef _DEBUG
	int i = 0;					// counter of steps
#endif

	//base = 9;					// to print spline for fixed base
	// progressive calculate PCC with unknown key points & summit
	for (; base; base--) {
	//for (int i=0; !i; i++) {	// to print spline for fixed base
		const fpair keypts = GetKeyPoints(base, summit0);
		SetPCC(type, keypts, dParams0, summit0.first);
#ifdef _DEBUG
		*_s << setw(4) << setfill(SPACE) << left << ++i;
		*_s << "base: " << setw(2) << base << "  summitX: " << keypts.first << "\tpcc: " << dParams0.PCC;
		if (dParams0 > dParams)	*_s << "\t>";
		*_s << LF;
		if (_fillSpline) { for (point p : _spline)	*_s << p.first << TAB << p.second << LF; _fillSpline = false; }
#endif
		if (dParams0 > dParams) {
			dParams = dParams0;
			summit = summit0;
			failCnt = 0;
		}
		else {
			if (dParams0.PCC > 0)	failCnt++;		// negative PCC is possible in rare cases
			else if (dParams0.IsUndefPcc()) {
				dParams.SetUndefPcc();
				break;
			}
			if (failCnt > failCntLim)	break;
		}
	}
	_allParams.SetParams(type, dParams);

#ifdef _DEBUG
	*_s << LF;
#endif
}

// Prints original distrib features
//	@s: print stream
//	@base: moving window half-length
//	@summit: X,Y coordinates of spliced (smoothed) summit
void LenFreq::PrintTraits(dostream& s, fraglen base, const LenFreq::point& summit)
{
	if (base == smoothBase)		s << Spec(eSpec::SMOOTH) << LF;
	if (summit.first - begin()->first < MA::Size(base)
		|| begin()->second / summit.second > 0.95)
		Err(Spec(eSpec::HCROP) + SepSCl + sParams + sInaccurate).Warning();
	else if (_spec == eSpec::MODUL)
		s << Spec(_spec) << LF;
	else if (begin()->second / summit.second > 0.5)
		Err(Spec(eSpec::CROP)).Warning();
	else {
		DParams dParams, dParams0;
		CalcPCC(_allParams.GetBestParams(dParams), dParams0, summit.first, false);	// sorts params
		float diffPCC = dParams0.PCC - dParams.PCC;
#ifdef _DEBUG
		s << "summit: " << summit.first << "\tPCCsummit: " << dParams.PCC << "\tdiff PCC: " << diffPCC << LF;
#endif
		if (diffPCC > 0.01)
			Err(Spec(eSpec::DEFECT) + SepSCl + sParams + sInaccurate).Warning();
		else if (diffPCC > 0.002)
			Err(Spec(eSpec::SDEFECT)).Warning();
	}
}

// Prints original distribution as a set of <frequency>-<size> pairs
void LenFreq::PrintSeq(dostream& s) const
{
	const chrlen maxLen = INT_MAX / 10;

	s << "\nOriginal " << sDistrib << COLON << "\nlength\tfrequency\n";
	for (const value_type& f : *this) {
		if (f.first > maxLen)	break;
		s << f.first << TAB << f.second << LF;
	}
}

// Constructor by pre-prepared frequency distribution file
//	@fname: name of pre-prepared frequency distribution file
LenFreq::LenFreq(const char* fName)
{
	TabFile file(fName, FT::eType::DIST);

	for(int x; file.GetNextLine();)
		if(x = file.IntField(0))	// IntField(0) returns 0 if zero field is not an integer
			(*this)[x] = file.IntField(1);
}

//#define _TIME
#ifdef _TIME
#include <chrono> 
using namespace std::chrono;
#endif

// Calculate and print dist params
//	@s: print stream
//	@ctype: combined type of distribution
//	@prDistr: if true then print distribution additionally
void LenFreq::Print(dostream& s, eCType ctype, bool prDistr)
{
	if(empty())		s << "empty " << sDistrib << LF;
	else {
		fraglen base = GetBase();	// initialized returned value
		//if (IsDegenerate())
		if (!base)
			s << "Degenerate " << sDistrib << " (only " << size() << " points)\n";
		else {
#ifdef _TIME
			auto start = high_resolution_clock::now();
			const int	tmCycleCnt = 1000;
#endif			
			// For optimization purposes, we can initialize base, keypts & summit at the first call of CallParams,
			// and use them on subsequent calls to avoid repeated PCC iterations.
			// However, the same base (and, as a consequence, keypts & summit) only works well for LNORM and GAMMA.
			// For the best NORM, base may be less, therefore, for simplicity and reliability, all parameters are always recalculated
#ifdef _DEBUG
			_s = &s;
			if(_fillSpline)	_spline.reserve(size() / 2);
#endif
#ifdef _TIME
			for(int i=0; i < tmCycleCnt; i++)
#endif
			point summit;				// returned value
			for (dtype i = 0; i < eCType::CNT; i++)
				if (IsType(ctype, i))
					CallParams(i, base, summit);
#ifdef _TIME
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			s << duration.count() / tmCycleCnt << " mcs\n";
#else
			// check for NORM if LNORM is defined
			if (IsType(ctype, eCType::LNORM) && !IsType(ctype, eCType::NORM)) {
				CallParams(GetDType(eCType::NORM), base, summit);
				_allParams.ClearNormDistBelowThreshold(1.02F);	// threshold 2%
			}
#endif
			PrintTraits(s, base, summit);
			_allParams.Print(s);
		}
		if (prDistr)	PrintSeq(s);
	}
	fflush(stdout);		// when called from a package
}

/************************ LenFreq: end ************************/
#endif	// _ISCHIP & _CALLDIST

#if defined _ISCHIP || defined _BSDEC
/************************ class AccumCover ************************/

// Adds fragment to accumulate the coverage
void AccumCover::AddRegion(const Region& frag)
{
	covmap::iterator it1 = find(frag.Start), it2;	// 'start', 'end' entries iterator

	// *** set up 'start' entry
	if (it1 == end()) {							// 'start' entry doesn't exist
		it2 = it1 = emplace(frag.Start, 1).first;
		if (it1 != begin())						// 'start' point is not the first one at all
			it1->second += (--it2)->second;		// correct val by prev point; keep it1 unchanged
	}
	else {
		it1->second++;							// incr val at existed 'start' entry
		if (--(it2 = it1) != end()				// decr it2; previous entry exists
			&& it2->second == it1->second)		// previous and current entries have the same value
			erase(it1), it1 = it2;				// remove current entry as duplicated
	}

	// *** set up 'end' entry
	it2 = find(frag.End);
	fraglen val = 0;							// 'end' entry value
	const bool newEnd = it2 == end();			// true if 'end' entry is new 

	if (newEnd) {								// 'end' entry is new
		it2 = emplace(frag.End, 0).first;
		val = next(it2) == end() ? 1 :			// 'end' entry is the last one at all
			prev(it2)->second;					// grab val by prev entry
	}

	// *** correct range between 'start' and 'end', set 'end' entry value
	for (it1++; /*it1 != end() &&*/ it1 != it2; it1++)		// correct values within range
		val = ++it1->second;					// increase value
	if ((--it1)->second == it2->second)			// is the last added entry a duplicate?
		erase(it2);								// remove duplicated entry
	else if (newEnd)
		it2->second = --val;					// set new 'end' entry value
}

/************************ class AccumCover: end ************************/
#endif	// _ISCHIP || _BSDEC