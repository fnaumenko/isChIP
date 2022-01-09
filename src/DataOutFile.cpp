/**********************************************************
DataOutFile.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides output data text files functionality
***********************************************************/

#include "DataOutFile.h"

const char* GM::title[] = {"test","control"};	// title: printed member's name

const char* Seq::modeTitles[] = { "single", "paired" };	// printed modes's name

/************************ ReadName ************************/

ReadName::tfAddRInfo ReadName::fAddInfo;
BYTE	ReadName::len = 0;		// Initialized in DataOutFile::Init();
								// if initialize by declaration, seg fault by compiling with gcc 4.1.2
bool	ReadName::MultiThread;	// true if program is executed in a multi thread

// Adds to Read name its position as string
//void ReadName::AddPos(const string& s)
//{
//	_len = BYTE(_headChrLen + s.length());
//	move(s.begin(), s.end(), _name + _headChrLen);
//}

// Adds to Read name its number followed be delimitor
//	@shift: starting position in buffer
void ReadName::AddNumb(const Region&)
{
	//char* buf = _name + _headChrLen;
	//*buf = Read::NmNumbDelimiter;
	//const string s = to_string(CountIncr());
	//move(s.begin(), s.end(), ++buf);
	//_len = BYTE(_headChrLen + s.size() + 1);

	_len = _headChrLen + BYTE(sprintf(_name + _headChrLen, "%c%llu", Read::NmNumbDelimiter, CountIncr()));
}

void ReadName::AddPosSE(const Region& frag)
{
	//ostringstream ss;
	//ss << frag.Start << Read::NmNumbDelimiter << CountIncr();
	//AddPos(ss.str());

	//#ifdef OS_Windows
//	AddPos(frag.Start);
//	AddNumb(_len);
//	_len = _headChrLen + PrintNumbToBuff(_name + _headChrLen, frag.Start);
//	_len = _len + PrintDelimNumbToBuff(_name + _len, Read::NmNumbDelimiter, CountIncr());
//#else
	_len = _headChrLen + BYTE(sprintf(_name + _headChrLen, "%u%c%llu", frag.Start, Read::NmNumbDelimiter, CountIncr()));
//#endif
}

void ReadName::AddPosPE	(const Region& frag)
{
	//ostringstream ss;
	//ss << frag.Start << Read::NmPos2Delimiter << frag.End << Read::NmNumbDelimiter << CountIncr();
	//AddPos(ss.str());

//#ifdef OS_Windows
//	AddPos(frag.Start);
//	_len += PrintDelimNumbToBuff(_name + _len, Read::NmPos2Delimiter, frag.End);
//	AddNumb(_len);
//#else
	_len = _headChrLen + BYTE(sprintf(_name + _headChrLen, "%u%c%u%c%llu", 
		frag.Start, Read::NmPos2Delimiter, frag.End, Read::NmNumbDelimiter, CountIncr()));
//#endif
}

void ReadName::Init()
{
	len =	BYTE(Product::Title.length()) + 1 +	// + delimiter
			Chrom::MaxAbbrNameLength + 1 +		// chrom's name + delimiter
			20 + DataOutFile::MateLen;			// number + Mate suffix

	if(Read::IsPosInName())
		if(Seq::IsPE())		fAddInfo = &ReadName::AddPosPE,	len += 2*CHRLEN_CAPAC + 1;
		else				fAddInfo = &ReadName::AddPosSE,	len += CHRLEN_CAPAC;
	else					fAddInfo = &ReadName::AddNumb;
}

// Initializes instance by constant part of Read name
//	@rCnt: external Read counter
ReadName::ReadName(ULLONG& rCnt) : _rCnt(rCnt), _headLen(BYTE(Product::Title.length()))
{
	_name = new char[len];
	memcpy(_name, Product::Title.c_str(), _headLen);
	_name[_headLen++] = Read::NmDelimiter;
	const BYTE cAbbrLen = BYTE(strlen(Chrom::Abbr));
	memcpy(_name + _headLen, Chrom::Abbr, cAbbrLen);
	_headLen += cAbbrLen;
}

// Copy constructor
//ReadName::ReadName(const ReadName& rName) : _rCnt(rName._rCnt), _headLen(rName._headLen)
//{
//	_name = new char[len];
//	memcpy(_name, rName._name, _headLen);
//	//cout << "ReadName copy constructor\n";
//}

// Sets current chrom's mark
void ReadName::SetChrom(const string&& cMark)
{
	_headChrLen = BYTE(cMark.length());
	move(cMark.begin(), cMark.end(), _name + _headLen);
	//memcpy(_name + _headLen, cMark.c_str(), _headChrLen = BYTE(cMark.length()));
	_headChrLen += _headLen;
	if (Read::IsPosInName())
		*(_name + _headChrLen++) = Read::NmPos1Delimiter;
}

/************************ ReadName: end ************************/

/************************ class ReadQualPattern ************************/

// Creates Read quality pattern buffer and fills it by first valid line from file.
//	@rqPattFName: name of valid file with a quality line
DataOutFile::ReadQualPattern::ReadQualPattern(const char* rqPattFName)
{
	const readlen rlen = DistrParams::IsRVL() ? Read::VarMinLen : Read::FixedLen;
	_rqPatt.reset(new char[rlen]);			// don't check allocation: small size
	Read::FillBySeqQual(_rqPatt.get(), rlen);	// fill buffer by default quality

	if (rqPattFName) {			// file exists
		TabFile file(rqPattFName);

		// variant with the last line in file
		//readlen lineLen = 0;
		//const char * ln, * line = NULL;
		//while (ln = file.GetLine())
		//	if (ln) { line = ln; lineLen = (readlen)file.LineLength(); }

		const char* line = file.GetNextLine();	// first line in file
		if (line) {
			readlen lineLen = (readlen)file.LineLength();
			memcpy(_rqPatt.get(), line, _defLen = lineLen > rlen ? rlen : lineLen);
		}
		else
			Err("has no meaningful lines; ignored", rqPattFName).Warning();
	}
}

// Fills external FQ|SAM template by variable Read quality pattern
//	@templ: pointer to external FQ|SAM Read quality template
void DataOutFile::ReadQualPattern::Fill(char* templ, readlen rlen) const
{
	Read::FillBySeqQual(templ, rlen);
	memcpy(templ, _rqPatt.get(), _defLen);
}

// Returns Read quality pattern buffer to print
const void DataOutFile::ReadQualPattern::Print() const
{
	cout << "pattern: ";
	if (_defLen)	cout << string(_rqPatt.get()).substr(0, _defLen);
	if (DistrParams::IsRVL() || _defLen < Read::FixedLen)
		Read::PrintSeqQuality();
}


/************************ end of class ReadQualPattern ************************/

/************************ class DataOutFile ************************/

const BYTE	DataOutFile::MateLen = BYTE(strlen(DataOutFile::Mate[0]));	// The length of Mate suffix
const char* DataOutFile::Mate[] = { "/1", "/2" };		// Array of Mate suffixes
const string* DataOutFile::CommLine;					// command line
const string DataOutFile::ReadLenTitle = " length=";
string DataOutFile::sReadConstLen;						// constant string "length=XX"

unique_ptr<DataOutFile::ReadQualPattern> DataOutFile::RqPattern;	// Read quality pattern

DataOutFile::tfAddReadName DataOutFile::fAddReadNames[] = {	// 'Add qualified Read name' methods
	&DataOutFile::AddReadNameEmpty,
	&DataOutFile::LineAddReadMate
};

// Copies block of chars before the current position in the line write buffer.
//	@src:  pointer to the block of chars
//	@len: number of chars
void DataOutFile::LineAddCharsBack(const char* src, size_t len)
{
	LineAddCharBack(_delim);		// add delimiter
	_lineBuffOffset -= rowlen(len);
	memcpy(_lineBuff + _lineBuffOffset, src, len);
}

// Copies default Read name and pos extention after current position in the line write buffer,
//	and increases current position.
//	@mate: mate number for PE Reads, or 0 for SE Read
void DataOutFile::LineAddReadName(BYTE mate)
{
	LineAddChars(_rName.Name(), _rName.Length(), !mate);
	(this->*fAddReadNames[bool(mate)])(mate);
}

// Copies qualified Read name started with '@' and read variable legth after current position
//	in the line write buffer, and increases current position
//	@len: Read length
//	Used by FqOutFile.
void DataOutFile::LineAddReadVarName(readlen len)
{
	LineAddChar(AT);
	LineAddReadName(false);
	LineAddStr(ReadLenTitle, false);
	LineAddInt(len);
}

// Copies qualified Read name started with '@', positions and Read constant length before current position
//	in the line write buffer, adds delimiter after Read Name and decreases current position.
//	Invoked in FqOutFile.
void DataOutFile::LineAddReadConstNameBack()
{
	LineAddCharsBack(_rName.Name(), _rName.Length() + sReadConstLen.size());
	memcpy(_lineBuff + _lineBuffOffset + _rName.Length(), sReadConstLen.c_str(), sReadConstLen.size());
	LineAddCharBack(AT);
}

// Fills line by Read variable quality pattern from the current position and increases current position
//	@rlen: Read's length
void DataOutFile::LineFillReadVarPatt(readlen rlen)
{
	RqPattern->Fill(_lineBuff + CurrBuffPos(), rlen);
	LineIncrOffset(rlen);
}

/************************ class DataOutFile: end ************************/

Seq::eMode	Seq::mode;		// sequence mode
ULLONG	Seq::maxFragCnt;	// up limit of saved fragments

// Prints sequencing modes
//	@signOut: output marker
void Seq::Print(const char* signOut)
{
	cout << signOut << "Sequencing: " << modeTitles[IsPE()] << "-end"
		<< SepSCl << FT::ItemTitle(FT::eType::ABED) << " limit = " << ReadsLimit() << LF;
}


/************************ class BedROutFile ************************/

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@rName: Read's name
BedROutFile::BedROutFile(const string& fName, const ReadName& rName)
	: DataOutFile(FT::eType::BED, fName, rName)
{
	CommLineToIOBuff(*DataOutFile::CommLine);
	if(ReadName::MultiThread)	Write();
	SetLineBuff(rowlen(
		Chrom::MaxAbbrNameLength +		// length of chrom name
		ReadName::MaxLength() + 		// length of Read name
		2 * CHRLEN_CAPAC +				// start + stop positions
		Output::MapQual.length() +		// score
		1 + 2 + 6));					// strand + HASH + SPACE + 5 TABs + LF
}

// Sets treated chrom's name to line write buffer
void BedROutFile::SetChrom(chrid cID)
{
	LineSetOffset();
	_offset = LineAddStr(Chrom::AbbrName(cID));
}

// Adds Read to the line's write buffer.
//	@pos: valid Read's start position
//	@len: valid Read's length
//	@reverse: if true then add complemented read
//	@mate: mate number for PE Reads, or 0 for SE Read
void BedROutFile::AddRead(const Read& read, bool reverse, BYTE mate)
{
	LineAddInts(read.Start(), read.End());		// start, end
	LineAddReadName(mate);						// Read name
	LineAddStr(Output::MapQual);				// score
	LineAddChar(Read::StrandMark(reverse));		// strand
	LineToIOBuff(_offset);
}

/************************ class BedROutFile: end ************************/

/************************ class FqOutFile ************************/

rowlen FqOutFile::ReadStartPos = 0;	// Read field constant start position 

FqOutFile::fAddRead FqOutFile::addRead;

// Adds Read with fixed length
//	@read: valid Read
//	@reverse: if true then add complemented read 
void FqOutFile::AddFLRead(const Read& read, bool reverse)
{
	LineSetOffset(ReadStartPos);
	read.Copy(LineCurrPosBuf(), reverse);
	LineAddReadConstNameBack();
	LineBackToBuffer();
}

// Adds Read with variable length
//	@read: valid Read
//	@reverse: if true then add complemented read 
void FqOutFile::AddVLRead(const Read& read, bool reverse)
{
	const readlen rlen = read.Length();

	LineAddReadVarName(rlen);
	read.Copy(LineCurrPosBuf(), reverse);
	LineIncrOffset(rlen);
	LineAddChars("\n+", 2);
	LineFillReadVarPatt(rlen);
	LineToIOBuff();
}

// Creates new instance for writing
//	@fName: file name without extention
//	@rName: Read's name
FqOutFile::FqOutFile(const string& fName, const ReadName& rName)
	: DataOutFile(FT::eType::FQ, fName, rName, LF) 
{
	const readlen ReadMaxStrLen = readlen(to_string(Read::VarMaxLen).length());

	if (DistrParams::IsRVL())
		SetLineBuff(
			ReadName::MaxLength() +
			rowlen(DataOutFile::ReadLenTitle.length()) +	// size of " length="
			ReadMaxStrLen +									// size of string representation of max Read len
			2 * Read::VarMaxLen +							// Read sequence + quality pattern
			5);												// AT + LF + 3 delimiters
	else {
		if (!ReadStartPos)							// if not initialized yet
			ReadStartPos = ReadName::MaxLength() + 2;		// Read name + AT + LF
		//== set line write buffer and writing start position
		const rowlen pos = ReadStartPos + Read::FixedLen + 3;	// + 3 delimiters
		SetLineBuff(pos + Read::FixedLen);
		LineSetOffset(pos - 2);
		LineAddChar(PLUS, true);
		LineFillReadConstPatt(pos);
	}
}

/************************ class FqOutFile: end ************************/

/************************ class SamOutFile ************************/

/*
SAM format:
https://samtools.github.io/hts-specs/SAMv1.pdf
Col	Field	Type	Regexp/Range Brief		description
----------------------------------------------------------------------------
1	QNAME	String	[!-?A-~]f1,255g			Query template NAME (read name)
2	FLAG	Int		[0,216-1]				bitwise FLAG
3	RNAME	String	\*|[!-()+-<>-~][!-~]*	Reference sequence NAME (chrom)
4	POS		Int		[0,231-1]				1-based leftmost mapping POSition
5	MAPQ	Int		[0,28-1]				MAPping Quality
6	CIGAR	String	\*|([0-9]+[MIDNSHPX=])+	CIGAR string (read length)
7	RNEXT	String	\*|=|[!-()+-<>-~][!-~]*	Ref. name of the mate/next read (SE:*; PE:=)
8	PNEXT	Int		[0,231-1]				Position of the mate/next read (SE:0)
9	TLEN	Int		[-231+1,231-1]			observed Template LENgth (SE:0; PE:frag length)
10	SEQ		String	\*|[A-Za-z=.]+			segment SEQuence
11	QUAL	String	[!-~]+					ASCII of Phred-scaled base QUALity+33
----------------------------------------------------------------------------

bitwise FLAG
Bit		Description
--------------------------------------------------------
1 0x1	template having multiple segments in sequencing
2 0x2	each segment properly aligned according to the aligner
3 0x4	segment unmapped
4 0x8	next segment in the template unmapped
5 0x10	SEQ being reverse complemented
6 0x20	SEQ of the next segment in the template being reversed
7 0x40	the rst segment in the template
8 0x80	the last segment in the template
9 0x100	secondary alignment
0 0x200	not passing quality controls
1 0x400	PCR or optical duplicate
2 0x800	supplementary alignment
--------------------------------------------------------

 0	   00000	SE +
16	   10000	SE -
99	01100011	PE +
147	10010011	PE -

*/

rowlen SamOutFile::ReadStartPos = 0;
string SamOutFile::Fld_5_6;		// combined value from 5 to 6 field: initialised in constructor
string SamOutFile::FLAG[2];		// FLAG value for SE/PE: nitialised in constructor

SamOutFile::tfAddRead SamOutFile::fAddRead;

// Adds Read with fixed length
//	@read: valid Read
//	@fld_7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
//	@flag: FLAG field value
void SamOutFile::AddFLRead(const Read& read, const string& fld_7_9, const string& flag)
{
	LineSetOffset(ReadStartPos);
	read.Copy(LineCurrPosBuf());					// 10: SEQ: Read
	LineAddStrBack(fld_7_9);						// 7-9: RNEXT + PNEXT + TLEN
	LineAddStrBack(Fld_5_6);						// 5-6: MAPQ + CIGAR
	LineAddStrBack(to_string(read.Start() + 1));	// 4: POS
	LineAddStrBack(_cName);							// 3: RNAME
	LineAddStrBack(flag);							// 2: FLAG
	LineAddReadNameBack();							// 1: QNAME: Read name

	LineBackToBuffer();
}

// Adds Read with variable length
//	@read: valid Read
//	@fld_7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
//	@flag: FLAG field value
void SamOutFile::AddVLRead(const Read& read, const string& fld_7_9, const string& flag)
{
	const readlen rlen = read.Length();

	LineAddReadName();							// 1: QNAME: Read name
	LineAddStr(flag);							// 2: FLAG
	LineAddStr(_cName);							// 3: RNAME
	LineAddStr(to_string(read.Start() + 1));	// 4: POS
	LineAddStr(Output::MapQual);				// 5: MAPQ
	LineAddStr(to_string(rlen), false);			// 6: CIGAR
	LineAddChar(CIGAR_M, true);
	LineAddStr(fld_7_9);						// 7-9: RNEXT + PNEXT + TLEN
	read.Copy(LineCurrPosBuf());				// 10: SEQ: Read
	LineIncrOffset(rlen);
	LineAddChar(TAB);
	LineFillReadVarPatt(rlen);

	LineToIOBuff();
}

// Creates new instance for writing, initializes line write buffer writes header.
//	@fName: file name without extention
//	@rName: Read's name
//	@cSizes: chrom sizes
SamOutFile::SamOutFile(const string& fName, const ReadName& rName, const ChromSizes& cSizes)
	: DataOutFile(FT::eType::SAM, fName, rName)
{
	if (Seq::IsPE())	FLAG[0] = "99", FLAG[1] = "147";
	else				FLAG[0] = "0",	FLAG[1] = "16";

	//== write header
	StrToIOBuff("@HD\tVN:1.0\tSO:unsorted");
	ostringstream oss;
	for (const auto& cs : cSizes) {
		oss << "@SQ\tSN:" << Chrom::Abbr << Chrom::Mark(cs.first) << "\tLN:" << cs.second.Data.Real;
		StrToIOBuff(oss.str());
		oss.str("");
	}
	oss << "@PG\tID:" << Product::Title << "\tPN:" << Product::Title
		<< "\tVN:" << Product::Version << "\tCL:\"" << *DataOutFile::CommLine << '\"';
	StrToIOBuff(oss.str());
	if(ReadName::MultiThread)	Write();

	//== maximal length of write line buffer without Read & Quality fields, with delimiters
	ReadStartPos =
		ReadName::MaxLength() +		//  1 QNAME: Read name
		3 +							//  2 FLAG: bitwise FLAG
		Chrom::MaxAbbrNameLength +	//  3 RNAME: AbbrChromName
		CHRLEN_CAPAC +				//  4 POS: 1-based start pos
		2 +							//  5 MAPQ: MAPping Quality
		4 +							//  6 CIGAR: Read length in 3 digits + letter 'M' or '='
		1 +							//  7 RNEXT: Ref. name of the mate/next read: SE '*', PE '='
		CHRLEN_CAPAC +				//  8 PNEXT: Position of the mate/next read
		3 + 						//  9 TLEN: observed Template LENgth
		10 + 1 + 1;					// 10 number of TABs + LF  + one for safety

	if (DistrParams::IsRVL())
		SetLineBuff(ReadStartPos + 2 * Read::VarMaxLen);
	else {
		//== set line buffer
		ReadStartPos += Read::FixedLen;
		SetLineBuff(ReadStartPos + Read::FixedLen);
		LineFillReadConstPatt(ReadStartPos);
		LineSetChar(--ReadStartPos, TAB);
		ReadStartPos -= Read::FixedLen;		// restore true value

		//== set Fields 5-6
		Fld_5_6 = Output::MapQual + TAB			// MAPQ
			+ to_string(Read::FixedLen) + CIGAR_M;	// CIGAR
		
		//=== set SE pattern
		if (!Seq::IsPE()) {
			LineSetOffset(rowlen(ReadStartPos - Fld_7_9.length() - Fld_5_6.length() - 2));
			LineAddStr(Fld_5_6);
			LineAddStr(Fld_7_9);
		}
	}
}

// Returns string containing initialized PE feilds 7-9
string GetPeFld_7_9(chrlen pos, int fLen)
{
	ostringstream ss;
	ss << "=\t" << ++pos << TAB << fLen;
	return ss.str();
}

// Adds two mate Reads to the line's write buffer.
//	@read1: valid first mate Read
//	@read2: valid second mate Read
//	@pos1: valid first mate Read's start position
//	@pos2: valid second mate Read's start position
//	@fLen: fragment's length
void SamOutFile::AddTwoReads(const Read& read1, const Read& read2, int fLen)
{
	(this->*fAddRead)(read1, GetPeFld_7_9(read2.Start(), fLen), FLAG[0]);
	(this->*fAddRead)(read2, GetPeFld_7_9(read1.Start(), -fLen), FLAG[1]);
}

/************************ class SamOutFile: end ************************/

/************************ class WigOutFile ************************/

const string WigOutFile::WigFormats[] = { "bedGraph", "wiggle_0" };

// Initializes line write buffer, adds command and definition lines
//	@ftype: BGRAPH or WIG_VAR
//	@fName: file name without extention
//	@declDescr: brief file description in declaration line
//	@strand: C-string described strand, or NULL
void WigOutFile::Init(FT::eType ftype, const string& fName, const char* declDescr, const char* strand)
{
	SetLineBuff(Chrom::MaxAbbrNameLength + 3 * CHRLEN_CAPAC);

	// comm line and decl line go directly to IO buff
	CommLineToIOBuff(*DataOutFile::CommLine);
	ostringstream oss;
	oss << "track type=" << WigFormats[int(ftype) - int(FT::eType::BGRAPH)]
		<< " name=\"" << fName << FT::Ext(ftype, DataOutFile::Zipped)
		<< "\" description=\"" << Product::Title << SPACE << declDescr;
	if (strand)		oss << SepCl << strand << " strand";
	oss << "\" color=50,130,190 autoScale=on";
	StrToIOBuff(oss.str());
}

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@ftype: BGRAP or WIG_VAR
//	@descr: brief file description in declaration line
//	@mtype: mutex locker type
//	@cSizes: chrom sizes
WigOutFile::WigOutFile(
	FT::eType ftype,
	const string& fName,
	const char* descr,
	Mutex::eType mtype,
	const ChromSizesExt& cSizes
)	: _mType(mtype), TxtOutFile(ftype, fName, TAB)
{
	Init(ftype, fName, descr);
	for (const auto& cs : cSizes)
		if (cs.second.Treated)
			AddElem(cs.first, DensCover());
}

// Close data container for given chrom for increment.
//	@cID: chrom ID
void WigOutFile::CloseChromData(chrid cID)
{
	// AccumCovers are created and filled in different threads independently.
	// To save them in chrom sorted order the total pool is examined 
	// each time the next coverage is completed.
	// If previous coveraged are filled without gups, they are recorded and removed from the pool.

	Mutex::Lock(_mType);		// different mutexes for different WIG files
	bool save = true;
	At(cID).Data.Closed = true;
	for (cIter it = Begin(); CID(it) < cID; it++)
		if (it->second.Data.Unsaved())
			if (save = it->second.Data.Closed)	// chrom cover is closed, so save it
				WriteChromData(CID(it));
			else 	break;						// some of chrom covers aren't saved
	//for (auto& c : Container())
	//	if (c.second.Data.Unsaved())
	//		if (save = c.second.Data.Closed)	// chrom cover is closed, so save it
	//			WriteChromData(c.first);
	//		else 	break;						// some of chrom covers aren't saved
	if (save)
		WriteChromData(cID);
	Mutex::Unlock(_mType);
}

/************************ class WigOutFile: end ************************/

/************************ class BedGrOutFile ************************/

const char* BedGrOutFiles::BedGrOutFile::DeclDescr = "actual coverage";	//	brief file description in declaration line
const char* BedGrOutFiles::BedGrOutFile::StrandTitles[] = { "positive", "negative" };

// Fill IO buffer by chrom data
void BedGrOutFiles::BedGrOutFile::WriteChromData(chrid cID)
{
	LineSetOffset();
	const rowlen offset = LineAddStr(Chrom::AbbrName(cID));

	auto& d = At(cID).Data;				// not const because of d.Clear() at the end
	auto it0 = d.cbegin(), it = it0;
	const auto end = d.cend();

	for (++it; it != end; it0 = it++)
		if (it0->second)
			LineAddInts(it0->first, it->first, it0->second, false),		// start, end, coverage
			LineToIOBuff(offset);

	//d.DoWith2Items([&](const auto& it0, const auto& it) {
	//	LineAddInts(it0->first, it->first, it0->second, false);		// start, end, coverage
	//	LineToIOBuff(offset);
	//	});
	d.Clear();
}

/************************ class BedGrOutFile: end ************************/

/************************ class BedGrOutFiles ************************/

bool BedGrOutFiles::IsStrands;		// true if wigs with different strands should be generated

	// Applies function fn to each of the item in _files
void BedGrOutFiles::DoForFiles(function<void(BedGrOutFiles::BedGrOutFile*)> fn)
{
	for (BYTE i = !IsStrands * (Count - 1); i < Count; i++)	// last item (basic) or all ones (plus strands)
		fn(_files[i]);
}

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@cSizes: chrom sizes
BedGrOutFiles::BedGrOutFiles(const string& fName, const ChromSizesExt& cSizes)
{
	_files[Count-1] = new BedGrOutFile(fName, cSizes);
	if(IsStrands)
		_files[0] = new BedGrOutFile(0, fName + "_pos", *_files[Count-1]),
		_files[1] = new BedGrOutFile(1, fName + "_neg", *_files[Count-1]);
}

// Adds SE fragment to total coverage, and strand coverage if set
//	@frag: added fragment
//	@reverse: true if read is reversed (neg strand)
void BedGrOutFiles::AddFrag(const Region& frag, bool reverse)
{
	_files[Count - 1]->AddFragCov(frag);
	if(_files[0])
		_files[reverse]->AddFragCov(frag);
}

// Prints output file names separated by comma
void BedGrOutFiles::PrintNames() const
{
	cout << _files[Count-1]->FileName();
	if (_files[0])
		for(BYTE i=0; i<Count-1; i++)
			cout << SepCm << _files[i]->FileName();
}

/************************ class BedGrOutFiles: end ************************/

/************************ class Wig0OutFile ************************/

// Fill IO buffer by chrom data
void Wig0OutFile::WriteChromData(chrid cID)
{
	// write declaration line
	LineSetOffset();
	StrToIOBuff("variableStep chrom=chr" + Chrom::Mark(cID) + " span=1");

	// write data lines
	auto& d = At(cID).Data;
	for (const auto& f : d) {					// frequency
		LineAddInts(f.first, f.second, false);	// pos, freq
		LineToIOBuff();
	}
	d.Clear();
}

/************************ class Wig0OutFile: end ************************/

/************************ class OutFile ************************/

Output::OutFile::tfAddRead	Output::OutFile::fAddRead = &Output::OutFile::AddReadSE;
float Output::OutFile::StrandErrProb;	// the probability of strand error

// Creates and initializes new instance for writing.
//	@fName: common file name without extention
//	@cSizes: chrom sizes
Output::OutFile::OutFile(const string& fName, const ChromSizesExt& cSizes)
{
	if (HasFormat(eFormat::FG))
		if (Seq::IsPE())
			_fqFile1 = new FqOutFile(fName + "_1", _rName),
			_fqFile2 = new FqOutFile(fName + "_2", _rName);
		else
			_fqFile1 = new FqOutFile(fName, _rName);
	if (HasFormat(eFormat::BED))	_bedFile = new BedROutFile(fName, _rName);
	if (HasFormat(eFormat::SAM))	_samFile = new SamOutFile(fName, _rName, cSizes);
	if (HasFormat(eFormat::BGR))	_bgFile = new BedGrOutFiles(fName, cSizes);
	if (HasFormat(eFormat::FDENS))	_coverFile[0] = new Wig0OutFile(fName, false, cSizes);
	if (HasFormat(eFormat::RDENS))	_coverFile[1] = new Wig0OutFile(fName, true, cSizes);
}

// Clone constructor for multithreading
//	@file: original instance
//	@threadNumb: number of thread
Output::OutFile::OutFile(const OutFile& file) : _primer(false)
{
	_rName.SetReadCounter(file._rCnt);

	if (file._fqFile1)	_fqFile1 = new FqOutFile(*file._fqFile1, _rName);
	if (file._fqFile2)	_fqFile2 = new FqOutFile(*file._fqFile2, _rName);
	if (file._bedFile)	_bedFile = new BedROutFile(*file._bedFile, _rName);
	if (file._samFile)	_samFile = new SamOutFile(*file._samFile, _rName);
	_bgFile = file._bgFile;						// primer files are common for the all clones
	_coverFile[0] = file._coverFile[0];			// primer file is common for the all clones
	_coverFile[1] = file._coverFile[1];			// primer file is common for the all clones
}

Output::OutFile::~OutFile()
{
	delete _fqFile1;
	delete _fqFile2;
	delete _bedFile;
	delete _samFile;
	if (_primer) {
		delete _bgFile;
		for (const auto& f : _coverFile)
			delete f;
	}
}

// Start recording chrom
void Output::OutFile::BeginWriteChrom(const RefSeq& seq)
{
	_seq = &seq;
	_rName.SetChrom(Chrom::Mark(seq.ID()));
	if (_bedFile)	_bedFile->SetChrom(seq.ID());	// set chrom's name for writing.
	if (_samFile)	_samFile->SetChrom(seq.ID());	// set chrom's name for writing.
	if (_bgFile)	_bgFile->SetChrom(seq.ID());	// set chrom's coverage as current
	for(auto& cf : _coverFile)
		if (cf)		cf->OpenChromData(seq.ID());	// set chrom's density as current
}

// Stop recording chrom
void Output::OutFile::EndWriteChrom() const
{
	if (_bgFile)	_bgFile->CloseChrom(_seq->ID());
	for (const auto& f : _coverFile)
		if(f)	f->CloseChromData(_seq->ID());
}

// Adds one SE Read
//	@frag: added fragment
//	@rLen: Read's length
//	@reverse: if true then add complemented read
//	return:	1: fragment is out of range (end of chrom)
//			0: Read is added successfully
//			-1: N limit is exceeded
int Output::OutFile::AddReadSE(const Region& frag, readlen rLen, bool reverse)
{
	const chrlen rPos = reverse ? frag.End - rLen : frag.Start;	// Read's position
	const Read read(_seq->Seq(rPos), rPos, rLen);
	int ret = read.CheckNLimit();
	if(ret)		return ret;
	/*
	if(RandomReverse && g==Gr::FG && _rng.Sample(OutFile::StrandErrProb) ) {
		reverse = !reverse;
		//short diff = ftr->Centre() - rPos - (Read::FixedLen>>1);
		//short halfDiffPeak = (200 - ftr->Length())>>1;
		//if(reverse)		diff += halfDiffPeak;
		//else			diff -= halfDiffPeak;

		//short diff = short( reverse ? (rPos + Read::FixedLen - ftr->Start) : ftr->End - rPos);
		//if(!_rng.Sample(float(diff) / fLen)) {
			//if(diff/10 < _fDist.Length()) _fDist[diff/10]++;
			//reverse = !reverse;			// imitate strand error
		//}
	}
	*/

	if (_bgFile)		_bgFile->AddFrag(frag, reverse);		// coverage
	if (_coverFile[0])	_coverFile[0]->AddFrag(frag);			// frag density
	if (_coverFile[1])	_coverFile[1]->AddRead(read, reverse);	// read density
	if (InclReadName()) {
		_rName.AddInfo(frag);
		if (_fqFile1)	_fqFile1->AddRead(read, reverse);
		if (_bedFile)	_bedFile->AddRead(read, reverse);
		if (_samFile)	_samFile->AddRead(read, reverse);
	}
	return 0;
}

// Adds two PE Reads
//	@frag: added fragment
//	@rLen: Read's length
///	@g: FG or BG; needs for strand error imitation; not used
//	@reverse: not used
//	return:	1: fragment is out of range (end of chrom)
//			0: Reads are added successfully
//			-1: N limit is exceeded
int Output::OutFile::AddReadPE(const Region& frag, readlen rLen, bool reverse)
{
	const Read read1(_seq->Seq(frag.Start), frag.Start, rLen);
	int ret = read1.CheckNLimit();
	if (ret)	return ret;
	chrlen pos2 = frag.End - rLen;
	const Read read2(_seq->Seq(pos2), pos2, rLen);
	ret = read2.CheckNLimit();
	if (ret)	return ret;

	if (_bgFile)		_bgFile->AddFrag(frag, reverse);		// coverage
	if (_coverFile[0])	_coverFile[0]->AddFrag(frag);			// frag density
	if (_coverFile[1])	_coverFile[1]->AddRead(read1, reverse),	// read density
						_coverFile[1]->AddRead(read2, reverse);
	if (InclReadName()) {
		_rName.AddInfo(frag);
		if (_fqFile1)	_fqFile1->AddRead(read1, false),
						_fqFile2->AddRead(read2, true);
		if (_bedFile)	_bedFile->AddRead(read1, false, 1),
						_bedFile->AddRead(read2, true, 2);
		if (_samFile)	_samFile->AddTwoReads(read1, read2, frag.Length());
	}


	return 0;
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
//	@predicate: 'output' marker
void Output::OutFile::PrintFormat(const char* signOut, const char* predicate) const
{
	if(HasFormat(eFormat::FG)) {
		cout << signOut << predicate << "sequence: " << _fqFile1->FileName();
		if(Seq::IsPE())		cout << SepCm << _fqFile2->FileName();
		cout << LF;
	}
	if(HasFormat(eFormat::BED, eFormat::SAM)) {
		cout << signOut << predicate << "alignment: ";
		if(HasFormat(eFormat::BED)) {
			cout << _bedFile->FileName();
			if(HasFormat(eFormat::SAM))	cout << SepCm;
		}
		if(HasFormat(eFormat::SAM))	cout << _samFile->FileName();
		cout << LF;
	}
	if(HasFormat(eFormat::BGR)) {
		cout << signOut << predicate << "coverage: ";
		_bgFile->PrintNames();
		cout << LF;
	}
	if (HasFormat(eFormat::FDENS, eFormat::RDENS)) {
		cout << signOut << predicate << "density" << SepDCl;
		for (BYTE i = 0; i < ND; i++) {
			if (HasContigFormat(eFormat::FDENS, i))
				cout << entityTitles[i] << SepCl << _coverFile[i]->FileName();
			if (i == 0 && HasBothFormats(eFormat::FDENS, eFormat::RDENS))
				cout << SepSCl;
		}
		cout << LF;
	}
}

/************************ class OutFile: end ************************/

/************************ class DistrFiles ************************/

const string Output::DistrFiles::fExt[] = { ".fdist", ".rdist" };
const char* Output::DistrFiles::entityAdjust[] = { " size",	" length" };

Output::DistrFiles::DistrFiles(const string& fName, bool isFragDist, bool isReadDist)
	: _fName(fName)
{
	if (isFragDist) {
		_dist[0] = new LenFreq();
		_fAddFrag = [&](fraglen flen) {_dist[0]->AddLen(flen); };
	}
	if (isReadDist) {
		_dist[1] = new LenFreq();
		_fAddRead = [&](readlen rlen) {_dist[1]->AddLen(rlen); };
	}
}

// Writes distributions to files and delete them
Output::DistrFiles::~DistrFiles()
{
	const char* sSet = "Set ";

	for(BYTE i=0; i<ND; ++i)
		if (_dist[i]) {
			LenFreq::eCType dtype = LenFreq::eCType::NORM;
			ofstream s;
			s.open(FileName(i));
			if (s.is_open()) {
				if (i)		// Reads
					if (DistrParams::IsRVL())	// variable reads
						DistrParams::PrintReadDistr(s, sSet, Read::title);
					else
						s << sSet << "constant " << Read::title << " length" << LF;
				else {		// fragments
					DistrParams::PrintFragDistr(s, sSet, false);
					if (!DistrParams::IsSS())	dtype = LenFreq::eCType::LNORM;
				}
				_dist[i]->Print(s, dtype);
				s.close();
			}
			else
				Err(Err::F_OPEN, FileName(i).c_str()).Throw(false);		// no exception from destructor
			delete _dist[i];
		}
}

// Adds frag/read length to statistics
void Output::DistrFiles::AddFrag(fraglen flen, readlen rlen)
{
	_fAddFrag(flen);		// fragments
	_fAddRead(rlen);		// reads
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
//	@predicate: output common title
void Output::DistrFiles::PrintFormat(const char* signOut, const char* predicate) const
{
	if (HasFormat(eFormat::FDIST, eFormat::RDIST)) {
		cout << signOut << predicate << sDistrib << SepDCl;
		for (BYTE i = 0; i < ND; i++) {
			if (HasContigFormat(eFormat::FDIST, i))
				cout << entityTitles[i] << entityAdjust[i] << SepCl << FileName(i);
			if (i == 0 && HasBothFormats(eFormat::FDIST, eFormat::RDIST))
				cout << SepSCl;
		}
		cout << LF;
	}
}

/************************ class DistrFiles: end ************************/

/************************ class Output ************************/

//bool	Output::RandomReverse = true;	// true if Read should be reversed randomly
string	Output::MapQual;				// the mapping quality
int		Output::Format;					// output formats as int
bool	Output::inclReadName;			// true if Read name is included into output data
const char* Output::entityTitles[] = { "fragment", Read::title };

// Initializes static members
//	@fFormat: types of output files
//	@mapQual: the mapping quality
//	@bgStrand: true if bedGraphs with different strands should be generated
//	@strandErrProb: the probability of strand error
//	@zipped: true if output files should be zipped
void Output::Init(int fFormat, BYTE mapQual, bool bgStrand, float strandErrProb, bool zipped)
{
	Format = int(eFormat(fFormat));
	inclReadName = HasFormat(eFormat::FG, eFormat::BED, eFormat::SAM);
	MapQual = to_string(mapQual);
	BedGrOutFiles::IsStrands = !Seq::IsPE() && bgStrand;
	TxtOutFile::Zipped = zipped;
	OutFile::Init(strandErrProb);
	ReadName::Init();
	FqOutFile::Init();
	SamOutFile::Init();
}

// Prints item title ("reads/fragments") according to output formats
void Output::PrintItemTitle()
{
	PrintItemsSummary(string(entityTitles[1]) + 's', string(entityTitles[0]) + 's');
	cout << COLON;
}

// Prints item title ("reads/fragments") according to output formats
//	@fCnt: number of fragments
void Output::PrintItemCount(ULLONG fCnt)
{
	if (Seq::IsPE())	PrintItemsSummary(2 * fCnt, fCnt);
	else				cout << fCnt;
}

// Creates new instance for writing.
//	@fName: common file name without extention
//	@control: if true, then control ('input') is generated
//	@cmLine: command line to add as a comment in the first file line
//	@cSizes: chrom sizes, or NULL
Output::Output(
	const string& fName, bool control, const string& cmLine, const ChromSizesExt& cSizes)
	: _dists(new DistrFiles(fName, HasFormat(eFormat::FDIST), HasFormat(eFormat::RDIST))),
	_gMode(BYTE(GM::eMode::Test))
{
	DataOutFile::CommLine = &cmLine;
	_oFiles[0].reset(new OutFile(fName, cSizes));
	if (control)	_oFiles[1].reset(new OutFile(fName + "_input", cSizes));
}

// Clone constructor for multithreading.
//	@file: original instance
Output::Output(const Output& file) : _dists(file._dists), _gMode(file._gMode)
{
	_oFiles[0].reset(new OutFile(*file._oFiles[0]));
	if (file._oFiles[1])	_oFiles[1].reset(new OutFile(*file._oFiles[1]));
}
	
// Starts recording chrom
void Output::BeginWriteChrom(const RefSeq& seq)
{
	_oFiles[0]->BeginWriteChrom(seq);
	if(_oFiles[1])	_oFiles[1]->BeginWriteChrom(seq);
}

// Stops recording chrom
void Output::EndWriteChrom()
{
	_oFiles[0]->EndWriteChrom();
	if(_oFiles[1])	_oFiles[1]->EndWriteChrom();
}

// Adds read(s) to output file
//	@pos: current fragment's position
//	@flen: length of current fragment
///	@g: FG or BG; needs for strand error imitation
//	@reverse: if true then add complemented read
//	return:	1: fragment is out of range (end chrom)
//			0: Read(s) is(are) added, or nothing (trial)
//			-1: N limit is exceeded; Read(s) is(are) not added
int Output::AddRead(chrlen pos, fraglen flen, /*Gr::eType g,*/ bool reverse)
{
	/*****
	 1) Fill distribution files doesn't need by the trial pass,
	 Acceptable, because
		a) trial pass is performed just once, and it's very unlikely that the format DIST
		will be set on the first pass
		b) even so, it almost doesn't matter for the distribution statistics
	 2) Generation Read variable length generation is not needed if one format BG is set.
	 Acceptable, because it's very unlikely that the only format BG and RVL are set at the same time
	*****/
	readlen rlen = Read::FixedLen;
	if (DistrParams::IsRVL()) {
		rlen = readlen(_rng.Normal() * DistrParams::rdSigma + DistrParams::rdMean);
		if (rlen > Read::VarMaxLen || rlen > flen)	rlen = flen;
	}
	_dists->AddFrag(flen, rlen);
	return _oFiles[_gMode]->AddRead(Region(pos, pos + flen), rlen, reverse);
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
void Output::PrintFormat(const char* signOut) const
{
	const char* output = "Output ";

	_oFiles[0]->PrintFormat(signOut, output);
	_dists->PrintFormat(signOut, output);
	if(_oFiles[1])	cout << signOut << output << "control supplied\n";
}

// Prints Read quality settins
//	@signOut: output marker
void Output::PrintReadQual (const char* signOut) const 
{
	if(!InclReadName())	return;
	bool prQualPatt = HasFormat(eFormat::SAM, eFormat::FG);
	cout << signOut << Read::Title << " quality" << SepDCl;
	if(prQualPatt)
		DataOutFile::PrintReadQualPatt();
	if (HasFormat(eFormat::SAM, eFormat::BED)) { 
		if (prQualPatt)	cout << SepSCl;
		cout << "mapping = " << MapQual;
	}
	cout << LF;
}

/************************ class Output: end ************************/