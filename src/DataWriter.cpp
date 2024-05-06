/**********************************************************
DataWriter.cpp
Last modified: 05/06/2024
***********************************************************/

#include "DataWriter.h"
#include <fstream>

const char* GM::title[] = { "test","control" };	// title: printed member's name

const char* SeqMode::modeTitles[] = { "single", "paired" };	// printed modes's name

class DataWriter;

/************************ ReadName ************************/

ReadName::tfAddRInfo ReadName::fAddNumber;
BYTE ReadName::Len = 0;		// Initialized in ReadWriter::Init();
							// if initialize by declaration, seg fault by compiling with gcc 4.1.2
BYTE ReadName::AddLen = Chrom::MaxMarkLength;	// + minimum one delimiter (DOT)


void ReadName::AddNumb(const Region&)
{
	const BYTE len = _constLen + AddLen + 1;
	_len = len + BYTE(sprintf(_name + len, "%ld", _rCnt));
}

void ReadName::AddPosSE(const Region& frag)
{
	const BYTE len = _constLen + AddLen + 1;
	_len = len + BYTE(sprintf(_name + len, "%u%c%ld", frag.Start, DOT, _rCnt));
}

void ReadName::AddPosPE(const Region& frag)
{
	const BYTE len = _constLen + AddLen + 1;
	_len = len + BYTE(sprintf(_name + len, "%u%c%u%c%ld",
		frag.Start, Read::NmPos2Delimiter, frag.End, DOT, _rCnt));
}

void ReadName::Init()
{
	Len = BYTE(Product::Title.length()) + 1 +	// + delimiter
		Chrom::MaxAbbrNameLength + 1 +			// chrom's name + delimiter
		20 + ReadWriter::MateLen;				// number + Mate suffix

	if (Read::IsPosInName())
		if (SeqMode::IsPE())	fAddNumber = &ReadName::AddPosPE, Len += 2 * CHRLEN_CAPAC + 1;
		else					fAddNumber = &ReadName::AddPosSE, Len += CHRLEN_CAPAC;
	else						fAddNumber = &ReadName::AddNumb;
}

ReadName::ReadName(LONG& rCnt) : _rCnt(rCnt), _constLen(BYTE(Product::Title.length()))
{
	assert(Len != 0);
	_name = new char[Len];
	memcpy(_name, Product::Title.c_str(), _constLen);
	_name[_constLen++] = Read::NmDelimiter;
	const BYTE cAbbrLen = BYTE(strlen(Chrom::Abbr));
	memcpy(_name + _constLen, Chrom::Abbr, cAbbrLen);
	_constLen += cAbbrLen;
	memset(_name + _constLen, '=', Chrom::MaxMarkLength);
}

void ReadName::SetChrom(const string& cMark)
{
	memcpy(_name + _constLen, cMark.c_str(), cMark.length());
	*(_name + _constLen + AddLen) = Read::IsPosInName() ? Read::NmPos1Delimiter : DOT;
}

/************************ ReadName: end ************************/

/************************ class ReadQualPattern ************************/

ReadWriter::ReadQualPattern::ReadQualPattern(const char* rqPattFName)
{
	const readlen rlen = DistrParams::IsRVL() ? Read::VarMinLen : Read::FixedLen;
	_rqPatt.reset(new char[rlen]);			// don't check allocation: small size
	Read::FillBySeqQual(_rqPatt.get(), rlen);	// fill buffer by default quality

	if (rqPattFName) {			// file exists
		TabReader file(rqPattFName);

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

void ReadWriter::ReadQualPattern::Fill(char* templ, readlen rlen) const
{
	Read::FillBySeqQual(templ, rlen);
	memcpy(templ, _rqPatt.get(), _defLen);
}

const void ReadWriter::ReadQualPattern::Print() const
{
	cout << "pattern: ";
	if (_defLen)	cout << string(_rqPatt.get()).substr(0, _defLen);
	if (DistrParams::IsRVL() || _defLen < Read::FixedLen)
		Read::PrintSeqQuality();
}


/************************ end of class ReadQualPattern ************************/

/************************ class ReadWriter ************************/

bool	ReadWriter::MultiThread;	// true if program is executed in a multi thread
const BYTE	ReadWriter::MateLen = BYTE(strlen(ReadWriter::Mate[0]));	// The length of Mate suffix
const char* ReadWriter::Mate[] = { "/1", "/2" };		// Array of Mate suffixes
const string ReadWriter::ReadLenTitle = " length=";
string ReadWriter::sReadConstLen;						// constant string "length=XX"

unique_ptr<ReadWriter::ReadQualPattern> ReadWriter::RqPattern;	// Read quality pattern

ReadWriter::tfAddReadMate ReadWriter::fAddReadMate[] = {
	&ReadWriter::LineNoAddReadMate,
	&ReadWriter::LineAddReadMate
};

void ReadWriter::LineAddCharsBack(const char* src, size_t len)
{
	LineAddCharBack(_delim);		// add delimiter
	_lineBuffOffset -= reclen(len);
	memcpy(_lineBuff + _lineBuffOffset, src, len);
}

void ReadWriter::LineAddReadName(BYTE mate)
{
	LineAddChars(_rName.Name(), _rName.Length(), !mate);
	(this->*fAddReadMate[bool(mate)])(mate);
}

void ReadWriter::LineAddReadVarName(readlen rLen)
{
	LineAddChar(AT);
	LineAddReadName(false);
	LineAddStr(ReadLenTitle, false);
	LineAddInt(rLen);
}

void ReadWriter::LineAddReadConstNameBack()
{
	LineAddCharsBack(_rName.Name(), _rName.Length() + sReadConstLen.size());
	memcpy(_lineBuff + _lineBuffOffset + _rName.Length(), sReadConstLen.c_str(), sReadConstLen.size());
	LineAddCharBack(AT);
}

void ReadWriter::LineFillReadVarPatt(readlen rlen)
{
	RqPattern->Fill(_lineBuff + CurrBuffPos(), rlen);
	LineIncrOffset(rlen);
}

/************************ class ReadWriter: end ************************/

SeqMode::eMode	SeqMode::mode;		// sequence mode
ULLONG	SeqMode::maxFragCnt;	// up limit of saved fragments

void SeqMode::Print(const char* signOut)
{
	cout << signOut << "Sequencing: " << modeTitles[IsPE()] << "-end"
		<< SepSCl << FT::ItemTitle(FT::ABED) << " limit = " << ReadsLimit() << LF;
}


/************************ class RBedWriter ************************/

RBedWriter::RBedWriter(const string& fName, ReadName& rName, const string* commLine)
	: ReadWriter(FT::BED, fName, rName)
{
	if (commLine)	CommLineToIOBuff(*DataWriter::CommLine());
	if (MultiThread)	Write();
	SetLineBuff(reclen(
		Chrom::MaxAbbrNameLength +		// length of chrom name
		ReadName::MaxLength() + 		// length of Read name
		2 * CHRLEN_CAPAC +				// start + stop positions
		DataWriter::MapQual.length() +		// score
		1 + 2 + 6));					// strand + HASH + SPACE + 5 TABs + LF
}

void RBedWriter::SetChrom(chrid cid)
{
	LineSetOffset();
	_offset = LineAddStr(Chrom::AbbrName(cid));
}

void RBedWriter::AddRead(const Read& read, bool reverse, BYTE mate)
{
	LineAddInts(read.Start, read.End);		// start, end
	LineAddReadName(mate);					// Read name
	LineAddStr(DataWriter::MapQual);		// score
	LineAddChar(Read::StrandMark(reverse));	// strand
	LineToIOBuff(_offset);
}

void RBedWriter::AddReads(const Read& read1, const Read& read2)
{
	AddRead(read1, false, 1);
	AddRead(read2, true, 2);
}


/************************ class RBedWriter: end ************************/

/************************ class FqWriter ************************/

reclen FqWriter::ReadStartPos = 0;	// Read field constant start position 

FqWriter::fAddRead FqWriter::addRead;

void FqWriter::AddFLRead(const Read& read, bool reverse)
{
	LineSetOffset(ReadStartPos);
	read.Copy(LineCurrPosBuf(), reverse);
	LineAddReadConstNameBack();
	LineBackToBuffer();
}

void FqWriter::AddVLRead(const Read& read, bool reverse)
{
	const readlen rlen = read.Length();

	LineAddReadVarName(rlen);
	read.Copy(LineCurrPosBuf(), reverse);
	LineIncrOffset(rlen);
	LineAddChars("\n+", 2);
	LineFillReadVarPatt(rlen);
	LineToIOBuff();
}

FqWriter::FqWriter(const string& fName, ReadName& rName)
	: ReadWriter(FT::FQ, fName, rName, LF)
{
	const readlen ReadMaxStrLen = readlen(to_string(Read::VarMaxLen).length());

	if (DistrParams::IsRVL())
		SetLineBuff(
			ReadName::MaxLength() +
			reclen(ReadWriter::ReadLenTitle.length()) +	// size of " length="
			ReadMaxStrLen +									// size of string representation of max Read len
			2 * Read::VarMaxLen +							// Read sequence + quality pattern
			5);												// AT + LF + 3 delimiters
	else {
		if (!ReadStartPos)							// if not initialized yet
			ReadStartPos = ReadName::MaxLength() + 2;		// Read name + AT + LF
		//== set line write buffer and writing start position
		const reclen pos = ReadStartPos + Read::FixedLen + 3;	// + 3 delimiters
		SetLineBuff(pos + Read::FixedLen);
		LineSetOffset(pos - 2);
		LineAddChar(PLUS, true);
		LineFillReadConstPatt(pos);
	}
}

/************************ class FqWriter: end ************************/

/************************ class SamWriter ************************/

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

reclen SamWriter::ReadStartPos = 0;
string SamWriter::Fld_5_6;		// combined value from 5 to 6 field: initialised in constructor
string SamWriter::FLAG[2];		// FLAG value for SE/PE: nitialised in constructor

SamWriter::tfAddRead SamWriter::fAddRead;

void SamWriter::AddFLRead(const Read& read, const string& fld_7_9, const string& flag)
{
	LineSetOffset(ReadStartPos);
	read.Copy(LineCurrPosBuf());				// 10: SEQ: Read
	LineAddStrBack(fld_7_9);					// 7-9: RNEXT + PNEXT + TLEN
	LineAddStrBack(Fld_5_6);					// 5-6: MAPQ + CIGAR
	LineAddStrBack(to_string(read.Start + 1));	// 4: POS
	LineAddStrBack(_cName);						// 3: RNAME
	LineAddStrBack(flag);						// 2: FLAG
	LineAddReadNameBack();						// 1: QNAME: Read name

	LineBackToBuffer();
}

void SamWriter::AddVLRead(const Read& read, const string& fld_7_9, const string& flag)
{
	const readlen rlen = read.Length();

	LineAddReadName();						// 1: QNAME: Read name
	LineAddStr(flag);						// 2: FLAG
	LineAddStr(_cName);						// 3: RNAME
	LineAddStr(to_string(read.Start + 1));	// 4: POS
	LineAddStr(DataWriter::MapQual);		// 5: MAPQ
	LineAddStr(to_string(rlen), false);		// 6: CIGAR
	LineAddChar(CIGAR_M, true);
	LineAddStr(fld_7_9);					// 7-9: RNEXT + PNEXT + TLEN
	read.Copy(LineCurrPosBuf());			// 10: SEQ: Read
	LineIncrOffset(rlen);
	LineAddChar(TAB);
	LineFillReadVarPatt(rlen);

	LineToIOBuff();
}

SamWriter::SamWriter(const string& fName, ReadName& rName, const ChromSizes& cSizes)
	: ReadWriter(FT::SAM, fName, rName)
{
	if (SeqMode::IsPE())	FLAG[0] = "99", FLAG[1] = "147";
	else					FLAG[0] = "0", FLAG[1] = "16";

	//== write header
	StrToIOBuff("@HD\tVN:1.0\tSO:unsorted");
	ostringstream oss;
	for (const auto& cs : cSizes) {
		oss << "@SQ\tSN:" << Chrom::Abbr << Chrom::Mark(cs.first) << "\tLN:" << cs.second.Data.Real;
		StrToIOBuff(oss.str());
		oss.str("");
	}
	oss << "@PG\tID:" << Product::Title << "\tPN:" << Product::Title
		<< "\tVN:" << Product::Version << "\tCL:\"" << *DataWriter::CommLine() << '\"';
	StrToIOBuff(oss.str());
	if (MultiThread)	Write();

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
		Fld_5_6 = DataWriter::MapQual + TAB			// MAPQ
			+ to_string(Read::FixedLen) + CIGAR_M;	// CIGAR

		//=== set SE pattern
		if (!SeqMode::IsPE()) {
			LineSetOffset(reclen(ReadStartPos - Fld_7_9.length() - Fld_5_6.length() - 2));
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

void SamWriter::AddRead(const Read& read, bool reverse)
{
	(this->*fAddRead)(read, Fld_7_9, FLAG[reverse]);
}

void SamWriter::AddReads(const Read& read1, const Read& read2, int fLen)
{
	(this->*fAddRead)(read1, GetPeFld_7_9(read2.Start, fLen), FLAG[0]);
	(this->*fAddRead)(read2, GetPeFld_7_9(read1.Start, -fLen), FLAG[1]);
}

/************************ class SamWriter: end ************************/

/************************ class BioWriters ************************/

DataWriter::BioWriters::tfAddRead	DataWriter::BioWriters::fAddRead = &DataWriter::BioWriters::AddReadSE;
float DataWriter::BioWriters::StrandErrProb;	// the probability of strand error

DataWriter::BioWriters::BioWriters(const string& fName, LONG* prCnt, const ChromSizesExt& cSizes)
	: _rName(*prCnt)
{
	_prCnt = prCnt;
	if (HasFormat(eFormat::BED))	_bedFile.reset(new RBedWriter(fName, _rName, CommLine()));
	if (HasFormat(eFormat::SAM))	_samFile.reset(new SamWriter(fName, _rName, cSizes));
	if (HasFormat(eFormat::FG))
		if (SeqMode::IsPE())
			_fqFile1.reset(new FqWriter(fName + "_1", _rName)),
			_fqFile2.reset(new FqWriter(fName + "_2", _rName));
		else
			_fqFile1.reset(new FqWriter(fName, _rName));
	if (HasFormat(eFormat::BED))	_bedFile.reset(new RBedWriter(fName, _rName, CommLine()));
	if (HasFormat(eFormat::SAM))	_samFile.reset(new SamWriter(fName, _rName, cSizes));
	if (HasFormat(eFormat::BGR))
		_bgFiles.reset(new OrderedCover(cSizes, isStrand ? 3 : 1, true, fName, "actual coverage", CommLine()));
	if (HasFormat(eFormat::FDENS)) {
		string name = fName + ".fdens";
		_fragWgFile.reset(new OrderedFreq(cSizes, 1, true, name, "fragment density", CommLine()));
	}
	if (HasFormat(eFormat::RDENS)) {
		string name = fName + ".rdens";
		_readWgFile.reset(new OrderedFreq(cSizes, 1, true, name, "read density", CommLine()));
	}
}

DataWriter::BioWriters::BioWriters(const BioWriters& primer)
	: _rName(*primer._prCnt)
{
	_prCnt = primer._prCnt;
	if (primer._fqFile1)	_fqFile1.reset(new FqWriter(*primer._fqFile1, _rName));
	if (primer._fqFile2)	_fqFile2.reset(new FqWriter(*primer._fqFile2, _rName));
	if (primer._bedFile)	_bedFile.reset(new RBedWriter(*primer._bedFile, _rName));
	if (primer._samFile)	_samFile.reset(new SamWriter (*primer._samFile, _rName));
	if (primer._bgFiles)	_bgFiles.reset(new OrderedCover(*primer._bgFiles));
	if (primer._fragWgFile)	_fragWgFile.reset(new OrderedFreq(*primer._fragWgFile));
	if (primer._readWgFile)	_readWgFile.reset(new OrderedFreq(*primer._readWgFile));
}

void DataWriter::BioWriters::BeginWriteChrom(const ChromSeq& seq)
{
	_seq = &seq;
	_rName.SetChrom(Chrom::Mark(seq.ID()));
	if (_bedFile)		_bedFile->SetChrom(seq.ID());
	if (_samFile)		_samFile->SetChrom(seq.ID());
	if (_bgFiles)		_bgFiles->SetChrom(seq.ID());
	if (_fragWgFile)	_fragWgFile->SetChrom(seq.ID());
	if (_readWgFile)	_readWgFile->SetChrom(seq.ID());
}

void DataWriter::BioWriters::EndWriteChrom() const
{
	if (_bgFiles)		_bgFiles->WriteChrom(_seq->ID());
	if (_fragWgFile)	_fragWgFile->WriteChrom(_seq->ID());
	if (_readWgFile)	_readWgFile->WriteChrom(_seq->ID());
}

int DataWriter::BioWriters::AddReadSE(const Region& frag, readlen rLen, bool reverse)
{
	const chrlen rPos = reverse ? frag.End - rLen : frag.Start;	// Read's position
	const Read read(_seq->Seq(rPos), rPos, rLen);
	int ret = read.CheckNLimit();
	if (ret)		return ret;

	//if(RandomReverse && g==Gr::FG && _rng.Sample(BioWriters::StrandErrProb) ) {
	//	reverse = !reverse;
	//	//short diff = ftr->Centre() - rPos - (Read::FixedLen>>1);
	//	//short halfDiffPeak = (200 - ftr->Length())>>1;
	//	//if(reverse)		diff += halfDiffPeak;
	//	//else			diff -= halfDiffPeak;

	//	//short diff = short( reverse ? (rPos + Read::FixedLen - ftr->Start) : ftr->End - rPos);
	//	//if(!_rng.Sample(float(diff) / fLen)) {
	//		//if(diff/10 < _fDist.Length()) _fDist[diff/10]++;
	//		//reverse = !reverse;			// imitate strand error
	//	//}
	//}

	if (_bgFiles)		_bgFiles->AddFrag(frag, reverse);		// coverage
	if (_fragWgFile)	_fragWgFile->AddFragDens(frag);			// frag frequency
	if (_readWgFile)	_readWgFile->AddReadDens(read, reverse);// read frequency
	if (InclReadName()) {
		InterlockedIncrement(_prCnt);
		_rName.AddNumber(frag);
		if (_fqFile1)	_fqFile1->AddRead(read, reverse);
		if (_bedFile)	_bedFile->AddRead(read, reverse);
		if (_samFile)	_samFile->AddRead(read, reverse);
	}
	return 0;
}

int DataWriter::BioWriters::AddReadPE(const Region& frag, readlen rLen, bool reverse)
{
	const Read read1(_seq->Seq(frag.Start), frag.Start, rLen);
	int ret = read1.CheckNLimit();
	if (ret)	return ret;
	chrlen pos2 = frag.End - rLen;
	const Read read2(_seq->Seq(pos2), pos2, rLen);
	ret = read2.CheckNLimit();
	if (ret)	return ret;

	if (_bgFiles)		_bgFiles->AddFrag(frag, reverse);			// coverage
	if (_fragWgFile)	_fragWgFile->AddFragDens(frag);				// frag frequency
	if (_readWgFile)	_readWgFile->AddReadDens(read1, reverse),	// read frequency
						_readWgFile->AddReadDens(read2, reverse);
	if (InclReadName()) {
		InterlockedIncrement(_prCnt);
		_rName.AddNumber(frag);
		if (_fqFile1)	_fqFile1->AddRead(read1, false),
						_fqFile2->AddRead(read2, true);
		if (_bedFile)	_bedFile->AddReads(read1, read2);
		if (_samFile)	_samFile->AddReads(read1, read2, frag.Length());
	}

	return 0;
}

void DataWriter::BioWriters::PrintFormat(const char* signOut, const char* predicate) const
{
	if (HasFormat(eFormat::FG)) {
		cout << signOut << predicate << "sequence: " << _fqFile1->FileName();
		if (SeqMode::IsPE())		cout << SepCm << _fqFile2->FileName();
		cout << LF;
	}
	if (HasFormat(eFormat::BED, eFormat::SAM)) {
		cout << signOut << predicate << "alignment: ";
		if (HasFormat(eFormat::BED)) {
			cout << _bedFile->FileName();
			if (HasFormat(eFormat::SAM))	cout << SepCm;
		}
		if (HasFormat(eFormat::SAM))	cout << _samFile->FileName();
		cout << LF;
	}
	if (HasFormat(eFormat::BGR)) {
		cout << signOut << predicate << "coverage: ";
		_bgFiles->PrintWritersName();
		cout << LF;
	}
	if (HasFormat(eFormat::FDENS, eFormat::RDENS)) {
		cout << signOut << predicate << "density" << SepDCl;
		const char* sep = strEmpty.c_str();
		if (_fragWgFile) {
			cout << entityTitles[0] << SepCl;
			_fragWgFile->PrintWritersName();
			sep = SepSCl;
		}
		if (_readWgFile) {
			cout << sep << entityTitles[1] << SepCl;
			_readWgFile->PrintWritersName();
		}
		cout << LF;
	}
}

/************************ class BioWriters: end ************************/

/************************ class DistrWriters ************************/

const string DataWriter::DistrWriters::fExt[] = { ".fdist", ".rdist" };
const char* DataWriter::DistrWriters::entityAdjust[] = { " size",	" length" };

DataWriter::DistrWriters::DistrWriters(const string& fName, bool isFragDist, bool isReadDist)
	: _fName(fName)
{
	if (isFragDist) {
		_dist[0] = new Distrib();
		_fAddFrag = [&](fraglen flen) {_dist[0]->AddVal(flen); };
	}
	if (isReadDist) {
		_dist[1] = new Distrib();
		_fAddRead = [&](readlen rlen) {_dist[1]->AddVal(rlen); };
	}
}

// Writes distributions to files and delete them
DataWriter::DistrWriters::~DistrWriters()
{
	const char* sSet = "Set ";

	for (BYTE i = 0; i < ND; ++i)
		if (_dist[i]) {
			Distrib::eCType dtype = Distrib::eCType::NORM;
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
					if (!DistrParams::IsSS())	dtype = Distrib::eCType::LNORM;
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
void DataWriter::DistrWriters::AddFrag(fraglen flen, readlen rlen)
{
	_fAddFrag(flen);		// fragments
	_fAddRead(rlen);		// reads
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
//	@predicate: output common title
void DataWriter::DistrWriters::PrintFormat(const char* signOut, const char* predicate) const
{
	if (HasFormat(eFormat::FDIST, eFormat::RDIST)) {
		cout << signOut << predicate << Distrib::sDistrib << SepDCl;
		for (BYTE i = 0; i < ND; i++) {
			if (HasContigFormat(eFormat::FDIST, i))
				cout << entityTitles[i] << entityAdjust[i] << SepCl << FileName(i);
			if (i == 0 && HasBothFormats(eFormat::FDIST, eFormat::RDIST))
				cout << SepSCl;
		}
		cout << LF;
	}
}

/************************ class DistrWriters: end ************************/

/************************ class DataWriter ************************/

//bool	DataWriter::RandomReverse = true;	// true if Read should be reversed randomly
string	DataWriter::MapQual;				// the mapping quality
int		DataWriter::Format;					// output formats as int
bool	DataWriter::inclReadName;			// true if Read name is included into output data
bool	DataWriter::isStrand;
const char* DataWriter::entityTitles[] = { "fragment", Read::title };
const string* DataWriter::commLine;			// command line

void DataWriter::Init(int fFormat, BYTE mapQual, bool bgStrand, float strandErrProb, bool zipped)
{
	Format = int(eFormat(fFormat));
	inclReadName = HasFormat(eFormat::FG, eFormat::BED, eFormat::SAM);
	MapQual = to_string(mapQual);
	isStrand = !SeqMode::IsPE() && bgStrand;
	TxtWriter::Zipped = zipped;
	BioWriters::Init(strandErrProb);
	ReadWriter::Init();
	FqWriter::Init();
	SamWriter::Init();
}

// Prints item title ("reads/fragments") according to output formats
void DataWriter::PrintItemTitle()
{
	PrintItemsSummary(string(entityTitles[1]) + 's', string(entityTitles[0]) + 's');
	cout << COLON;
}

void DataWriter::PrintItemCount(size_t fCnt)
{
	if (SeqMode::IsPE())	PrintItemsSummary(2 * fCnt, fCnt);
	else					cout << fCnt;
}

DataWriter::DataWriter(
	const string& fName, bool control, const string& cmLine, const ChromSizesExt& cSizes)
	: _dists(new DistrWriters(
		fName,
		HasFormat(eFormat::FDIST), HasFormat(eFormat::RDIST))),
		_gMode(BYTE(GM::eMode::Test)
		)
{
	commLine = &cmLine;
	_writers[0].reset(new BioWriters(fName, &_rCnt[0], cSizes));
	if (control)
		_writers[1].reset(new BioWriters(fName + "_input", &_rCnt[1], cSizes));
}

DataWriter::DataWriter(const DataWriter& file) : _dists(file._dists), _gMode(file._gMode)
{
	_writers[0].reset(new BioWriters(*file._writers[0]));
	if (file._writers[1])
		_writers[1].reset(new BioWriters(*file._writers[1]));
}

void DataWriter::BeginWriteChrom(const ChromSeq& seq)
{
	_writers[0]->BeginWriteChrom(seq);
	if (_writers[1])	_writers[1]->BeginWriteChrom(seq);
}

void DataWriter::EndWriteChrom()
{
	_writers[0]->EndWriteChrom();
	if (_writers[1])
		_writers[1]->EndWriteChrom();
}

int DataWriter::AddRead(chrlen pos, fraglen flen, bool reverse)
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
	return _writers[_gMode]->AddRead(Region(pos, pos + flen), rlen, reverse);
}

void DataWriter::PrintFormat(const char* signOut) const
{
	const char* output = "Output ";

	_writers[0]->PrintFormat(signOut, output);
	_dists->PrintFormat(signOut, output);
	if (_writers[1])	cout << signOut << output << "control supplied\n";
}

void DataWriter::PrintReadQual(const char* signOut) const
{
	if (!InclReadName())	return;
	bool prQualPatt = HasFormat(eFormat::SAM, eFormat::FG);
	cout << signOut << Read::Title << " quality" << SepDCl;
	if (prQualPatt)
		ReadWriter::PrintReadQualPatt();
	if (HasFormat(eFormat::SAM, eFormat::BED)) {
		if (prQualPatt)	cout << SepSCl;
		cout << "mapping = " << MapQual;
	}
	cout << LF;
}

/************************ class DataWriter: end ************************/