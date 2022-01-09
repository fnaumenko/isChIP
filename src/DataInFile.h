/**********************************************************
DataInFile.h (c) 2021 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 28.12.2021
-------------------------
Provides read|write text file functionality
***********************************************************/
#pragma once

//#ifndef _DATAINFILE_H
//#define _DATAINFILE_H

#include "TxtFile.h"
#include <map>
#include <functional>
//#include <variant>

// 'eOInfo' defines types of outputted oinfo
enum class eOInfo {	
	NONE,	// nothing printed: it is never pointed in command line
	LAC,	// laconic:		print file name if needs
	NM,		// name:		print file name
	STD,	// standard:	print file name and items number
	STAT,	// statistics:	print file name, items number and statistics
};

// 'DataInFile' is a common program interface (PI) of bed/bam input files
class DataInFile
{
protected:
	int	_cID = Chrom::UnID;			// current readed chrom ID; int for BAM PI compatibility

	// Sets the next chromosome as the current one if they are different
	//	@cID: next chrom ID
	//	@return: true, if new chromosome is set as current one
	bool SetNextChrom(chrid cID);

public:
	// Returns estimated number of items
	virtual ULONG EstItemCount() const = 0;

	// Sets the next chromosome as the current one if they are different
	// @cID: returned next chrom ID
	//	@return: true, if new chromosome is set as current one
	virtual bool GetNextChrom(chrid& cID) = 0;

	// Retrieves next item's record
	virtual bool GetNextItem() = 0;

	// Returns current item start position
	virtual chrlen ItemStart() const = 0;
	
	// Returns current item end position
	virtual chrlen ItemEnd() const = 0;
	
	// Returns current item length
	virtual readlen ItemLength() const = 0;

	// Returns true if alignment part of paired-end read
	virtual bool IsPairedItem() const  = 0;
	
	// Returns current item value (score)
	virtual float ItemValue() const  = 0;

	// Returns current item name
	virtual const char* ItemName() const = 0;
	
	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	virtual const string LineNumbToStr(Err::eCode code = Err::EMPTY) const = 0;

	// Throws exception with message included current reading line number
	//	@msg: exception message
	virtual void ThrowExceptWithLineNumb(const string& msg) const = 0;

	// Gets conditional file name: name if it's printable, otherwise empty string.
	virtual const string CondFileName() const = 0;

	// Returns true if item contains the strand sign
	//	Is invoked in the Feature constructor only.
	//virtual bool IsItemHoldStrand() const = 0;

	// Returns current item strand: true - positive, false - negative
	virtual bool ItemStrand() const = 0;
};

// 'BedInFile' represents unified PI for reading bed file
class BedInFile : public DataInFile, public TabFile
{
	const BYTE	NameFieldInd = 3;	// inner index of name field
	const BYTE	StrandFieldInd = 5;	// inner index of strand field

	BYTE _scoreInd;		// 0-based index of 'score' filed (used for FBED and all WIGs)
	BYTE _chrMarkPos;	// chrom's mark position in line (BED, BedGraph) or definition line (wiggle_0)
	function<bool()> _getStrand;	// returns current item strand; different for BED and ABED

	// Reset WIG type, score index, chrom mark position offset and estimated number of lines
	void ResetWigType(FT::eType type, BYTE scoreInd, size_t cMarkPosOffset);

	// Inserts '0' after chrom in current line and returns point to the next decl parameter if exists
	//const char* SplitLineOnChrom();

	// Checks for Fixed or Variable step wiggle type and corrects it if found
	//	@line: possible declaration line
	//	return: true if Fixed or Variable step type is specified
	bool DefineWigType(const char* line);

public:
	// Creates new instance for reading and open file
	//	@fName: name of file
	//	@type: file type
	//	@scoreNumb: number of 'score' filed (0 by default for ABED and BAM)
	//	@msgFName: true if file name should be printed in exception's message
	//	@abortInval: true if invalid instance should be completed by throwing exception
	BedInFile(const char* fName, FT::eType type, BYTE scoreNumb, bool msgFName, bool abortInval);

	// Gets pointer to the chrom mark in current line without check up
	inline const char* ChromMark() const { return GetLine() + _chrMarkPos; }

	// Returns estimated number of items
	inline ULONG EstItemCount() const { return EstLineCount(); }

	// Gets file bioinfo type
	inline FT::eType Type() const { return TabFile::Type(); }

	// Sets the next chromosome as the current one if they are different
	//	@cID: returned next chrom ID
	//	@str: null-terminated string started with short chrom name
	//	@return: true, if new chromosome is set as current one
	//	Used in Calc.cpp
	bool GetNextChrom(chrid& cID, const char* str) {
		return SetNextChrom(cID = Chrom::ValidateID(str, strlen(Chrom::Chrom::Abbr)));
	}

	// Sets the next chromosome as the current one if they are different
	//	@cID: returned next chrom ID
	//	@return: true, if new chromosome is set as current one
	//	To implement DataInFile virtual GetNextChrom(chrid& cID)
	bool GetNextChrom(chrid& cID) { 
		return SetNextChrom(cID = Chrom::ValidateID(ChromMark()));
	}

	// Retrieves next item's record
	inline bool GetNextItem()	{ return TabFile::GetNextLine(); }

	// Returns current item start position
	inline chrlen ItemStart()	const { return LongField(1); }

	// Returns current item end position
	inline chrlen ItemEnd()		const { return LongField(2); }

	// Returns current item length
	inline readlen ItemLength()	const { return readlen(ItemEnd() - ItemStart()); }

	// Returns true if alignment part of paired-end read
	inline bool IsPairedItem()	const { return strchr(ItemName() + 1, '/'); }

	// Returns current item value (score)
	inline float ItemValue()	const {
		return FloatFieldValid(_scoreInd);
	}

	// Returns current item name
	inline const char* ItemName() const { return StrFieldValid(NameFieldInd); }

	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	inline const string LineNumbToStr(Err::eCode code) const {
		return TxtInFile::LineNumbToStr(code);
	}

	// Throws exception with message included current reading line number
	//	@msg: exception message
	inline void ThrowExceptWithLineNumb(const string& msg) const { return TxtInFile::ThrowExceptWithLineNumb(msg); }

	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string CondFileName() const { return TxtFile::CondFileName(); }

	// Returns current item strand: true - positive, false - negative
	inline bool ItemStrand() const { return _getStrand(); }
};

class ChromSizes;	// Data.h

#ifdef _BAM
#ifdef _BIOSTAT		// defined in bioStat makefile
#include "../bam/BamReader.h"	// path in bioStat package
#else
#include "bam/BamReader.h"
#endif	// _BIOSTAT
#include <algorithm>
using namespace BamTools;

// 'BamInFile' represents unified PI for reading bam file
class BamInFile : public DataInFile
{
	// BamTools: http://pezmaster31.github.io/bamtools/struct_bam_tools_1_1_bam_alignment.html

	bool			_prFName;
	BamReader		_reader;
	BamAlignment	_read;
	mutable string	_rName;
	ULONG _estItemCnt = vUNDEF;	// estimated number of items

	// Returns SAM header data
	inline const string GetHeaderText() const { return _reader.GetHeaderText(); }

public:
	// Creates new instance for reading and open file
	//	@fName: name of file
	//	@cSizes: chrom sizes to be initialized or NULL
	//	@prName: true if file name should be printed in exception's message
	BamInFile(const char* fName, ChromSizes* cSizes, bool prName);

	// Returns estimated number of items
	inline ULONG EstItemCount() const { return _estItemCnt; }

	// returns chroms count
	inline chrid ChromCount() const { return _reader.GetReferenceCount(); }

	// Sets the next chromosome as the current one if they are different
	// @cID: returned next chrom ID
	//	@return: true, if new chromosome is set as current one
	bool GetNextChrom(chrid& cID) { return SetNextChrom(cID = Chrom::ValidateID(_read.RefID)); }

	// Retrieves next item's record
	bool GetNextItem()	{ 
		return _reader.GetNextAlignmentCore(_read)
			&& _read.Position >= 0;	// additional check because of bag: GetNextAlignment() doesn't
									// return false after last read while reading the whole genome
	}

	// Returns current item start position
	inline chrlen ItemStart()	const { return _read.Position; }

	// Returns current item end position
	inline chrlen ItemEnd()		const { return _read.Position + _read.Length; }

	// Returns current item length
	inline readlen ItemLength()	const { return readlen(_read.Length); }

	// Returns true if alignment part of paired-end read
	inline bool IsPairedItem()	const { return _read.IsPaired(); }

	// Returns current item value (score)
	inline float ItemValue()	const { return _read.MapQuality; }

	// Returns current item name
	inline const char* ItemName() const { return (_rName = _read.Name).c_str(); }

	// Gets string containing file name and current line number.
	inline const string LineNumbToStr(Err::eCode) const { return strEmpty; }

	// Throws exception with message included current reading line number
	inline void ThrowExceptWithLineNumb(const string& msg) const { Err(msg).Throw(); }

	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string CondFileName() const { return _prFName ? _reader.GetFilename() : strEmpty; }

	// DataInFile method empty implementation.
	//inline bool IsItemHoldStrand() const { return true; }

	// Returns current item strand: true - positive, false - negative
	inline bool  ItemStrand() const	{ return !_read.IsReverseStrand(); }

};
#endif	// _BAM

class UniBedInFile
{
public:
	// Action types
	enum eAction {
		ACCEPT,	// accepted
		TRUNC,	// truncated
		JOIN,	// joined
		OMIT,	// omitted
		ABORT,	// interruption at the first issue
	};
	
	// Ñonsolidated issue information; public becauseod use in CallDist (class FragDist)
	struct Issue {
		size_t	Cnt = 0;				// total number of issue cases
		const char*	Title;				// issue description
		const char* Ext = nullptr;		// extension of treatment description
		eAction Action = eAction::OMIT;	// issue treatment

		Issue(const char* title) : Title(title) {}
	};

private:
	FT::eType _type;			// should be const, but can be edited (from BEDGRAPF to WIGGLE)
	const char	_MaxDuplCnt;	// max allowed number of duplicates; -1 if keep all
	const bool	_abortInv;		// true if invalid instance should be completed by throwing exception
	const eOInfo _oinfo;		// level of output stat info

	chrid	_cCnt = 0;			// number of readed chroms
	Region	_rgn0{ 0,0 };		// previous item's region
	Region	_rgn{ 0,0 };		// current item's region
	size_t	_cDuplCnt = 0;		// number of duplicates per chrom; the first 'originals' are not counted
	char	_duplCnt = 0;		// current number of duplicates
	bool	_strand = true;		// current item's strand
	bool	_strand0 = true;	// previous item's strand; first sorted read is always negative
	DataInFile* _file = nullptr;// data file; unique_ptr is useless because of different type in constructor/destructor
	//variant<BedInFile, BamInFile> file;
	const ChromSizes* _cSizes;
	//mutable char _calledPaired = -1;	// initialized by first call of _file->IsPairedItem()

	// Resets the current accounting of items
	void ResetChrom();

	// Validate item
	//	@cLen: current chrom length or 0 if _cSizes is undefined
	//	return: true if item is valid
	bool CheckItem(chrlen cLen);

	// Validate item by final class
	//	return: true if item is valid
	virtual bool ChildCheckItem() { return true; }

	// Returns overlapping action
	virtual eAction GetOverlAction() const { return eAction::ACCEPT; }

	// Prints items statistics
	//	@cnt: total count of items
	void PrintStats(ULONG cnt);

	// Returns chrom size
	//	Defined in cpp because of call in template function (otherwise ''ChromSize' is no defined')
	chrlen ChromSize(chrid cID) const;

protected:
	// Item essue types
	enum eIssue {
		DUPL,		// duplicates
		OVERL,		// overlapping
		STARTOUT,	// starting outside the chromosome
		ENDOUT		// ending outside the chromosome
	};

	vector<Issue> _issues = {	// all possible issue info collection
		"duplicates",
		"overlapping" ,
		"starting outside the chromosome",
		"ending outside the chromosome"
};
	map<readlen, ULONG>	_lenFreq;		// item length frequency

	// Returns true if adjacent items overlap
	inline bool IsOverlap() const { return _rgn.Start <= _rgn0.End;  }

public:
	static bool IsTimer;	// if true then manage timer by Timer::Enabled, otherwise no timer

	// Prints count of items
	//	@cnt: total count of items
	//	@title: item title
	static void PrintItemCount(ULONG cnt, const string& title);

	// Prints items statistics
	//	@cnt: total count of items
	//	@issCnt: count of item issues
	//	@issues: issue info collection
	//	@prStat: it TRUE then print issue statsistics
	static void PrintStats(ULONG cnt, ULONG issCnt, const vector<Issue>& issues, bool prStat);

	// Creates new instance for reading and open file
	//	@fName: file name
	//	@type: file type
	//	@cSizes: chrom sizes
	//	@scoreNumb: number of 'score' filed (0 by default for ABED and BAM)
	//	@dupLevel: number of additional duplicates allowed; -1 - keep all additional duplicates
	//	@oinfo: output stat info level
	//	@prName: true if file name should be printed unconditionally
	//	@abortInval: true if invalid instance should be completed by throwing exception
	UniBedInFile(const char* fName, const FT::eType type, ChromSizes* cSizes,
		BYTE scoreNumb, char dupLevel, eOInfo oinfo, bool prName, bool abortInval);

	// explicit destructor
	~UniBedInFile();

	// pass through records
	template<typename Functor>
	void Pass(Functor& f)
	{
		const bool setCustom = Chrom::CustomID() != Chrom::UnID;	// 	chrom is specified by user
		size_t	cItemCnt = 0;					// count of chrom entries
		ULONG	tItemCnt = 0;					// total count of entries
		chrid cID = Chrom::UnID, nextcID = cID;	// current, next chrom
		chrlen	cLen = 0;						// current chrom length
		bool skipChrom = false;
		bool userChromInProc = false;
		Timer timer(IsTimer);

		while (_file->GetNextItem()) {
			if (_file->GetNextChrom(nextcID)) {			// the next chrom
				if (setCustom) {
					if (userChromInProc)		break;
					if (skipChrom = nextcID != Chrom::CustomID()) continue;
					userChromInProc = true;
				}
				f(cID, cLen, cItemCnt, nextcID);					// close current chrom, open next one
				ResetChrom();
				cID = nextcID;
				cItemCnt = 0;				
				if (_cSizes)	cLen = ChromSize(cID);
			}
			else if (skipChrom)		continue;
			_rgn.Set(_file->ItemStart(), _file->ItemEnd());	// the single invoke file->ItemStart(), ItemEnd()
			if (CheckItem(cLen)) {
				cItemCnt += f(); 							// treat entry
				_rgn0 = _rgn;
			}
			tItemCnt++;
		}
		f(cID, cLen, cItemCnt, tItemCnt);				// close last chrom
		
		if (_oinfo >= eOInfo::STD)	PrintStats(tItemCnt);
		timer.Stop(1, true, _oinfo > eOInfo::NM);
	}

	inline DataInFile& BaseFile() const { return *_file; }

	// Returns estimated number of items
	ULONG EstItemCount() const;

	// Gets file bioinfo type
	inline FT::eType Type() const { return _type; }

	// Gets count of chromosomes read
	inline chrid ReadedChromCount() const { return _cCnt; }

	// Returns current item region
	inline const Region& ItemRegion() const { return _rgn; }

	// Returns current item start position
	inline chrlen ItemStart() const { return _rgn.Start; }

	// Returns current item end position
	inline chrlen ItemEnd()	const	{ return _rgn.End; }

	// Returns previous accepted item end position
	inline chrlen PrevItemEnd() const { return _rgn0.End; }

	// Returns current item length
	inline readlen ItemLength()	const { return _rgn.Length(); }

	// Returns current item strand: true - positive, false - negative
	inline bool ItemStrand() const { return _strand; }

	// Returns current item value (score)
	inline float ItemValue() const { return _file->ItemValue(); }

	// Returns current item name
	inline const char* ItemName() const { return _file->ItemName(); }

	// Returns true if item contains the strand sign
	//	Is invoked in the Feature constructor only.
	//inline bool IsItemHoldStrand() const { return _file->IsItemHoldStrand(); }

	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	inline const string LineNumbToStr(Err::eCode code) const { return _file->LineNumbToStr(code); }

	// Throws exception with message included current reading line number
	//	@msg: exception message
	inline void ThrowExceptWithLineNumb(const string& msg) const { return _file->ThrowExceptWithLineNumb(msg); }

	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string CondFileName() const { return _file->CondFileName(); }

	// Returns number of duplicated items in current chrom (without the first 'original')
	inline size_t DuplCount() const { return _cDuplCnt; }

	// Returns total number of duplicate items (without the first 'original')
	size_t DuplTotalCount() const { 
		return _issues[DUPL].Cnt + _cDuplCnt;	// _cDuplCnt separate because the last chrom is not counted in a loop
	}
};

#ifdef _READS

class RBedInFile : public UniBedInFile
{
	mutable readlen _rLen = 0;		// most frequent (common) read length

public:
	// Creates new instance for reading and open file
	//	@fName: file name
	//	@cSizes: chrom sizes
	//	@dupLevel: number of additional duplicates allowed; -1 - keep all additional duplicates
	//	@oinfo: verbose level
	//	@prName: true if file name should be printed unconditionally
	//	@abortInval: true if invalid instance should be completed by throwing exception
	inline RBedInFile(const char* fName, ChromSizes* cSizes,
		char dupLevel, eOInfo oinfo, bool prName, bool abortInval = true) :
			UniBedInFile(fName, FT::GetType(fName, true), cSizes, 0, dupLevel, oinfo, prName, abortInval)
	{}

	// Returns the most frequent Read length
	// last in _rfreq because typically reads can be shorter, not longer
	readlen ReadLength() const { return _rLen ? _rLen : (_rLen = prev(_lenFreq.cend())->first); }

	// Returns true if alignment part of paired-end read
	inline bool IsPairedRead() const { return BaseFile().IsPairedItem(); }
};

#endif	// _READS
#ifdef _FEATURES

class FBedInFile : public UniBedInFile
{
private:
	const eAction _overlAction;	// set overlapping action
	const bool _isJoin;			// true if overlapping features should be joined
	bool _isOverlap;			// true if features are overlapping
	function<bool()> _action;	// overlapping features treatment

	// Returns: true if feature is valid
	bool ChildCheckItem();

	// Returns proposed overlapping action
	inline eAction GetOverlAction() const { return _overlAction; }

public:
	// Creates new instance for reading and open file
	//	@fName: file name
	//	@cSizes: chrom sizes
	//	@scoreNmb: number of 'score' filed
	//	@action: action for overlapping features
	//	@prName: true if file name should be printed unconditionally
	//	@abortInval: true if invalid instance should be completed by throwing exception
	FBedInFile(const char* fName, ChromSizes* cSizes,
		BYTE scoreNmb, eAction action, eOInfo oinfo, bool prName, bool abortInval);

	// If true then join overlapping feature
	bool IsJoined() const { return _isJoin &&_isOverlap; }

	// Returns true if features length distribution is degenerate
	bool NarrowLenDistr() const;
};
#endif	// _FEATURES

#if !defined _WIGREG && !defined _FQSTATN

// 'Read' represents Read (with name and score in case of _VALIGN) as item
class Read
{
public:
	static const readlen VarMinLen = 20;	// minimum Read length in variable Read mode
	static const readlen VarMaxLen = 3000;	// maximum Read length in variable Read mode
	static readlen	FixedLen;				// fixed length of Read
private:
#if defined _ISCHIP || defined _VALIGN || defined _PE_READ
	static const char	Strands[2];			// strand markers: [0] - positive, [1] - negative
public:
	static const char	NmDelimiter = ':';		// delimiter between progTitle and chrom
	static const char	NmNumbDelimiter = DOT;	// delimiter between prog title and number
	static const char	NmPos1Delimiter = ':';	// delimiter before first recorded position
	static const char	NmPos2Delimiter = '-';	// delimiter between two recorded positions in pair
#endif

#ifdef _ISCHIP

private:
	static char		SeqQuality;			// the quality values for the sequence (ASCII)
	static readlen	LimitN;				// maximal permitted number of 'N' in Read or vUNDEF if all
	static bool		PosInName;			// true if Read name includes a position
	static const char Complements[];	// template for complementing Read

	typedef void (Read::* pCopyRead)(char*) const;
	static pCopyRead CopyRead[2];

	const char* _seq;
	Region	_rgn;

	// Copies complemented Read into dst
	void CopyComplement(char* dst) const;

public:
	static const char* Title;
	static const char* title;

	// Initializes static members
	//	@len: length of Read
	//	@posInName: true if Read position is included in Read name
	//	@seqQual: quality values for the sequence
	//	@limN: maximal permitted number of 'N'
	static void Init(readlen len, bool posInName, char seqQual, short limN);

	inline static char StrandMark(bool reverse) { return Strands[int(reverse)]; }

	inline static bool IsPosInName() { return PosInName; }

	// Fills external buffer by quality values for the sequence
	inline static void FillBySeqQual(char* dst, readlen rlen) { memset(dst, SeqQuality, rlen); }

	// Constructor by sequence, start position and length
	Read(const char* seq, chrlen pos, readlen len) : _seq(seq) { _rgn.Set(pos, pos + len); }

	// Gets Read's region
	//inline const Region& Rgn() const { return _rgn; }

	// Gets Read's length
	inline readlen Length() const { return _rgn.Length(); }

	// Gets Read's start position
	inline chrlen Start() const { return _rgn.Start; }

	// Gets Read's end position
	inline chrlen End() const { return _rgn.End; }

	// Gets Read's sequence
	inline const char* Seq() const { return _seq; }

	// Copies Read into dst
	inline void Copy(char* dst) const { memcpy(dst, _seq, Length()); }

	// Copies initial or complemented Read into dst
	inline void Copy(char* dst, bool reverse) const { (this->*CopyRead[reverse])(dst); }

	// Checks Read for number of 'N'
	//	return:	1: NULL Read; 0: success; -1: N limit is exceeded
	int CheckNLimit() const;

	// Prints quality values for the sequence
	inline static void PrintSeqQuality() { cout << '[' << SeqQuality << ']'; }

	// Prints Read values - parameters
	//	@signOut: output marker
	//	@isRVL: true if Read variable length is set
	static void PrintParams(const char* signOut, bool isRVL);

#else
public:
	chrlen	Pos;		// Read's actual start position
	readlen Len;		// Read's length
	bool	Strand;		// true if strand is positive

	void InitBase(const RBedInFile& file);

	chrlen Centre() const { return Pos + (Len >> 1); }
#ifdef _PE_READ
	ULONG	Numb;		// read number keeped in name

	// PE Read constructor
	Read(const RBedInFile& file);

	inline chrlen End() const { return Pos + Len; }

	// Returns frag length
	//	@r: second read in a pair
	inline fraglen FragLen(const Read& r) const { return Strand ? r.End() - Pos : End() - r.Pos; }

	void Print() const { dout << Pos << TAB << Numb << TAB << Strand << LF; }

#elif defined _VALIGN
	chrid	RecCID;		// recorded (true) chrom - owner
	chrlen	RecPos;		// recorded (true) Read start position
	float	Score;		// Read's score

	// Extended (with saved chrom & position in name) Read constructor
	Read(const RBedInFile& file);

	void Print() const;
#endif

	// Compares two Reads by position. For sorting a container.
	inline static bool CompareByStartPos(const Read& r1, const Read& r2) { return r1.Pos < r2.Pos; }

	//inline static bool CompareByNum(const Read& r1, const Read& r2) {	return r1.Num < r2.Num; }
#endif
};
#endif

//#endif	//_DATAINFILE_H