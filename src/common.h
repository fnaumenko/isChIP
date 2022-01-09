/**********************************************************
common.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides common functionality
***********************************************************/

#pragma once
#ifndef _COMMON_H
#define _COMMON_H

#include "def.h"
#include <stdint.h>		// uint32_t
#include <string>
#include <iostream>	
#include <sstream>		// for NSTR( x )
#include <iomanip>		// setprecision(), setw()
#include <vector>
#include <sys/stat.h>	// struct stat
#include <limits>       // std::numeric_limits
#ifdef _MULTITHREAD
#include <mutex>
#endif

#ifdef __unix__
	#include <unistd.h>
	#include <stdlib.h>
	#include <memory>
	#include <limits.h>
	//#include <math.h>
	#include <dirent.h>
	#include <stdio.h>
#ifdef _MULTITHREAD
	//#include <pthread.h>
	#define InterlockedExchangeAdd	__sync_fetch_and_add
	#define InterlockedIncrement(p)	__sync_add_and_fetch(p, 1)
#endif
	#include <string.h>		// strerror_r()
	#include <stdexcept>	// throw std exceptions

	typedef unsigned char	BYTE;
	typedef	unsigned short	USHORT;
	typedef	unsigned int	UINT;
	typedef unsigned long	ULONG;
	typedef unsigned long long ULLONG;
	typedef long long LLONG;
	typedef struct stat struct_stat64;
	//typedef void*	thrRetValType;

	#define _stricmp strcasecmp	// case-sensitive comparison
	#define _fseeki64 fseeko64
	#define _ftelli64 ftello64
	#define _fileno	fileno
	#define _gcvt	gcvt
	#define _stat64 stat
	#define fopen	fopen64
	#define TRUE	1
	#define FALSE	0
	#define isNaN(x)	(x!=x)	// will not work by using --ffast-math compiler flag
	#define	LOCALE_ENG	"en_US.utf8"
	//#define SLASH '/'	// standard Linux path separator
	//#define F_READ_MODE O_RDONLY
#elif defined _WIN32
	#define OS_Windows
	#include <windows.h>
//#ifdef _MULTITHREAD
	//#define pthread_mutex_t CRITICAL_SECTION
//#endif	
	typedef unsigned __int64 ULLONG;
	typedef __int64 LLONG;
	typedef struct __stat64 struct_stat64;

	#define isNaN isnan
	#define	LOCALE_ENG	"English"
#endif		// OS Windows

#ifndef _NO_ZLIB
	#ifdef OS_Windows
		#define ZLIB_WINAPI
	#endif
	#include "zlib.h"
	#if ZLIB_VERNUM >= 0x1240
		#define ZLIB_NEW
	#else	
		#define ZLIB_OLD
	#endif
#endif	// _NO_ZLIB

// specific types
typedef BYTE		thrid;		// type number of thread
typedef BYTE		chrid;		// type number of chromosome
typedef uint16_t	readlen;	// type length of Read
typedef uint32_t	chrlen;		// type length of chromosome
typedef chrlen		fraglen;	// type length of fragment
typedef uint64_t	genlen;		// type length of genome

#define	CHRLEN_UNDEF	-1	// undefined length of chromosome
#define	CHRLEN_MAX		-1	// max length of chromosome
//#define CHRLEN_CAPAC	10	// capacity of max chrom length;
							// may be count by DigitsCount() every time,
							// but is the same if chrlen defined as int or long,
							// so is defined as static value

#define cNULL	'\0'
#define cN		'N'
#define HPH		'-'
#define USCORE	'_'
#define SPACE	' '
#define sSPACE	" "
#define QUOT	'\''
#define DOT		'.'
#define COLON	':'
#define COMMA	','
#define PERS	'%'
#define AT		'@'
#define PLUS	'+'
#define HASH	'#'
#define TAB		'\t'
#define LF		'\n'	// 10
#define CR		'\r'	// 13
//#define sBACK	"\b"

using namespace std;

static const chrlen CHRLEN_CAPAC = numeric_limits<chrlen>::digits10 + 1;

static const std::string ZipFileExt = ".gz";
static const string strEmpty = "";
static const char* SepCl = ": ";		// colon separator
static const char* SepSCl = "; ";		// semicolon separator
#ifdef _CONST_MARKS
static const char* SepCm = ", ";		// comma separator
static const char* SepDCl = ":: ";		// double colon separator
//static const char* SepClTab = ":\t";	// colon + tab separator
static const char* Equel = " = ";
#endif
// common help string
static const char* sGen = "gen";		// isChIP, readDens
static const char* sOutput = "out";
static const char* sTime = "time";
static const char* sVers = "version";
static const char* sHelp = "help";

//static const char* Done = " done";
static const char* sUnitDens = "r/kbp";	// unit of density measure; isChIP, readDens
static const char* sTotal = "total";
//static const string sTotal = "total";

//#ifdef _DUP_OUTPUT
static const char* sSumm = "summ";		// to invoke app from bioStat
static const string sFileDuplBegin = "duplicate standard output to specified file\nor to ";
static const string sFileDuplEnd = " if file is not specified";
//#endif
#if defined _ISCHIP || defined _CALLDIST
static const char* SepGroup = ";  ";
static const char* sMean = "mean";
static const char* sSigma = "sigma";
static const char* sDistrib = "distribution";
#endif	// _ISCHIP || _CALLDIST
#if !defined _WIGREG && !defined _FQSTATN
static const char* sTemplate = "template";
#endif

/*** COMMON MACROS & FUNCTION ***/

#define vUNDEF	-1	// undefined value
#define NO_VAL	-1	// option value is prohibited
#define NO_DEF	-1	// do not print option default value

/*************************************************************
* time testing integer to char[]
* printing string <numb><char><numb>[<char><numb>] and returing number of charachters printed
* 
* func		string			len	time, ss.mcs
*								Windows	Linux
* ------------------------------------------------------------
itoa    123456789-123456789	19	04.79
to_str  123456789-123456789	19	02.46	03.99
sprintf 123456789-123456789	19	03.59	01.47	common the best
oss     123456789-123456789	19	23.36	04.98
one oss 123456789-123456789	19	18.76	01.67
**************************************************************/

// Prints number to buffer without checkup
//	@buf: buffer to print
//	@numb: number to print
//	return: total number of characters written
//	NOT USED in Linux (sprintf() is used instead due to performance) 
//template <typename T>
//BYTE	PrintNumbToBuff(char* const buf, T numb)
//{
//	const string s = to_string(numb);
//	move(s.begin(), s.end(), buf);
//	return BYTE(s.size());
//}

// Prints the number following the delimiter to buffer without checkup
//	@buf: buffer to print
//	@delim: delimiter to print
//	@numb: number to print
//	return: total number of characters written
//	NOT USED in Linux (sprintf() is used instead due to performance) 
//template <typename T>
//BYTE	PrintDelimNumbToBuff(char* const buf, char delim, T numb) {
//	*buf = delim;
//	return PrintNumbToBuff(buf + 1, numb) + 1;
//}

// Gets number of members in static array
#define ArrCnt(arr)	sizeof(arr)/sizeof(arr[0])

// Digital to STRing
// Returns value's string representation. http://rootdirectory.de/wiki/NSTR()
// Instead of std::to_string(x) [C++11] because of compatibility
//#define NSTR(x) std::to_string(x)

// Returns number of ones in an bynary integer
int OnesCount(int n);

// Returns right position of right one in an bynary integer
int RightOnePos(int n);

// Gets number of digist in a integral value
//	@val: integral value
//	@isLocale: if true then adds number of '1000' separators
//	return: number of digist without minus symbol or 0 if value is 0
int DigitsCount (LLONG val, bool isLocale = false);

// Returns percent of @part relatively @total
inline float Percent(ULLONG part, ULLONG total) { 
	return total ? 100.f * part / total : 0.f;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of fractional digits; 
//	if count of value's mapped digits is more then that, printed "<X%", or exactly by default
//	@fieldWith: displayed width of value and '%' or '<' sign (excluding parentheses), or exactly if 0;
//	@parentheses: if true then parenthesize the value (not considering fieldWith)
string sPercent(float percent, BYTE precision=0, BYTE fieldWith=0, bool parentheses=false);

// Returns string represents the percent of part relatively total
//	@part: value represents desired %
//	@total: value represents 100%
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default;
//	should include a '%' mark
//	@parentheses: if true parenthesize the value
inline string sPercent(ULONG part, ULONG total,
	BYTE precision=0, BYTE fieldWith=0, bool parentheses=false) {
		return sPercent(Percent(part, total), precision, fieldWith, parentheses);
}

// Gets linear density (density per 1000 bs) of some elements
//	@cnt: number of elements
//	@len: length on which the density is determined
//	Used in isChIP (Imitator.cpp) and readDens (readDens.h).
inline float LinearDens(ULLONG cnt, chrlen len) { return len ? 1000.f * cnt / len : 0; }

// Prints horizontal line
//	@w: width of line
void PrintHorLine(int w);

#if defined _WIGREG || defined _BIOCC

// Align position to the up or down resoluation level
// f.e. by resoluation==5 pos 102 -> 100+relative, pos 104 -> 105++relative
//	@pos: chromosome's position
//	@res: resoluation
//	@relative: 0 or 1
//	1 used for 1-relative position (the first base is 1, WIG)
//	0 used for 0-relative position (BED)
//	return: aligned position
chrlen AlignPos(chrlen pos, BYTE res, BYTE relative);

#endif

//#ifdef OS_Windows
//string	Wchar_tToString(const wchar_t* wchar);
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar);
//#endif

/*** end of COMMON FUNCTION ***/

#ifdef _DUP_OUTPUT
#include <fstream>
// 'dostream' duplicates outstream to stdout & file
class dostream : public std::ostream
{
	std::ofstream file;
public:
	//inline dostream(std::ostream& s) : std::ostream(cout.rdbuf()) {}

	inline dostream() : std::ostream(cout.rdbuf()) {}

	inline ~dostream() { if (file.is_open())	file.close(); }		// in case of exception or holding execution by user

	// Open output file with given name
	//	return: true if file is open
	bool OpenFile(const string fname);

	template <typename T> dostream& operator<< (T val) {
		cout << val;
		file << val;
		return *this;
	}

	template <typename T> 
	inline void ToFile(T val) { file << val; }
};

//class dostream
//{
//    std::ostream &first, &second;
//public:
//    dostream(std::ostream &f, std::ostream &s) : first(f), second(s) {}
//	
//	template <typename T> dostream& operator<< (T val) {
//		first << val;
//		second << val;
//		return *this;
//	}
//};

extern dostream dout;		// stream's duplicator
#elif defined _WIGREG
	#define dout	cerr
#else
	#define dout	cout
	#define dostream ostream
#endif	// _DUP_OUTPUT

// 'Product' keeps Product oinfo
static struct Product
{
	static const string	Title;
	static const string	Version;
	static const string	Descr;

	// Gets length of product name (title plus version).
	//inline static BYTE MarkLength() { return BYTE(Title.length()); // + strlen(Version) + 1);
	//}

	// Gets title plus version.
	//inline static const string& Name() { return Title + string(1, HPH) + string(Version); }
} product;

typedef pair<float,float> pairVal;	// pair of float values

// 'Options' implements main() options and parameters treatment.
static class Options
{
private:
	// types of option values
	// do not forget to support the correlation with Options::TypeNames []
	// tVERS option can be absent, so this 'special' group should be started by tHELP
	enum valType {	// using 'enum class' in this case is inconvenient
		tUNDEF,		// not used
		tNAME,		// string value
		tCHAR,		// character value
		tINT,		// integer value
		tFLOAT,		// float value
		tLONG,		// long value
		tENUM,		// enum value
		tCOMB,		// combo value
		tPR_INT,	// tPAIR_INT pair of integers value
		tPR_FL,		// tPAIR_FL pair of floats value
		tHELP,		// help value
		tVERS,		// version value
		tSUMM,		// special value used to call a program from another program
	};
	
	enum eFlag {		// using 'enum class' in this case is inconvenient
		fNone = 0,		// no flags
		fOblig = 0x01,	// true if option is obligatory
		fOptnal = 0x02,	// true if option's value is optional
		fAllow0 = 0x04,	// true if option's value allowed be 0 while checking min value
		fHidden = 0x08,	// hidden option (not printed in help)
		fTrimmed = 0x10,// true if option has been processed already
		fWord = 0x20,	// true if option is stated as multi-chars
	};

	// Option signes
	struct Signs {
	private:
		BYTE	signs;
	public:
		inline Signs(int x)	{ signs = (BYTE)x; }	// to initialize obligatory in main()
		// Returns true if given sign is set
		inline bool Is(eFlag mark) const	{ return (signs & BYTE(mark)) != 0; }	// '!= 0' to avoid warning
		// Sets given sign
		inline void MarkAs(eFlag mark)	{ signs |= BYTE(mark); }

#ifdef DEBUG
		void Print() const {
			cout << int(Is(fWord)) << int(Is(fTrimmed)) << int(Is(fHidden)) << int(Is(fOblig));
		}
#endif
	};

	// structure 'Option' keeps options attributes.
	struct Option {
		const char	Char;		// option - character
		const char*	Str;		// option - string
			  Signs	Sign;		// option's signs
		const valType ValType;	// type of value
		const BYTE	OptGroup;	// option's category
			  float NVal;		// default or established numeric|enum value
		const float MinNVal;	// minimal permissible value;
								// for enum should be a first defined value
		const float MaxNVal;	// maximal permissible value or count of enum values
		const char*	SVal;		// string value or pointer to enum values array;
								// in the last case should be cast to char**.
								// If enum option hase ValRequired==false,
								// this pointer dousn't use and may be NULL
		const char*	Descr;		// tip string
		const char*	AddDescr;	// additional tip string

		// Return true if option value is escapable (not obligatory)
		inline bool IsValEsc() const { return Sign.Is(eFlag::fOptnal); }

		// Sets option value.
		//	@opt: option
		//	@isword: true if option is a word, false if option is a char
		//	@val: value of option
		//	@nextItem: next token after opt and val, or NULL
		//	@argInd: the current index in argc; increased by 1 if value is accepted
		//	return: 0 if success, -1 if not found, 1 if option or value is wrong
		int	SetVal(const char* opt, bool isword, char* val, char* nextItem, int& argInd);

		// Check option for obligatory.
		//	return: -1 if option is obligatory but not stated, otherwise 1
		int CheckOblig() const;

		// Prints option if it's obligatory
		inline void PrintOblig() const		{ if(Sign.Is(fOblig)) Print(false); }

		// Prints option if it belongs to a group g
		inline void PrintGroup(BYTE g) const { if(OptGroup == g) Print(true); }

		// Returns option name and value
		string ToStr(bool prVal) const;
		
		// Prints option in full or short way.
		//	@descr: if true, prints in full way: 
		//	signature, description (marks as Required if needed), default value, 
		//	otherwise signature only
		void Print(bool descr) const;

#ifdef DEBUG
		void Print() const;
#endif
	private:
		static const char EnumDelims[];	// specifies delimiter for enum [0] and combi [1] values

		// Checks digital value representation. Prints 'wrong val' message in case of failure
		//	@str: defined (no NULL) string  representing digital value
		//	return: true if digital value representation is correct
		bool IsValidFloat(const char *str, bool isInt, bool isPair = false);

		// Returns option's signature as a string
		//	@asPointed: true if returns signature as it was stated by user
		string NameToStr (bool asPointed) const;

		// Returns string represented pair of value's separated by delimiter.
		const string PairValsToStr(const pairVal* vals) const;

		// returns true if option value is required
		bool inline ValRequired() const { return MinNVal != vUNDEF; }

		// Checks limits and set numerical value
		//	@val: numerical value
		//	return: 1 if limits are exceeded, otherwise 0
		int SetTriedFloat(float val, float min, float max);

		// Checks and sets enum option value. Prints 'wrong val' message in case of failure
		//	@val: input value as C string or NULL if value is optional and not defined
		//	return: 0 if success, 1 if wrong value
		int SetEnum(const char* val);

		// Checks and sets enum option value. Prints 'wrong val' message in case of failure
		//	@val: input value as C string
		//	return: 0 if success, 1 if wrong value
		int SetComb(char* val);

		// Checks and sets pair option value
		//	@vals: pair of values as C string or NULL if value isn't set
		//	return: 0 if success, 1 if wrong values
		int SetPair(const char* vals, bool isInt);

		// Prints enum or combi values
		//	return: number of printed symbols
		BYTE PrintEnumVals() const;
		
		// Performs a case-insensitive search of given string value among enum values.
		//	@val: input value as string
		//	return: index of finded value in enum,
		//	or -1 if the value is not present in enum,
		//	or -2 if the wrong delimiter is encountered
		int GetEnumInd (const char* val);

		// Ouptuts option with error message to cerr
		//	@val: value or NULL
		//	@msg: error message about value
		//	@return: always 1
		int PrintWrong(const char* val, const string& msg=strEmpty) const;
	};

	// structure 'Usage' is used to output some Usage variants in PrintUsage()
	struct Usage {
		const int	Opt;		// option that should be printed in Usage as obligatory
		const char*	Par;		// prog parameter that should be printed in Usage as obligatory
		const bool	IsParOblig;	// true if parameter is obligatory, otherwise printed in square brackets
		const char* ParDescr;	// description of prog parameter that should be printed in Usage

		// Prints Usage params
		void Print(Option* opts) const;
	};

	static const char* sPrSummary;		// summary string printed in help
	static const char* sPrTime;			// run time string printed in help
	static const char* sPrUsage;		// usage string printed in help
	static const char* sPrVersion;		// version string printed in help
	static const char*	Booleans [];	// boolean values
	static const char*	TypeNames[];	// names of option value types in help
	static const char*	OptGroups[];	// names of option groups in help
	static const Usage	Usages	 [];	// content of 'Usage' variants in help
	static const BYTE	OptCount,		// count of options
						GroupCount,		// count of option groups in help
						UsageCount;		// count of 'Usage' variants in help
	static 	Option		List[];			// list of options. Option 'help' always shuld be
										// the last one, option 'version' - before last.
	
	// Check obligatory options and output message about first absent obligatory option.
	//	return: -1 if some of obligatory options does not exists, otherwise 1
	static int	CheckObligs();
	
	// Set option [with value] or splitted short options
	//	@opt: option without HYPHEN
	//	@val: option's value
	//	@nextItem: next token after opt and val, or NULL
	//	@argInd: the current index in argc; increased by 1 if value is accepted
	//	Return: 0 - success, 1 - false
	static int	SetOption(char* opt, char* val, char* nextItem, int& argInd);

	// Returns true if long option opt is defined
	static bool Find(const char* opt);

	// Ouptuts ambiguous option with error message to cerr
	//	@opt: option
	//	@isWord: true if option is a word
	//	@headMsg: message at the beginning
	//	@inOpt: initial option (in case of ambiguous composite)
	//	return: always 1
	static int PrintAmbigOpt(const char* opt, bool isWord,
		const char* headMsg, const char* inOpt = NULL);

	// Prints version
	//	return: always 1
	static int PrintVersion();

	static int PrintSummary(bool prTitle);

public:
	// 'PairVals' holds pair of values and their min and max limits
	class PairVals
	{
		pairVal vals[3];

	public:
		enum Type { SET=0, MIN=1, MAX=2 };

		PairVals(float val1, float val2, float min1, float min2, float max1, float max2) {
			vals[SET] = make_pair(val1, val2);
			vals[MIN] = make_pair(min1, min2);
			vals[MAX] = make_pair(max1, max2);
		}

		// Gets values
		inline const pairVal& Values(Type t = SET) const { return vals[t]; }
	};

	// Prints 'usage' information
	//	@title: if true prints title before information
	//	return: 1 if title is settinf to true, 0 otherwise
	static int PrintUsage (bool title);

	// Returns option name [and value]
	inline static string OptionToStr(int opt, bool prVal = false) { return List[opt].ToStr(prVal); }
	
	// Returns command line.
	//	@argc: count of main() parameters
	//	@argv: array of main() parameters
	static const string CommandLine(int argc, char* argv[]);

	// Reset int option value to 0
	inline static void ResetIntVal(int opt) {  List[opt].NVal = 0; }

	// Parses and checks main() parameters and their values.
	//	Output message if some of them is wrong.
	//	@argc: count of main() pearmeters
	//	@argv: array of main() pearmeters
	//	@obligPar: name of required application parameter or NULL if not required
	//	return: index of first parameter (not option) in argv[],
	//	or argc if it is absent,
	//	or negative if tokenize complets wrong
	static int Parse(int argc, char* argv[], const char* obligPar=NULL);
	
	// Get float value by index
	inline static float GetFVal	(int opt)	{ return List[opt].NVal; }
	// Get string value by index
	inline static const char* GetSVal(int opt)	{ return List[opt].SVal; }
	// Get booling value by index
	inline static bool GetBVal	(int opt)	{ return List[opt].NVal != 0; }
	// Get UINT value by index
	inline static UINT GetUIVal	(int opt)	{ return UINT(List[opt].NVal); }
	// Get int value by index
	inline static int GetIVal	(int opt)	{ return int(List[opt].NVal); }
#ifdef _READS
	// Get read duplicates level
	inline static char GetRDuplLevel(int opt) { return GetBVal(opt) ? vUNDEF : 0; }
#endif

	// Returns true if the option value is assigned by user
	static bool Assigned	(int opt)	{ return List[opt].Sign.Is(fTrimmed); }
	// True if maximal enum value is setting
	static bool IsMaxEnum	(int opt)	{ return List[opt].NVal == List[opt].MaxNVal - 1; }
	// Get maximal permissible numeric value by index
	inline static UINT GetMaxIVal(int opt)	{ return UINT(List[opt].MaxNVal); }
	// Returns C string 'ON' or 'OFF'
	inline static const char* BoolToStr(int opt)	{ return Booleans[GetBVal(opt)]; }
	inline static const char* BoolToStr(bool val)	{ return Booleans[val]; }
	// Gets pointer to the C string contained 'ON' or 'OFF'
	//inline static const char* BoolToStr(bool val)	{ return Booleans[int(val)]; }

	// Returns file name by index: if value is not oblig and is not specified, then defName with given extention
	static const string GetFileName(int opt, const char* defName, const string& ext="_out.txt");

#ifdef DEBUG
	static void Print();
#endif
} options;

#define ErrWARNING	Err::NONE

// 'Error' implements error's (warning's) messages and treatment.
class Err
{
public:
	enum eCode {
		NONE,		// warning message
		MISSED,		// something missed
		F_NONE,		// no file
		D_NONE,		// no directory
		FD_NONE,	// nor file nor directory 
		F_MEM,		// TxtFile(): memory exceeded
		F_OPEN,		// TxtFile(): file open error
		F_CLOSE,	// TxtFile(): file close error
		F_READ,		// TxtFile(): file read error
		F_EMPTY,
		F_BIGLINE,	// TxtFile(): too big line
		FZ_MEM,		// TxtFile(): not enough gzip buffer
		FZ_OPEN,	// TxtFile(): gzip open error
		FZ_BUILD,	// TxtFile(): doen't support gzip
		F_WRITE,	// file write error
		//F_FORMAT,
#ifndef _FQSTATN
		TF_FIELD,	// TabFile: number of fields is less than expected
		TF_EMPTY,	// TabFile: no records
#endif
		EMPTY
	};

private:
	static const char* _msgs[];
	enum eCode	_code;			// error code
	char * _outText;			// output message
	//unique_ptr<char[]> _outText;			// output message

	// Initializes _outText by C string contained message kind of
	// "<sender>: <text> <specifyText>".
	void set_message(const char* sender, const char* text, const char* specifyText=NULL);

	//inline void set_message(const string& sender, const char* text, const char* specifyText=NULL) {
	//	set_message(sender==strEmpty ? NULL : sender.c_str(), text, specifyText);
	//}

public:
	static const char* FailOpenOFile;

	// Returns string containing file name and issue number.
	//	@issName: name of issue
	//	@issNumb: number of issue
	//	@fName: file name
	//static const string IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName);

	// Gets message "no @fileName.@fileExt[.gz] files in this directory"
	static const string MsgNoFiles (const string & fileName, const string fileExt);

	// Code-attached constructor.
	//	@code: exception/warning message as code
	//	@sender: name of object who has generated exception/warning, or NULL if no sender
	//	@specifyText: aditional text to specify exception/warning message
	inline Err(eCode code, const char* sender, const char* specifyText=NULL): _code(code) {
		set_message(sender, _msgs[code], specifyText);
	}

	// Code-attached constructor.
	//	@code: exception/warning message as code
	//	@sender: name of object who has generated exception/warning, or NULL if no sender
	//	@specifyText: aditional text to specify exception/warning message
	inline Err(eCode code, const char* sender, const string& specifyText) : _code(code) {
		set_message(sender, _msgs[code], specifyText.c_str());
	}

	// C-string-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const char* text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text);
	}

	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const char* text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text);
	}

	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const string& text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text.c_str());
	}

	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const string& text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text.c_str());
	}

	// copy constructor
	Err(const Err& src);

	~Err() { delete[] _outText; _outText = NULL; }

	inline const char* what() const /*throw()*/ { return _outText; }

	inline eCode Code() const		{ return _code; }

	// Returns point to the aditional text to specify exception/warning message or NULL
	//inline const char* SpecifyText() const	{ return _specifyText; }

	// Throws exception or outputs Err message.
	//	@throwExc: if true then throws exception, otherwise outputs Err message
	//	@eol: if true then carriage return after Err message
	void Throw(bool throwExc = true, bool eol = true);
	
	// Outputs warning
	//	@prefix: output ": " before "WARNING"
	//	@eol: if true then carriage return after Err message
	void Warning(bool eol = true, bool prefix = false);
};

// 'FileSystem' implements common file system routines
static class FS
{
	// Returns true if file system's object exists
	//	@name: object's name
	//	@st_mode: object's system mode
	static bool IsExist(const char* name, int st_mode);

	// Checks if file system's object doesn't exist
	//	@name: object's name
	//	@st_mode: object's system mode
	//	@throwExcept: if true throws exception,
	//	@ecode: error's code
	//	otherwise outputs Err message as warning without LF
	//	return: true if file or directory doesn't exist
	static bool CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode);

	// Searches through a file name for the any extention ('/' or '\' insensible).
	//	return: the index of the DOT matched extention; otherwise npos
	static size_t GetLastExtPos	(const string &fname);

	// Returns true if file name has specified extension ignoring zip extension. Case insensitive
	static bool HasCaseInsExt(const string &fname, const string &ext, bool knownZip, bool composite = true);

	// Returns file name without extention ('/' or '\' insensible)
	inline static string const FileNameWithoutLastExt (const string& fname) {
		return fname.substr(0, GetLastExtPos(fname));
	}


public:
	// === file size

	// Gets size of file or -1 if file doesn't exist
	static LLONG Size 	(const char*);

	// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
	static LLONG UncomressSize	(const char*);

	// === check dir/file existing

	// Returns true if file exists
	inline static bool IsFileExist	 (const char* name) { return IsExist(name, S_IFREG); }
	
	// Returns true if directory exists
	inline static bool IsDirExist	 (const char* name) { return IsExist(name, S_IFDIR); }
	
	// Returns true if file or directory exists
	inline static bool IsFileDirExist(const char* name) { return IsExist(name, S_IFDIR|S_IFREG); }
	
	// Checks if file doesn't exist
	//	@name: name of file
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without LF
	//	return: true if file doesn't exist
	inline static bool CheckFileExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFREG, throwExcept, Err::F_NONE);
	}

	// Checks if directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without LF
	//	return: true if file or directory doesn't exist
	inline static bool CheckDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR, throwExcept, Err::D_NONE);
	}

	// Checks if file or directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without LF
	//	return: true if file or directory doesn't exist
	inline static bool CheckFileDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR|S_IFREG, throwExcept, Err::FD_NONE);
	}

	// Throws exsception if file or directory doesn't exist
	//	@name: name of file or directory
	//	@ext: file extention; if set, check for file first
	//	@throwExcept: if true throws exception,
	//	otherwise outputs Err message as warning without LF
	//	return: true if file or directory doesn't exist
	static bool CheckFileDirExist(const char* name, const string & ext, bool throwExcept);

	// === check dir/file name

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optsVal: Options char* value
	//	return: pointer to the checked file name
	static const char* CheckedFileDirName	(const char* name);

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optVal: Options value
	//	return: pointer to the checked file name
	inline static const char* CheckedFileDirName	(int optVal) {
		return CheckedFileDirName(Options::GetSVal(optVal));
	}

	// Returns a pointer to the file name checked if it exist, otherwise throws exception
	//	@name: pointer to the file name
	//	return: pointer to the checked file name
	static const char* CheckedFileName	(const char* name);

	// Returns a pointer to the path checked if it exist, otherwise throws exception
	//	@opt: Options value
	//	return: pointer to the checked path
	static const char* CheckedDirName	(int opt);

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@opt: Options value
	//	return: pointer to the checked file name
	//	return: pointer to the checked file name
	inline static const char* CheckedFileName	(int opt) {
		return CheckedFileName(Options::GetSVal(opt));
	}

	// Returns true if directory is open to writing
	inline static bool IsDirWritable(const char* name) { 
#ifdef OS_Windows
		return true;
#else
		return !access(name, W_OK); 
#endif
	}

	// === check file extension

	// Returns true if file has any extension.
	//	@fname: file name
	inline static bool HasExt	(const string &fname) { return GetLastExtPos(fname) != string::npos; }

	// Returns true if file has a specified  extension. Case insensitive search
	//	@ext: extension includes dot symbol
	//	@composite: true if extension is (strictly last)
	inline static bool HasExt	(const string& fname, const string& ext, bool composite = true) {
		return HasCaseInsExt(fname, ext, HasGzipExt(fname), composite);
	}

	// Returns true if file has '.gz' extension. Case insensitive search
	inline static bool HasGzipExt(const string& fname) { return HasCaseInsExt(fname, ZipFileExt, false); }

	// Returns string containing real file extension (without gzip).
	//	@fname: pointer to the file name
	//	return: string containing real file extension or empty string if no real extention
	static string const GetExt(const char* fname);

	// === dir/file name transform

	// Returns file name without extention ('/' or '\' insensible)
	static string const FileNameWithoutExt (const string& fname) {
		return FileNameWithoutLastExt( HasGzipExt(fname) ? FileNameWithoutLastExt(fname) : fname );
	}

	// Returns true if file name is short (without path)
	static bool IsShortFileName(const string& fname);

	// Returns short file name by long one
	//	@fname: long file name
	static string const ShortFileName (const string& fname);

	// Returns directory name by long file name
	//	@fname: long file name
	//	@addSlash: true if slash sould be added at the end
	static string const DirName (const string& fname, bool addSlash = false);

	// Returns the name of last subdirectory by long file name
	//	@name: long dir name
	static string const LastSubDirName (const string& name);

	// Returns the name of last subdirectory
	//	@name: long dir name
	static string const LastDirName (const string& fname);

	// Returns the name ended by slash without checking
	static string const MakePath(const string& name);

	// === files in dir

#if !defined _WIGREG && !defined _FQSTATN
	// Fills external vector of strings by file's names found in given directory
	// Implementation depends of OS.
	//	@files: external vector of strings that should be filled by file's names
	//	@dirName: name of directory
	//	@ext: file's extention as a choosing filter
	//	@all: true if all files with given extention should be placed into external vector,
	//	otherwise only one (any)
	//	return: true if files with given extention are found
	static bool GetFiles (vector<string>& files, const string& dirName, const string& ext, bool all = true);
#endif	// _WIGREG, _FQSTATN

//	inline static void	Delete		(const char* fname) {
//#ifdef OS_Windows
//		DeleteFile(fname);
//#else
//		unlink(fname);
//#endif
//	}

} fs;

// Basic class for wall time measurement
class TimerBasic
{
protected:
	// Prints elapsed wall time interval
	//	@elapsed: elapsed time in seconds
	//	@title: string printed before time output
	//	@parentheses: if true then output time in parentheses
	//	@isLF: if true then ended output by LF
	static void Print(long elapsed, const char *title, bool parentheses, bool isLF);

	mutable time_t	_startTime;
	mutable bool	_enabled;	// True if local timing is enabled

	// Creates a new TimerBasic
	//	@enabled: if true then set according total timing enabling
	TimerBasic(bool enabled = true) : _startTime(0) { _enabled = enabled ? Enabled : false;	}
	
	// Stops timer and return elapsed wall time in seconds
	long GetElapsed() const;
	
public:
	// True if total timing is enabled
	static bool		Enabled;

	// True if instance timing is enabled
	inline bool IsEnabled() const { return _enabled; }

	// Starts timer
	inline void Start()		{ if(_enabled) time( &_startTime ); }
};

// 'Timer' measures the single wall time interval
class Timer : public TimerBasic
{
private:
	static clock_t	_StartCPUClock;

public:
	// Starts enabled CPU timer, if it is enabled
	inline static void StartCPU()	{ if( Enabled ) _StartCPUClock = clock(); }
	
	// Stops enabled CPU timer and print elapsed time
	//	@isLF: if true then ended output by LF
	static void StopCPU(bool isLF=true) {
		if(Enabled)	Print((clock()-_StartCPUClock)/CLOCKS_PER_SEC, "CPU: ", false, isLF);
	}

	// Creates a new Timer and starts it if timing is enabled
	//	@enabled: if true then set according total timing enabling
	inline Timer(bool enabled = true) : TimerBasic(enabled) { Start(); }
	
	// Stops enabled timer and prints elapsed time with title
	//	@title: string printed before time output
	//	@parentheses: if true then output time in parentheses
	//	@isLF: if true then ended output by LF
	void Stop(const char *title, bool parentheses, bool isLF) {
		if(_enabled)	Print(GetElapsed(), title, parentheses, isLF);
	}

	// Stops enabled timer and prints elapsed time
	//	@offset: space before time output
	//	@parentheses: if true then output time in parentheses
	//	@isLF: if true then ended output by LF
	void Stop(int offset = 0, bool parentheses = false, bool isLF = false);

	// Stops enabled timer and prints elapsed time
	//	@parentheses: if true then output time in parentheses
	//	@isLF: if true then ended output by LF
	//inline void Stop(bool parentheses = false, bool isLF = true)	{
	//	Stop(NULL, parentheses, isLF); }
};

#ifdef _TEST

// 'Stopwatch' measures the sum of wall time intervals
class Stopwatch : public TimerBasic
{
	private:
		mutable long	_sumTime;
		mutable bool	_isStarted;		// true if Start() was called even ones

	public:
		inline Stopwatch() : _sumTime(0), _isStarted(false), TimerBasic() {}

		// Starts Stopwatch
		void Start() const		{ ((TimerBasic*)this)->Start(); _isStarted = true; }

		// Stops Stopwatch
		//	@title: if not empty, and if instance was launched, output sum wall time with title
		//	'const' to apply to constant objects
		void Stop(const string title = strEmpty) const;
};

#endif	// _TEST

// 'Stopwatch' measures the sum of CPU time (clocks) intervals
class StopwatchCPU
{
	private:
		clock_t	_clock;
		clock_t	_sumclock;

	public:
		inline StopwatchCPU() : _sumclock(_clock=0) {}

		inline void Start(bool reset=false)	{ _clock = clock(); if(reset) _sumclock = 0; }

		// Stops StopwatchCPU
		//	@title: string printed before time output
		//	@print: if true time should be printed
		//	@isLF: if true then ended output by LF
		void Stop(const char* title, bool print = false, bool isLF = false);
};

static class Mutex
{
public:
	enum class eType { OUTPUT, INCR_SUM, WR_BED, WR_SAM, WR_FQ, 
#ifdef _ISCHIP
		WR_RDENS, WR_BG, WR_BGPOS, WR_BGNEG,
#endif
		NONE };
#ifdef _MULTITHREAD
private:
	static bool	_active;	// true if the mutex really should work
	//static pthread_mutex_t	_mutexes[];
	static mutex	_mutexes[];
public:
	// Returns true if Mutex type is not NONE
	inline static bool IsReal(eType mtype) { return mtype != eType::NONE;	}
	inline static void Init(bool active) { _active = active; }
	//inline static void Finalize() {};
	static void Lock(eType type) { if (_active) _mutexes[int(type)].lock();	}
	static void Unlock(eType type) { if (_active) _mutexes[int(type)].unlock(); }
#endif
} myMutex;

// 'Chrom' establishes correspondence between chromosome's ID and it's name.
static class Chrom
/**********************************************************************************

===== TERMS ======
MARK: human view of chrom's number:	1, ..., Y
ABBREVIATION CHROM's NAME:			chr1, ..., chrY
SHORT CHROM's NAME:					chrom 1, ..., chrom Y
TITLE CHROM's NAME:					chromosome X, ..., chromosome Y
PREFIX: substring before mark:		chr | chrom | chromosome.

===== CHROM DATA DEPENDENCIES =====
User selectable chromosome (customization)		needs:	absolute | relative chrom ID discipline
Control of the item belonging to the chromosome	needs:	absolute | relative chrom ID discipline
Chrom size control								needs:	absolute | relative chrom ID discipline + ChromSizes
Using BAM										needs:	relative chrom ID discipline

===== CHROM's ID NUMBERING DISCIPLINES =====
RELATIVE DISCIPLINE:
The autosomes' ID corresponds to the chrom number reduced by 1 (0-based numbering)
The heterosomes' ID starts with the last autosomes' ID plus 1, in the next order: X, Y, M
This discipline is coming from BamTools API.

ABSOLUTE DISCIPLINE:
The autosomes' ID numering is the same as relative.
The heterosomes' ID is a code of mark character.

Chrom class distinguishes these disciplines by the value of the 'firstHeteroID' variable:
a zero value indicates an absolute discipline, a non-zero value indicates a relative one.

***********************************************************************************/
{
public:
	static const char*	Abbr;				// Chromosome abbreviation
	static const BYTE	MaxMarkLength = 2;	// Maximal length of chrom's mark
	static const BYTE	MaxAbbrNameLength;	// Maximal length of abbreviation chrom's name
#ifndef _FQSTATN
	static const string	Short;				// Chromosome shortening; do not convert to string in run-time
	static const chrid	UnID = -1;			// Undefined ID (255)
	static const chrid	Count = 24;			// Count of chromosomes by default (for container reserving)
	static const BYTE	MaxShortNameLength;	// Maximal length of short chrom's name
	static const BYTE	MaxNamedPosLength;	// Maximal length of named chrom's position 'chrX:12345'

private:
	static const string	sTitle;		// Chromosome title; do not convert to string in run-time
	static		 BYTE	CustomOpt;	// user chrom option number; used to set custom cID
	static const char*	Marks;		// heterosome marks
	static const string	UndefName;	// string not to convert in run-time

	static chrid cID;				// user-defined chrom ID
	static chrid firstHeteroID;		// first heterosome (X,Y) ID

	// Returns true if relative ID numbering discipline is used
	inline static bool IsRelativeID() { return firstHeteroID; }

	// Gets heterosome ID by mark without control, or undefined ID
	static chrid HeteroID	(const char cMark);

	// Gets heterosome ID by case insensitive mark, or undefined ID
	inline static chrid CaseInsHeteroID(const char cMark) { return HeteroID( toupper(cMark) ); }

	// Gets chrom ID by case insensitive mark
	//	firstHeteroID should be initialized!
	static chrid CaseInsID	(const char* cMark);

public:
	//*** common methods

	// Returns true if chrom is autosome, false for somatic; for relative ID discipline only
	inline static bool IsAutosome(chrid cid) { return cid < firstHeteroID; }

	// Gets the length of prefix, or -1 if name is not finded
	static short PrefixLength(const char* cName);

	// Sets relative ID numbering
	//static void SetRelativeID() { if(!IsRelativeID()) firstHeteroID = 1; }
	
	//*** ID getters

	// Gets chrom's ID by name without control of case insensitivity and undefined ID
	//	@cName: chrom's name
	//  @prefixLen: length of name prefix
	static chrid ID(const char* cName, size_t prefixLen=0);
	
	// Gets chrom ID by mark without control of case insensitivity and undefined ID
	inline static chrid ID(const string& cMark) { return ID(cMark.c_str()); }

	// Gets chrom ID by abbreviation name without control of case insensitivity and undefined ID
	inline static chrid IDbyAbbrName(const char* cAbbrName) { return ID(cAbbrName, strlen(Abbr)); }

	// Gets chrom's ID by long name without control of case insensitivity and undefined ID
	inline static chrid IDbyLongName(const char* cName) { return ID(cName, PrefixLength(cName)); }

	//*** validation methods

	// Checks chrom ID and returns this ID or undefined ID; used in BamInFile::GetNextChrom() only
	inline static chrid ValidateID(chrid cID) { return cID < firstHeteroID + strlen(Marks) ? cID : UnID; }

	// Validates chrom name and returns chrom ID
	//	@cName: string of arbitrary length containing chrom's name
	//  @prefixLen: length of prefix before chrom's mark
	static chrid ValidateID(const char* cName, size_t prefixLen = 0);

	// Validates chrom mark and returns chrom ID
	//	@cMark: string of arbitrary length, starting with chrom's mark
	inline static chrid ValidateID(const string& cMark)	{ return ValidateID(cMark.c_str()); }

	// Validates chrom abbreviation name and returns chrom ID
	//	@cName: string of arbitrary length, starting with chrom's name
	inline static chrid ValidateIDbyAbbrName(const char* cName) { return ValidateID(cName, strlen(Abbr)); }

	// Validates all chrom ID by SAM header data, establishes relative ID discipline, and sets custom ID
	template<typename Functor>
	static void ValidateIDs(const string& samHeader, Functor f)
	{
		for (const char* header = samHeader.c_str();
			header = strstr(header, Abbr);
			header = strchr(header, LF) + strlen("\n@SQ\tSN:") )
		{
			chrid cID = ValidateIDbyAbbrName(header);
			header = strchr(header, TAB) + strlen("\tLN:");
			f(cID, header);
		}	
		SetCustomID(true);
	}

	//*** custom ID getters, setters

	// Gets custom chrom's ID
	inline static chrid CustomID()	{ return cID; }

	// Returns true if chrom is set by user, or if no chroms have been set by user (by default)
	//	@cid: chrom id specified by user
	inline static bool IsCustom(chrid cid = UnID) { return cID == UnID || cID == cid; }

	// Sets custom chrom ID with control
	//	@prColon: if true then print ": " before exception message
	//	exception: wrong chrom
	static void SetCustomID(bool prColon = false);

	// Sets number of 'custom chrom' progr option
	//	@absIDNumb: true if absolute ID numbering discipline is applied
	static void SetCustomOption(int opt/*, bool absIDNumb = false*/);

	//*** work with name

	static const string Title(bool pl = false) { return pl ? (sTitle + 's') : sTitle;	}

	// Gets mark length by ID; used in WigMap (bioCC)
	inline static size_t MarkLength(chrid cID) {
		return cID == UnID ? UndefName.length() : (cID > firstHeteroID || cID < 9 ? 1 : 2);
	}

	// Returns mark by ID
	static const string Mark(chrid cid);

	// Locate chrom's mark in string.
	//	@str: string checked for chrom number
	//	return: pointer to the chrom number in str, or a null pointer if Chrom::Abbr is not part of str
	static const char* FindMark(const char* str);

	// Gets chrom's abbreviation name'chrX'
	//	@numbSep: if true then separate chrom's number
	static string AbbrName(chrid cid, bool numbSep = false);

	// Gets short name 'chrom X'
	inline static string ShortName(chrid cid)	{ return Short + SPACE + Mark(cid); }

	// Gets title name 'chromosome X' or 'chromosomes'
	//	@cid: chromosome's ID or UnID if plural
	inline static string TitleName(chrid cid = UnID) { return Short + (cid==UnID ? "s" : Mark(cid)); }

	inline static const string Absent(chrid cid, const string & what) {
		return AbbrName(cid) + " is absent in " + what + " file: skipped";
	}

#endif	// _FQSTATN
} chrom;

// 'Region' represents a simple region within nucleotides array (chromosome).
struct Region
{
	chrlen Start;	// start position of the region in standard chromosomal coordinates
	chrlen End;		// end position of the region in standard chromosomal coordinates

	inline Region(chrlen start = 0, chrlen end = 0) : Start(start), End(end) {}

	// Gets length of region.
	// The End is not included in the bases https://genome.ucsc.edu/FAQ/FAQformat.html#format1
	inline chrlen Length()	const { return End - Start; }

	inline bool Empty()		const { return !End; }

	inline chrlen Centre()	const { return Start + (Length() >> 1); }

	// Initializes instance
	inline void Set(chrlen start = 0, chrlen end = 0) { Start = start; End = end; }

	inline bool operator==(const Region& r) const { return Start == r.Start && End == r.End; }
	inline bool operator!=(const Region& r) const { return Start != r.Start || End != r.End; }

	// Returns true if this instance is invalid
	inline bool Invalid() const { return Start >= End; }

	// Returns true if Region r is covered by this instance.
	//	@r: tested Region; should be sorted by start position
	bool PlainCover(const Region& r) const { return r.End <= End && r.Start >= Start; }

	// Returns true if Region r is adjoined with this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool Adjoin(const Region& r) const { return r.Start == End; }

	// Returns true if Region r is crossed with this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool Cross(const Region& r) const { return r.Start < End&& r.End > Start; }

	// Compares two Regions by start position. For sorting a container.
	static inline bool CompareByStartPos(const Region& r1, const Region& r2) {
		return r1.Start < r2.Start;
	}

	// Extends Region with chrom length control.
	// If extended Region starts from negative, or ends after chrom length, it is fitted.
	//	@extLen: extension length in both directions
	//	@cLen: chrom length; if 0 then no check
	void Extend(chrlen extLen, chrlen cLen);

#ifdef _DEBUG
	inline void Print() const { cout << Start << TAB << End << LF; }
#endif
};

// 'Regions' represents a container of defined regions within chromosome
// Defined in TxtFile.h since it's used in Fa class
class Regions
{
protected:
	vector<Region> _regions;

public:
	typedef vector<Region>::const_iterator Iter;

	// Iterator to region Begin position
	inline const Iter Begin()	const { return _regions.begin(); }

	// Iterator to region End position
	inline const Iter End()		const { return _regions.end(); }

	// Default (empty) constructor to form Chroms collection
	inline Regions() {}

	// Single region constructor
	inline Regions(chrlen start, chrlen end) { _regions.emplace_back(start, end); }

	// Copying constructor
	//inline Regions(const Regions& rgns) { _regions = rgns._regions;	}

	// Gets total length of regions.
	//chrlen Length() const;

	// Gets count of regions.
	inline chrlen Count()		const { return chrlen(_regions.size()); }

	// Gets first start position.
	inline chrlen FirstStart()	const { return _regions.front().Start; }

	// Gets last end position.
	inline chrlen LastEnd()		const { return _regions.back().End; }

	// Gets conditionally defined length: distance between first start and last end
	inline chrlen DefLength()	const { return LastEnd() - FirstStart(); }

	//Regions& operator=(const Regions& rgn);

	inline const Region& operator[](chrlen ind) const { return _regions[ind]; }

	// Reserves container's capacity.
	//	@count: reserved number of regions. The real number may be differ.
	inline void Reserve(chrlen count) { _regions.reserve(count); }

	// Clears container.
	inline void Clear() { _regions.clear(); }

	inline void Add(const Region& rgn) { _regions.push_back(rgn); }

	inline void Add(chrlen start, chrlen end) { _regions.emplace_back(start, end); }

	// Copies subregions
	inline void Copy(const vector<Region>& source, chrlen start, chrlen stop) {
		_regions = vector<Region>(source.begin() + start, source.begin() + stop + 1);
	}

#if defined _READDENS || defined _BIOCC

	// Returns an iterator referring to the past-the-end element, where end is external
	//	@curr_it: region's const iterator, from which the search is started
	//	@end: external pre-defined end coordinate
	Iter ExtEnd(Iter curr_it, chrlen end) const;

	// Initializes this instance by intersection of two Regions.
	void FillOverlap(const Regions& regn1, const Regions& regn2);

	// Initializes this instance by inverted external Regions.
	//	@regn: external Regions
	//	@masEnd: the maximum possible end-coordinate of region:
	//	the chromosome length in case of nucleotides sequance.
	void FillInvert(const Regions& regn, chrlen maxEnd);

	// Initializes this instance by external Regions.
	inline void Copy(const Regions& rgns) { _regions = rgns._regions; }

#endif	// _READDENS, _BIOCC
#ifdef DEBUG
	void Print() const;
#endif
};


// 'MemStatus' outputs RAM remainder  
//static class MemStatus
//{
//	static bool	_enable;
//	static LLONG _startVolume;
//
//	static LLONG getAvailMemory();
//public:
//	static void StartObserve(bool enable);
//	static void StopObserve();
//} ramControl;

#endif	// _COMMON_H