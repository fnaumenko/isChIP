/**********************************************************
common.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 22.07.2019
-------------------------
Provides common functionality
***********************************************************/

#pragma once
#include "def.h"
#include <stdint.h>		// uint32_t
#include <string>
#include <iostream>	
#include <sstream>		// for NSTR( x )
#include <iomanip>		// setprecision(), setw()
#include <vector>
#include <sys/stat.h>	// struct stat

#ifdef __unix__
	#include <unistd.h>
	#include <stdlib.h>
	#include <memory.h>
	#include <limits.h>
	#include <math.h>
	#include <dirent.h>
	#include <stdio.h>
#ifdef _MULTITHREAD
	#include <pthread.h>
	#define InterlockedExchangeAdd	__sync_fetch_and_add
	#define InterlockedIncrement(p)	__sync_add_and_fetch(p, 1)
#endif
	#include <string.h>		// strerror_r()
	#include <stdexcept>	// throw std exceptions
	//#include <errno.h>	// pthread_join errcodes

	typedef unsigned char	BYTE;
	typedef	unsigned short	USHORT;
	typedef	unsigned int	UINT;
	typedef unsigned long	ULONG;
	typedef unsigned long long ULLONG;
	typedef long long LLONG;
	typedef struct stat struct_stat64;
	//typedef off64_t __int64;
	typedef void*	thrRetValType;
	//#define retThreadValTrue	(void*)1

	#define _stricmp strcasecmp	// case-sensitive comparison
	#define _fseeki64 fseeko64
	#define _ftelli64 ftello64
	#define _fileno	fileno
	#define _gcvt	gcvt
	#define _stat64 stat
	//#define _itoa	itoa
	#define fopen	fopen64
	#define TRUE	1
	#define FALSE	0
	#define	LOCALE_ENG	"en_US.utf8"
	//#define SLASH '/'	// standard Linux path separator
	//#define F_READ_MODE O_RDONLY
#elif defined _WIN32
	#define OS_Windows
	#include <windows.h>
#ifdef _MULTITHREAD
	#include <process.h>	    // _beginthread, _endthread
	#define pthread_t HANDLE
	#define pthread_mutex_t CRITICAL_SECTION
#endif		// OS Windows
	typedef unsigned __int64 ULLONG;
	typedef __int64 LLONG;
	typedef struct __stat64 struct_stat64;
	typedef UINT	thrRetValType;
	//#define retThreadValTrue	1

	#define atol _atoi64
	#define isnan _isnan
	#define	LOCALE_ENG	"English"
	//#define SLASH '\\'		// standard Windows path separator
	//#define REAL_SLASH '/'	// is permitted in Windows too
#endif

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
typedef USHORT		readlen;	// type length of Read
typedef uint32_t	chrlen;		// type length of chromosome
typedef chrlen		fraglen;	// type length of fragment
typedef ULONG		genlen;		// type length of genome

#define	CHRLEN_UNDEF	-1	// undefined length of chromosome
#define	CHRLEN_MAX		-1	// max length of chromosome
#define thrRetValFalse	0
#define CHRLEN_CAPAC	10	// capacity of max chrom length;
							// may be count by DigitsCount() every time,
							// but is the same if chrlen defined as int or long,
							// so is defined as static value

#define cNULL	'\0'
#define cN		'N'
#define HPH		'-'
#define USCORE	'_'
#define BLANK	' '
#define sBLANK	" "
#define QUOT	'\''
#define DOT		'.'
#define COLON	':'
#define PERS	'%'
#define AT		'@'
#define PLUS	'+'
#define HASH	'#'
#define TAB		'\t'
#define EOL		'\n'	// LF, 10
#define CR		'\r'	// 13
//#define sBACK	"\b"

using namespace std;

#define SepCm		", "	// comma separator
#define SepCl		": "	// colon separator
#define SepDCl		":: "	// double colon separator
#define SepSCl		"; "	// semicolon separator
#define SepClTab	":\t"	// colon + tab separator
#define Equel		" = "

static const string ZipFileExt = ".gz";
static const string strEmpty = "";

static const char* Done = " done";
static const char* Notice = "NOTICE: ";
static const string Total = "total";
static const char* Version = "version";
static const char* UnitDens = "r/kbp";	// unit of density measure
#ifdef _ISCHIP
#define SepGroup	";  "
#endif	// _ISCHIP
#if !defined _WIGREG && !defined _FQSTATN
static const char* sTemplate = "template";
#endif

/*** COMMON MACROS & FUNCTION ***/

#define vUNDEF	-1	// undefined value
#define NO_VAL	-1	// option value is prohibited
#define NO_DEF	-1	// do not print option default value

// Digital to STRing
// Returns value's string representation. http://rootdirectory.de/wiki/NSTR()
// Instead of std::to_string(x) [C++11] because of compatibility
#define NSTR(x) static_cast<ostringstream&>( ostringstream() << dec << (x) ).str()

// Returns two int value's separated by delimiter string representation.
#define NNSTR(x, delim, y) static_cast<ostringstream&>(ostringstream()<<dec<<int(x)<<delim<<int(y)).str()

// Byte to STRing
#define BSTR(x) static_cast<ostringstream & >( ostringstream() << dec << int(x) ).str()

// Float to STRing
// Returns float's string representation.
//	@p: precision (number of digits after separator)
//#define FSTR( x, p ) static_cast<ostringstream & >( ostringstream() << dec << setprecision(p)<< x ).str()

// Returns number of value's symbols, including minus & separator
//#define SCNT(x)	NSTR(x).length()

// Returns number of float's symbols, including minus and dot
//	@p: precision (number of digits after separator)
//#define SFCNT(x, p)	FSTR(x, p).length()

// Gets number of digist in a integral value
//	@val: integral value
//	@isLocale: if true then adds number of '1000' separators
//	return: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCount (LLONG val, bool isLocale = false);

// Returns percent of @part relatively @total
inline float Percent(ULLONG part, ULLONG total) { 
	return total ? 100.f * part / total : 0.f;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that, printed "<X%", or exactly by default
//	@fieldWith: the width of the display field insine parentheses, or exactly by default;
//	should include '%' and possible '<' marks
//	@parentheses: if true then parenthesize the value (not considering fieldWith)
string	sPercent(float percent, BYTE precision=0, BYTE fieldWith=0, bool parentheses=false);

// Returns string represents the percent of part relatively total
//	@part: value represents desired %
//	@total: value represents 100%
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default;
//	should include a '%' mark
//	@parentheses: if true parenthesize the value
inline string	sPercent(ULLONG part, ULLONG total,
	BYTE precision=0, BYTE fieldWith=0, bool parentheses=false) {
		return sPercent(Percent(part, total), precision, fieldWith, parentheses);
}

// Gets Read sensity per 1000 bs
//	@rCnt: number of Reads
//	@densLen: length on which the density is determined
inline float ReadDens(ULLONG rCnt, chrlen densLen) {
	return rCnt ? 1000.f * rCnt / densLen : 0;
}

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

#ifdef _NO_DOUT
	#define dout	cout
	#define dostream ostream
#elif defined _WIGREG
	#define dout	cerr
#else
#include <fstream>
// 'dostream' duplicates outstream to stdout & file
class dostream : public std::ostream
{
    //std::ostream& dupl;
    std::ofstream file;
public:
    //inline dostream(std::ostream& s) : dupl(s), std::ostream(cout.rdbuf()) {}
    inline dostream() : std::ostream(cout.rdbuf()) {}

	inline ~dostream() { if(file.is_open())	file.close(); }		// in case of holding execution by user

	// Open output file with given name
	//	return: true if file is open
	bool OpenFile(const string fname);
	
	template <typename T> dostream& operator<< (T val) {
		cout << val;
		file << val;
		return *this;
	}
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

#endif	// _NO_DOUT

// 'Gr' defines ground types
static class Gr
{
	static const char* title[];
public:
	enum Type {
		FG = 0,	// foreground
		BG = 1	// background
	};
	// Count of grounds
	static const int Cnt = 2;

	// Maximum title length
	static const BYTE TitleLength = 2;

	// Gets ground title
	inline static const char* Title(Type g) { return title[g]; }
} ground;

// 'Product' keeps Product info
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

typedef pair<float,float> pairVal;	// pair of values

// 'Options' implements main() options and parameters treatment.
static class Options
{
private:
	// types of option values
	// do not forget to support the correlation with Options::TypeNames []
	// tVERS option can be absent, so this 'special' group should be started by tHELP
	enum valType { tUNDEF, tNAME, tCHAR, tINT, tFLOAT, tLONG, tENUM, tCOMB,
		tPAIR_INT, tPAIR_FL, tHELP, tVERS, tSUMM, };
	
	struct Signs {
	private:
		BYTE	signs;
	public:
		static const BYTE Oblig		= 0x01;	// true if option is obligatory
		static const BYTE EscVal	= 0x02;	// true if opt value is escapable (can be missed)
		static const BYTE Allow0	= 0x04;	// true if allowed 0 while checking min value
		static const BYTE Hidden	= 0x08;	// hidden option (not printed in help)
		static const BYTE Trimmed	= 0x10;	// true if option has been processed
		static const BYTE Word		= 0x20;	// true if option is stated as multi-chars

		inline Signs(int x)	{ signs = (BYTE)x; }	// to initialize obligatory in main()
		// Returns true if given sign is set
		inline bool Is(BYTE mark) const	{ return (signs & mark) != 0; }	// '!= 0' to avoid warning
		// Sets given sign
		inline void MarkAs(BYTE mark)	{ signs |= mark; }

#ifdef DEBUG
		void Print() const {
			cout << int(Is(Word)) << int(Is(Trimmed)) << int(Is(Hidden)) << int(Is(Oblig));
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
		inline bool IsValEsc() const { return Sign.Is(Signs::EscVal); }

		// Sets option value.
		//	@opt: option
		//	@isword: true if option is a word, false if option is a char
		//	@val: value of option
		//	@isValLastToken: true if value is the last token
		//	@argInd: the current index in argc; increased by 1 if value is accepted
		//	return: 0 if success, -1 if not found, 1 if option or value is wrong
		int	SetVal(const char* opt, bool isword, char* val, bool isValLastToken, int* argInd);

		// Check option for obligatory.
		//	return: -1 if option is obligatory but not stated, otherwise 1
		int CheckOblig() const;

		// Prints option if it's obligatory
		inline void PrintOblig() const		{ if(Sign.Is(Signs::Oblig)) Print(false); }

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

		// Checks is string represents digital value.  Prints 'wrong val' message in case of failure
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
		//	@val: input value as C string
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
	//	@isValLastToken: true if value is the last token
	//	@argInd: the current index in argc; increased by 1 if value is accepted
	//	Return: 0 - success, 1 - false
	static int	SetOption(char* opt, char* val, bool isValLastToken, int *argInd);

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

	// Returns true if the option value is assigned by user
	inline static bool Assigned	(int opt)	{ return List[opt].Sign.Is(Signs::Trimmed); }
	// True if maximal enum value is setting
	inline static bool IsMaxEnum(int opt)	{ return List[opt].NVal == List[opt].MaxNVal - 1; }
	// Get maximal permissible numeric value by index
	inline static UINT GetMaxIVal(int opt)	{ return UINT(List[opt].MaxNVal); }
	// Returns C string 'ON' or 'OFF'
	inline static const char* BoolToStr(int opt)	{ return Booleans[GetBVal(opt)]; }
	inline static const char* BoolToStr(bool val)	{ return Booleans[val]; }
	// Gets pointer to the C string contained 'ON' or 'OFF'
	//inline static const char* BoolToStr(bool val)	{ return Booleans[int(val)]; }

	// Return string value by index: if value is not oblig and is not specified, then defName with given extention
	static const string GetSVal(int opt, const char* defName, const string& ext="_out.txt");

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
		F_BIGLINE,	// TxtFile(): too big line
		FZ_MEM,		// TxtFile(): not enough gzip buffer
		FZ_OPEN,	// TxtFile(): gzip open error
		FZ_BUILD,	// TxtFile(): doen't support gzip
		F_WRITE,	// file write error
#ifndef _FQSTATN
		TF_FIELD,	// TabFile: number of fields is less than expected
		TF_EMPTY,	// TabFile: none item (should be specified)
		B_BADEND,	// bed: start is is equal or more than end
		B_NEGPOS,	// bed: negative position
		ARR_OUTRANGE,
		SUM_EXCEED,
#endif
		EMPTY
	};

private:
	static const char* _msgs[];
	enum eCode	_code;			// error code
	char * _outText;			// output message
	
	// Initializes _outText by C string contained message kind of
	// "<sender>: <text> <specifyText>".
	void set_message(const char* sender, const char* text, const char* specifyText=NULL);

	//inline void set_message(const string& sender, const char* text, const char* specifyText=NULL) {
	//	set_message(sender==strEmpty ? NULL : sender.c_str(), text, specifyText);
	//}

public:
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

	// empty constructor for silent quit
	//inline Err() { set_message(NULL, strEmpty.c_str()); }

	// copy constructor
	Err(const Err& src);

	inline ~Err() { if( _outText) delete [] _outText; }

	inline const char* what() const /*throw()*/ { return _outText;	}

	//inline bool IsEmpty() const	{ return strlen(_outText)==0; }

	inline eCode Code() const		{ return _code; }

	// Returns point to the aditional text to specify exception/warning message or NULL
	//inline const char* SpecifyText() const	{ return _specifyText; }

	// Throws exception or outputs Err message.
	//	@throwExc: if true then throws exception, otherwise outputs Err message
	//	@eol: if true then carriage should be return while output Err message
	void	Throw(bool throwExc = true, bool eol = true);
	
	// Outputs warning with prefix "WARNING" and additional text, if it is setting.
	void	Warning	(string const & addText = strEmpty);
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
	//	otherwise outputs Err message as warning without EOL
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
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file doesn't exist
	inline static bool CheckFileExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFREG, throwExcept, Err::F_NONE);
	}

	// Checks if directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	inline static bool CheckDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR, throwExcept, Err::D_NONE);
	}

	// Checks if file or directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	inline static bool CheckFileDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR|S_IFREG, throwExcept, Err::FD_NONE);
	}

	// Throws exsception if file or directory doesn't exist
	//	@name: name of file or directory
	//	@ext: file extention; if set, check for file first
	//	@throwExcept: if true throws exception,
	//	otherwise outputs Err message as warning without EOL
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
	//	@isEOL: if true then ended output by EOL
	static void Print(long elapsed, const char *title, bool parentheses, bool isEOL);

	mutable time_t	_startTime;
	mutable bool	_enabled;	// True if local timing is enabled

	// Creates a new TimerBasic
	//	@enabled: if true then set according total timing enabling
	TimerBasic(bool enabled = true)	{ _enabled = enabled ? Enabled : false;	}
	
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
	//	@isEOL: if true then ended output by EOL
	static void StopCPU(bool isEOL=true) {
		if(Enabled)	Print((clock()-_StartCPUClock)/CLOCKS_PER_SEC, "CPU: ", false, isEOL);
	}

	// Creates a new Timer and starts it if timing is enabled
	//	@enabled: if true then set according total timing enabling
	inline Timer(bool enabled = true) : TimerBasic(enabled) { Start(); }
	
	// Stops enabled timer and prints elapsed time with title
	//	@title: string printed before time output
	//	@parentheses: if true then output time in parentheses
	//	@isEOL: if true then ended output by EOL
	void Stop(const char *title, bool parentheses, bool isEOL) {
		if(_enabled)	Print(GetElapsed(), title, parentheses, isEOL);
	}

	// Stops enabled timer and prints elapsed time
	//	@offset: space before time output
	//	@isEOL: if true then ended output by EOL
	void Stop(BYTE offset, bool isEOL = false);

	// Stops enabled timer and prints elapsed time
	//	@parentheses: if true then output time in parentheses
	//	@isEOL: if true then ended output by EOL
	inline void Stop(bool parentheses = false, bool isEOL = true)	{
		Stop(NULL, parentheses, isEOL); }
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
		inline StopwatchCPU() : _sumclock(0) {}

		inline void Start(bool reset=false)	{ _clock = clock(); if(reset) _sumclock = 0; }

		// Stops StopwatchCPU
		//	@title: string printed before time output
		//	@print: if true time should be printed
		//	@isEOL: if true then ended output by EOL
		void Stop(const char* title, bool print = false, bool isEOL = false);
};

#ifdef _MULTITHREAD

static class Mutex
{
private:
	static pthread_mutex_t	_mutexes[];
	static bool	_active;	// true if the mutex really should work
	static const BYTE Count = 6;
public:
	enum eType { OUTPUT, WR_FILE, INCR_SUM, OTHER1, OTHER2, OTHER3 };
	static void Init(bool);
	static void Finalize();
	static void Lock(const eType type);
	static void Unlock(const eType type);
} mutex;

class Thread
{
private:
	pthread_t _thread;

public:
	inline Thread(thrRetValType(
		#ifdef OS_Windows
		__stdcall
		#endif
		*proc)( void*), void *arglist)
	{
	#ifdef OS_Windows
		//_thread = (HANDLE) _beginthread(proc, 0, arglist);	// gets unstable call WaitForSingleObject
		_thread = (HANDLE) _beginthreadex(NULL, 0, proc, arglist, 0, NULL);
	#else
		pthread_create(&_thread, NULL, proc, arglist);
		//int code = pthread_create(&_thread, NULL, proc, arglist);
		//if (code) { 
		//	char buf[256]; 
		//	strerror_r(code, buf, sizeof buf);
		//	cout << "Thread constructor: " << buf << endl;
		//}
	#endif
	}

	inline ~Thread()
	{
	#ifdef OS_Windows
		CloseHandle(_thread);
	#endif
	}

	// Waits for a thread to terminate and detaches the thread.
	inline void WaitFor()
	{
#ifdef OS_Windows
		//DWORD res = WaitForSingleObject(_thread, INFINITE);
		//switch(res) {
		//	case WAIT_ABANDONED: cout << " WAIT_ABANDONED\n"; break;
		//	case WAIT_OBJECT_0:	cout << " WAIT_OBJECT_0\n"; break;
		//	case WAIT_TIMEOUT: cout << " WAIT_TIMEOUT\n"; break;
		//	case WAIT_FAILED: cout << " WAIT_FAILED\n"; break;
		//}
		WaitForSingleObject(_thread, INFINITE);
#else
		//int res = pthread_join(_thread, NULL);
		//switch(res) {
		//	case 0: cout << " WELL\n"; break;
		//	case EINVAL: cout << " EINVAL\n"; break;
		//	case ESRCH:	cout << " ESRCH\n"; break;
		//	case EDEADLK: cout << " EDEADLK\n"; break;
		//}
		pthread_join(_thread, NULL);
#endif
	}

// Terminates the calling thread. Does not used because of uninvoke destructors under Windows
//	void Thread::EndThread()
//	{
//#ifdef OS_Windows
//		_endthread();
//#else
//		pthread_exit(NULL);
//#endif
//	}
};

#endif	// _MULTITHREAD

template <typename T> class Array
{
protected:
	T*	_data;

private:
	ULONG	_len;

	//void Dump() { if(_len) { delete [] _data; _len = 0; } }

	// Prepaires memory array to initialization or copying
	 ULONG MakeData(ULONG len) {	
		 if(len) _data = new T[_len = len];
		 return len;
	 }

public:
	// Creates an instance with capacity, initialized by 0.
	//	@len: capacity
	Array(ULONG len=0) : _len(0), _data(NULL)	// inline is forbidden because of T& Chroms::AddEmptyElem()
	{ Reserve(len); }

	//inline Array(const Array& arr) { Copy(arr); }

	inline ~Array()		// only the data created in the constructor or MakeData() is freed
	{ if(_len) delete [] _data; _len = 0; }

	inline bool Empty() const	{ return !_len; }

	inline ULONG Length() const	{ return _len; }

	// Access data
	//	returns: a direct pointer to the memory array used internally
	inline T* Data() const	{ return _data; }

	inline T& operator[](const ULONG i)				{ return _data[i]; }

	inline const T& operator[](const ULONG i) const	{ return _data[i]; }

	inline Array& operator=(const Array& arr) { Copy(arr); return *this; }

	inline void Copy (const Array& arr)
	{ if(MakeData(arr._len))	copy(arr._data, arr._data + _len, _data); }

	// Deletes previous data and reserves array capacity, initialized by 0.
	//	@len: capacity
	inline void Reserve(long len)
	{ if(MakeData(len)) memset(_data, 0, _len*sizeof(T)); }	// no control of existing data

	// CLears an instance: set all members to zero
	inline void Clear() { if(_len)	memset(_data, 0, _len*sizeof(T)); }

#ifndef _FQSTATN
	//---------- extended methods

	// Adds array to this instance.
	//	@src: added array
	//	@start: position from which array should be added
	//	@startSrc: start position in added array; start of added array by default
	//	@size: length of added region; whole added array by default
	void Concat(const Array& src, ULONG start, ULONG startSrc=0, ULONG size=0) {
		if( !size )	size = src._len - startSrc;
		if( start + size > _len )
			Err(Err::ARR_OUTRANGE, "Array.Concat()", NNSTR(start + size, " > ", _len)).Throw();
		memcpy(_data + start, src._data + startSrc, size*sizeof(T));
	}

	// Fills subarray by value
	//	@begin: first position of subarray
	//	@end: last position of subarray
	//	@val: value to fill
	void Fill(long begin, long end, T val) {
		if( end > long(_len) )
			Err(Err::ARR_OUTRANGE, "Array.Fill()", NSTR(_len)).Throw();
		for(long i=begin; i<=end; i++)
			_data[i] = val;
	}
#endif	// _FQSTATN
#ifdef DEBUG
	void Print(long begin=0, long end=0) const {
		if( !end )	end = _len;
		for (long i=begin; i<end; i++)
			dout << i << TAB << _data[i] << EOL;
	}
#endif	// _DEBUG
};

// 'Chrom' establishes correspondence between chromosome's ID and it's name.
static class Chrom
/**********************************************************************************
TERMS.
MARK: human view of chromo's number:	1, ..., Y
ABBREVIATION CHROM's NAME:				chr1, ..., chrY
SHORT CHROM's NAME:						chrom 1, ..., chrom Y
TITLE CHROM's NAME:						chromosome X, ..., chromosome Y
PREFIX: substring before chrom's mark:	chr, chromosome.

CHROM's ID NUMBERING.
There are two ID numbering discipline: absolute and relative.

Relative discipline is applied when there is an extern list of chroms.
This list can be a reference  genome library or chrom.sizes file (explicit or embedded in the SAM format).
The autosomes' ID correspond to the chrom human number reduced by 1 (0-based numbering)
The heterosomes' ID start with the last autosomes' ID plus 1, in the next order: X, Y, M
This discipline is maintained by BamTools API.
It allows to control the presence in the genome of an user-specified ID.

Absolute discipline is applied when there is no extern list of chroms.
The autosomes' ID numering is the same as relative.
The heterosomes' ID is just a code of mark character.
It doesn't allow to control the presence in the genome of an user-specified ID.

Chrom class distinguishes these disciplines by the value of the 'firstHeteroID' variable:
a zero value indicates an absolute discipline, a non-zero value indicates a relative one.
***********************************************************************************/
{
public:
	static const char*	Abbr;				// Chromosome abbreviation
	static const string	Short;				// Chromosome shortening; do not convert to string in run-time
#ifndef _FQSTATN
	static const string	Title;				// Chromosome title; do not convert to string in run-time
	static const chrid	UnID = -1;			// Undefined ID (255)
	static const chrid	Count = 24;			// Count of chromosomes by default (for container reserving)
	static const BYTE	MaxMarkLength = 2;	// Maximal length of chrom's mark
	static const BYTE	MaxAbbrNameLength;	// Maximal length of abbreviation chrom's name
	static const BYTE	MaxShortNameLength;	// Maximal length of short chrom's name
	static const BYTE	MaxNamedPosLength;	// Maximal length of named chrom's position 'chrX:12345'

private:
	static		 BYTE	CustomOpt;	// user chrom option number; used to set custom cID
	static const char*	Marks;		// heterosome marks
	static const string	UndefName;	// string not to convert in run-time

	static chrid cID;				// user-defined chrom ID
	static chrid firstHeteroID;		// first heterosome (X,Y) ID

	// Gets heterosome ID by mark without control, or undefined ID
	static chrid HeteroID	(const char cMark);

	// Gets heterosome ID by case insensitive mark, or undefined ID
	inline static chrid CaseInsHeteroID(const char cMark) { return HeteroID( toupper(cMark) ); }

	// Gets chrom ID by case insensitive mark
	static chrid CaseInsID	(const char* cMark);

public:
	//*** common methods

	// Returns true if chrom is autosome, false for somatic
	inline static bool IsAutosome(chrid cid) { return cid < firstHeteroID; }

	// Gets the length of prefix, or -1 if name is not finded
	static short PrefixLength(const char* cName);

	// Sets relative ID numbering
	inline static void SetRelativeID() { if(!firstHeteroID)	firstHeteroID = 1; }
	
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

	// Checks chrom ID and returns this ID or undefined ID; used in BamInFile
	inline static chrid ValidateID(chrid cID) { return cID < firstHeteroID + 3 ? cID : UnID; }

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

#ifdef _FRAGDIST
	// Validates all chrom ID by SAM header data and sets custom ID
	static void Validate(const string& samHeader);
#endif

	//*** custom ID getters, setters

	// Gets custom chrom's ID
	inline static chrid CustomID()	{ return cID; }

	// Returns true if no custom chrom at all
	inline static bool NoCustom()	{ return cID == UnID; }

	// Returns true if no custom chrom chrom is set by user
	inline static bool IsCustom(chrid cid) { return cID == UnID || cID == cid; }

	// Sets custom chrom without control; used in bioCC
	inline static void SetCustomID(chrid cid) { cID = cid; }

	// Sets custom chrom ID with control
	//	exception: wrong chrom
	static void SetCustomID(bool prColon = false);

	// Sets number of 'custom chrom' progr option
	//	@absID: true if absolute discipline is applied
	static void SetCustomOption(int opt, bool absID = false);

	//*** work with name

	// Gets mark length by ID; used in WigMap (bioCC)
	inline static BYTE MarkLength(chrid cid) {
		return cID == UnID ? UndefName.length() : (cid > firstHeteroID || cid < 9 ? 1 : 2);
	}

	// Returns mark by ID
	static const string Mark(chrid cid);

	// Locate chrom's mark in string.
	//	@str: string checked for chrom number
	//	return: pointer to the chrom number in str, or a null pointer if Chrom::Abbr is not part of str
	static const char* FindMark(const char* str);

	// Gets chrom's abbreviation name by ID
	//	@numbSep: if true then separate chrom's number
	static string AbbrName(chrid cid, bool numbSep = false);

	// Gets short name 'chrom X'
	inline static string ShortName(chrid cid)	{ return Short + BLANK + Mark(cid); }

	// Gets title name 'chromosome X' or 'chromosomes'
	//	@cid: chromosome's ID or UnID if plural
	inline static string TitleName(chrid cid = UnID) { return Short + (cid==UnID ? "s" : Mark(cid)); }

	inline static const string Absent(chrid cid, const string & what) {
		return AbbrName(cid) + " is absent in " + what + " file: skipped";
	}

	//static char* LongToShortName(char* longName) ;
#endif	// _FQSTATN
} chrom;

#if !defined _WIGREG && !defined _FQSTATN
// 'Read' represents Read (with name and score in case of _VALIGN) as item
struct Read
{
#if defined _ISCHIP || defined _VALIGN
	static const char	NmNumbDelimiter;	// delimiter between prog title and number
	static const char	NmPos1Delimiter;	// delimiter before first position
	static const char	NmPos2Delimiter;	// delimiter between two positions in pair
	static const char	Strands[2];			// strand markers: [0] - positive, [1] - negative
#endif
	static readlen	Len;					// length of Read

#ifdef _ISCHIP

private:
	static char		SeqQuality;			// the quality values for the sequence (ASCII)
	static short	LimitN;				// maximal permitted number of 'N' in Read or vUNDEF if all
	static bool		PosInName;			// true if Read name includes a position
	static const char ToUp;				// shift beween lowercase and uppercase characters
	static const char Complements[];	// template for complementing Read

public:
	// Initializes static members
	//	@len: length of Read
	//	@posInName: true if Read position is included in Read name
	//	@seqQual: quality values for the sequence
	//	@limN: maximal permitted number of 'N'
	static void Init(readlen len, bool posInName, char seqQual, short limN);

	inline static bool IsPosInName	() { return PosInName; }

	// Fills external buffer by quality values for the sequence
	inline static void FillBySeqQual(char* dst) { memset(dst, SeqQuality, Len); }
	
	// Copies complemented Read.
	static void CopyComplement(char* dst, const char* src);

	// Checks Read for number of 'N'
	//	@read: checked Read
	//	return:	1: NULL Read; 0: success; -1: N limit is exceeded
	static int CheckNLimit(const char* read);

	// Prints quality values for the sequence
	static void PrintSeqQuality() { cout << '[' << SeqQuality << ']'; }

	// Prints Read values - parameters.
	static void Print();

#else
public:
	chrlen	Pos;		// Read's actual start position
	ULLONG	Numb;		// read number keeped in name
	bool	Strand;		// true if strand is positive

	chrlen Centre() const { return Pos + (Len>>1); }

#ifdef _VALIGN
	chrid	InitCID;	// initial chrom - owner
	float	Score;		// Read's score

	// Cobstructs extended Read
	//	@cid: chrom ID
	//	@num: uniq number within chrom or initial position
	Read(chrlen pos, chrid cid, ULLONG numb, char strand, float score)
		: Pos(pos), InitCID(cid), Numb(numb), Strand(strand), Score(score) {}
#elif defined _FRAGDIST

	Read(chrlen pos, ULLONG numb, bool strand) : Pos(pos), Numb(numb), Strand(strand) {}
#else
	inline Read(chrlen pos) : Pos(pos) {}
#endif
	// Compares two Reads by position. For sorting a container.
	inline static bool CompareByStartPos(const Read& r1, const Read& r2) { return r1.Pos < r2.Pos; }

	//inline static bool CompareByNum(const Read& r1, const Read& r2) {	return r1.Num < r2.Num; }

	void inline Print() const {
		cout << Pos  << TAB << Numb
#ifdef _VALIGN
			 << TAB << int(Score)
#endif
			 << EOL;
	}
#endif
};
#endif

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

