#pragma once

#include "def.h"
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
	typedef void*	retThreadValType;
	#define retThreadValTrue	(void*)1

	#define rand()	rand_r(&_seed)
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
	//#define SLASH '/'	// standard Linux path separator
	//#define F_READ_MODE O_RDONLY
#elif defined _WIN32
	#define OS_Windows
	#include <windows.h>
#ifdef _MULTITHREAD
	#include <process.h>	    // _beginthread, _endthread
	#define pthread_t HANDLE
	#define pthread_mutex_t CRITICAL_SECTION
#endif
	typedef unsigned __int64 ULLONG;
	typedef __int64 LLONG;
	typedef struct __stat64 struct_stat64;
	typedef UINT	retThreadValType;
	#define retThreadValTrue	1

	#define atol _atoi64
	#define isnan _isnan

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
typedef char threadnumb;	// type number of thread
typedef BYTE chrid;			// type number of chromosome
typedef USHORT readlen;		// type length of Read
typedef UINT chrlen;		// type length of chromosome
typedef ULONG genlen;		// type length of genome
typedef float readscr;		// type score of Read

#define CHRLEN_CAPAC	10	// capacity of max chrom length;
							// may be count by DigitsCount() every time,
							// but is the same if chrlen defined as int or long,
							// so is defined as static value
#define	CHRLEN_UNDEF	UINT_MAX	// undefined length of chromosome
#define retThreadValFalse	0

#define cN	'N'
#define HPH	'-'
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
#define SepSCl		"; "	// semicolon separator
#define SepClTab	":\t"	// colon + tab separator
#define Equel		" = "

static const string ZipFileExt = ".gz";
static const string strEmpty = "";

static const char* MsgDone = " done\t";
static const char* Notice = "NOTICE: ";
static const char* Total = "total";
static const char* Version = "version";
#ifdef _ISCHIP
#define SepGroup	";  "
#endif	// _ISCHIP
#if !defined _WIGREG && !defined _FQSTATN
static const char* Template = "template";
#endif

/*** COMMON MACROS & FUNCTION ***/

#define vUNDEF	-1	// undefined value

// Digital to STRing
// Returns value's string representation. http://rootdirectory.de/wiki/NSTR()
// Instead of std::to_string(x) [C++11] because of compatibility
#define NSTR( x ) static_cast<ostringstream & >( ostringstream() << dec << (x) ).str()

// Returns two int value's separated by delimiter string representation.
#define NNSTR( x, delim, y ) static_cast<ostringstream & >( ostringstream() << dec << int(x) << delim << int(y) ).str()

// Byte to STRing
#define BSTR(x) static_cast<ostringstream & >( ostringstream() << dec << int(x) ).str()

// Float to STRing
// Returns float's string representation.
//	@p: precision (number of digits after separator)
#define FSTR( x, p ) static_cast<ostringstream & >( ostringstream() << dec << setprecision(p)<< x ).str()

// Returns number of value's symbols, including minus & separator
#define SCNT( x ) NSTR(x).length()

// Returns number of float's symbols, including minus and dot
//	@p: precision (number of digits after separator)
#define SFCNT( x, p ) FSTR(x, p).length()

// Gets number of digist in a integral value
//	@val: integral value
//	return: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCount (LLONG val);

// Returns percent of @part relatively @total
inline float Percent(ULLONG part, ULLONG total) { 
	return total ? 100.f * part / total : 0.f;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
string	sPercent(float percent, BYTE precision=0, BYTE fieldWith=0, bool parentheses=true);

// Returns string represents the percent of part relatively total
//	@part: value represents desired %
//	@total: value represents 100%
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
inline string	sPercent(ULLONG part, ULLONG total,
	BYTE precision=0, BYTE fieldWith=0, bool parentheses=true) {
		return sPercent(Percent(part, total), precision, fieldWith, parentheses);
}

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

// Gets available system memory
//size_t	getAvailSystemMemory();

//#ifdef OS_Windows
//string	Wchar_tToString(const wchar_t* wchar);
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar);
//#endif

/*** end of COMMON FUNCTION ***/

#ifdef _ISCHIP
	#define dout	cout
#elif defined _WIGREG
	#define dout	cerr
#else	// _DENPRO || _BIOCC || _FQSTATN
// 'dostream' duplicates outstream to stdout & file
class dostream
{
    std::ostream &first, &second;
public:
    dostream(std::ostream &f, std::ostream &s) : first(f), second(s) {}
	template <typename T> dostream& operator<< (T val) {
		first << val;
		second << val;
		return *this;
	}
};

extern ofstream outfile;	// file ostream duplicated stdout; inizialised by file in code
extern dostream dout;		// stream's duplicator

#endif	// _DENPRO || _BIOCC || _FQSTATN

// 'Product' keeps Product info
static struct Product
{
	static const string	Title;
	static const string	Version;
	static const string	Descr;

	// Gets length of product name (title plus version).
	//static inline BYTE NameLength() { return BYTE(Title.length()); // + strlen(Version) + 1);
	//}

	// Gets title plus version.
	//static inline const string& Name() { return Title + string(1, HPH) + string(Version); }
} product;

// 'Options' implements main() options and parameters treatment.
static class Options
{
public:
	static const char* Booleans [];		// boolean values

	// Gets pointer to the C string contained 'ON' or 'OFF'
	//	@val: true of false
	static inline const char* GetBoolean(bool val) { return Options::Booleans[int(val)]; }

	// Gets pointer to the C string contained 'ON' or 'OFF'
	//	@i: option index
	static inline const char* GetBoolean(int i) { return Options::Booleans[Options::GetBVal(i)]; }

	// Prints 'usage' information
	//	@title: if true prints title before information
	//	return: 1 if title is settinf to true, 0 otherwise
	static int PrintUsage (bool title);
	
	// Returns command line.
	//	@argc: count of main() parameters
	//	@argv: array of main() parameters
	static const string CommandLine(int argc, char* argv[]);

	// Parses and checks main() parameters and their values.
	//	Output message if some of them is wrong.
	//	@argc: count of main() pearmeters
	//	@argv: array of main() pearmeters
	//	@obligPar: name of required application parameter or NULL if not required
	//	return: index of first parameter (not option) in argv[], argc if it is absent,
	//	negative if  tokenize complets wrong
	static int Tokenize(int argc, char* argv[], const char* obligPar=NULL);
	
	// Get double value by index
	static inline double GetDVal(int i)	{ return _Options[i].NVal; }
	// Get float value by index
	static inline float GetFVal	(int i)	{ return float(_Options[i].NVal); }
	// Get string value by index
	static inline char* GetSVal (int i)	{ return _Options[i].SVal; }
	//Get booling value by index
	static inline bool GetBVal	(int i)	{ return _Options[i].NVal != 0; }
	// Get UINT value by index
	static inline UINT GetIVal	(int i)	{ return UINT(_Options[i].NVal); }
	// True if maximal enunm value is setting
	static inline bool IsMaxEnum(int i)	{ 
		return _Options[i].NVal == _Options[i].MaxNVal - 1;
	}
	// Get maximal permissible numeric value by index
	static inline UINT GetMaxIVal(int i){ return UINT(_Options[i].MaxNVal); }
	// Get descriptor by index
	//static inline const char* GetDescr	(int i) { return _Options[i].Descr; }

private:
	// types of option values
	// do not forget to support the correlation with Options::_TypeNames []
	// tVERS option can be absent, so this 'special' group should be started by tHELP
	enum eValType	{ tUNDEF, tNAME, tCHAR, tINT, tFLOAT, tLONG, tENUM, tCOMB, tHELP, tVERS };
	
	struct Signs {
	private:
		BYTE	signs;

	public:
		static const BYTE Oblig		= 0x1;	// obligatory option sign
		static const BYTE Trimmed	= 0x2;	// trimmed option sign: true if option has already been processed
		//static const BYTE Printed	= 0x4;	// printed option mask

		inline Signs(int x)	{ signs = (BYTE)x; }	// to initialize obligatory in main()
		// Returns true if given sign is set
		bool inline Is(BYTE mark)		{ return (signs & mark) != 0; }	// '!= 0' to avoid warning
		// Sets given sign
		void inline MarkAs(BYTE mark)	{ signs |= mark; }
	};

	// structure 'Option' keeps options attributes.
	struct Option {
		const char	Char;		// option - character
		const char*	Str;		// option - string
			  Signs	Sign;		// initialize obligatory in main():
								// 1 if option is obligatory, 0 otherwise;
		//const bool	_ValRequired;// true if option's value is required
		const eValType ValType;	// type of value
		const BYTE	OptGroup;	// option's category
			  double NVal;		// default or established numeric or enum value
		const double MinNVal;	// minimal permissible value;
								// for enum should be a first defined value
		const double MaxNVal;	// maximal permissible value or count of enum values
			  char*	SVal;		// string value or pointer to enum values array;
								// in last case should be cast to char**.
								// If enum option hase ValRequired==false,
								// this pointer dousn't use and may be NULL
		const char*	Descr;		// tip string
		const char*	AddDescr;	// additional tip string

		static const char EnumDelims[];	// specifies delimiter for enum [0] and combi [1] values

		// returns true if option value is required
		bool inline ValRequired() { return MinNVal != vUNDEF; }

		// Sets option value.
		//	@isWord: true if option is a word, false if option is a char
		//	@opt: option without HYPHENs
		//	@val: value of option
		//	@isNextOpt: true if next parameter is option
		//	return: 0 if success, -1 if not found, 1 1 if option or value is wrong
		int	SetVal(bool isWord, char* opt, char* val, bool isNextOpt);

		// Check option for obligatory.
		//	return: -1 if option is obligatory but not stated, otherwise 1
		int CheckOblig();

		// Prints option in full or short way.
		//	@descr: if true, prints in full way: 
		//	signature, description (marks as Required if needed), default value, 
		//	otherwise signature only
		void Print(bool descr);

		void PrintCharDefValue();

		void (*PrintExtDefValue)();


	private:
		// Checks limits and set numerical value
		//	@val: numerical value
		//	@isWord: true if option enters by long name
		//	return: 1 if limits are exceeded, otherwise 0
		int SetTriedDigit(double val, bool isWord);

		// Checks and sets enum option value.
		//	@val: input value as C string
		//	return: true if success
		bool SetEnum(char* val);

		// Checks and sets enum option value.
		//	@val: input value as C string
		//	return: true if success
		bool SetComb(char* val);

		// Prints enum or combi values
		//	return: number of printed symbols
		BYTE PrintEnumVals();
		
		// Performs a case-insensitive search of given string value among enum values.
		//	@val: input value as string
		//	return: index of finded value in enum, or -1 if the value is not present in enum
		int GetEnumInd (char* val);
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
//public:
//	static const char*	TypeExt;
	static const char*	_TypeNames [];	// names of option value types in help
private:
	static const char*	_OptGroups [];	// names of option groups in help
	static	Option		_Options[];		// options. Option 'help' always shuld be the last one,
										// option 'version' - before last.
	static const Usage	_Usages[];		// content of 'Usage' variants in help
	static const BYTE	_OptCount,		// count of options
						_GroupCount,	// count of option groups in help
						_UsageCount;	// count of 'Usage' variants in help
	
	static inline string OptToStr	(bool isWord, const char* opt) {
		return string(1 + int(isWord), HPH) + opt;
	}
	
	// Check obligatory options and output message about first absent obligatory option.
	//	return: -1 if some of obligatory options does not exists, otherwise 1
	static int	CheckObligs();
	
	// Set value of option 
	//	@isWord: true if option is a word, false if option is a char
	//	@opt: option with HYPHENs
	//	@val: option's value
	//	@isNextOpt: true if next parameter is option
	//	@argIndex: the current index in argc; is increasing by 1 in case of required value
	//	Return: 0 if success, 1 otherwise
	static int	SetOption(bool isWord, char* opt, char* val, bool isNextOpt, int *argIndex);

	// Prints version
	//	return: always 1
	static int PrintVersion();

	// Ouptuts option with error message to cerr
	// Ouptuts option with error message to cerr
	//	@isWord: true if option is long
	//	@opt: option
	//	@val: value or NULL
	//	@msg: error message about value
	static int PrintWrongOpt(bool isWord, const char* opt, const char* val,
		const string msg = strEmpty);

	// Ouptuts ambiguous option with error message to cerr
	//	@isWord: true if option is long
	//	@opt: option
	//	@headMsg: message at the beginning
	//	@tailMsg: message at the end or NULL
	static int PrintAmbigOpt(bool isWord, const char* opt, const char* headMsg, const char* tailMsg);

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
		//TF_SPEC,
		TF_EMPTY,	// TabFile: none item (should be specified)
		BP_BADEND,	// bed: start is is equal or more than end
		BP_NEGPOS,	// bed: negative position
#ifdef _BEDR_EXT
		BR_RNAME,	// bed read: wrong read name format
#endif
		//FA_LONGLEN,
#endif
#if defined _ISCHIP || defined _FQSTATN
		FQ_HEADER,
		FQ_HEADER2,
#elif defined _DENPRO || defined _BIOCC
		ARR_OUTRANGE,
		SUM_EXCEED,
#endif
		EMPTY
	};

	//static const char* TREAT_BED_EXT;		// clarifying message in the bed stretch operation
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
	static const string IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName)
	{
		string res = fName;
		if(fName != strEmpty)	res += SepSCl;
		return res + issName + BLANK + NSTR(issNumb);
	}

	// Gets message "no @fileName.@fileExt[.gz] files in this directory"
	static const string MsgNoFiles (const string & fileName, const string fileExt)
	{
		return string("no " + fileName + fileExt +
			"[" + ZipFileExt + "] files in this directory");
	}

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

	inline const char* what() const { return _outText;	}
	//inline const char* what() const throw() { return _outText;	}

	//inline bool IsEmpty() const	{ return strlen(_outText)==0; }

	inline eCode Code() const		{ return _code; }

	// Returns point to the aditional text to specify exception/warning message or NULL
	//inline const char* SpecifyText() const	{ return _specifyText; }

	// Throws exception or outputs warning.
	//	@throwException: if true throws exception, otherwise outputs Err message
	//	@endOfLine: if true print EOL after Err message
	//	return: false
	void	Throw(bool throwException = true, bool endOfLine = true);
	
	// Outputs warning with prefix "WARNING" and additional text, if it is setting.
	void	Warning	(string const & addText = strEmpty);
};

// 'FileSystem' implements common file system routines
static class FS
{
private:
	// Returns true if file system's object exists
	//	@name: object's name
	//	@st_mode: object's system mode
	static bool IsExist(const char* name, int st_mode);

	// Checks if file system's object doesn't exist
	//	@name: object's name
	//	@st_mode: object's system mode
	//	@throwExcept: if true throws excwption,
	//	@ecode: error's code
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static bool CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode);

	// Searches through a file name for the any extention (slash|back-slash insensible).
	//	@fname: file name
	//	return: the index of the DOT mathed extention; otherwise npos
	static size_t GetExtPos	(const string &fname);

public:
	// Gets size of file or -1 if file doesn't exist
	static LLONG Size 	(const char*);

	// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
	static LLONG UncomressSize	(const char*);

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

	// Checks if file or directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static inline bool CheckFileDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR|S_IFREG, throwExcept, Err::FD_NONE);
	}

	// Throws exsception if file or directory doesn't exist
	//	@name: name of file or directory
	//	@ext: file extention; if set, check for file first
	//	@throwExcept: if true throws exception,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static bool CheckFileDirExist(const char* name, const string & ext, bool throwExcept) 
	{
		return HasExt(name, ext) ?
			CheckFileExist(name, throwExcept) :
			CheckFileDirExist(name, throwExcept);
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optsVal: Options char* value
	//	return: pointer to the checked file name
	static const char* CheckedFileDirName	(const char* name) {
		if( !IsFileDirExist(name) )	Err(Err::FD_NONE, name).Throw();
		return name;
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optVal: Options value
	//	return: pointer to the checked file name
	static inline const char* CheckedFileDirName	(int optVal) {
		return CheckedFileDirName(Options::GetSVal(optVal));
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@name: pointer to the file name
	//	return: pointer to the checked file name
	static const char* CheckedFileName	(const char* name) {
		if( name && !IsFileExist(name) )	Err(Err::F_NONE, name).Throw();
		return name;
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optVal: Options value
	//	return: pointer to the checked file name
	static inline const char* CheckedFileName	(int optVal) {
		return CheckedFileName(Options::GetSVal(optVal));
	}

	// Returns true if file has any extension.
	//	@fname: file name
	static inline bool HasExt	(const string &fname) { 
		return GetExtPos(fname) != string::npos;
	}

	// Returns true if file has a specified  extension.
	//	@fname: file name
	//	@ext: extension includes dot symbol and can be composite
	static bool HasExt	(const string& fname, const string& ext);

	// Returns true if file has '.gz' extension
	static inline bool HasGzipExt(const string& fname) { return HasExt(fname,ZipFileExt); }

	// Returns string containing real file extension (without zip extention).
	//	@fname: pointer to the file name
	//	return: string containing real file extension or empty string if no real extention
	static string const GetExt(const char* fname);

	// Returns file name without extentiom (slash|back-slash insensible)
	static string const FileNameWithoutExt (const string& fname);

	// Returns short file name by long one
	//	@fname: long file name
	static string const ShortFileName (const string& fname);

	// Returns directory name by long file name
	//	@fname: long file name
	//	@addSlash: true if slash sould be added at the end
	static string const DirName (const string& fname, bool addSlash = false);

	// Returns the name of last subdirectory by long file name
	//	@fname: long file name
	static string const LastSubDirName (const string& fname);

	// Returns the name ended by slash without checking the name
	static string const MakePath(const string& name);

#ifndef _WIGREG
	// Fills external vector of strings by file's names found in given directory
	// Implementation depends of OS.
	//	@files: external vector of strings that should be filled by file's names
	//	@dirName: name of directory
	//	@fileExt: file's extention as a choosing filter
	//	@all: true if all files with given extention should be placed into external vector,
	//	otherwise only one (any)
	//	return: true if files with given extention are found
	static bool GetFiles	(vector<string>& files, const string& dirName, const string& fileExt, bool all = true);
#endif	// _WIGREG

//	static inline void	Delete		(const char* fname) {
//#ifdef OS_Windows
//		DeleteFile(fname);
//#else
//		unlink(fname);
//#endif
//	}

} fs;

// 'TabFilePar' keeps basic parameters for TabFile
struct TabFilePar {
	const BYTE	MinFieldCnt;	// minimum number of feilds in file line
	const BYTE	MaxFieldCnt;	// maximum used number of feilds in data line
	const char Comment;			// char indicates that line is comment
	const char* LineSpec;		// substring on which each data line is beginning

	inline TabFilePar(BYTE	minTabCnt, BYTE maxTabCnt, char comm, const char* lSpec)
		: MinFieldCnt(minTabCnt), MaxFieldCnt(maxTabCnt), Comment(comm), LineSpec(lSpec) {}
};

// 'File Type' implements bioinformatics file type routines 
static class FT
{
private:
	struct fType {
		const char* Extens;			// file extension
		const string Item;			// item title
		const string ItemPl;		// item title in plural
		TabFilePar FileParam;		// TabFile parameters, defined feilds
	};
	static const fType Types[];
	static const BYTE	Count = 7;

	// Validates file format
	//	@fName: file name (with case insensitive extension and [.gz])
	//	@t: file type
	//	return: true if file extension correspondes to file type
	//static bool CheckType(const char* fName, eTypes t) { 
	//	return GetType(fName) == (t == ABED ? BED : t); 
	//}

public:
	enum eTypes {
		UNDEF,	// undefined type
		ABED,	// alignment bed
		BED,	// ordinary bed
		WIG,	// wiggle
		SAM,	// sam
		FQ,		// fastQ
		FA		// fasta
	};

	// Gets file format
	//	@fName: file name (with case insensitive extension)
	static eTypes GetType(const char* fName);

	// Validates file extension
	//	@fName: file name (with case insensitive extension and [.gz])
	//	@t: file type
	//	@printfName: true if file name should be ptinted
	//	@throwExc: true if throw exception, otherwise throw warning
	//	return: true if file extension correspondes to file type
	static bool CheckType(const char* fName, eTypes t, bool printfName=true, bool throwExc=true);

	// Gets file extension
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string Ext(eTypes t, bool isZip = false)	{ return Types[t].Extens; }
#ifdef _ISCHIP
	// Gets file extension, beginning at DOT and adding .gz if needed
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string RealExt(eTypes t, bool isZip = false);
#endif
	// Gets an item's title
	//	@t: file type
	//	@pl: true if plural form
	static const string& ItemTitle(eTypes t, bool pl = false)
	{ return pl ? Types[t].ItemPl : Types[t].Item; }

	// Gets TabFile params
	//	@t: file type
	static inline TabFilePar FileParams(eTypes t) { return Types[t].FileParam; }

} fformat;

class Timer
{
private:
	time_t	_startTime;
	bool	_enabled;	// True if local timing is enabled
	static clock_t	_StartCPUClock;

	// Prints elapsed time interval
	//	@title: string printed before time output
	//	@elapsed: elapsed time interval
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	static void PrintElapsed(const char *title, long elapsed, bool parentheses, bool isCarriageReturn);

public:
	// True if total timing is enabled
	static bool		Enabled;

	// Starts enabled CPU timer, if it is enabled
	static inline void StartCPU()	{ if( Enabled ) _StartCPUClock = clock(); }
	
	// Stops enabled CPU timer and print elapsed time
	//	@isCarrgReturn: if true then ended output by EOL
	static void StopCPU(bool isCarrgReturn=true);

	// Creates a new Timer and starts it if timing is enabled
	//	@enabled: if true then set according total timing enabling
	Timer(bool enabled = true)	{ _enabled = enabled ? Enabled : false; Start(); }
	
	// Restarts timer, if timing is enabled
	inline void Start()				{ if( _enabled ) time( &_startTime ); }

	// Stops enabled timer and print elapsed time with title
	//	@title: string printed before time output
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	void Stop(const char *title, bool parentheses, bool isCarrgReturn);
	
	// Stops enabled timer and print elapsed time
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	inline void Stop(bool parentheses = false, bool isCarrgReturn = true)	{
		Stop(NULL, parentheses, isCarrgReturn); }

	// True if instance timing is enabled
	inline bool IsEnabled() const { return _enabled; }
};

#ifdef _MULTITHREAD
static class Mutex
{
private:
	static const BYTE Count = 2;
	static pthread_mutex_t	_mutexes[];
public:
	enum eType { OUTPUT, WR_FILE };
	static void Init();
	static void Finalize();
	static void Lock(const eType type);
	static void Unlock(const eType type);
} mutex;

class Thread
{
private:
	pthread_t _thread;

public:
	Thread(retThreadValType(
		#ifdef OS_Windows
		__stdcall
		#endif
		*proc)( void*), void *arglist);

	inline ~Thread()
	{
	#ifdef OS_Windows
		CloseHandle(_thread);
	#endif
		//;
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
	Array(ULONG len=0) : _len(0), _data(NULL)	// inline is forbidden because of T& Chroms::AddEmptyClass()
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

#if defined _DENPRO || defined _BIOCC
	// Adds array to this instance.
	//	@src: added array
	//	@start: position from which array should be added
	//	@startSrc: start position in added array; start of added array by default
	//	@size: length of added region; whole added array by default
	void Concat(const Array& src, ULONG start, ULONG startSrc=0, ULONG size=0) {
		if( !size )	size = src._len - startSrc;
		if( start + size > _len )
			Err(Err::ARR_OUTRANGE, "Array.Concat()", NSTR(start + size) + " > " + NSTR(_len)).Throw();
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
#endif	// _DENPRO, _BIOCC
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
/*
 * Terms:
 * Short chromosome's name: 1, ..., Y.
 * Abbreviation chromosome's name: chr1, ..., chrY
 * Title chromosome's name: chromosome X, ..., chromosome Y
 * Prefix: substring before short chromosome's name: chr, chromosome.
 * Chromosome's ID is an integer that corresponds to short name for digital names,
 * and is a letter's ASCII code for literal names (X, Y).
 */
{
public:
	static const char*	Abbr;				// Chromosome abbreviation
#ifndef _FQSTATN
	static const char*	Short;				// Chromosome shortening
	static const string	Title;				// Chromosome title
	static const char	M = 'M';
	static const chrid	UnID = 0;				// Undefined ID
	static const chrid	Count = 24;				// Count of chromosomes by default
	static const BYTE	MaxShortNameLength = 2;	// Maximal length of short chrom's name
	static const BYTE	MaxAbbrNameLength;		// Maximal length of abbreviation chrom's name
	static const BYTE	MaxNamedPosLength;	// Maximal length of named chrom's position 'chrX:12345'

private:
		static const char	X = 'X';
	static const char	Y = 'Y';
	static const char*	UndefName;

	static chrid _cID;	// user-defined chrom ID
public:

	// Gets chromosome's ID stated by user
	static inline chrid StatedID()	{ return _cID; }

	// Returns true if user stated all chromosomes
	static inline bool StatedAll()	{ return _cID == UnID; }

	// Sets chromosome's ID with validation.
	//	@cID: chromosome's ID
	//	return: true if the check has passed, otherwise do not set and print message to cerr
	static bool SetStatedID(chrid cID);

	// Sets chromosome's ID stated by user with validation.
	//	@cName: chromosome's name
	//	return: true if the check has passed, otherwise do not set and print message to cerr
	static inline bool SetStatedID(const char* cName) { return SetStatedID(ID(cName)); }

	// Gets chromosome's ID by name
	//	@cName: any name ('file.fa' etc) or short name by default
	//  @prefixLen: length of prefix or 0 by default
	static chrid ID(const char* cName, size_t prefixLen=0);

	//// Checks chrom abbreviation
	////	@cName: chrom abbreviation name
	////	return: true if name is wrong
	//inline static bool CheckAbbr(const char* cName) { 
	//	return strncmp(cName, Abbr, strlen(Abbr)) != 0;
	//}

	// Locate chrom number in string.
	//	@str: string checked for chrom number
	//	return: pointer to the chrom number in str,
	//	or a null pointer if Chrom::Abbr is not part of str.
	static const char* FindNumb(const char* str);

	// Gets chromosome's short name by ID.
	inline static string Name(chrid cID) {
		return cID != UnID ? (cID < M ? (BSTR(cID)) : string(1, cID)) : UndefName;
	}

	// Gets chromosome's abbreviation name 'chrX' by its ID.
	inline static string AbbrName(chrid cID)	{ return Abbr + Name(cID); }

	// Gets chromosome's shortening name 'chrom X'
	//	@cID: chromosome's ID
	inline static string ShortName(chrid cID)	{ return Short + Name(cID); }

	// Gets chromosome's title name 'chromosome X' or 'chromosomes'
	//	@cID: chromosome's ID or UnID if plural
	inline static string TitleName(chrid cID = UnID)	{ 
		return Title + (cID==UnID ? "s" : (sBLANK + Name(cID)));
	}

	// Gets chromosome's ID by abbreviation name (with prefix 'chr')
	inline static chrid IDbyAbbrName(const char* cAbbrName) { return ID(cAbbrName, strlen(Abbr)); }

	// Gets chromosome's ID by string short (without prefix) name
	inline static chrid ID(const string & cName){ return ID(cName.c_str(), 0); }

	// Gets chromosome's ID by long name
	inline static chrid IDbyLongName(const char* cLongName) { 
		return ID(cLongName, PrefixLength(cLongName)); }

	// Gets the length of short chromosome's name by ID
	static inline BYTE NameLength(chrid cID) {
		return cID == UnID ? strlen(UndefName) : (cID >= 10 && cID < M ? 2 : 1);
		//return cID != UnID ? (cID < M ? (cID < 10 ? 1 : 2) : 1) : strlen(UndefName);
	}

	//// Gets the length of abbreviation chromosome's name by ID
	//static inline BYTE AbbrNameLength(chrid cID) {
	//	return NameLength(cID) + strlen(Abbr);
	//}

	// Gets the length of prefix (substring before short chromosome's name)
	//	@cLongName: long chromosome's name
	//	return: length of substring before short chromosome's name
	//	or -1 if short name is not finded
	static short PrefixLength(const char* cLongName);

	inline static const string Absent(chrid cID, const string & what) {
		return AbbrName(cID) + " is absent in " + what + " file: skipped";
	}

	//static char* LongToShortName(char* longName) ;
#endif	// _FQSTATN
} chrom;

#if !defined _WIGREG && !defined _FQSTATN
// 'Read' represents Read (with name and score in case of _BEDR_EXT) as item
struct Read
{
#if defined _ISCHIP || defined _BEDR_EXT
	enum rNameType {	// defines types of Read's name
		nmUndef = 0,	// undefined name: never stated by user
		nmNumb	= 1,	// a unique number stated as name
		nmPos	= 2		// a start position stated as name
	};

	static char	SeqQuality;						// the quality values for the sequence (ASCII)
	static const char	NmPosDelimiter;			// delimiter between two positions in pair
	static const char*	NmNumbDelimiter;		// delimiter between chrom name and number
	static const BYTE	NmDelimiterShift = 2;	// shift to pass ':N' for nmPos type
	static const string	NmSuffMate1,			// suffix of Read's name on pair-end first mate
						NmSuffMate2;			// suffix of Read's name on pair-end second mate
	static BYTE	OutNameLength;			// Maximum length of Read name in output file
#endif	// _ISCHIP || _BEDR_EXT
	static	readlen	Len;						// length of Read
	static const char	Strand[2];				// strand markers: [0] - positive, [1] - negative

#ifdef _ISCHIP

private:
	static short	LimitN;				// maximal permitted number of 'N' in Read or vUNDEF if all
	static ULONG	Count;				// counter of total writed Reads
	static rNameType	NameType;		// type of name of Read in output files
	static const char ToUp;				// shift beween lowercase and uppercase characters
	static const char Complements[];	// template for complementing Read

public:
	static ULONG	MaxCount;	// up limit of writes Reads
	static const char*	NmDelimiter;	// delimiter between chrom name & value: ":N" or ":"

	static void Init(readlen rLen, rNameType name, char seqQual, 
		//BYTE mapQual, 
		short limN, ULONG maxCnt);

	// Gets true if start position is stated as name
	//static inline bool IsPositionName () { return NameType == nmPos; }
	// Gets true if Read name keeps its number
	static inline bool IsNameAsNumber () { return NameType == nmNumb; }

	// Gets the common part of Read name in output files':'
	static inline const string& Name () { return Product::Title; }

	// Increments counter of total writed Reads thread-safely.
	//	return: true if limit is exceeded.
	static inline bool IncrementCounter() {	return InterlockedIncrement(&Count) >= MaxCount; }

	// Copies complemented Read.
	static void CopyComplement(char* dst, const char* src);

	// Checks Read for number of 'N'
	//	@read: checked Read
	//	return: -1 if Read is NULL, 0 if N limit is exceeded, 1 if success
	static int CheckNLimit(const char* frag);

	// Prints Read values - parameters.
	static void Print();

#else
	chrlen	Pos;		// Read's actual start position

#ifdef _BEDR_EXT
	chrid	InitCID;	// initial chrom - owner
	size_t	Num;		// Read's number or initial start position
	readscr	Score;		// Read's score

	inline Read(chrlen pos, chrid cid, size_t num, readscr score)
		: Pos(pos), InitCID(cid), Num(num), Score(score) {}
#else
	inline Read(chrlen pos) : Pos(pos) {}

#endif
	// Compares two Reads by position. For sorting a container.
	static inline bool CompareByStartPos(const Read& r1, const Read& r2) { return r1.Pos < r2.Pos; }

	//static inline bool CompareByNum(const Read& r1, const Read& r2) {
	//	return r1.Num < r2.Num;
	//}

	void inline Print() const {
		cout << Pos 
#ifdef _BEDR_EXT
			 << TAB << Num << TAB << int(Score)
#endif
			 << EOL;
	}
#endif
};
#endif
