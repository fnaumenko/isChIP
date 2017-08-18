#pragma once
#include "bioTxtFile.h"

static const char* SignPar = "# ";	// Marker of output parameter
static const char* SignDbg = "## ";	// Marker of output debug info

//typedef enum {	// three states boolean 
//   False	= 0,
//   True		= 1,
//   Neutral	= 2
//} triad;

//#define RAND_STD			// rand() (Windows) or rand_r(int *seed) (Linux)
//#ifdef OS_Windows
	// Xorshift by George Marsaglia http://en.wikipedia.org/wiki/Xorshift
	#define RAND_XORSHIFT
//#else
	// Mersenne Twister by Agner Fog, 2008-11-16 http://www.agner.org/random/
	//#define RAND_MT
//#endif

using namespace std;

// Task modes
enum eMode { TEST, CONTROL, REGULAR	};

#define	TestMode	(Imitator::Mode==TEST)
#define	ControlMode (Imitator::Mode==CONTROL)
#define	RegularMode (Imitator::Mode==REGULAR)

#define NoAmplification	(Amplification::Coefficient==1)

class Average
/*
 * Class 'Average' encapsulates count of average.
 */
{
private:
	bool	_keep;		// true if statistic is accumulated
	ULONG	_count;
	ULLONG	_summator;
public:
	inline Average() : _summator(0), _count(0), _keep(true)	{}
	inline float Value()	const { return (float)_summator/_count; }
	inline ULLONG Sum()		const { return _summator; }
	inline ULONG  Count()	const { return _count; }
	void operator+=(UINT val);
};

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
  // Compilers supporting C99 or C++0x have stdint.h defining these integer types
  #include <stdint.h>
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
  // 16 bit systems use long int for 32 bit integer.
  typedef   signed long int int32_t;
  typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Older Microsoft compilers have their own definition
  typedef   signed __int32  int32_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int64  int64_t;
  typedef unsigned __int64 uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
  // This works with most compilers
  typedef signed int          int32_t;
  typedef unsigned int       uint32_t;
  typedef long long           int64_t;
  typedef unsigned long long uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif

class Random
/*
 * Class 'Random' encapsulates random number generator.
 */
{
#ifdef RAND_MT
// Version of Mersenne Twister: MT11213A or MT19937
#if 0	// Constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else   // Constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
#endif

private:
	static int Seed;

public:
	// Sets and returns seed
	//	@random: if true, random seed
	static	int SetSeed(bool random);

	Random();
	// Returns random integer within interval [1, max]
	int	Range(int max);
	
	// Returns random true ot false with probability 0.5
	bool Boolean();
	
	// Returns true with given likelihood
	//	@sample: probability of returning true; from 0.0. to 1.0
	bool RequestSample(float sample);
	
	// Normal distribution
	double Normal();

private:
#ifdef RAND_STD
	int _seed;
#elif defined(RAND_MT)
	int mti;                            // Index into mt
	uint32_t y, mt[MERS_N];             // State vector
	
	void Init0(int seed);               // Basic initialization procedure
	uint32_t BRandom();                 // Generates random bits
#elif defined(RAND_XORSHIFT)
	uint32_t x, y, z, w;

	// Generates 32 random bits
	uint32_t IRand();
#endif

	double normal_x2;  
	//double _V1, _V2, _S;
	short _phase;

	// Generates random double number in the interval 0 <= x < 1
	double	 DRand();
};

class LognormDistribution : public Random
/*
 * 'LognormDistribution' encapsulates lognormal generator and it's average counter.
 */
{
private:
	static Average _Average;
	static float _RelSigma;		// Sigma/LnFactor
	static float _RelMean;		// Mean/LnFactor + LnTerm
	static float _szSelSigma;	// size selection normal distribution sigma
	static int k;	

	const float _relSigma, _relMean;	/// non-static for hight throughput in different threads.

	// Random number distribution that produces floating-point values according to a lognormal distribution,
	// About 1.5 times faster then std::lognormal_distribution (<random>)
	inline fraglen Next()	{ 
		return fraglen(exp( Normal() * _relSigma + _relMean ));
	// canonical form:  exp( (Normal()*Sigma + Mean) / LnFactor + LnTerm )
	}	

public:
	//inline LognormDistribution() { LognormDistribution(false); }
	bool  _saveAverage;		// true if Average should be calculated
	inline LognormDistribution()
		: _relSigma(_RelSigma), _relMean(_RelMean), _saveAverage(false) {}
	//LognormDistribution(bool recordAverage) : _phase(0), _relSigma(_RelSigma), _relMean(_RelMean), _saveAverage(!recordAverage) {}

	inline fraglen NormalNext() { return fraglen(Normal() * _szSelSigma);	}

	// Random number distribution that produces floating-point values according to a lognormal distribution,
	// with or without output accumulation to calculate average
	// About 1.5 times faster then std::lognormal_distribution (<random>)
	fraglen NextWithAccum();

	// Initializes distribution values
	//	@mean: expectation of the based normal distribution
	//	@sigma: standart deviation of the normal distribution
	//	@lnFactor: power multiplication factor in lognormal distribution
	//	@lnTerm: power summand in lognormal distribution
	static inline void Init(float mean, float sigma, float lnFactor, float lnTerm, float szSelSigma)	{
		_RelSigma = sigma / lnFactor;
		_RelMean = mean / lnFactor + lnTerm;
		_szSelSigma = szSelSigma;
	}
	//static inline ULLONG Sum()	{ return _Average.Sum(); }
	static inline ULONG CallsCnt()	{ return _Average.Count(); }
	static inline float Mean()		{ return _Average.Value(); }
	static inline float SavedMean(ULONG savedCnt)	{ return (float)_Average.Sum()/savedCnt; }
};

class Amplification
/*
 * Class 'Amplification' implements Multiple Displacement Amplification (MDA).
 */
{
private:
	struct Fraction {
		short length;
		short shift;

		inline void Set(short len, short shft) { length = len; shift = shft; }
	};
	static Average _Average;// the average of Amplification

	Fraction*_fractions;	// array of Fraction structures
	short	_count;			// the real number of amplifications
	short	_fragLen;		// the length of fragment
	short	_initCnt;		// the initial number of selected fragments in _fractions
	short	_currCnt;		// the current number of added selected fragments in _fractions
	BYTE	_memFactor;		// the coefficient of memory allocation

	void	Split	(bool firstCall);
	inline UINT	FractionsSize() { return sizeof(struct Fraction)*_memFactor*Coefficient; }

public:
	LognormDistribution* random;

	bool	calcAverage;	// true if Average should be calculate
	
	// The coefficient (formal number) of amplifications
	static short Coefficient;
	
	// Returns arithmetic mean of number of amplifications
	static inline double SimpleMean ()	{ return _Average.Value(); }

	//inline Amplification() : _memFactor(3), calcAverage(false), _fractions(Coefficient > 1 ? new Fraction[_memFactor*Coefficient] : NULL) {}
	Amplification();
	inline ~Amplification() { if( _fractions ) delete [] _fractions; }
	// Prepare instance to the new generate cycle
	void  Reset (int fragLen);
	// Gets next fraction in a cycle. shift: shift of fraction; return: length of fraction or 0 if fractions are ended
	short GetFraction (short* shift);
};

class ChromsThreads
/*
 * Class 'ChromsThreads' is a container of chromosome's threads.
 * Chromosome's thread containes thread number,
 * summary treated sizes of chromosomes in thread and cromosomes IDs.
 * Chromosomes are distributing among threads possibly according equally runtime,
 * mining that runtime is in proportion to chroms treated length.
 * First thread with number 1 is the main.
 */
{
public:
	struct ChromsThread
	{
		threadnumb		Numb;		// thread number from 1
		ULONG			sumSize;	// summary treated sizes of chromosomes in thread
		vector<chrid>	chrIDs;		// croms ID container

		// Creates an empty instance: (without chromosomes) for sampling & regular cutting
		inline ChromsThread(bool writable) : Numb(threadnumb(writable)), sumSize(0) {}

		// Creates a 'real' instance for imitation
		//	@maxcCntInThread: reserved capacity of chroms container
		inline ChromsThread(chrid maxcCntInThread) : Numb(0), sumSize(0) { 
			chrIDs.reserve(maxcCntInThread);
		}

		// for sumSize ascending sorting by default
		inline bool operator < (const ChromsThread& thrPool) const {
			return (sumSize < thrPool.sumSize);
		}

		// Returns true if this thread is trial (without writing output files).
		inline bool IsTrial() const	{ return Numb == 0; }

		// Returns true if this thread is slave.
		inline bool IsSlave() const	{ return Numb > 1; }
	};
private:
	vector<ChromsThread>	_threads;
public:
	// Creates chroms distributions among threads possibly according equally runtime,
	// mining runtime is in proportion to chroms treated length
	ChromsThreads(threadnumb thrCnt, const ChromFiles& chrFiles);

	// Gets chromosomes IDs by index of thread
	ChromsThread& operator[](threadnumb thrInd) { return _threads[thrInd]; }

	void Print();
};

class AvrFragLengths
/*
 * Class 'AvrFragLengths' reads and writes average legths of fragment into plain text file.
 * File stores next info, each on separate line:
 *   Value on the first line with AmplCoeff == -1 is a size factor;
 *   Value on the second line with AmplCoeff == 0 is a common average;
 *   Value on the each next line is an average for given AmplCoeff.
 */
{
private:
	struct AvrFragLen {
		short AmplCoeff;
		float Value;

		inline AvrFragLen(short coeff, float average) : AmplCoeff(coeff), Value(average) {}
	};
	short	_cntAverages;			// number of features (Feature structure)
	bool	_isChanged;
	string	_fileName;
	vector<AvrFragLen> _avrLens;

public:
	// Creates a new instance of AvrFragLengths from given path;
	// a file name will be constructed.
	// If file does not exist, the instance is empty.
	//  Exception: Err.
	AvrFragLengths	(const string path);
	
	// Writes AvrFragLengths to file if it is changed.
	//  Exception: Err.
	~AvrFragLengths();
	
	// Gets size Factor and averages for given amplification Coefficient; commonAverage is always the same.
	//  Return: true if given amplification Coefficient is defined
	bool Get(short amplCoeff, float *sizeFactor, float *commonAverage, float *savedAverage) const;
	
	// Adds size Factor and savedAverage for given amplCoeff and averages commonAverage.
	void Add(short amplCoeff, float sizeFactor, float commonAverage, float savedAverage);

#ifdef DEBUG
	void Print() const;
#endif
	//friend TxtFile;
};

#define	GR_CNT 2	// count of grounds, or count of eGround elements

class Imitator
/*
 * Class 'Imitator' implements main algorithm of simulation.
 */
{
private:
	enum eGround {
		FG,	// foreground
		BG	// background
	};	
	class ChromCutter
	/*
	 * 'ChromCutter' encapsulates thread context and methods to cut chromosomes.
	 * Each thread represents new instance of class and its execution.
	 */
	{
	private:
		string	_chrName;				// abbr name of current chrom; needs for output files
		bool	_isTerminated;			// true if thread is cancelled by exception
		ULONG	_selReadsCnt[GR_CNT];	// local array of counts of all selected Reads
										// for current chromosome: [0] - fg, [1] - bg
		ULONG	_wrReadsCnt	[GR_CNT];	// local array of counts of writed Reads
										// for current chromosome: [0] - fg, [1] - bg
		const ChromFiles& _chrFiles;	// input genome library
		OutFile	*_partoFile;			// partial output file
		ChromsThreads::ChromsThread& _thread;
		Amplification _ampl;
		LognormDistribution _lnDist;

		// Creates instance
		//	@imitator: the owner
		//	@csThread: thread contained treated chromosomes
		//	@calcAverage: true if averages should be calculated
		ChromCutter(const Imitator* imitator, ChromsThreads::ChromsThread* csThread,
			bool calcAverage);
		
		~ChromCutter ();
		
		// Returns random true or false with probability set by index i
		inline bool RequestSample(BYTE i)	{
			return _lnDist.RequestSample( Imitator::Samples[i] );
		}
		// Returns random true or false with adjusted sample probability
		inline bool RequestAdjSample()		{
			return _lnDist.RequestSample( Imitator::AdjSample );
		}
		// Returns random fragment's length within interval [1, max]
		inline fraglen Range(fraglen max)	{ return fraglen(_lnDist.Range(max)); }
		// Outputs count and percent of writes Reads
		//  gr: fore/background
		//  title: string printed at first
		void OutputReadCnt(Imitator::eGround gr, const string title);
		
		// Outputs chromosome's name and treatment info
		//  @nts: current chromosome
		//	@timer: current timer to thread-saves time output 
		//	@exceedLimit: true if limit is exceeded
		void OutputChromInfo (const Nts& nts, Timer& timer, bool xceedLimit);

		// Clears all imitation counters
		void ClearCounters()	{
			fill(_selReadsCnt, _selReadsCnt+GR_CNT, 0);
			fill(_wrReadsCnt, _wrReadsCnt+GR_CNT, 0);
		}

		// Sets terminate's sign and output message
		void Terminate(const char*msg);
		
		// Treats chromosomes given for current thread
		//	@singleThread: true if single thread execution: just for print chrom name
		void Execute(bool singleThread);
		
		// Cuts chromosome 
		//	@nts: cutted chromosome
		//	@currPos: cutting start position
		//	@feature: current treated feature
		//	@fgInFeature: if true accept foreground keeps insinde feature
		//	return: 0 if success,
		//		1 if end chromosome is reached (continue treatment),
		//		-1 if limit is achieved (cancel treatment)
		int	CutChrom	(
			const Nts& nts,
			chrlen* const currPos,
			const Featr& feature,
			bool fgInFeature
		);

		// Adds read(s) to output file
		//	@nts: cutted chromosome
		//	@currPos: current cutting position
		//	@fragLen: length of current fragment
		//	@revrsSample: the probability of generating reversing Read: from -1 to 1
		//		if < 0 than increasing to the end of "left" fragment;
		//		if 0 than fifty-fifty;
		//		if > 0 than increasing to the end of "right" fragment;
		//	return: -1 if fragment is NULL, 0 if limitN is exceeded,
		//		1 if Read(s) is(are) added, 2 if output file is NULL
		int AddRead	(const Nts& nts, chrlen currPos, short fragLen, bool isReverse);
			//float revrsSample);

		friend class Imitator;
	};

	static ULONG	TotalWrReadsCnts[];		// total counts of writed Reads: [0] - fg, [1] - bg
	static ULLONG	TotalSelReadsCnts[];	// total counts of selected Reads: [0] - fg, [1] - bg
	static ULLONG	TotalSlaveWrReadsCnt;	// total count of all writed Reads in slaved threads
	static ULONG	CellsCnt;	// count of cells
	static float	AdjSample;	// Adjusted Sample to stay in limit
	static float	Samples[];	// User samples: [0] - fg, [1] - bg
	static readlen	FlatLen;	// Boundary flattening length
	static fraglen	FragLenMax;	// maximal length of selected fragments:
								// established by --frag-dev or
								// by SHRT_MAX if size filter is OFF
	static BYTE	Verb;			// verbose level
	static BYTE	DigitsCnt[];	// maximal counts of Reads digits: [0] - fg, [1] - bg
	static bool	TreatOutFtrs;	// true if out_of_features areas are treated
	static bool	LetN;			// true if 'N' nucleotides should be counted
	static bool	UniformScore;	// true if template features scores are ignored
	static bool	StrandAdmix;	// true if opposite strand should be admixed at the bound of the binding site
	static Imitator	*Imit;		// singletone instance: to call threads only
	static const BedF *Bed;		// template bed-file (test mode) or NULL (control mode)
	
	const ChromFiles& _chrFiles;// input genome library
	OutFile& _oFile;			// output file

	// Prints chromosome's name and treatment info
	//	@cID: chromosomes ID
	//	@isOutput: true if chromosomes name should be printed
	static void OutputChromName(chrid cID, bool isOutput=true) {
		if( isOutput && Verbose(V_RT) ) {
			cout << Chrom::TitleName(cID) << MSGSEP_TAB;
			fflush(stdout);
		}
	}
	static void	OutputReadCnt(eGround gr, const string title);
	// Sets maximal count of digits for given ground if val is maximal
	static void	SetMaxDigitCnt(eGround gr, ULLONG val);
	
	// Sets adjusted Samples and clear all counter and means.
	// Samples are needed for the control of BF&FG levels by percent (given by user),
	// and to prorate number of written reads for each chromosome depending on reads limit.
	void	SetSample	();

	bool	CutGenome	();
	
	bool	CutRegular	();

	inline retThreadValType CutChrom	(void* arg, bool singleThread)	{
		ChromCutter(this, (ChromsThreads::ChromsThread*)arg, false).Execute(singleThread);
		return retThreadValFalse; 
	}

	static inline retThreadValType 
		#ifdef OS_Windows
		__stdcall 
		#endif
		StatCutChrom(void* arg)		{ return Imit->CutChrom(arg, false); }
	
public:
	static eMode Mode;			// current task mode
	// minimal length of selected fragments: established by --frag-dev or Read::Len if size filter is OFF
	static fraglen FragLenMin;
	// true if total genome is treated.
	// Set to false in Test mode only if single chrom is defined and BG_ALL is false.
	static bool All;
	
	static inline bool	Verbose(BYTE level)	{ return Verb >= level; }

	// Initializes fragment lengths.
	static void	InitFragLen(fraglen fragLen, fraglen fragDev, bool isSizeSelect) {
		FragLenMin = Read::Len;
		if( isSizeSelect ) {
			if(fragLen > fragDev + Read::Len)
				FragLenMin = fragLen - fragDev;
			FragLenMax = fragLen + fragDev;
		}
		else
			FragLenMax = SHRT_MAX;
	}
	// Initializes static values.
	static void	Init(
		eMode mode,
		ULONG cellsCnt,
		UINT verb,
		bool allBg,
		bool letN,
		bool uniformScore,
		bool strandAdmix,
		readlen flatLen
	) {
		Mode = mode;
		CellsCnt = cellsCnt;
		Verb = verb;
		All = (mode == CONTROL) || allBg;
		LetN = letN;
		UniformScore = uniformScore;
		StrandAdmix = strandAdmix;
		FlatLen = flatLen;
	}

	// Creates singleton instance.
	//  @cFiles: list of chromosomes as fa-files
	//	@oFile: output files
	//	@templ: input template or NULL
	inline Imitator(const ChromFiles& cFiles, OutFile& oFile, BedF* templ)
		: _chrFiles(cFiles), _oFile(oFile)
	{
		Bed = templ;
		Imit = this;
	}

	// Executes task in current mode and write result to output fq file
	void Execute();
};


