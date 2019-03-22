/**********************************************************
Imitator.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
Last modified: 21.03.2019
Provides chip-seq imitation functionality
***********************************************************/
#pragma once
#include "OutTxtFile.h"
#include "effPartition.h"

using namespace std;

typedef	short	a_cycle;	// type amplification cycles
typedef	BYTE	a_coeff;	// type coefficient of amplification

// Task modes
enum eMode { TEST, CONTROL };

#define SignPar	"# "	// Marker of output parameter in Usage
#define	TestMode	(Imitator::TMode==TEST)
#define	ControlMode (Imitator::TMode==CONTROL)

// 'Average' encapsulates average counting
class Average
{
private:
	bool	_keep;		// true if statistic is accumulated
	ULONG	_count;
	ULLONG	_summator;
public:
	inline Average() : _summator(0), _count(0), _keep(true)	{}
	// Gets arithmetic average.
	inline float Mean()	 const { return (float)_summator/_count; }
	inline ULLONG Sum()	 const { return _summator; }
	inline ULONG Count() const { return _count; }
	// Returns parameter val
	UINT operator+=(UINT val);
};

//static struct DistrParams
//{
//	static float lnMean;	// lognorm: Mean/LnFactor + LnTerm
//	static float lnSigma;	// lognorm: Sigma/LnFactor
//	static float ssMean;	// mean of size selection normal distribution
//	static float ssSigma;	// sigma of size selection normal distribution
//
//	// Initializes lognorm fragment and size selection distribution values;
//	//	@mean: expectation of frag normal distribution
//	//	@sigma: standard deviation of frag normal distribution
//	//	@ssMean: expectation of size sel normal distribution
//	//	@ssSigma: standard deviation of size sel normal distribution
//	//	@ssActive: true if size selection mode is ON
//	static void Init(float lnM, float lnS, float ssM, float ssS, bool ssActive) {
//		lnMean = lnM;  lnSigma = lnS;  ssMean = ssActive ? ssM : 0; ssSigma = ssS;
//	}
//
//	inline bool IsSS() { return ssMean; }
//} distrParams;
//

// 'Imitator' implements main algorithm of simulation.
class Imitator
{
	// 'FragLenStat' keeps statistics of selected frag length; used in SetSample()
	struct FragLenStat
	{
		fraglen Min;	// minimum selected frag length
		fraglen	Max;	// maximum selected frag length
		Average	SelAvr;	// average selected frag length 

		inline FragLenStat() : Min(FRAG_MAX), Max(0) {}

		// Accepts frag length in purpose of statistics
		//	@len: current selected frag length
		void TakeFragLen(fraglen len) {
			SelAvr += len;
			if(len > Max)		Max = len;
			else if(len < Min)	Min = len;
		}
	};

	// 'FragCnt' - Fragment's Counter - keeps statistics for selected and recorded Reads
	struct FragCnt {
	private:
		typedef ULLONG	(FragCnt::*pRecIncr)(bool);
		typedef void	(FragCnt::*pSelAdd)	(ULLONG);

		// pointer to the thread-saved 'recorded frag's number increment' method
		static pRecIncr	pRecIncrSaved;
		// pointer to the thread-saved 'selected frag's number adding' method
		static pSelAdd	pSelAddSaved;

		ULLONG sel;		// number of selected fragments
		ULLONG rec[2];	// number of recorded fragments: 0: derived (amplified), 1: primer (initial)

		// Thread-savely increments number of recorded frags
		//	@primer: if true then increment counter of primer (not amplified) frag
		//	return: number of recorded frags after increment
		inline ULLONG RecInterlIncr(bool primer) { return InterlockedIncrement(&rec[primer]); }

		// Thread-savely adds value to number of selected frags
		inline void SelAddIncr(ULLONG val) { InterlockedExchangeAdd(&sel, val);	}

		// Adds value to number of selected frags
		inline void SelAdd(ULLONG val) { sel += val; }

	public:
		inline static void Init(bool singleThread) {
			pRecIncrSaved = singleThread ? &Imitator::FragCnt::RecIncr : &Imitator::FragCnt::RecInterlIncr;
			pSelAddSaved =	singleThread ? &Imitator::FragCnt::SelAdd : &Imitator::FragCnt::SelAddIncr;
		}

		inline void Clear() { memset(rec, 0, 2*sizeof(ULLONG));	}

		// Gets total number of selected frags
		inline ULLONG SelCnt()	const { return sel; }

		// Gets total number of recorded frags
		inline ULLONG RecCnt()	const { return rec[0] + rec[1]; }

		// Gets percentage of recorded reads in relation to the selected
		inline float Sample()	const { return Percent(rec[1], sel); }
		
		// Gets ratio primer/derived frag's numbers
		inline float RealAmplCoeff() const { return rec[1] ? float(rec[0]) / rec[1] : 0; }

		// Increments number of recorded frags
		//	@primer: if true then increment counter of primer (not amplified) frag
		//	return: number of recorded frags after increment
		inline void SelIncr() { sel++; }

		// Increments number of recorded frags
		//	@primer: if true then increment counter of primer (not amplified) frag
		//	return: number of recorded frags after increment
		inline ULLONG RecIncr(bool primer) { return rec[primer]++; }

		// Thread-savely increments number of recorded frags
		//	@primer: true if increment counter of primer (not amplified) frag
		//	return: number recorded frags after increment
		inline ULLONG RecIncrSaved(bool primer) {
			return (this->*pRecIncrSaved)(primer), RecCnt();
		}

		// Thread-savely adds value to number of selected frags
		inline void SelAddSaved(ULLONG val) { (this->*pSelAddSaved)(val); }
	};

	// 'FragCnts' keeps statistics for selected and recorded Reads, for both Test and MakeControl modes
	struct FragCnts {
	private:
		FragCnt	fCnts[2][Gr::Cnt];
		GM::Mode	gMode;

	public:
		inline void SetGMode(GM::Mode mode) { gMode = mode; }

		inline const FragCnt* GetFragCnts() const { return fCnts[gMode]; }

		inline FragCnt& operator[](BYTE g)  { return fCnts[gMode][g]; }

		inline const FragCnt& operator[](BYTE g) const { return fCnts[gMode][g]; }

		// Clears all fragment counters
		inline void Clear() { memset(fCnts, 0, 2*Gr::Cnt*sizeof(FragCnt)); }
	};

	// 'ChromCutter' encapsulates thread context and methods to cut chromosomes.
	class ChromCutter
	{
	public:
		// 'FragDistr' encapsulates normal and lognormal random generator and it's average counter.
		class FragDistr : public Random
		{
		private:
			//typedef fraglen	(FragDistr::*tpLognormNext)();

			//static tpLognormNext pLognormNext;
			//static float ssDVar;	// Doubled variance 2*sigma^2 in the size sel normal distribution
			//static float ssFactor;	// Factor sigma*sqr(2PI) in the size sel normal distribution
			//static float ssAlRatio; // Size sel normal distribution aligning up ratio
			static float ssFactor0;		// factor sigma*sqrt(2) in the size sel norm distr
			static float ssFactor1;		// factor 2.5/sqrt(2PI) in the size sel norm distr

			double	_normal_x2;		// second random coordinate (for normal RNG)
			short	_phase;			// phase (for normal RNG)
			Average*_avr;			// average to count the frags (for statistics), or NULL

			// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
			//	return:  value with Gaussian likelihood between about -5 and +5
			//	(from -6 to 6 in 1000000000 cycles) 
			double Normal();

			// Random number distrib that produces values according to a lognormal distrib.
			// canonical form:  exp( (Normal()*Sigma + Mean) / LnFactor + LnTerm )
			// About 1.5 times faster then std::lognormal_distribution (<random>)
			inline fraglen LognormNextDirect() {
				return fraglen(exp(Normal() * DistrParams::lnSigma + DistrParams::lnMean));
			}

			//inline fraglen LognormNextWithStat() { return *_avr += LognormNextDirect();	}

		public:
			//	Initializes size sel factors
			static void Init();

			// Creates instance
			//	@avr: generated frags average (for trial mode), or NULL (for working mode)
			// _normal_x2 is initialized by random double to avoid
			//		undesirable 'out of range' random initialization
			FragDistr(Average* avr) : _avr(avr), _normal_x2(DRand()) {
				//pLognormNext = avr ? 
				//	&Imitator::ChromCutter::FragDistr::LognormNextWithStat:
				//	&Imitator::ChromCutter::FragDistr::LognormNextDirect;
			}
			
			// Returns next random frag length according to a lognormal distribution,
			// with or without output accumulation to calculate average
			// About 1.5 times faster then std::lognormal_distribution (<random>)
			inline fraglen LognormNext() {
				return _avr ? *_avr += LognormNextDirect() : LognormNextDirect();
				//return (this->*pLognormNext)(); 
			}

			//// Returns next random size selection length according to a normal distribution
			//inline fraglen SzSelNext() { 
			//	return fraglen(Normal() * DistrParams::ssSigma + DistrParams::ssMean);
			//}

			// Returns next random size sel limits
			//	@min: random min limit
			//	@max: random max limit
			inline void SizeSelLimNext(fraglen& min, fraglen& max)
			{
				if(DistrParams::IsSS()) {
					float ssDev = ssFactor0 * (float)sqrt(log(ssFactor1 / DRand()));
					min = fraglen(DistrParams::ssMean - ssDev);
					if(min < Read::Len)		min = Read::Len;
					// below doesn't compiled by gcc
					//min = max(fraglen(DistrParams::ssMean - ssDev), fraglen(Read::Len));
					max = fraglen(DistrParams::ssMean + ssDev);
				}
			}

			// Returns true if size selection is failure
			//bool SzSelFailure(fraglen fLen) {
			//	return ssON ? 
			//	DRand() > ssAlRatio * exp(-pow(fLen - ssMean, 2) / ssDVar) / ssFactor : false;
			//}

			// Returns random fragment's length within 0 and mean frag length
			// should to round ssMean??
			inline fraglen RandFragLen() { return Range(int(DistrParams::ssMean)); }
		};

	private:
		// 'MDA' implements Multiple Displacement Amplification (MDA).
		class MDA
		{
		private:
			struct Fraction {
				fraglen Length;
				fraglen Shift;

				// Saves fraction
				//	@frac: fraction's length
				//	@shift: fraction's shift
				inline void Set(fraglen frac, fraglen shift) { Length = frac; Shift = shift; }
			};

			Fraction*_fracts;	// array of Fraction structures
			a_cycle	_iterCnt;	// number of iteration, i.e. the real number of amplifications
			fraglen	_fragLen;	// length of fragment
			fraglen	_initCnt;	// initial number of selected fragments in _fracts
			fraglen	_currCnt;	// current number of added selected fragments in _fracts
			BYTE	_memFactor;	// coeff of memory allocation, i.e. ratio possible_needed/allocated memory
			Random& _rng;		// used to invoke Range() only

			// Recursively fills array _fracts by splitted fragments (segments).
			//  @firstCall: true if it is a first call
			//	@fragLenMin: current minimal fragment length
			void	Split	(bool firstCall, fraglen fragLenMin);

			inline UINT	FractionsSize() { return sizeof(struct Fraction)*_memFactor*Coeff; }

		public:
			// The coefficient (formal number) of amplifications
			static a_cycle Coeff;
	
			inline MDA(Random& rng) : _memFactor(3), _rng(rng) {
				_fracts = Coeff > 1 ? new Fraction[_memFactor*Coeff] : NULL;
			}

			inline ~MDA() { if( _fracts ) delete [] _fracts; }

			// Prepare instance to the new generate cycle
			//	@fragLen: current fragment length
			//	@fragLenMin: current minimal fragment length
			void  Reset (fraglen fragLen, fraglen fragLenMin);
	
			// Gets next fraction in a cycle.
			//	@shift: shift of fraction
			//	return: length of fraction exceeded min frag length, established in Reset(),
			//	or 0 if fractions are ended
			fraglen GetFraction (fraglen* shift);
		};

		static fraglen	_SsDev;			// deviation of frag size selection
		static a_cycle	_PCRdcycles;	// PCR cycles: read doubling cycles
		static a_coeff	_amplCoeff;		// user-stated amplification coefficient
		
		bool		_slave;		// if true then this instance is slave
		GM::Mode	_gMode;		// generating mode: 0 - Test, 1 - Control
		OutFiles*	_output;	// partial output files
		FragCnts	_fragCnt;	// numbers of selected/recorded fragments for FG & BG, for both Teat & Input
		FragDistr	_fragDistr;	// normal & lognormal random number generator
		MDA			_ampl;
		const ChromFiles& _cFiles;	// reference genome

		// Sets global mode
		ULONG SetGMode(GM::Mode gmode);

		// Increments counters of local and total recorded fragments thread-safely
		//	@g: ground
		//	@primer: true if increment derived (amplified) frag's counter
		//	return: true if Reads limit is exceeded.
		bool IncrRecFragCount(Gr::Type g, bool primer);

		// Increments counter of total selected fragments thread-safely
		void IncrTotalSelFragCount();

		// Returns random true or false per ground sample
		inline bool PerSample(Gr::Type g) { return _fragDistr.Sample(Sample(_gMode, g)); }

		// Returns random true or false per auto sample
		inline bool PerAutoSample()	{ return _fragDistr.Sample(AutoSample); }

		// Prints thread-safe info about treated chroms and stops timer
		//	@cID: chrom ID
		//  @nts: current chromosome
		//	@enRegLen: length of all enriched regions
		//	@timer: current timer to thread-saves time output 
		//	@excLimit: true if limit is exceeded
		void PrintChrom (chrid cID, const Nts& nts, chrlen enRegLen, Timer& timer, bool excLimit);

		// Sets terminate's output message
		//	@tID: thread ID
		//	@msg: message
		void Terminate(thrid tID, const char*msg) {
			cerr << sThread << int(tID) << SepCl << msg << endl;
		}

	public:
		// Returns amplification coefficient
		static inline a_coeff AmplCoeff()	{ return _amplCoeff; }

		// Returns true if MDA amplification is established
		static inline bool IsMDA()	{ return MDA::Coeff; }

		// Returns true if PCR amplification is established
		static inline bool IsPCR()	{ return _PCRdcycles > 1; }

		// Set amplification
		//	@isPCR: true if PCR set
		//	@coeff: amplification coefficient; 0 if no amplification
		static void SetAmpl(bool isPCR, a_coeff coeff);

		// Creates instance
		//	@imitator: the owner
		//	@avr: recorded frags average (trial for sampling), or NULL
		//	@slave: if true then this instance is slave
		ChromCutter(const Imitator* imitator, Average* avr, bool slave);
		
		inline ~ChromCutter () { if(_slave)	delete _output; }

		// Returns number of FG recorded frags
		inline ULLONG RecFgFragCnt() const { return _fragCnt[Gr::FG].RecCnt(); }
		
		// Treats chromosomes given for current thread
		//	@cSubset: pointer to ChrSubset - set of chrom IDs treated in this thread
		thrRetValType Execute(const effPartition::Subset& cSubset);
		
		// Cuts chromosome until reaching end position of current treated feature
		//	@nts: cutted chromosome
		//	@currPos: cutting start position
		//	@feature: current treated feature
		//	@bg: if true then generate background (swap foreground and background)
		//	@fragStat: statistics of selected fragments average counter, or NULL under work mode
		//	return: 0 if success,
		//		1 if end chromosome is reached (continue treatment),
		//		-1 if limit is achieved (cancel treatment)
		int	CutChrom (const Nts& nts, chrlen* const currPos, const Featr& feature,
			bool bg, FragLenStat* fragStat = NULL);
	};

	// 'ChromView' provides template and methods for viewing chrom's treatment results
	struct ChromView
	{
	/*
	* output view:
	* A: always aligned                    A                                A
	* chrom       reads  sample ampl r/kbp|        reads sample ampl  r/kbp|   N   N_excl  mm:ss
	* ------------------------------------|--------------------------------|--------------------
	* chr 17: FG: 111111  100%   1.4  11.3|   BG: 222222   10%  0.622  0.59|  3.5%  3.1%   00:03
	* -------^---^-------^-----^------^------^----^-------^-----^------^-----^------^------^-----
	* chrNam      CountW  SmpW  AmplW  DnsW        CountW  SmpW  AmplW  DnsW   time
    * margN: 1   2       3     4      5      6    2       3     4      5     7      8      9
	*/
	private:
		static const BYTE margC_	= 1;	// length of margin 1 after chrom name
		static const BYTE marg_R	= 1;	// length of margin 2 before Reads count
		static const BYTE margR_S	= 3;	// length of margin 3 between Reads count and sample
		static const BYTE margS_	= 3;	// length of margin 4 after sample
		static const BYTE marg_D	= 3;	// length of margin 5 before density
		static const BYTE margD_	= 3;	// length of margin 6 after density
		//static const BYTE margD_N	= 3;	// length of margin 7 between density and N
		static const BYTE margN_Nex	= 2;	// length of margin 8 between N and N excleded
		static const BYTE marg_T	= 2;	// length of margin 9 before time
		static const BYTE AmplW		= 4;	// width of ampl coefficient field
		static const BYTE NW		= 4;	// length of 'N' field
		static const BYTE AmplPr	= 1;	// precision of ampl coefficient field
		static const BYTE SamplePr	= 2;	// precision of sample field

		// Prints empty margin
		//	@width: margin width
		static inline void PrintMarg(BYTE width) { cout << setw(width) << BLANK; }

		BYTE	CountW;		// width of 'number of Reads' field
		BYTE	SampleW;	// width of 'sample' field
		BYTE	DensW;		// width of 'density' field
		BYTE	DensPr;		// precision of density
		float	DensMax;	// maximal density
		Gr::Type GrType;	// FG|BG

	public:
		static const char* sMargTime;	// margin string before time output

		// Gets the maximum length of chrom name field with blank after
		static const BYTE ChromNameW() { 
			// 2 = 1 for blank inside chrom name + 1 for COLON
			return Chrom::MaxAbbrNameLength + 2 + margC_;
		}

		// returns ground title width or 0 if not TEST mode
		static BYTE GrTitleW() { return TestMode * (Gr::TitleLength + strlen(SepCl)); }

		inline ChromView(Gr::Type grType) : GrType(grType) {}

		// Prints chroms 'N' statistics
		static void PrintN(const Nts& nts);

		// Prints chrom's header: Reads
		//	@FgHeader: true for foreground header
		//	return: header width
		int PrintHeader(bool FgHeader);

		// Prints chrom's header: 'N' statistics
		static int PrintHeaderN();

		// Print info about recorded Reads (count, sample, ampl, density) for ground type of this nstance
		//  @fragCnt: Reads counters
		//	@gMode: gen mode for which info is printed
		//	@densLen: length on which the density is determined
		void PrintReads(const FragCnt fragCnt, GM::Mode gMode, chrlen densLen);

		// Set task mode
		//static void SetMode() { margC_R -= (gMode = BYTE(TestMode)); }

		// Initializes viewing data
		void Init(ULONG maxCnt, float sample, float maxDens);
#ifdef DEBUG
		void Print(ULONG maxCnt);
#endif
	};

	struct Context {
		UINT	CellCnt;			// count of cells
		float	Sample[Gr::Cnt];	// user-defined samples: [0] - fg, [1] - bg
		FragCnt	fCnts[Gr::Cnt];		// fragment counters: [0] - fg, [1] - bg

		// Gets total number of recorded frags
		inline ULLONG RecCnt()	const{ return fCnts[Gr::FG].RecCnt() + fCnts[Gr::BG].RecCnt(); }

		// Sets sample for both grounds
		inline void SetSample(float sample) { Sample[Gr::FG] = Sample[Gr::BG] = sample; }

		inline void ClearFragCounters() { fCnts[Gr::FG].Clear(); fCnts[Gr::BG].Clear(); }
		
		// Returns exact number of 'background' cells according to background sample
		inline float GetExactBGCellCnt() const { return Sample[Gr::BG] * CellCnt; }

		// Sets Control number of cells and sample according to exact count of Test 'background' cells
		//	@exactCellCnt: exact count of Test 'background' cells
		void SetControlSample(float exactCellCnt) {
			Sample[Gr::BG] = exactCellCnt/(CellCnt = UINT(ceil(exactCellCnt)));
		}

		// Thread-safely increments counters of local and total recorded fragments
		//	@g: ground
		//	@primer: true if increment amplified (second) frag's counter
		//	return: true if Reads limit is exceeded.
		bool IncrRecFragCount(Gr::Type g, bool primer) {
			return fCnts[g].RecIncrSaved(primer) + fCnts[!g].RecCnt()
				>= Seq::FragsLimit();
		}
	};

	static Context	GlobContext[];	// global generation context
	static ChromView ChrView[];	// chrom views per ground
	static ULONG	TotalLen[];	// FG, BF genome length; needed to define total density
	static float	AutoSample;	// adjusted FG sample to stay in limit
	static float	SelFragAvr;	// mean length of selected fragments
	static BYTE	Verb;			// verbose level
	static bool	MakeControl;	// true if control file (input) should be produced
	static bool	UniformScore;	// true if template features scores are ignored
	static Imitator	*Imit;		// singletone instance: to call threads only
	static const BedF *Bed;		// template bed-file (test mode) or NULL (control mode)
	//static readlen	BindLen;	// binding length

	const ChromFiles& _cFiles;	// ref genome library
	OutFiles& _oFile;			// output file

	// Returns stated count of cells
	//	@gm: generated mode
	inline static UINT CellCnt(GM::Mode gm) { return GlobContext[gm].CellCnt; }

	// Returns stated sample
	//	@gm: generated mode
	//	@g: ground 
	inline static float Sample(GM::Mode gm, Gr::Type g) { return GlobContext[gm].Sample[g]; }

	// Prints chromosome's name
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@gm: generation mode Test|Control
	//	@print: true if chromosomes name should be printed
	static void PrintChromName(chrid cID, GM::Mode gm, bool print);

	// Prints chromosome's name and full statistics
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@fCnts: array of fragment's counters, for FG and BG
	//	@regLens: array of region's lengths to print FG|BG density
	//	@printChrName: true if chromosome's name should be printed
	static void PrintChromInfo(chrid cID, GM::Mode gMode,
		const FragCnt* fCnts, const ULONG regLens[], bool printChrName=true);

	// Prints header (FG and BG, or single) and solid line
	//	@header: if false, then print solid line only
	static void PrintHeader(bool header);

	// Increments grounds total length.
	static void IncrementTotalLength(chrlen defRegLen, chrlen enrRegLen);
	
	//// Increments counter of total recorded fragments thread-safely.
	////	@g: ground
	////	@primer: if true then increment amplified (second) frag's counter
	////	return: true if limit is exceeded.
	//static bool IncrementRecFragCount(Gr::Type g, bool primer) {
	//	ULLONG cnt = ChromCutter::IsAmpl() ?
	//		TotalFragCnt[g].InterlockedIncrRec(primer):
	//		OutFiles::SingleThread ?
	//			TotalFragCnt[g].PrimerRec++:
	//			InterlockedIncrement(&(TotalFragCnt[g].PrimerRec));
	//	return cnt + TotalFragCnt[!g].Rec() >= Seq::FragsLimit();
	//}

	// Gets Reads estimated count per chrom and corrects maximum count and density
	//	@g: ground
	//	@densLen: length on which density is defined
	//	@factor: CellsCnt / recAvrLen
	//	@numeric: 1 for diploid chromosomes, 0 for haploid ones
	//	@maxCnt: FG|BG maximum counters to fill
	//	@maxDens: FG|BG maximum densities to fill
	//	return: estimated number of Reads
	static ULONG GetReadsCnt(Gr::Type g, chrlen densLen, float factor, 
		BYTE numeric, ULONG maxCnt[], float maxDens[]);

	// Initializes Reads view
	//	@g: ground
	//	@maxCnt: FG|BG maximum count of reads
	//	@maxDens: FG|BG maximum density
	static void InitReadsView(Gr::Type g, ULONG maxCnt[], float maxDens[])
	{
		ChrView[g].Init(
			ULONG(maxCnt[g] * AutoSample),
			// max() works only in case MakeControl==true and g==Gr::BG,
			// to set maximum BG sample among Test and Control.
			// In other cases Sample(GM::Control, g) always returns 0.
			//max(Sample(GM::Test, g), Sample(GM::Control, g)) * AutoSample,
			Sample(GM::Test, g) * AutoSample,
			maxDens[g] * AutoSample
		);
	}

	// Curs genome into fragments and generate output
	void	CutGenome	();
	
	// Cuts genome into fragments and generate output
	//	@arg: pointer to ChrSubset - set of chrom IDs treated in this thread
	//	@slave: if true then this instance it slave
	inline thrRetValType CutChrom(const effPartition::Subset& arg, bool slave)
	{ 	return ChromCutter(this, NULL, slave).Execute(arg);	}

	// Starts cutting chrom in a separate thread
	//	@arg: pointer to ChrSubset - set of chrom IDs treated in this thread
	static inline thrRetValType 
		#ifdef OS_Windows
		__stdcall 
		#endif
		StatCutChrom(void* arg)	// arg should be void*
		{	return Imit->CutChrom(*(effPartition::Subset*)arg, true); }

public:
	static BYTE	ThrCnt;			// actual number of threads
	static eMode TMode;			// current task mode
	// minimal length of selected fragments: established by --frag-dev or Read::Len if size filter is OFF
	static fraglen FragLenMin;
	static short	FlatLen;	// Boundary flattening length
	// true if total genome is treated.
	// Set to false in Test mode only if single chrom is defined and BG_ALL is false.
	static bool All;
	
	// Print amplification info
	static void PrintAmpl(const char* ampls[]);

	// Returns true if control is set
	static inline bool IsControl() { return MakeControl; }

	// Returns true if single thread is set
	static inline bool IsSingleThread() { return ThrCnt == 1; }

	// Returns true if given verbose level is active
	static inline bool	Verbose(eVerb level)	{ return Verb >= level; }

	// Set number of threads
	static inline void SetThreadNumb(BYTE numb) { 
		FragCnt::Init((ThrCnt=numb) == 1);
		OutFiles::MultiThread = numb > 1;
	}

	// Initializes static values
	//	@tmode: task mode
	//	@input: if true, then 'input' (control) is generated
	static void	Init(
		eMode	tmode,
		bool	input,
		ULONG	cellsCnt,
		bool	ampl,		// false - MDA, true - PCR
		a_coeff	amplCoeff,
		UINT	verb,
		bool	allBg,
		bool	uniformScore,
		//readlen bindLen,
		pairVal& flattens
	) {
		TMode = tmode;
		MakeControl = TestMode ? input : false;
		GlobContext[GM::Test].CellCnt = cellsCnt;
		ChromCutter::SetAmpl(ampl, amplCoeff);
		ChromCutter::FragDistr::Init();
		Verb = verb;
		All = (tmode == CONTROL) || allBg;
		UniformScore = uniformScore;
		//BindLen = bindLen;
		FlatLen = flattens.first + flattens.second;
	}

	// Creates singleton instance.
	//  @cFiles: list of chromosomes as fa-files
	//	@oFile: output files
	inline Imitator(const ChromFiles& cFiles, OutFiles& oFile)
		: _cFiles(cFiles), _oFile(oFile) { Imit = this; }

	// Runs task in current mode and write result to output files
	//	@templ: input template or NULL
	void Execute(BedF* templ);

	// Sets adjusted Sample and clear all counter and means.
	// Sample are needed for the control of BF&FG levels by percent (given by user),
	// and to prorate number of written reads for each chromosome depending on reads limit.
	void	SetSample	();
};
