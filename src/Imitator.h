/**********************************************************
Imitator.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 11/12/2023
-------------------------
Provides chip-seq imitation functionality
***********************************************************/
#pragma once
#include "DataOutFile.h"
#include "effPartition.h"
//#include <math.h>       /* log */

using namespace std;

typedef	short	a_cycle;	// type amplification cycles
typedef	BYTE	a_coeff;	// type coefficient of amplification

#define FRAG_MAX	UINT32_MAX		// maximum fragment length

extern const Options::PairVals grounds;

// Task modes
enum eMode { TEST, CONTROL };	// using 'enum class' in this case is inconvenient

enum class eVerb {	// verbose level
	CRIT,	// print critical messages only
	RES,	// print vCRIT + results
	RT,		// print vRES + runtime info
	PAR,	// print vRT + parameters
	DBG		// print vPAR + additional info
};


#define SignPar	"# "	// Marker of output parameter in Usage
#define	TestMode	(Imitator::TMode==TEST)
#define	ControlMode (Imitator::TMode==CONTROL)

// 'Gr' defines ground types
static class Gr
{
	static const char* title[];
public:
	enum eType {		// using 'class' in this case is inconvenient
		FG = 0,	// foreground
		BG = 1	// background
	};
	// Count of grounds
	static const int Cnt = 2;

	// Maximum title length
	static const BYTE TitleLength = 2;

	// Gets ground title
	inline static const char* Title(eType g) { return title[g]; }
} ground;



// 'Imitator' implements main algorithm of simulation.
class Imitator
{
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
		inline ULLONG Sum()	 const { 
			return _summator; }
		inline ULONG Count() const { return _count; }
		// Returns parameter val
		UINT operator+=(UINT val);
	};

//public:
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

//private:
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
		inline void SelIncr() { sel += 2; }	// because of double fragment (forward & backward) for each selection

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
		GM::eMode	gMode;

	public:
		inline void SetGMode(GM::eMode mode) { gMode = mode; }

		inline const FragCnt* GetFragCnts() const { return fCnts[int(gMode)]; }

		inline FragCnt& operator[](BYTE g)  { return fCnts[int(gMode)][g]; }

		inline const FragCnt& operator[](BYTE g) const { return fCnts[int(gMode)][g]; }

		// Clears all fragment counters
		inline void Clear() { memset(fCnts, 0, 2*Gr::Cnt*sizeof(FragCnt)); }
	};

	// 'GenomeSizes' keeps total, defined and gaps length in purpose to print average gap statistics
	struct GenomeSizes
	{
	private:
		ULONG	Total;		// sum of chroms length
		ULONG	Defined;	// sum if chroms defined (effective) length
		ULONG	Gaps;		// sum of chroms gaps length

	public:
		inline GenomeSizes() { Total = Defined = Gaps = 0; }

		// constructor by single chrom
		inline GenomeSizes(const RefSeq& seq) : 
			Total(seq.Length()), Defined(seq.DefRegion().Length()), Gaps(seq.GapLen()) {}

		// Thread-safety increment sizes by RefSeq
		void IncrSizes(const RefSeq& seq);

		inline float GapsInPers() const { return 100.f * Gaps / Total; }

		inline float UndefInPers() const { return 100.f * (Total - Defined) / Total; }
	};

	//public:
	// 'ChromCutter' encapsulates thread context and methods to cut chromosomes.
	class ChromCutter
	{
	public:
		// 'FragDistr' encapsulates normal and lognormal random generator and it's average counter.
		class FragDistr : public Random
		{
		private:
			//typedef fraglen	(FragDistr::*tpLognormNext)();

			Average*_avr;			// average to count the frags (for statistics), or NULL
			//tpLognormNext _pLognormNext;

			//fraglen LognormNextWithStat() { return *_avr += Lognormal(); }

		public:
			//	Initializes size sel factors
			inline static void Init() { ssFactor0 = float(DistrParams::ssSigma * sqrt(2.f)); }

			// Creates instance
			//	@avr: generated frags average (for trial mode), or NULL (for working mode)
			// _normal_x2 is initialized by random double to avoid
			//		undesirable 'out of range' random initialization
			FragDistr(Average* avr) : _avr(avr)
			{
				//_pLognormNext = avr ? 
				//	&Imitator::ChromCutter::FragDistr::LognormNextWithStat:
				//	&Random::Lognormal;
			}
			
			// Returns next random frag length according to a lognormal distribution,
			// with or without output accumulation to calculate average
			// About 1.5 times faster then std::lognormal_distribution (<random>)
			inline fraglen LognormNext() {
				return _avr ? *_avr += Lognormal() : Lognormal();
				//return (this->*_pLognormNext)(); 
			}

			// Returns next random size sel limits
			//	@min: random min limit
			//	@max: random max limit
			void SizeSelLimits(fraglen& min, fraglen& max)
			{
				const float ssDev = ssFactor0 * float(sqrt(log(ssFactor1 / DRand())));
				
				min = fraglen(DistrParams::ssMean - ssDev);
				if(min < Read::FixedLen)		min = Read::FixedLen;
				max = fraglen(DistrParams::ssMean + ssDev);
			}

			// Returns random fragment's length within 0 and mean frag length
			// should to round ssMean??
			inline fraglen RandFragLen() { return Range(int(DistrParams::ssMean)); }

			//bool ExpoProbability(const readlen len) {
			//	//return -log(DRand()) / len;
			//	//return DRand() <= (2 ^ Range(len) - 1) / (2 ^ len - 1);
			//	return DRand() <= (2 ^ Range(len)) / (2 ^ len);
			//}

			//readlen ExpoProbability(const readlen len) {
			//	return readlen(pow(len * DRand(), 3) / (len * len));
			//}
		};

	private:
		typedef pair<fraglen, fraglen> Fraction;	// first: frac shift, second: frac length

		// 'MDA' implements Multiple Displacement Amplification
		class MDA : public vector<Fraction>
		{
			Random&	_rng;		// used to invoke Range() only

			void Split(fraglen shift, fraglen len, fraglen minLen);
		public:
			MDA(Random& rng) :_rng(rng)
			{ reserve(IsMDA ? size_t(4 * SelFragAvr / Read::FixedLen) : 1); }

			// Generats amplified fragment collection
			//	@fLen: current fragment length
			//	@fLenMin: current minimal fragment length
			void  Generate (fraglen fLen, fraglen fLenMin);
		};

		static fraglen	_SsDev;			// deviation of frag size selection
		static a_cycle	_PCRdcycles;	// PCR cycles: read doubling cycles
		
		bool		_master;	// if true then this instance is master
		GM::eMode	_gMode;		// generating mode: 0 - Test, 1 - Control
		Output*		_output;	// partial output files
		FragCnts	_fragCnt;	// numbers of selected/recorded fragments for FG & BG, for both Teat & Input
		FragDistr	_fragDistr;	// normal & lognormal random number generator
		MDA			_ampl;
		const ChromSizesExt& _cSizes;	// reference genome

		// Sets global mode
		void SetGMode(GM::eMode gmode);

		// Increments counters of local and total recorded fragments thread-safely
		//	@g: ground
		//	@primer: true if increment derived (amplified) frag's counter
		//	return: true if Reads limit is exceeded.
		bool IncrRecFragCount(Gr::eType g, bool primer);

		// Increments counter of total selected fragments thread-safely
		void IncrTotalSelFragCount();

		// Returns local (instance-defined) random true or false according to rate
		inline bool Sample(float rate) { return _fragDistr.Sample(rate); }

		// Prints thread-safe info about treated chroms and stops timer
		//  @seq: current reference chromosome
		//	@enRegLen: length of all enriched regions
		//	@timer: current timer to thread-saves time output 
		//	@excLimit: true if limit is exceeded
		void PrintChrom (const RefSeq& seq, chrlen enRegLen, Timer& timer, bool excLimit);

		// Sets terminate's output message
		//	@tID: thread ID
		//	@msg: message
		void Terminate(thrid tID, const char*msg) {
			cerr << sThread << int(tID) << SepCl << msg << endl;
		}

		// Returns sample of Flattening of binding site suburb
		//	@fStart: tested frag's start
		//	@fEnd: tested frag's end
		//	@ft: binding site
		//	@sample: corrected sample
		void GetFlattSample(chrlen fStart, chrlen fEnd, const Featr& ft, bool& sample);

	public:
		// Returns true if PCR amplification is established
		inline static bool IsPCR()	{ return _PCRdcycles > 1; }

		// Set amplification
		inline static void SetAmpl() { _PCRdcycles = 1<<PCRCoeff; }

		// Creates instance
		//	@imitator: the owner
		//	@avr: recorded frags average (trial for sampling), or NULL
		//	@master: if true then this instance is master
		ChromCutter(const Imitator* imitator, Average* avr, bool master);
		
		inline ~ChromCutter () { if(!_master)	delete _output; }

		// Returns number of FG recorded frags
		inline ULLONG RecFgFragCnt() const { return _fragCnt[Gr::FG].RecCnt(); }
		
		// Set mode, print chrom name, start timer, return number of cells
		ULONG PrepareCutting(GM::eMode gm, chrid cID, Timer& timer);

		// Treats chromosomes given for current thread
		//	@cSubset: pointer to ChrSubset - set of chrom IDs treated in this thread
		void Execute(const effPartition::Subset& cSubset);
		
		// Cuts chromosome until reaching end position of current treated feature
		//	@cLen: chromosome's 'end' position
		//	@fStart: fragment start position
		//	@ft: current treated feature
		//	@scores: FG (in-feature, first) and BG (out-feature, second) scores
		//	@bg: if true then generate background (swap foreground and background)
		//	@fStat: statistics of selected fragments average counter, or NULL under work mode
		//	return: 0 if success,
		//		1 if end chromosome is reached (continue treatment),
		//		-1 if limit is achieved (cancel treatment)
		int	CutChrom (chrlen cLen, chrlen& fStart, const Featr& ft, 
			float scores[], bool bg, FragLenStat* fStat = NULL);
	};

	// 'ChromView' provides template and methods for viewing chrom's treatment results
	struct ChromView
	{
	/*
	* output view:
	* A: always aligned                    A                                A
	* chrom       reads  sample ampl r/kbp|        reads sample ampl  r/kbp|  gaps g_excl  mm:ss
	* ------------------------------------|--------------------------------|--------------------
	* chr 17: FG: 111111  100%   1.4  11.3|   BG: 222222   10%  0.622  0.59|  3.5%  3.1%   00:03
	* -------^---^-------^-----^------^------^----^-------^-----^------^-----^------^------^-----
	* chrNam      CountW  SmpW  AmplW  DnsW        CountW  SmpW  AmplW  DnsW  GapsW  GapsW  time
    * margN: 1   2       3     4      5      6    2       3     4      5     7      8      9
	*/
	public:
		static const BYTE marg_T	= 2;	// length of margin 9 before time
	private:
		static const BYTE margC_	= 1;	// length of margin 1 after chrom name
		static const BYTE marg_R	= 1;	// length of margin 2 before Reads count
		static const BYTE marg_S	= 2;	// length of margin 3 before sample
		static const BYTE marg_A	= 3;	// length of margin 4 before ampl coeff
		static const BYTE marg_D	= 1;	// length of margin 5 before density
		static const BYTE margD_	= 2;	// length of margin 6 after density
		static const BYTE margGexcl	= 1;	// length of margin 8 between gaps and gaps excleded
		static const BYTE SampleW	= 6;	// width of 'sample' field	[strlen("sample");]
		static const BYTE DensW		= 6;	// width of 'density' field [strlen(UnitDens)]
		static const BYTE AmplW		= 4;	// width of ampl coefficient field
		static const BYTE GapsW		= 4;	// width of gaps field
		static const BYTE AmplPr	= 1;	// precision of ampl coefficient
		static const BYTE SamplePr	= 2;	// precision of sample
		static const BYTE GapsPr	= 2;	// precision of gaps

		static const char* tGapsExcl;		// title of excluded gaps
		static const char* tTime;			// title of time output

		static BYTE GapsWexcl;		// width of excluded gaps field

		// Prints empty margin
		//	@width: margin width
		static inline void PrintMarg(BYTE width) { cout << setw(width) << SPACE; }

		BYTE CountW;		// width of 'number of Reads' field
		BYTE SampleWRP;		// 'sample' field right padding: space to the right of the value (0|1)
		BYTE DensWRP;		// 'density' field right padding: space to the right of the value (0|1)
		BYTE DensPr;		// precision of density
		Gr::eType GrType;	// FG|BG

	public:
		// Gets the maximum length of chrom name field with blank after
		static const BYTE ChromNameW() { 
			// 2 = 1 for blank inside chrom name + 1 for COLON
			return Chrom::MaxAbbrNameLength + 2 + margC_;
		}

		// returns ground title width or 0 if not TEST mode
		static BYTE GrTitleW() { return BYTE(TestMode * (Gr::TitleLength + strlen(SepCl))); }

		// Prints chroms gaps statistics
		static void PrintGaps(const GenomeSizes& s);

		static inline void PrintGapsMarg() { PrintMarg(margD_ + GapsW + margGexcl + GapsWexcl); }

		inline ChromView(Gr::eType grType) : GrType(grType) {}

		// Prints chrom's header: Reads
		//	@FgHeader: true for foreground header
		//	return: header width
		int PrintHeader(bool FgHeader);

		// Prints chrom's header: gaps statistics
		static int PrintHeaderGaps();

		// Print info about recorded Reads (count, sample, ampl, density) for ground type of this nstance
		//	@gMode: gen mode for which info is printed
		//  @fragCnt: Reads counters
		//	@densLen: length on which the density is determined
		void PrintReads(GM::eMode gMode, const FragCnt fragCnt, chrlen densLen);

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
		bool IncrRecFragCount(Gr::eType g, bool primer) {
			return fCnts[g].RecIncrSaved(primer) + fCnts[!g].RecCnt()
				>= Seq::FragsLimit();
		}
	};

	static GenomeSizes gSizes;		// genome total sizes; needed to define average % of gaps, excl gaps
	static ChromView ChrView[];		// FG, BF chrom views
	static ULONG	TreatedLen[];	// FG, BF genome treated length; needed to define total density
	static Context	GlobContext[];	// TM, CM global generation context
	static float	AutoSample;		// adjusted FG sample to stay in limit
	static float	SelFragAvr;		// mean length of selected fragments
	static a_coeff	PCRCoeff;		// user-stated amplification coefficient
	static eVerb	Verb;				// verbose level
	static bool	MakeControl;		// true if control file (input) should be produced
	static Imitator	*Imit;			// singletone instance: to call threads only
	static const Features *Templ;	// template features or NULL (control mode)
	//static readlen	BindLen;	// binding length

	const ChromSizesExt& _cSizes;	// ref genome library
	Output& _oFile;				// output file

	// Returns stated count of cells
	//	@gm: generated mode
	inline static UINT CellCnt(GM::eMode gm) { return GlobContext[int(gm)].CellCnt; }

	// Gets Reads estimated count per chrom and corrects maximum count and density
	//	@g: ground
	//	@densLen: length on which density is defined
	//	@factor: CellsCnt / recAvrLen
	//	@numeric: 1 for diploid chromosomes, 0 for haploid ones
	//	@maxCnt: FG|BG maximum counters to fill
	//	@maxDens: FG|BG maximum densities to fill
	//	return: estimated number of Reads
	static ULONG GetReadsCnt(Gr::eType g, chrlen densLen, float factor, 
		BYTE numeric, ULONG maxCnt[], float maxDens[]);

	// Returns stated sample
	//	@gm: generated mode
	//	@g: ground 
	inline static float Sample(GM::eMode gm, Gr::eType g) { return GlobContext[int(gm)].Sample[g]; }

	// Prints chromosome's name
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@gm: generation mode Test|Control
	//	@print: true if chromosomes name should be printed
	static void PrintChromName(chrid cID, GM::eMode gm, bool print);

	inline static void PrintReadInfo(
		Gr::eType gr, GM::eMode gMode, const FragCnt fCnts[], const ULONG rgnLens[])
	{ ChrView[gr].PrintReads(gMode, fCnts[gr], rgnLens[gr]); }

	// Prints chromosome's name and full statistics
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@gMode: Test|Control
	//	@fCnts: array of fragment's counters, for FG and BG
	//	@rgnLens: array of region's lengths to print FG|BG density
	//	@prChrName: true if chromosome's name should be printed
	static void PrintChromInfo(chrid cID, GM::eMode gMode,
		const FragCnt fCnts[], const ULONG rgnLens[], bool prChrName=true);

	// Prints header (FG and BG, or single) and solid line
	//	@header: if false, then print solid line only
	static void PrintHeader(bool header);

	// Prints total outcome
	static void PrintTotal();

	// Increments grounds total length
	static void IncrementTotalLength(const RefSeq& seq, chrlen enrRgnLen);
	
	//// Increments counter of total recorded fragments thread-safely.
	////	@g: ground
	////	@primer: if true then increment amplified (second) frag's counter
	////	return: true if limit is exceeded.
	//static bool IncrementRecFragCount(Gr::eType g, bool primer) {
	//	ULLONG cnt = ChromCutter::IsAmpl() ?
	//		TotalFragCnt[g].InterlockedIncrRec(primer):
	//		Output::SingleThread ?
	//			TotalFragCnt[g].PrimerRec++:
	//			InterlockedIncrement(&(TotalFragCnt[g].PrimerRec));
	//	return cnt + TotalFragCnt[!g].Rec() >= Seq::FragsLimit();
	//}

	// Initializes Reads view
	//	@g: ground
	//	@maxCnt: FG|BG maximum count of reads
	//	@maxDens: FG|BG maximum density
	static void InitReadsView(Gr::eType g, ULONG maxCnt[], float maxDens[])
	{
		ChrView[int(g)].Init(
			ULONG(maxCnt[int(g)] * AutoSample),
			// max() works only in case MakeControl==true and g==Gr::BG,
			// to set maximum BG sample among Test and Control.
			// In other cases Sample(GM::Control, g) always returns 0.
			//max(Sample(GM::Test, g), Sample(GM::Control, g)) * AutoSample,
			Sample(GM::eMode::Test, g) * AutoSample,
			maxDens[int(g)] * AutoSample
		);
	}

	// Curs genome into fragments and generate output
	void	CutGenome	();
	
	// Cuts genome into fragments and generate output
	//	@arg: pointer to ChrSubset - set of chrom IDs treated in this thread
	//	@master: if true then this instance it master
	inline void CutChrom(const effPartition::Subset& arg, bool master)
	{ ChromCutter(this, NULL, master).Execute(arg); }

	// Imitates cutting chromosome to reach statistics
	//	return: estimated total number of Reads
	ULLONG CutForSample(Average& genFrAvr, FragLenStat* frLenStat);

	// Sets adjusted Sample and clear all counter and means.
	// Sample are needed for the control of BF&FG levels by percent (given by user),
	// and to prorate number of written reads for each chromosome depending on reads limit.
	void SetSample	();

public:
	static bool	UniScore;	// true if template features scores are ignored
	static BYTE	ThrCnt;			// actual number of threads
	static bool	IsExo;
	static bool	IsMDA;
	static eMode TMode;			// current task mode
	static short FlatLen;		// BS edge flattening length
	static fraglen FragLenMin;	// minimal length of selected fragments:
								// established by --frag-dev or Read::FixedLen if size filter is OFF
	static bool All;	// true if total genome is treated
						// false only in Test mode if custom chrom and BG_ALL is false.
	
	// Print amplification info
	//	@signOut: output marker
	static void PrintAmpl(const char* signOut);

	// Returns true if control is set
	static inline bool IsControl()			{ return MakeControl; }

	// Returns true if single thread is set
	static inline bool IsSingleThread()		{ return ThrCnt == 1; }

	// Returns true if given verbose level is active
	static inline bool Verbose(eVerb level)	{ return Imitator::Verb >= level; }

	// Set number of threads
	static void SetThreadNumb(BYTE numb) { 
		FragCnt::Init((ThrCnt=numb) == 1);
		ReadOutFile::MultiThread = numb > 1;
	}

	// Initializes static values
	static void	Init(
		eMode	tmode,				// task mode
		bool	input,				// true if control should be generated as well
		ULONG	cellsCnt,			// count of cells
		bool	isExo,				// true if EXO mode
		bool	isReadLenAssigned,	// true if Read length is assigned by user
		bool	isMDA,				// true if MDA is assigned
		a_coeff	amplCoeff,			// coefficient of PCR
		UINT	verb,				// verbosity level
		bool	allBg,				// true if all background mode is assigned
		//bool	uniformScore,		// true if uniform template score is assigned
		//readlen bindLen,
		//const pairVal& flattens
		UINT	unstBindLen			// unstable binding length
	);

	// Creates singleton instance.
	//  @cFiles: list of chromosomes as fa-files
	//	@ocSizes: chrom sizes
	inline Imitator(const ChromSizesExt& cSizes, Output& oFile)
		: _cSizes(cSizes), _oFile(oFile) { Imit = this; }

	// Runs task in current mode and write result to output files
	//	@templ: input template or NULL
	void Execute(Features* templ);
};
