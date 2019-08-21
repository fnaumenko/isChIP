/**********************************************************
Imitator.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 21.08.2019
-------------------------
Provides chip-seq imitation functionality
***********************************************************/
#pragma once
#include "OutTxtFile.h"
#include "effPartition.h"

using namespace std;

typedef	short	a_cycle;	// type amplification cycles
typedef	BYTE	a_coeff;	// type coefficient of amplification

#define FRAG_MAX	UINT32_MAX		// maximum fragment length

extern const Options::PairVals grounds;

// Task modes
enum eMode { TEST, CONTROL };

enum eVerb {	// verbose level
	vCRIT,		// print critical messages only
	vRES,		// print vCRIT + results
	vRT,		// print vRES + runtime info
	vPAR,		// print vRT + parameters
	vDBG		// print vPAR + additional info
};


#define SignPar	"# "	// Marker of output parameter in Usage
#define	TestMode	(Imitator::TMode==TEST)
#define	ControlMode (Imitator::TMode==CONTROL)

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

				  static float ssFactor0;	// factor sigma*sqrt(2) in the size sel norm distr
			const static float ssFactor1;	// factor 2.5/sqrt(2PI) in the size sel norm distr

			double	_normal_x2;		// second random coordinate (for normal RNG)
			short	_phase;			// phase (for normal RNG)
			Average*_avr;			// average to count the frags (for statistics), or NULL
			//tpLognormNext _pLognormNext;

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

			//fraglen LognormNextWithStat() { return *_avr += LognormNextDirect(); }

		public:
			//	Initializes size sel factors
			inline static void Init() { ssFactor0 = float(DistrParams::ssSigma * sqrt(2.f)); }

			// Creates instance
			//	@avr: generated frags average (for trial mode), or NULL (for working mode)
			// _normal_x2 is initialized by random double to avoid
			//		undesirable 'out of range' random initialization
			FragDistr(Average* avr) : _avr(avr), _normal_x2(DRand())
			{
				//_pLognormNext = avr ? 
				//	&Imitator::ChromCutter::FragDistr::LognormNextWithStat:
				//	&Imitator::ChromCutter::FragDistr::LognormNextDirect;
			}
			
			// Returns next random frag length according to a lognormal distribution,
			// with or without output accumulation to calculate average
			// About 1.5 times faster then std::lognormal_distribution (<random>)
			inline fraglen LognormNext() {
				return _avr ? *_avr += LognormNextDirect() : LognormNextDirect();
				//return (this->*_pLognormNext)(); 
			}

			// Returns next random size sel limits
			//	@min: random min limit
			//	@max: random max limit
			inline void SizeSelLimNext(fraglen& min, fraglen& max)
			{
				float ssDev = ssFactor0 * sqrt(log(ssFactor1 / DRand()));
				min = fraglen(DistrParams::ssMean - ssDev);
				if(min < Read::Len)		min = Read::Len;
				max = fraglen(DistrParams::ssMean + ssDev);
			}

			// Returns random fragment's length within 0 and mean frag length
			// should to round ssMean??
			inline fraglen RandFragLen() { return Range(int(DistrParams::ssMean)); }
		};

	private:
		// 'MDA' implements Multiple Displacement Amplification
		class MDA
		{
		private:
			typedef pair<fraglen,fraglen> Fraction;

			vector<Fraction> _fracs;
			fraglen		 _minLen;
			Random&		_rng;		// used to invoke Range() only
			//a_cycle	_iterCnt;	// number of iteration (displacement reaction)
			vector<Fraction>::const_iterator _it;

			void	Split	(fraglen shift, fraglen len);

		public:
			inline MDA(Random& rng) :_rng(rng) {
				_fracs.reserve(Imitator::IsMDA ? 2*Imitator::SelFragAvr / Read::Len : 1);
			}

			// Prepare instance to the new generating cycle
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
		
		bool		_slave;		// if true then this instance is slave
		GM::Mode	_gMode;		// generating mode: 0 - Test, 1 - Control
		Output*	_output;	// partial output files
		FragCnts	_fragCnt;	// numbers of selected/recorded fragments for FG & BG, for both Teat & Input
		FragDistr	_fragDistr;	// normal & lognormal random number generator
		MDA			_ampl;
		const ChromSizesExt& _cSizes;	// reference genome

		// Sets global mode
		void SetGMode(GM::Mode gmode);

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
		//  @seq: current reference chromosome
		//	@enRegLen: length of all enriched regions
		//	@timer: current timer to thread-saves time output 
		//	@excLimit: true if limit is exceeded
		void PrintChrom (chrid cID, const RefSeq& seq, chrlen enRegLen, Timer& timer, bool excLimit);

		// Sets terminate's output message
		//	@tID: thread ID
		//	@msg: message
		void Terminate(thrid tID, const char*msg) {
			cerr << sThread << int(tID) << SepCl << msg << endl;
		}

	public:
		// Returns true if PCR amplification is established
		inline static bool IsPCR()	{ return _PCRdcycles > 1; }

		// Set amplification
		inline static void SetAmpl() { _PCRdcycles = 1<<AmplCoeff; }

		// Creates instance
		//	@imitator: the owner
		//	@avr: recorded frags average (trial for sampling), or NULL
		//	@slave: if true then this instance is slave
		ChromCutter(const Imitator* imitator, Average* avr, bool slave);
		
		inline ~ChromCutter () { if(_slave)	delete _output; }

		// Returns number of FG recorded frags
		inline ULLONG RecFgFragCnt() const { return _fragCnt[Gr::FG].RecCnt(); }
		
		// Set mode, print chrom name, start timer, return number of cells
		ULLONG PrepareCutting(GM::Mode gm, chrid cID, Timer& timer);

		// Treats chromosomes given for current thread
		//	@cSubset: pointer to ChrSubset - set of chrom IDs treated in this thread
		thrRetValType Execute(const effPartition::Subset& cSubset);
		
		// Cuts chromosome until reaching end position of current treated feature
		//	@seq: cutted reference chromosome
		//	@fragStart: fragment start position
		//	@feature: current treated feature
		//	@bg: if true then generate background (swap foreground and background)
		//	@fragStat: statistics of selected fragments average counter, or NULL under work mode
		//	return: 0 if success,
		//		1 if end chromosome is reached (continue treatment),
		//		-1 if limit is achieved (cancel treatment)
		int	CutChrom (const RefSeq& seq, chrlen* const fragStart, const Featr& feature,
			bool bg, FragLenStat* fragStat = NULL);
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
		static inline void PrintMarg(BYTE width) { cout << setw(width) << BLANK; }

		BYTE CountW;		// width of 'number of Reads' field
		BYTE SampleWRP;		// 'sample' field right padding: space to the right of the value (0|1)
		BYTE DensWRP;		// 'density' field right padding: space to the right of the value (0|1)
		BYTE DensPr;		// precision of density
		Gr::Type GrType;	// FG|BG

	public:
		// Gets the maximum length of chrom name field with blank after
		static const BYTE ChromNameW() { 
			// 2 = 1 for blank inside chrom name + 1 for COLON
			return Chrom::MaxAbbrNameLength + 2 + margC_;
		}

		// returns ground title width or 0 if not TEST mode
		static BYTE GrTitleW() { return TestMode * (Gr::TitleLength + strlen(SepCl)); }

		// Prints chroms gaps statistics
		static void PrintGaps(const GenomeSizes& s);

		static inline void PrintGapsMarg() { PrintMarg(margD_ + GapsW + margGexcl + GapsWexcl); }

		inline ChromView(Gr::Type grType) : GrType(grType) {}

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
		void PrintReads(GM::Mode gMode, const FragCnt fragCnt, chrlen densLen);

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

	static GenomeSizes gSizes;		// genome total sizes; needed to define average % of gaps, excl gaps
	static ChromView ChrView[];		// FG, BF chrom views
	static ULONG	TreatedLen[];	// FG, BF genome treated length; needed to define total density
	static Context	GlobContext[];	// TM, CM global generation context
	static float	AutoSample;		// adjusted FG sample to stay in limit
	static float	SelFragAvr;		// mean length of selected fragments
	static a_coeff	AmplCoeff;		// user-stated amplification coefficient
	static BYTE	Verb;				// verbose level
	static bool	MakeControl;		// true if control file (input) should be produced
	static Imitator	*Imit;			// singletone instance: to call threads only
	static const Features *Templ;	// template features or NULL (control mode)
	//static readlen	BindLen;	// binding length

	const ChromSizesExt& _cSizes;	// ref genome library
	Output& _oFile;				// output file

	// Returns stated count of cells
	//	@gm: generated mode
	inline static UINT CellCnt(GM::Mode gm) { return GlobContext[gm].CellCnt; }

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

	// Returns stated sample
	//	@gm: generated mode
	//	@g: ground 
	inline static float Sample(GM::Mode gm, Gr::Type g) { return GlobContext[gm].Sample[g]; }

	// Prints chromosome's name
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@gm: generation mode Test|Control
	//	@print: true if chromosomes name should be printed
	static void PrintChromName(chrid cID, GM::Mode gm, bool print);

	inline static void PrintReadInfo(
		Gr::Type gr, GM::Mode gMode, const FragCnt fCnts[], const ULONG rgnLens[])
	{
		ChrView[gr].PrintReads(gMode, fCnts[gr], rgnLens[gr]);
	}

	// Prints chromosome's name and full statistics
	//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@gMode: Test|Control
	//	@fCnts: array of fragment's counters, for FG and BG
	//	@rgnLens: array of region's lengths to print FG|BG density
	//	@prChrName: true if chromosome's name should be printed
	static void PrintChromInfo(chrid cID, GM::Mode gMode,
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
	//static bool IncrementRecFragCount(Gr::Type g, bool primer) {
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

	// Imitates cutting chromosome to reach statistics
	ULLONG CutForSample(Average& genFrAvr, FragLenStat* frLenStat);

	// Sets adjusted Sample and clear all counter and means.
	// Sample are needed for the control of BF&FG levels by percent (given by user),
	// and to prorate number of written reads for each chromosome depending on reads limit.
	void SetSample	();

public:
	static bool	UniformScore;	// true if template features scores are ignored
	static BYTE	ThrCnt;			// actual number of threads
	static bool	IsMDA;
	static eMode TMode;			// current task mode
	static short FlatLen;		// BS edge flattening length
	static fraglen FragLenMin;	// minimal length of selected fragments:
								// established by --frag-dev or Read::Len if size filter is OFF
	static bool All;	// true if total genome is treated
						// false only in Test mode if custom chrom and BG_ALL is false.
	
	// Print amplification info
	static void PrintAmpl();

	// Returns true if control is set
	static inline bool IsControl()			{ return MakeControl; }

	// Returns true if single thread is set
	static inline bool IsSingleThread()		{ return ThrCnt == 1; }

	// Returns true if given verbose level is active
	static inline bool Verbose(eVerb level)	{ return Imitator::Verb >= level; }

	// Set number of threads
	static void SetThreadNumb(BYTE numb) { 
		FragCnt::Init((ThrCnt=numb) == 1);
		ReadName::MultiThread = numb > 1;
	}

	// Initializes static values
	//	@tmode: task mode
	//	@input: if true, then 'input' (control) is generated
	static void	Init(
		eMode	tmode,
		bool	input,
		ULONG	cellsCnt,
		bool	isMDA,
		a_coeff	amplCoeff,
		UINT	verb,
		bool	allBg,
		bool	uniformScore,
		//readlen bindLen,
		//const pairVal& flattens
		UINT	unstBindLen		// unstable binding length
	) {
		TMode = tmode;
		MakeControl = TestMode ? input : false;
		GlobContext[GM::Test].CellCnt = cellsCnt;
		IsMDA = isMDA;
		AmplCoeff = amplCoeff;	// the actual ChromCutter ampl coeff will be set in Sample()
		ChromCutter::FragDistr::Init();
		Verb = verb;
		All = (tmode == CONTROL) || allBg;
		UniformScore = uniformScore;
		//BindLen = bindLen;
		//FlatLen = flattens.first + flattens.second;
		FlatLen = unstBindLen;
		SelFragAvr = DistrParams::LnMean();	// for the Sample(), before set actual ln mean
	}

	// Creates singleton instance.
	//  @cFiles: list of chromosomes as fa-files
	//	@ocSizes: chrom sizes
	inline Imitator(const ChromSizesExt& cSizes, Output& oFile)
		: _cSizes(cSizes), _oFile(oFile) { Imit = this; }

	// Runs task in current mode and write result to output files
	//	@templ: input template or NULL
	void Execute(Features* templ);
};
