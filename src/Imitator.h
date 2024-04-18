/**********************************************************
Imitator.h
Provides chip-seq imitation functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 04/18/2024
***********************************************************/
#pragma once

#include "ChromSeq.h"
#include "DataWriter.h"
#include "effPartition.h"
#include "Features.h"

using cells_cnt = uint16_t;	// type number of cells
using a_cycle	= int16_t;	// type PCR amplification cycles
using a_coeff	= BYTE;		// type coefficient of PCR amplification

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
	
	static const BYTE Cnt = 2;			// Count of grounds
	static const BYTE TitleLength = 2;	// Maximum title length

	// Gets ground title
	static const char* Title(eType g) { return title[g]; }
} ground;


// 'Imitator' implements main algorithm of simulation.
class Imitator
{
	// 'Average' encapsulates average counting
	class Average
	{
		bool	 _keep = true;		// true if statistic is accumulated
		uint32_t _count = 0;
		size_t	 _summator = 0;

	public:
		// Gets arithmetic average.
		float	 Mean()	 const { return (float)_summator/_count; }
		size_t	 Sum()	 const { return _summator; }
		uint32_t Count() const { return _count; }
		// Returns parameter val
		uint32_t operator += (uint32_t val);
	};

	// 'FragLenStat' keeps statistics of selected frag length; used in SetSample()
	struct FragLenStat
	{
		fraglen Min;	// minimum selected frag length
		fraglen	Max;	// maximum selected frag length
		Average	SelAvr;	// average selected frag length 

		FragLenStat() : Min(FRAG_MAX), Max(0) {}

		// Accepts frag length in purpose of statistics
		//	@param len: current selected frag length
		//	invoked once
		void AddFragLen(fraglen len) {
			SelAvr += len;
			if(len > Max)		Max = len;
			else if(len < Min)	Min = len;
		}
	};

	// Fragment's Counter keeps statistics for selected and recorded Reads
	class FragCnt {
		typedef size_t	(FragCnt::*pRecIncr)(bool);
		typedef void	(FragCnt::*pSelAdd)	(size_t);

		// pointer to the thread-saved 'recorded frag's number increment' method
		static pRecIncr	pRecIncrSaved;
		// pointer to the thread-saved 'selected frag's number adding' method
		static pSelAdd	pSelAddSaved;

		size_t sel;		// number of selected fragments
		size_t rec[2];	// number of recorded fragments: 0: derived (amplified), 1: primer (initial)

		// Thread-savely increments number of recorded frags
		//	@param primer: if true then increment counter of primer (not amplified) frag
		//	@returns: number of recorded frags after increment
		size_t RecInterlIncr(bool primer) { return InterlockedIncrement(&rec[primer]); }

		// Thread-savely adds value to number of selected frags
		void SelAddIncr(size_t val) { InterlockedExchangeAdd(&sel, val);	}

		// Adds value to number of selected frags
		void SelAdd(size_t val) { sel += val; }

	public:
		static void Init(bool singleThread) {
			pRecIncrSaved = singleThread ? &Imitator::FragCnt::RecIncr : &Imitator::FragCnt::RecInterlIncr;
			pSelAddSaved =	singleThread ? &Imitator::FragCnt::SelAdd : &Imitator::FragCnt::SelAddIncr;
		}

		void Clear() { memset(rec, 0, 2*sizeof(size_t));	}

		// Gets total number of selected frags
		size_t SelCnt()	const { return sel; }

		// Gets total number of recorded frags
		size_t RecCnt()	const { return rec[0] + rec[1]; }

		// Gets percentage of recorded reads in relation to the selected
		float Sample()	const { return Percent(rec[1], sel); }
		
		// Gets ratio primer/derived frag's numbers
		float RealAmplCoeff() const { return rec[1] ? float(rec[0]) / rec[1] : 0; }

		// Increments number of recorded frags
		//	@param primer: if true then increment counter of primer (not amplified) frag
		//	@returns: number of recorded frags after increment
		void SelIncr() { sel += 2; }	// because of double fragment (forward & backward) for each selection

		// Increments number of recorded frags
		//	@param primer: if true then increment counter of primer (not amplified) frag
		//	@returns: number of recorded frags after increment
		size_t RecIncr(bool primer) { return rec[primer]++; }

		// Thread-savely increments number of recorded frags
		//	@param primer: true if increment counter of primer (not amplified) frag
		//	@returns: number recorded frags after increment
		size_t RecIncrSaved(bool primer) {
			return (this->*pRecIncrSaved)(primer), RecCnt();
		}

		// Thread-savely adds value to number of selected frags
		void SelAddSaved(size_t val) { (this->*pSelAddSaved)(val); }
	};

	// Fragment's Counters keeps statistics for selected and recorded Reads, for both Test and MakeControl modes
	class FragCnts {
		FragCnt	fCnts[2][Gr::Cnt];
		GM::eMode	gMode;

	public:
		void SetGMode(GM::eMode mode) { gMode = mode; }

		const FragCnt* GetFragCnts() const { return fCnts[int(gMode)]; }

		FragCnt& operator[](BYTE g)  { return fCnts[int(gMode)][g]; }

		const FragCnt& operator[](BYTE g) const { return fCnts[int(gMode)][g]; }

		// Clears all fragment counters
		void Clear() { memset(fCnts, 0, 2*Gr::Cnt*sizeof(FragCnt)); }
	};

	// Genome Sizes keeps total, defined and gaps length in purpose to print average gap statistics
	class GenomeSizes
	{
		genlen	Total = 0;		// sum of chroms length
		genlen	Defined = 0;	// sum of chroms defined (effective) length
		genlen	Gaps = 0;		// sum of chroms gaps length

	public:
		GenomeSizes() {}

		// constructor by single chrom
		GenomeSizes(const ChromSeq& seq) : 
			Total(seq.Length()), Defined(seq.DefRegion().Length()), Gaps(seq.GapLen()) {}

		// Thread-safety increment sizes by ChromSeq
		void IncrSizes(const ChromSeq& seq);

		float GapsInPers() const { return 100.f * Gaps / Total; }

		float UndefInPers() const { return 100.f * (Total - Defined) / Total; }
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

			Average*_avr;		// average to count the frags (for statistics), or NULL
			//tpLognormNext _pLognormNext;
			//fraglen LognormNextWithStat() { return *_avr += Lognormal(); }

		public:
			//	Initializes size sel factors
			static void Init() { ssFactor0 = float(DistrParams::ssSigma * sqrt(2.f)); }

			// Creates instance
			//	@param avr: generated frags average (for trial mode), or NULL (for working mode)
			// _normal_x2 is initialized by random double to avoid undesirable 'out of range' random initialization
			FragDistr(Average* avr) : _avr(avr)
			{
				//_pLognormNext = avr ? 
				//	&Imitator::ChromCutter::FragDistr::LognormNextWithStat:
				//	&Random::Lognormal;
			}
			
			// Returns next random frag length according to a lognormal distribution,
			// with or without output accumulation to calculate average
			// About 1.5 times faster then std::lognormal_distribution (<random>)
			fraglen LognormNext() {
				return _avr ? *_avr += Lognormal() : Lognormal();
				//return (this->*_pLognormNext)(); 
			}

			// Sets next random size sel limits
			//	@param min[out]: next set random min limit
			//	@param max[out]: next set random max limit
			void SizeSelLimits(fraglen& min, fraglen& max);

			// Returns random fragment's length within 0 and mean frag length
			// should to round ssMean??
			fraglen RandFragLen() { return Range(int(DistrParams::ssMean)); }

			//bool ExpoProbability(const readlen len) {
			//	//return -log(DRand()) / len;
			//	return DRand() <= (2 ^ Range(len)) / (2 ^ len);
			//}

			//readlen ExpoProbability(const readlen len) {
			//	return readlen(pow(len * DRand(), 3) / (len * len));
			//}
		};

		typedef pair<fraglen, fraglen> Fraction;	// first: fraction shift, second: fraction length

		// 'MDA' implements Multiple Displacement Amplification
		class MDA : public vector<Fraction>
		{
			Random&	_rng;			// used to invoke method Range() only
			fraglen _fLenMin = 0;	// keeps min frag length during a single amplification process

			void Split(fraglen shift, fraglen len);

		public:
			MDA(Random& rng) :_rng(rng)
			{ reserve(IsMDA ? size_t(4 * FragMean / Read::FixedLen) : 1); }

			// Generates amplified fragments collection
			//	@param fLen: initial fragment length
			//	@param fLenMin: minimum allowable fragment length
			void  Generate (fraglen fLen, fraglen fLenMin);
		};
	private:

		static fraglen	_SsDev;			// deviation of frag size selection
		static a_cycle	_PCRdcycles;	// PCR cycles: read doubling cycles
		
		bool		_master;	// if true then this instance is master
		GM::eMode	_gMode;		// generating mode: 0 - Test, 1 - Control
		DataWriter*	_writer;
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
		bool Sample(float rate) { return _fragDistr.Sample(rate); }

		// Prints thread-safe info about treated chroms and stops timer
		//  @seq: current reference chromosome
		//	@enRegLen: length of all enriched regions
		//	@timer: current timer to thread-saves time output 
		//	@excLimit: true if limit is exceeded
		void PrintChrom (const ChromSeq& seq, chrlen enRegLen, Timer& timer, bool excLimit);

		// Sets terminate's output message
		//	@param tID: thread ID
		//	@param msg: message to output to cerr
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
		static bool IsPCR()	{ return _PCRdcycles > 1; }

		// Set amplification
		static void SetAmpl() { _PCRdcycles = 1<<PCRCoeff; }

		// Constructor
		//	@param imitator: the owner
		//	@param avr: recorded frags average (trial for sampling), or NULL
		//	@param master: if true then this instance is master
		ChromCutter(const Imitator* imitator, Average* avr, bool master);
		
		~ChromCutter () { if(!_master)	delete _writer; }

		// Set mode, print chrom name, start timer
		//	@param gm: generation mode
		//	@param cID: chrom ID
		//	@param timer: timer to start
		//	@returns: number of cells
		cells_cnt Init(GM::eMode gm, chrid cID, Timer& timer);

		// Treats chromosomes given for current thread
		//	@param cSubset: pointer to ChrSubset - set of chrom IDs treated in this thread
		void Execute(const effPartition::Subset& cSubset);
		
		// Cuts chromosome until reaching end position of current treated feature
		//	@param cLen: chromosome's 'end' position
		//	@param fStart: fragment start position
		//	@param ft: current treated feature
		//	@param scores: FG (in-feature, first) and BG (out-feature, second) scores
		//	@param bg: if true then generate background (swap foreground and background)
		//	@param fStat: statistics of selected fragments average counter, or NULL under work mode
		//	@returns: 0 if success,
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
		//	@param width: margin width
		static void PrintMarg(BYTE width) { cout << setw(width) << SPACE; }

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

		static void PrintGapsMarg() { PrintMarg(margD_ + GapsW + margGexcl + GapsWexcl); }

		ChromView(Gr::eType grType) : GrType(grType) {}

		// Prints chrom's header: Reads
		//	@param FgHeader: true for foreground header
		//	@returns: header width
		int PrintHeader(bool FgHeader);

		// Prints chrom's header: gaps statistics
		static int PrintHeaderGaps();

		// Print info about recorded Reads (count, sample, ampl, density) for ground type of this nstance
		//	@param gMode: gen mode for which info is printed
		//  @param fragCnt: Reads counters
		//	@param densLen: length on which the density is determined
		void PrintReads(GM::eMode gMode, const FragCnt fragCnt, chrlen densLen);

		// Set task mode
		//static void SetMode() { margC_R -= (gMode = BYTE(TestMode)); }

		// Initializes viewing data
		//	@param maxCnt: estimated maximum number of generated reads (FG or BG)
		//	@param sample: sampling value
		//	@param maxDens: estimated maximum read density (FG or BG)
		void Init(size_t maxCnt, float sample, float maxDens);

#ifdef DEBUG
		void Print(ULONG maxCnt);
#endif
	};

	struct Context {
		cells_cnt	CellCnt;			// count of cells
		float		Sample[Gr::Cnt];	// user-defined samples: [0] - fg, [1] - bg
		FragCnt		fCnts[Gr::Cnt];		// fragment counters: [0] - fg, [1] - bg

		// Gets total number of recorded frags
		size_t RecCnt()	const{ return fCnts[Gr::FG].RecCnt() + fCnts[Gr::BG].RecCnt(); }

		// Sets sample for both grounds
		void SetSample(float sample) { Sample[Gr::FG] = Sample[Gr::BG] = sample; }

		// Returns exact number of 'background' cells according to background sample
		float GetExactBGCellCnt() const { return Sample[Gr::BG] * CellCnt; }

		// Sets Control number of cells and sample according to exact count of Test 'background' cells
		//	@param exactCellCnt: exact count of Test 'background' cells
		void SetControlSample(float exactCellCnt) {
			Sample[Gr::BG] = exactCellCnt/(CellCnt = cells_cnt(ceil(exactCellCnt)));
		}

		// Thread-safely increments counters of local and total recorded fragments
		//	@g: ground
		//	@primer: true if increment amplified (second) frag's counter
		//	@returns: true if Reads limit is exceeded.
		bool IncrRecFragCount(Gr::eType g, bool primer) {
			return fCnts[g].RecIncrSaved(primer) + fCnts[!g].RecCnt() >= SeqMode::FragsLimit();
		}
	};

	static GenomeSizes gSizes;		// genome total sizes; needed to define average % of gaps, excl gaps
	static ChromView ChrView[];		// FG, BF chrom views
	static ULONG	TreatedLen[];	// FG, BF genome treated length; needed to define total density
	static Context	GlobContext[];	// TestMode, ControlMode global generation context
	static float	AutoSample;		// adjusted FG sample to stay in limit
	static float	FragMean;		// mean length of selected fragments
	static a_coeff	PCRCoeff;		// user-stated amplification coefficient
	static eVerb	Verb;			// verbose level
	static bool		MakeControl;	// true if control file (input) should be produced
	static Imitator	*Imit;			// singletone instance: to call threads only
	static const Features *Templ;	// template features or NULL (control mode)
	//static readlen	BindLen;	// binding length

	const ChromSizesExt& _cSizes;	// chrom sizes
	DataWriter& _writer;			// data writer

	// Returns stated count of cells
	//	@param gm: generated mode
	static cells_cnt CellCnt(GM::eMode gm) { return GlobContext[int(gm)].CellCnt; }

	// Prints chromosome's name
	//	@param cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@param gm: generation mode Test|Control
	//	@param print: true if chromosomes name should be printed
	static void PrintChromName(chrid cID, GM::eMode gm, bool print);

	static void PrintReadInfo(
		Gr::eType gr, GM::eMode gMode, const FragCnt fCnts[], const ULONG rgnLens[])
	{ ChrView[gr].PrintReads(gMode, fCnts[gr], rgnLens[gr]); }

	// Prints chromosome's name and full statistics
	//	@param cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
	//	@param gMode: Test|Control
	//	@param fCnts: array of fragment's counters, for FG and BG
	//	@param rgnLens: array of region's lengths to print FG|BG density
	//	@param prChrName: true if chromosome's name should be printed
	static void PrintChromInfo(chrid cID, GM::eMode gMode,
		const FragCnt fCnts[], const ULONG rgnLens[], bool prChrName=true);

	// Prints header (FG and BG, or single) and solid line
	//	@param header: if false, then print solid line only
	static void PrintHeader(bool header);

	// Prints total outcome
	static void PrintTotal();

	// Returns the increasing MDA factor (for the purpose of estimating the number of reads during amplification)
	//	@param getMinGragLen: lambda returning minimum fragmemt length
	static float GetMDAfactor(function<readlen(ChromCutter::FragDistr&)> getMinFragLen);

	// Gets estimated number of reads generated per chrom and corrects maximum count and density.
	// Including MDA, excluding PCR
	//	@param gr: ground
	//	@param effLen: actual length on which density is defined
	//	@param factor[in,out]: CellsCnt / recAvrLen
	//	@param numeric: 1 for diploid chromosomes, 0 for haploid ones
	//	@param maxCnt[in,out]: FG|BG reads maximum counters to fill
	//	@param maxDens[in,out]: FG|BG maximum read densities to fill
	//	@returns: estimated number of reads generated per chrom
	static size_t GetEstimReadsCnt(Gr::eType gr, chrlen effLen,
		float factor, BYTE numeric, size_t maxCnt[], float maxDens[]);

	// Increments grounds total length
	//	@param seq: target sequence
	//	@param enrRgnLen: total enriged regions length
	static void IncrementTotalLength(const ChromSeq& seq, chrlen enrRgnLen);
	
	// Returns stated sample for given ground
	//	@param gm: generated mode
	//	@param gr: ground
	static float Sample(GM::eMode gm, Gr::eType gr) { return GlobContext[int(gm)].Sample[gr]; }

	//// Increments counter of total recorded fragments thread-safely.
	////	@g: ground
	////	@primer: if true then increment amplified (second) frag's counter
	////	return: true if limit is exceeded.
	//static bool IncrementRecFragCount(Gr::eType g, bool primer) {
	//	ULLONG cnt = ChromCutter::IsAmpl() ?
	//		TotalFragCnt[g].InterlockedIncrRec(primer):
	//		DataWriter::SingleThread ?
	//			TotalFragCnt[g].PrimerRec++:
	//			InterlockedIncrement(&(TotalFragCnt[g].PrimerRec));
	//	return cnt + TotalFragCnt[!g].Rec() >= SeqMode::FragsLimit();
	//}

	// Initializes Reads view
	//	@param g: ground
	//	@param maxCnt: FG|BG estimated maximum number of generated reads
	//	@param maxDens: FG|BG estimated maximum read density
	static void InitReadsView(Gr::eType g, size_t maxCnt[], float maxDens[]);

	// Curs genome into fragments and generate output
	void	CutGenome	();
	
	// Cuts genome into fragments and generate output
	//	@param arg: pointer to ChrSubset - set of chrom IDs treated in this thread
	//	@param master: if true then this a master instance
	void CutChrom(const effPartition::Subset& arg, bool master)
	{ ChromCutter(this, NULL, master).Execute(arg); }

	// Sets adjusted Sample and clear all counter and means.
	// Sample are needed for the control of BF&FG levels by percent (given by user),
	// and to prorate number of written reads for each chromosome depending on reads limit.
	void SetSample	();

public:
	static bool		UniScore;	// true if template features scores are ignored
	static BYTE		ThrCnt;		// actual number of threads
	static bool		IsExo;
	static bool		IsMDA;
	static eMode	TMode;		// current task mode
	static short	FlatLen;	// BS edge flattening length
	static fraglen	FragLenMin;	// minimal length of selected fragments:
								// established by --frag-dev or Read::FixedLen if size filter is OFF
	static bool All;	// true if total genome is treated
						// false only in Test mode if custom chrom and BG_ALL is false.
	
	// Print amplification info
	//	@param signOut: output marker
	static void PrintAmpl(const char* signOut);

	// Returns true if control is set
	static bool IsControl()			{ return MakeControl; }

	// Returns true if single thread is set
	static bool IsSingleThread()	{ return ThrCnt == 1; }

	// Returns true if given verbose level is active
	static bool Verbose(eVerb level){ return Imitator::Verb >= level; }

	// Set number of threads
	static void SetThreadNumb(BYTE thrNumb) { 
		FragCnt::Init((ThrCnt= thrNumb) == 1);
		ReadWriter::MultiThread = thrNumb > 1;
	}

	// Initializes static values
	//	@param tmode: task mode
	//	@param input: true if control should be generated as well
	//	@param cellsCnt: count of cells
	//	@param isExo: true if EXO mode
	//	@param isReadLenAssigned: true if Read length is assigned by user
	//	@param isMDA: true if MDA is assigned
	//	@param amplCoeff: coefficient of PCR
	//	@param verb: verbosity level
	//	@param allBg: true if all background mode is assigned
	//	@param unstBindLen: unstable binding length
	static void	Init(
		eMode		tmode,
		bool		input,
		cells_cnt	cellsCnt,
		bool		isExo,
		bool		isReadLenAssigned,
		bool		isMDA,
		a_coeff		amplCoeff,
		UINT		verb,
		bool		allBg,
		//bool	uniformScore,		// true if uniform template score is assigned
		//readlen bindLen,
		//const pairVal& flattens
		UINT		unstBindLen			// unstable binding length
	);

	// Creates singleton instance.
	//	@param ocSizes: chrom sizes
	//  @param writer: output data writer
	Imitator(const ChromSizesExt& cSizes, DataWriter& writer)
		: _cSizes(cSizes), _writer(writer) { Imit = this; }

	// Runs task in current mode and write result to output files
	//	@param templ: input template or NULL
	void Execute(Features* templ);
};
