/**********************************************************
Imitator.cpp
Provides chip-seq imitation functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 01/04/2025
***********************************************************/

#include "isChIP.h"
#include "Imitator.h"
#include "Options.h"
#include <algorithm>    // std::sort
#include <cfloat>		// FLT_MAX
#include <cwchar>		// long '-'
#include <random>		// std::exponential_distribution
#include <thread>
#include <fstream>

const char* Gr::title[] = { "FG","BG" };	// if change, correct TitleLength

/************************  class Average ************************/
UINT Imitator::Average::operator+=(UINT val)
{
	static const char* cAverage = "Average";
	if (_keep)
		if (_summator < (ULLONG_MAX - val)) {
			_summator += val;
			_count++;
			_keep = _count < ULONG_MAX;		//  cut off treatment due to counter overflow
			if (!_keep && Imitator::Verbose(eVerb::DBG))
				Err("block count", cAverage).Throw();
		}
		else {
			_keep = false;					//  cut off treatment due to summator overflow
			if (Imitator::Verbose(eVerb::DBG))
				Err("block summation", cAverage).Throw();
		}
	return val;
}
/************************  end of class Average ************************/

/************************ FragDistr ************************/

void Imitator::ChromCutter::FragDistr::SizeSelLimits(fraglen& min, fraglen& max)
{
	const float ssDev = ssFactor0 * float(sqrt(log(ssFactor1 / DRand())));

	min = fraglen(DistrParams::ssMean - ssDev);
	if (min < Read::FixedLen)	min = Read::FixedLen;
	max = fraglen(DistrParams::ssMean + ssDev);
}

/************************ end of FragDistr ************************/

/************************  Amplification ************************/

void Imitator::ChromCutter::MDA::Split(fraglen shift, fraglen len)
{
	if (len < _fLenMin)		return;
	emplace_back(shift, len);
	auto len1 = fraglen(_rng.Range(len));	// new right fraction relative position
	Split(shift, len1 - 1);					// split left fraction
	Split(shift + len1, len - len1);		// split right fraction
}

void Imitator::ChromCutter::MDA::Generate(fraglen fLen, fraglen fLenMin)
{
	clear();
	if (IsMDA) {
		_fLenMin = fLenMin;
		Split(0, fLen);
	}
	else if (fLen >= fLenMin)
		emplace_back(0, fLen);		// a single member: original fragment

	//printf("MDA %d  min %d\n", size(), fLenMin);
	//for (Fraction f : *this)	printf("%d\t%d\n", f.first, f.second);
}

/************************  Amplification: end ************************/

void Imitator::GenomeSizes::IncrSizes(const ChromSeq& seq)
{
	Mutex::Lock(Mutex::eType::INCR_SUM);
	Total += seq.Length();
	Defined += seq.DefRegion().Length();
	Gaps += seq.GapLen();
	Mutex::Unlock(Mutex::eType::INCR_SUM);
}

/************************ ChromView ************************/

/*
	difference between cout & printf precision:

	float v1 = 0.1234, v2 = 0.001234;
	int prec = 2;
	cout << setprecision(prec) << v1 << TAB << v2 << LF;
	printf("%.*f\t%.*f\n", prec, v1, prec, v2);

	== output
	0.12    0.0012
	0.12    0.00
*/

#define cFIL_VAL SPACE
const char* Imitator::ChromView::tGapsExcl = "g_excl";
const char* Imitator::ChromView::tTime = "time";	// "mm:ss"
BYTE Imitator::ChromView::GapsWexcl;

// Prints str on given field and return field width
//	@width: field width
int PrFittedStr(const char* str, int width)
{
	cout << setw(width) << str;
	return width;
}

// Prints unzero float val or nothing on given field
//	@width: field width
void PrFittedInt(ULONG val, int width)
{
	cout << right << setw(width);
	if (val)	cout << val;
	else	cout << cFIL_VAL;
}

// Prints unzero float value or nothing on given field
//	@percent: if true then print percent
//	@width: field width
//	@padd: right padding - space to the right of the value
void PrFittedFloat(bool percent, float val, int precision, int width, int padd = 0)
{
	if (val) {
		if (percent)	printf("%s", PercentToStr(val, precision, width - padd).c_str());
		else		printf("%*.*f", width - padd, precision, val);
		if (padd)	printf("%c", cFIL_VAL);
	}
	else	printf("%*c", width, cFIL_VAL);
}

int Imitator::ChromView::PrintHeader(bool FgHeader)
{
	int lineW = 0;
	// shift title to right in case of little reads count
	const string& rTitle = FT::ItemTitle(FT::ABED);
	const bool w = CountW < rTitle.length();

	cout << right;
	// if little reads count then extend "reads" field to 1, otherwhise extend "sample" field to 1
	lineW += PrFittedStr(rTitle.c_str(),
		GrTitleW() + CountW - w + marg_R * ControlMode);	// "reads"
	lineW += PrFittedStr("sample", marg_S + SampleW);		// "sample"
	// now we should compensate extention to 1:
	// if MDA is set then shrink "MDAc" to 1, otherwhise shrink "unit densuty" field to 1
	if (Imitator::IsMDA || PCRCoeff)
		lineW += PrFittedStr("ampl", marg_A + AmplW);		// amplification coefficient 
	lineW += PrFittedStr("r/kbp", marg_D + DensW);			// unit density
	if (FgHeader)	PrintMarg(margD_), lineW += margD_;		// empty space

	return lineW;
}

int Imitator::ChromView::PrintHeaderGaps()
{
	int lineW = 0;

	if (Verbose(eVerb::PAR)) {
		lineW = PrFittedStr("gaps", margD_ + GapsW);
		if (!ChromSeq::LetGaps)
			lineW += PrFittedStr(tGapsExcl, margGexcl + GapsWexcl);
	}
	if (Timer::Enabled)
		lineW += PrFittedStr(tTime, marg_T + int(strlen(tTime)));
	return lineW;
}

void Imitator::ChromView::PrintReads(GM::eMode gMode, const FragCnt fragCnt, chrlen densLen)
{
	ULONG rCnt = ULONG(fragCnt.RecCnt() << SeqMode::Mode());	// count of reads

	if (TestMode)								// Ground title
		if (gMode == GM::eMode::Control)		// can be true only if Imitator::MakeControl is true
			PrintMarg(Gr::TitleLength + 1);		// print blanks instead of Gr title
		else
			cout << right << Gr::Title(GrType) << COLON;	// print Gr title
	PrFittedInt(rCnt, marg_R + CountW);												// count of Reads
	PrFittedFloat(true, fragCnt.Sample(), SamplePr, marg_S + SampleW, SampleWRP);	// sample
	if (Imitator::IsMDA || PCRCoeff)
		PrFittedFloat(false, fragCnt.RealAmplCoeff(), AmplPr, marg_A + AmplW);			// MDA coef
	PrFittedFloat(false, LinearDens(rCnt, densLen), DensPr, marg_D + DensW, DensWRP);	// density
	if (GrType == Gr::FG)
		PrintMarg(margD_);
}

void Imitator::ChromView::PrintGaps(const GenomeSizes& s)
{
	PrFittedFloat(true, s.GapsInPers(), GapsPr, margD_ + GapsW);
	if (ChromSeq::LetGaps)		return;
	PrFittedFloat(true, s.UndefInPers(), GapsPr, margGexcl + int(strlen(tGapsExcl)));
}

void Imitator::ChromView::Init(size_t maxCnt, float sample, float maxDens)
{
	CountW = DigitsCountLocale(maxCnt, Options::GetBVal(oLOCALE)) + 1;	// +1 for total
	if (CountW < FT::ItemTitle(FT::ABED).length())
		CountW = BYTE(FT::ItemTitle(FT::ABED).length());

	//if(sample <= 0.1)		SampleW = 6;	// 0.099% or <0.01%
	//else if(sample <= 1.1)	SampleW = 5;	// 0.99%
	//else if(sample < 10.1)	SampleW = 4;	// 9.9%
	//else if(sample < 100)	SampleW = 3;	// 99%
	//else					SampleW = 4;	// 100%
	SampleWRP = BYTE(sample * 100 > 9.9);		// '9.9% '
	DensWRP = 0;
	if (maxDens < 0.015)		DensPr = 3;				// 0.0012
	else if (maxDens < 0.15)	DensPr = 3;				// 0.012
	else if (maxDens < 1.5)	DensPr = 2;				// 0.123
	else if (maxDens < 15)	DensPr = 1, DensWRP = 1;// 1.23
	else if (maxDens < 150)	DensPr = 0;				// 123.4
	else					DensPr = 0;				// 1234.5

	GapsWexcl = BYTE(strlen(tGapsExcl));
}

#ifdef DEBUG
void Imitator::ChromView::Print(ULONG maxCnt)
{
	cout << Gr::Title(GrType) << SepCl;
	cout << "MaxCnt = " << setw(6) << setfill(SPACE) << maxCnt << SepCm;
	cout << "CountW = " << int(CountW) << SepCm;
	cout << "SampleW = " << int(SampleW) << SepCm;
	cout << "DensPr = " << int(DensPr) << SepCm;
	cout << "DensW = " << int(DensW) << SepCm;
}
#endif

/************************ ChromView: end ************************/

// pointer to the 'Read counter increment' method
Imitator::FragCnt::pRecIncr	Imitator::FragCnt::pRecIncrSaved;
// pointer to the thread-saved 'selected frag's number adding' method
Imitator::FragCnt::pSelAdd	Imitator::FragCnt::pSelAddSaved;

/************************ class ChromCutter ************************/

fraglen	Imitator::ChromCutter::_SsDev;			// deviation of frag size selection
a_cycle	Imitator::ChromCutter::_PCRdcycles = 1;	// PCR cycles: read doubling cycles

// Sets global mode
void Imitator::ChromCutter::SetGMode(GM::eMode gmode)
{
	_gMode = gmode;
	_fragCnt.SetGMode(gmode);
	_writer->SetGMode(gmode);
}

// Increments counters of local and total recorded fragments thread-safely
//	@g: ground
//	@primer: true if increment derived (amplified) frag's counter
//	return: true if Reads limit is exceeded.
bool Imitator::ChromCutter::IncrRecFragCount(Gr::eType g, bool primer)
{
	_fragCnt[g].RecIncr(primer);	// incr of local recorded Reads
	return GlobContext[int(_gMode)].IncrRecFragCount(g, primer);
}

// Increments counter of total selected fragments thread-safely
void Imitator::ChromCutter::IncrTotalSelFragCount()
{
	for (BYTE i = 0; i < Gr::Cnt; i++)
		GlobContext[int(_gMode)].fCnts[i].SelAddSaved(_fragCnt[i].SelCnt());
}

Imitator::ChromCutter::ChromCutter(const Imitator* imitator, Average* avr, bool master) :
	_cSizes(imitator->_cSizes),
	_ampl(_fragDistr),
	_fragDistr(avr),
	_master(master),
	_gMode(GM::eMode::Test)
{
	DataWriter::SetSeqMode(/*avr*/);
	_fragCnt.Clear();
	_fragCnt.SetGMode(GM::eMode::Test);
	_writer = master ? &(imitator->_writer) : new DataWriter(imitator->_writer);
}

// Prints thread-safe info about treated chroms and stops timer
//  @seq: current reference chromosome
//	@enrRegLen: length of all enriched regions
//	@timer: current timer to thread-saves time output or NULL
//	@exceedLimit: true if limit is exceeded
void Imitator::ChromCutter::PrintChrom(
	const ChromSeq& seq, chrlen enrRegLen, Timer& timer, bool excLimit)
{
	if (!Verbose(eVerb::RT))	return;

	Mutex::Lock(Mutex::eType::OUTPUT);

	const ULONG rgnLens[] = { enrRegLen, seq.DefRegion().Length() - enrRegLen }; // FG, BG region's lengths
	PrintChromInfo(seq.ID(), _gMode, _fragCnt.GetFragCnts(), rgnLens, !IsSingleThread());
	if (Verbose(eVerb::PAR))
		if (_gMode == GM::eMode::Test) {
			const GenomeSizes s(seq);
			ChromView::PrintGaps(s);
		}
		else if (Timer::Enabled)	ChromView::PrintGapsMarg();
	timer.Stop(ChromView::marg_T, false, false);		// print time
	if (excLimit) {
		cout << " exceeded limit";
		if (!Verbose(eVerb::PAR))
			cout << " of " << SeqMode::ReadsLimit() << SPACE << FT::ItemTitle(FT::ABED, true);
	}
	cout << endl;

	Mutex::Unlock(Mutex::eType::OUTPUT);
}

cells_cnt Imitator::ChromCutter::Init(GM::eMode gm, chrid cID, Timer& timer)
{
	SetGMode(gm);
	PrintChromName(cID, gm, IsSingleThread());			// print chrom name before cutting
	timer.Start();
	return CellCnt(gm) << cells_cnt(Chrom::IsAutosome(cID));	// multiply twice for autosomes
}

void Imitator::ChromCutter::Execute(const effPartition::Subset& cIDSet)
{
	try {
		Timer	timer(Verbose(eVerb::RT));	// print local time on Verbose 'runtime info and above'

		for (chrid cID : cIDSet.NumbIDs()) {	// loop through chroms
			const ChromSeq seq(cID, _cSizes);
			const chrlen cLen = seq.End();		// chrom 'end' position
			const auto cellCnt = Init(GM::eMode::Test, cID, timer);
			float scores[]{ 1,1 };
			Features::cIter	cit;				// template chrom's iterator
			chrlen	fCnt = 0, enrRegLen = 0;	// count of features, length of enriched regions
			int		res = 0;					// result of cutting

			if (Templ && (cit = Templ->GetIter(cID)) != Templ->cEnd()) {
				fCnt = chrlen(Templ->ItemsCount(cID));
				enrRegLen = Templ->EnrRegnLength(cit, 0, FragMean);
			}
			_fragCnt.Clear();
			_writer->BeginWriteChrom(seq);
			for (auto n = 0; n < cellCnt; n++) {
				chrlen currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
				for (chrlen k = 0; k < fCnt; k++)
					if (res = CutChrom(cLen, currPos, Templ->Feature(cit, k), scores, false))
						goto A;			// achievement of Reads limit
				// add background after last 'end' position
				if ((res = CutChrom(cLen, currPos, seq.DefRegion(), scores, true)) < 0)
					break;				// achievement of Reads limit
			}
		A:	PrintChrom(seq, enrRegLen, timer, res < 0);		// timer stops and printed in here
			IncrTotalSelFragCount();
			// collect total enriched regions length to calculate total density
			IncrementTotalLength(seq, enrRegLen);
			if (MakeControl) {
				const auto cellCnt = Init(GM::eMode::Control, seq.ID(), timer);
				for (auto n = 0; n < cellCnt; n++) {
					chrlen currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
					CutChrom(cLen, currPos, seq.DefRegion(), scores, true);
				}
				PrintChrom(seq, enrRegLen, timer, false);		// timer stops in here
				IncrTotalSelFragCount();
			}
			_writer->EndWriteChrom();
			if (res < 0)		break;			// achievement of Reads limit
		}
	}
	catch (const Err& e) { Terminate(cIDSet.ID(), e.what()); }
	catch (const exception& e) { Terminate(cIDSet.ID(), e.what()); }
	//catch(...)					{ Terminate(cIDSet.ID(), "Unregistered error in thread"); }
	if (!IsSingleThread() && Verbose(eVerb::DBG)) {
		Mutex::Lock(Mutex::eType::OUTPUT);
		cout << SignDbg << sThread << int(cIDSet.ID()) << ":  end" << endl;
		Mutex::Unlock(Mutex::eType::OUTPUT);
	}
}

typedef pair<fraglen, fraglen>	frag;

// Returns sample of Flattening of binding site suburb
//	@fStart: tested frag's start
//	@fEnd: tested frag's end
//	@ft: binding site
//	@sample: corrected sample
void Imitator::ChromCutter::GetFlattSample(chrlen fStart, chrlen fEnd, const Featr& ft, bool& sample)
{
	short uZone = short(fEnd - ft.Start);	// ustable zone, 
											// with increasing probability of frag binding
	if (uZone > FlatLen)					// is frags end in unstable zone?
		uZone = short(ft.End - fStart);		// is frags start in unstable zone?
	if (uZone <= FlatLen)					// is frags start or end in unstable zone?
		sample = _fragDistr.Sample(float(uZone) / FlatLen);	// new unstable select
}

int Imitator::ChromCutter::CutChrom(
	chrlen cLen,
	chrlen& fStart,
	const Featr& ft,
	float scores[],
	bool bg,
	FragLenStat* fStat
)
{
	fraglen	fLenMin = Read::FixedLen;	// minimal fragment's length after size selection
	fraglen	fLenMax = FRAG_MAX;			// maximal fragment's length after size selection

	scores[0] = ft.Value;
	for (fraglen fEnd = 0; fStart <= ft.End; fStart = fEnd + 1) {	// ChIP: control right mark
		chrlen fLen, fTreatedEnd;

		if (DistrParams::IsLn()) {
			fLen = _fragDistr.LognormNext();
			fTreatedEnd = fEnd = fStart + fLen;							// fragment's end position
		}
		else {
			fLen = chrlen(DistrParams::lnMean);
			fTreatedEnd = fStart + fLen;
			fEnd = fTreatedEnd + _fragDistr.Range(fLen);
		}

		if (fEnd > cLen)
			if (fStart >= cLen - Read::FixedLen)	return 1;	// end of chrom
			else	fLen = (fTreatedEnd = fEnd = cLen) - fStart;				// cut last fragment

		if (DistrParams::IsSS())
			_fragDistr.SizeSelLimits(fLenMin, fLenMax);			// get next size selection limits

		/*
		 * Since the lower limit of the fragment length after size selection remains unchanged
		 * after amplification, in order to increase the efficiency, it is cut off immediately.
		 * The upper limit is checked after amplification, since the fragment length can decrease with MDA.
		 */

		 //== size selection check 1: skip short fragment
		if (fLen < fLenMin)	continue;
		/*
		 * control left mark:
		 * TestMode: foreground (g==0) is inside and background (g==1) outside template features;
		 * for other chromosomes background (g==0) is inside feature==chrom's length.
		 * ControlMode: foreground (g==0) is always inside feature==chrom's length
		 */
		//Gr::eType g = Gr::eType(!bg ^ (fEnd >= ft.Start));
		Gr::eType g = Gr::eType(!bg ^ (fTreatedEnd >= ft.Start));

		//== EXO processing
		// frags[0],frags[1] - sheared fragments of forward (0) and backward (1) strand DNA,
		// frags[2],frags[3] - additional fragments generated by EXO
		//Region frags[4]{ {fStart,fEnd},{fStart,fEnd},{0,0},{0,0} };		// set frags[0] and frags[1] identical
		Region frags[4]{ {fStart,fTreatedEnd},{fStart,fTreatedEnd},{0,0},{0,0} };		// set frags[0] and frags[1] identical

		if (IsExo) {
			int diff = ft.Start - fStart - _fragDistr.Expo();	// left difference
			if (diff > 0)	frags[0].Start += diff;

			//diff = fEnd - ft.End - _fragDistr.Expo();			// right difference
			diff = fTreatedEnd - ft.End - _fragDistr.Expo();			// right difference
			if (diff > 0)	frags[1].End -= diff;

			frags[2] = frags[1];
			frags[3] = frags[0];
		}

		//== flattening
		auto select = bool(g);				// selection by corrected bounds
		if (!select							// foreground; always true for BG
			//&& (select = fEnd >= ft.Start)	// does the fragment cover the feature?
			&& (select = fTreatedEnd >= ft.Start)	// does the fragment cover the feature?
			&& FlatLen)						// is flattening ON?
			//GetFlattSample(fStart, fEnd, ft, select);
			GetFlattSample(fStart, fTreatedEnd, ft, select);

		//== sequencing
		if (select && Sample(scores[g]))	// selection by bounds and feature score
			for (int x = (2 << int(IsExo)) - 1; x >= 0; x--)	// loop through the fragments from frags[]: 2 (no EXO) or 4 (EXO)
				if (Sample(Imitator::Sample(_gMode, g))
					&& frags[x].Length() > Read::FixedLen) {		// FG/BG loss and not short fragment after EXO
					// ** MDA amplification
					bool primer = true;
					_ampl.Generate(frags[x].Length(), fLenMin);
					for (const Fraction& frac : _ampl)
						if (Sample(AutoSample)				// adjusted limits sample
							&& frac.second <= fLenMax)		// ** size selection 2: skip long fragments
							// ** PCR amplification: _PCRdcycles is number of read doubling cycles
							for (a_cycle i = 0; i < _PCRdcycles; primer = false, i++)
								// while BG or without MDA frac.first is always 0
								if (!_writer->AddRead(frags[x].Start + frac.first, frac.second, x % 2)	// split identical fragments into forward and reverse stranded 
									&& IncrRecFragCount(g, primer))		// recorded reads
									return -1;							// Reads limit is exceeded
				}
		//== size selection check 2: statistics record
		if (fLen <= fLenMax) {
			_fragCnt[g].SelIncr();
			if (fStat)	fStat->AddFragLen(fLen);
		}
	}
	return 0;
}

/************************ end of class ChromCutter ************************/

/************************  class Imitator ************************/

Imitator::GenomeSizes Imitator::gSizes;
Imitator::Context	Imitator::GlobContext[2];	// global generation context
Imitator::ChromView Imitator::ChrView[] = { Gr::FG,Gr::BG };
ULONG	Imitator::TreatedLen[] = { 0,0 };	// FG, BF genome treated length
float	Imitator::AutoSample = 1;		// adjusted FG sample to stay in limit
//readlen	Imitator::BindLen;
short	Imitator::FlatLen = 0;		// BS edge flattening length
a_coeff	Imitator::PCRCoeff = 0;		// user-stated amplification coefficient
eVerb	Imitator::Verb;
BYTE	Imitator::ThrCnt;			// actual number of threads
bool	Imitator::IsExo;
bool	Imitator::IsMDA;
bool	Imitator::MakeControl;		// true if control file (input) should be produced
bool	Imitator::UniScore;
bool	Imitator::All;
eMode	Imitator::TMode;			// Current task mode
float	Imitator::FragMean;			// mean length of selected fragments
Imitator* Imitator::Imit = NULL;
const Features* Imitator::Templ = NULL;

void Imitator::PrintChromName(chrid cID, GM::eMode gm, bool print)
{
	if (Verbose(eVerb::RT) && print) {
		if (MakeControl) 	cout << *GM::Title(gm) << SPACE;
		cout << setw(ChromView::ChromNameW()) << left << setfill(SPACE)
			<< (cID == Chrom::UnID ? sTotal : Chrom::AbbrName(cID, true)) + COLON;
		fflush(stdout);		// including reset to default right and setfill
	}
}

void Imitator::PrintChromInfo(
	chrid cID, GM::eMode gMode, const FragCnt fCnts[], const ULONG rgnLens[], bool prChrName)
{
	PrintChromName(cID, gMode, prChrName);
	if (TestMode)
		PrintReadInfo(Gr::FG, gMode, fCnts, rgnLens);
	PrintReadInfo(Gr::BG, gMode, fCnts, rgnLens);
}

void Imitator::PrintHeader(bool header)
{
	if (!Verbose(eVerb::RT))	return;
	static int w = 0;	// width

	if (header) {
		w = 0;
		if (Verbose(eVerb::PAR))	cout << LF;
		if (MakeControl)
			w += PrFittedStr(sSPACE, 2);	// "t " or "c "
		cout << setfill(SPACE) << left;
		w += PrFittedStr(Chrom::Short.c_str(), ChromView::ChromNameW());
		if (TestMode)
			w += ChrView[Gr::FG].PrintHeader(true);
		w += ChrView[Gr::BG].PrintHeader(false);
		w += ChromView::PrintHeaderGaps();
		cout << LF;
	}
	PrintSolidLine(w + 1);
}

void Imitator::PrintTotal()
{
	cout << sTotal << " recorded "; DataWriter::PrintItemTitle();
	if (MakeControl)		// add "test:"
		cout << SepCl << GM::Title(GM::eMode::Test) << COLON;
	cout << SPACE;
	DataWriter::PrintItemCount(GlobContext[int(GM::eMode::Test)].RecCnt());
	if (MakeControl) {	// add "control:"
		cout << SepCm << GM::Title(GM::eMode::Control) << SepCl;
		DataWriter::PrintItemCount(GlobContext[int(GM::eMode::Control)].RecCnt());
	}
	cout << endl;		// flash cout buffer
}

float Imitator::GetMDAfactor(function<readlen(ChromCutter::FragDistr&)> getMinFragLen)
{
	const uint16_t cycles = 100;
	const fraglen lnMean = fraglen(DistrParams::LnMean());
	ChromCutter::FragDistr	fragDistr(nullptr);
	ChromCutter::MDA		ampl(fragDistr);
	size_t amplCnt = 0;

	for (uint16_t i = 0; i < cycles; i++) {
		ampl.Generate(lnMean, getMinFragLen(fragDistr));
		amplCnt += ampl.size();
	}
	//float coeff = float(amplCnt) / cycles;
	//printf(">>>DistrParams::lnMean: %d  amplCnt: %d  coeff: %f\n", lnMean, amplCnt, coeff);
	//return coeff;
	return float(amplCnt) / cycles;
}

size_t Imitator::GetEstimReadsCnt(Gr::eType gr, chrlen effLen, float factor,
	BYTE numeric, size_t maxCnt[], float maxDens[])
{
	// Reads count
	size_t cnt = 2 * size_t(Sample(GM::eMode::Test, gr) * factor * effLen) << SeqMode::Mode();
	
	//if(DistrParams::IsLn())
		if (DistrParams::IsSS()) {
			/*
			Size selection normal distribution reduces the reads number in inverse proportion to the increase in the Y-value
			at a point equal to the SS normal median, relative to the lognormal distribution.
			This is the ratio by which we need to multiply the 'nested' SS normal distribution
			so that it fits into the initial 'overarching' lognormal one.
			*/
			float ratio =	// from about 0.1 to 0.5
				float(Distrib::GetVal(Distrib::eCType::LNORM, DistrParams::lnMean, DistrParams::lnSigma, fraglen(DistrParams::ssMean))) /
				float(Distrib::GetVal(Distrib::eCType::NORM, DistrParams::ssMean, float(DistrParams::ssSigma), fraglen(DistrParams::ssMean)));

			if (IsMDA)
				ratio *= GetMDAfactor(
					[](ChromCutter::FragDistr& fragDistr) {
						fraglen	fLenMin, fLenMax;	// minimum, maximum fragment's length after size selection

						fragDistr.SizeSelLimits(fLenMin, fLenMax);
						/*
						This is a rough empirical estimate.
						Decreasing the fragment minimum length by read size somehow compensates 
						for the lack of long fragments after the lognormal distribution.
						Table of values of the estimated amplification factor compared to the actually measured one
						for different size selection distributions (by default lognormal parameters):
						mean sigma	est/real
						201	  30	146%	(mean = default lognorm Mean)
						230	  20	101%
						160	  15	122%
						130	  12	82%
						*/
						return fLenMin > Read::FixedLen ? fLenMin - Read::FixedLen : Read::FixedLen;
					}
				);
			cnt = size_t(cnt * ratio);
		}
		else
			if (IsMDA)
				cnt = size_t(cnt * GetMDAfactor([](ChromCutter::FragDistr&) { return Read::FixedLen; }));	//  99% prediction accuracy

	if (maxCnt[gr] < cnt)		maxCnt[gr] = cnt;
	
	// Reads density
	factor = LinearDens(cnt, effLen >> numeric);
	if (maxDens[gr] < factor)	maxDens[gr] = factor;
	
	return cnt;
}

void Imitator::IncrementTotalLength(const ChromSeq& seq, chrlen enrRgnLen)
{
	if (enrRgnLen)	InterlockedExchangeAdd(&(TreatedLen[Gr::FG]), enrRgnLen);
	InterlockedExchangeAdd(&(TreatedLen[Gr::BG]), seq.DefRegion().Length() - enrRgnLen);
	if (Verbose(eVerb::PAR))	gSizes.IncrSizes(seq);
}

void	Imitator::Init(
	eMode		tmode,				// task mode
	bool		input,				// true if control should be generated as well
	cells_cnt	cellsCnt,			// count of cells
	bool		isExo,				// true if EXO mode
	bool		isReadLenAssigned,	// true if Read length is assigned by user
	bool		isMDA,				// true if MDA is assigned
	a_coeff		amplCoeff,			// coefficient of PCR
	UINT		verb,				// verbosity level
	bool		allBg,				// true if all background mode is assigned
	//bool	uniScore,			// true if uniform template score is assigned
	//readlen bindLen,
	//const pairVal& flattens
	fraglen		unstBindLen,			// unstable binding length
	const pairVal&	samples
) {
	TMode = tmode;
	MakeControl = TestMode ? input : false;
	GlobContext[int(GM::eMode::Test)].CellCnt = cellsCnt;
	IsExo = isExo;
	IsMDA = isMDA;
	PCRCoeff = amplCoeff;	// the actual ChromCutter ampl coeff will be set in Sample()
	ChromCutter::FragDistr::Init();
	Verb = eVerb(verb);
	All = (tmode == CONTROL) || allBg;
	//UniScore = uniScore;
	//BindLen = bindLen;
	//FlatLen = flattens.first + flattens.second;
	FlatLen = unstBindLen;
	FragMean = DistrParams::LnMean();	// for the Sample(), before set actual ln mean
	if (DistrParams::IsRVL())			// Read variable length mode is set
		if (!isReadLenAssigned)
			Read::FixedLen = Read::VarMinLen;
		else if (Read::FixedLen > Read::VarMinLen && Verb >= eVerb::CRIT)
			Err(string("assigned read length exceeds the default value of ")
				+ to_string(Read::VarMinLen) + " when read variable mode is activated!").Warning();

	// *** Set Test samples and Control number of cells and sample
	GlobContext[int(GM::eMode::Test)].SetSample(samples.first / 100);
	if (TestMode)
		if (!(GlobContext[int(GM::eMode::Test)].Sample[Gr::BG] *= samples.second / 100))	All = false;
	if (MakeControl)	// can be true only in TEST task mode
		GlobContext[int(GM::eMode::Control)].SetControlSample(GlobContext[int(GM::eMode::Test)].GetExactBGCellCnt());
}

void Imitator::InitReadsView(Gr::eType g, size_t maxCnt[], float maxDens[])
{
	ChrView[g].Init(
		size_t(maxCnt[g] * AutoSample),
		// max() works only in case MakeControl==true and g==Gr::BG,
		// to set maximum BG sample among Test and Control.
		// In other cases Sample(GM::Control, g) always returns 0.
		//max(Sample(GM::Test, g), Sample(GM::Control, g)) * AutoSample,
		Sample(GM::eMode::Test, g) * AutoSample,
		maxDens[g] * AutoSample
	);
}

void Imitator::Execute(Features* templ)
{
	Templ = templ;
	CutGenome();

	// print statistics
	if (Verb == eVerb::RES)		PrintTotal();
	else if (Verbose(eVerb::RT)) {
		if (_cSizes.TreatedCount() > 1) {	// print summary test statistics?
			PrintChromInfo(Chrom::UnID, GM::eMode::Test, GlobContext[int(GM::eMode::Test)].fCnts, TreatedLen);
			if (Verbose(eVerb::PAR))	ChrView[Gr::BG].PrintGaps(gSizes);	// ground doesn't matter
			cout << LF;
		}
		if (TestMode)	PrintTotal();
	}
}

void Imitator::CutGenome()
{
	effPartition cSets(_cSizes, thrid(ThrCnt));

	if (ThrCnt > 1 && Verbose(eVerb::DBG))	cSets.Print();
	if (FlatLen < 0)		FlatLen = -FlatLen;
	SetSample();
	PrintHeader(true);

	bool master = true;
	vector<thread> threads;
	threads.reserve(ThrCnt);
	for (BYTE i = 0; i < ThrCnt; master = false, i++)
		threads.emplace_back(&Imitator::CutChrom, this, cSets[i], master);
	for (thread& t : threads)	t.join();

	PrintHeader(false);
}

void Imitator::SetSample()
{
	// *** Determine the estimated numbers of reads generated 
	size_t	rTotalCnt = 0;			// total estimated number of reads generated
	size_t	rMaxCnt[] = { 0,0 };	// FG|BG maximum estimated number of reads generated per chrom
	float	rMaxDens[] = { 0,0 };	// FG|BG maximum estimated read densities
	const float	countFactor = float(GlobContext[int(GM::eMode::Test)].CellCnt) / FragMean;

	for (ChromSizes::cIter it = _cSizes.cBegin(); it != _cSizes.cEnd(); it++) {
		if (!_cSizes.IsTreated(it))	continue;
		
		chrlen	enrRgnLen = 0;				// length of enriched regions
		// count is estimated according to diploid (numerical) sign,
		// but density is not, because basic length is single!
		if (Templ && Templ->FindChrom(CID(it))) {
			enrRgnLen = Templ->EnrRegnLength(CID(it), 0, FragMean);
			rTotalCnt += GetEstimReadsCnt(Gr::FG, enrRgnLen, countFactor, 0, rMaxCnt, rMaxDens);
		}
		rTotalCnt += GetEstimReadsCnt(
			Gr::BG,
			_cSizes.DefEffLength(it) - enrRgnLen,
			countFactor,
			Chrom::IsAutosome(CID(it)),
			rMaxCnt,
			rMaxDens
		);
	}
	if (PCRCoeff) {
		ChromCutter::SetAmpl();
		rTotalCnt *= size_t(pow(2.f, int(PCRCoeff)));
	}
	// *** Estimate adjusted Sample
	if (rTotalCnt > SeqMode::ReadsLimit())
		AutoSample = SeqMode::ReadsLimit() / rTotalCnt;
	// *** print debug info
	if (Verbose(eVerb::PAR))
		cout << SignPar << "Actual fragments Mean" << SepCl
		<< (DistrParams::IsSS() ? DistrParams::ssMean : FragMean) << endl;
	if (Verbose(eVerb::DBG)) {
		if (MakeControl)
			cout << SignDbg << "Generated " << GM::Title(GM::eMode::Control) << SepDCl
			<< "Count of cells" << Equel << CellCnt(GM::eMode::Control) << SepSCl
			<< "sample" << Equel << PercentToStr(Sample(GM::eMode::Control, Gr::BG) * 100, 2) << LF;
		cout << SignDbg << "Total recorded reads number estimate" << SepCl << rTotalCnt << endl;
	}
	if (AutoSample < 1 && Verbose(eVerb::RES))
		cout << "Added recovery sample = " << setprecision(3) << (AutoSample * 100)
		<< "% due to reads limit of " << SeqMode::ReadsLimit() << endl;

	// *** set Reads statistics params
	if (TestMode)
		InitReadsView(Gr::FG, rMaxCnt, rMaxDens);
	InitReadsView(Gr::BG, rMaxCnt, rMaxDens);
}

void Imitator::PrintAmpl(const char* signOut)
{
	cout << signOut << "Amplification" << SepCl;
	if (IsMDA)	cout << "MDA";
	if (PCRCoeff) {
		if (IsMDA)	cout << SepSCl;
		cout << "PCR cycles" << Equel << int(PCRCoeff);
	}
	else if (!IsMDA)	cout << Booleans[false];
	cout << LF;
}

/************************  end of Imitator ************************/
