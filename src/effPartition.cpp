/**********************************************************
effPartition.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
Last modified: 21.03.2019

Effective Partition Problem (optimization version) solution
Effectively distributes a set of input positive integers into k subsets,
such that the difference between the subset sums is minimized.
***********************************************************/

#include "effPartition.h"
#include <iostream>
#include <iomanip>      // std::setfill, std::setw

using namespace std;

typedef unsigned long long	ULLONG;

#define SUM_MAX	ULONG_MAX		// maximum sum

typedef std::vector<effPartition::Subset>::iterator ss_it;

// Gets count of digist in an integer value
//	return: count of digist
BYTE DigitsCnt (ss_sum val)
{
	BYTE res = 0;
	for(; val; val/=10, res++);
	return res;
}

// ********** IdNumbers **********

// Constructor by chrom sizes
effPartition::IdNumbers::IdNumbers(const ChromSizesExt& cSizes)
{
	reserve(cSizes.TreatedCount());
	for(ChromSizes::cIter it=cSizes.cBegin(); it!=cSizes.cEnd(); it++)
		if(cSizes.IsTreated(it))
			emplace_back(CID(it), cSizes.DefEffLength(it));
}

		// Initializes instance by numbers
void effPartition::IdNumbers::Init(const IdNumbers& numbers)
{
	reserve(numbers.size());
	for (auto it = numbers.begin(); it != numbers.end(); emplace_back(*it++));
}

// Returns minimum accuracy with which the available space of the average sum should be measured
float effPartition::IdNumbers::GetMinUp() const
{
	float res = float(back().Val)/minTol;
	return res > 1 ? res : 1;
}

// Returns average values sum
//	@ssCnt: count of subsets
float effPartition::IdNumbers::AvrSum(ss_id ssCnt) const
{
	float sum = 0;
	for(numb_cit it=begin(); it!=end(); sum += it++->Val);
	return sum/ssCnt;
}

// returns maximum numbers ID
numb_id effPartition::IdNumbers::MaxID() const
{
	numb_id res = 0;
	for(numb_cit it=begin(); it!=end(); it++)
		if(it->Id > res)	res = it->Id;
	return res;
}

// Marks all numbers as unallocated
void effPartition::IdNumbers::Reset() { for(numb_it it=begin(); it!=end(); it++->BinIInd = 0); }

// ********** Subset **********

// Adds number
//	@n: added Id number
void effPartition::Subset::AddNumb(const IdNumber& n)
{
	sumVal += n.Val;
	numbIDs.push_back(n.Id);
}


// Takes into account (charges) number's value
//	@it: number's iterator
//	@binInd: bin's inverse index
inline void effPartition::Subset::TakeVal(const numb_it& it, ss_id binInd)
{
	sumVal += it->Val;
	it->BinIInd = binInd;
}

// Eliminates (discharges) number's value
//	@it: number's iterator
inline void effPartition::Subset::ElimVal(const numb_it& it)
{
	sumVal -= it->Val;
	it->BinIInd = 0;
}

// Prints subset
//	@sumW: maximum sum width (count of digits) to align sums left
//	@valW: maximum number value width (count of digits) to align vals right
//	@prNumbsCnt: maximum count of printed number IDs or 0 if all (by default)
void effPartition::Subset::Print(BYTE sumW, BYTE valW, size_t prNumbsCnt) const
{
	cout << SignDbg << sThread << int(id) << ":\t"
		 //<< "sum " << setw(sumW) << sumVal << "  "
		 << Chrom::Abbr;
	if(prNumbsCnt && prNumbsCnt < numbIDs.size()) {
		for(size_t i=0; i<prNumbsCnt; i++)
			cout << SPACE << setw(valW) << setfill(SPACE) << Chrom::Mark(numbIDs[i]);
		cout << "... (" << numbIDs.size() << ')';
	}
	else
		for(numb_id_cit it=numbIDs.begin(); it!=numbIDs.end(); it++)
			cout << SPACE << setw(valW) << setfill(SPACE) << Chrom::Mark(*it);
	cout << LF;
}

// ********** Result *********

// Initializes instance: creates empty bins with reserved capacity
//	@nCnt: count of numbers
//	@ssCnt: count of subsets
void effPartition::Result::Init(size_t nCnt, ss_id ssCnt)
{
	SumDiff = SUM_MAX;
	Bins.reserve(ssCnt);
	for(ss_id i=0; i<ssCnt; i++)
		Bins.emplace_back(rowlen(nCnt) / ssCnt + 2);
}

// Sorts subsets in descending order and sets subsets ID starting with 1
void effPartition::Result::SetIDs()
{
	// set subset's ID
	ss_id i = 1;		// first subset ID
	for(ss_it it=Bins.begin(); it!=Bins.end(); it++->id = i++);
}

// Clears instance: removes all subsets number IDs and resets the sum of the number values
void effPartition::Result::Clear()
{
	SumDiff = SUM_MAX;
	for(ss_it it=Bins.begin(); it!=Bins.end(); it++->Clear());
}

// Outfits this result and sets subsets ID starting with 1
//	@numbers: numbers to be distributed
//	@ind: index of the first subset: 1 for SGreedy() or 0 for DSTree()
//	@diff: difference between the subset sums
void effPartition::Result::Fill(effPartition::IdNumbers& numbers, char ind, numb_val diff)
{
	ss_id ssCnt = SubsetCount() - ind;	// count of subsets minus first index
	Clear();
	for(const auto& n : numbers)
		Bins[ssCnt - n.BinIInd].AddNumb(n);
	SumDiff = diff;
	SetIDs();
}

// Sets the minimum/maximum subset sums and their difference for this instance
//	@minSum: minimum subset sum to be retrieved
//	@maxSum: maximum subset sum to be retrieved
//	@return: true if sum difference is updated
bool effPartition::Result::SetSumDiff(ss_sum& minSum, ss_sum& maxSum)
{
	ss_cit it = Bins.begin();

	minSum = maxSum = it->sumVal;
	for(it++; it!=Bins.end(); it++)
		if(it->sumVal > maxSum) 		maxSum = it->sumVal;
		else if(it->sumVal < minSum)	minSum = it->sumVal;
	ss_sum newSumDiff = maxSum - minSum;
	if(newSumDiff < SumDiff) {
		SumDiff = newSumDiff;
		return true;
	}
	return false;
}

// Sorts subsets by their sum
//	@ascend: if true then in in ascending order, otherwise in descending order
void effPartition::Result::Sort(bool ascend)
{ 
	sort(Bins.begin(), Bins.end(), ascend ? Subset::SortByAscend : Subset::SortByDescend);
}

// Prints sorted in descending order instance
//	@valW: maximum number value width (count of digits) to align vals right
//	@prNumbsCnt: maximum count of printed number IDs or 0 if all (by default)
void effPartition::Result::Print(BYTE valW, size_t prNumbsCnt) const
{
	if(!Bins.size())	return;
	BYTE sumW = DigitsCnt(Bins[0].sumVal);	// max sum width (count of digits); first sum is the biggest
	for(const auto& b : Bins)
		b.Print(sumW, valW, prNumbsCnt);
}

// ********** Partition *********
typedef void (Partition::*pMethod)(ss_id);

// Encapsulates partition methods
class Partition
{
	static const UINT minCallLimit = 1000000;	// the minimum DSTree() invokes limit
	static const BYTE limFlag = 0x1;			// flag of DSTree() completion by limit
	static const BYTE perfFlag = 0x2;			// flag of DSTree() completion by 'perfect' result

	const float	avrSum;			// average sum value among subsets
	const bool	isAvrFract;		// true if avrSum is not an integer (has a nonzero fractional part)
	const ULLONG callLimit;		// current DSTree() invokes limit
	ULLONG	callCnt;			// counter of SGreedy() iterations or DSTree() invokes
	ss_sum	minSum;				// current minimum sum among subsets
	ss_sum	maxSum;				// current maximum sum among subsets
	ss_sum	sumDiff;			// current maximum difference between subset sums (inaccuracy)
	ss_sum	lastSumDiff;		// last best inaccuracy: used in DSTree()
	ss_sum	standbySumDiff;		// temporarily stored inaccuracy: used in DSTree()
	BYTE complete;				// holder of DSTree() completion flags (signs)
	effPartition::Result&	finalResult;	// final result
	effPartition::Result	currResult;		// current result
	effPartition::IdNumbers& numbers;		// input numbers
	effPartition::IdNumbers	standbyNumbers;	// numbers with the best unsaved maximum sum difference:
											// used in DSTree()
	static pMethod methods[];	// pointers to partition method

    // Raises completion flag
    //	@flag: flag to be raised
	inline void RaiseComplFlag(BYTE flag) { complete |= flag; }

	// Returns true if 'perfect' completion flag is raised
	//	!= 0 to avoid warning C4800
	inline bool IsCompleteByPerfect() const { return (complete & perfFlag) != 0; }	

	// Returns true if current result is the best possible
	inline bool IsResultPerfect() const { return !sumDiff || (isAvrFract && sumDiff == 1); }

	// Sets the best minimum/maximum subset sum and their difference
	//	@res: current result that delivers the sums as candidates for the best
	//	return: true if the more narrow range is achieved
	bool SetRange(effPartition::Result& res) {
		ss_sum resMinSum, resMaxSum;

		if( res.SetSumDiff(resMinSum, resMaxSum) &&	res.SumDiff < sumDiff) {
			sumDiff = res.SumDiff;
			minSum = resMinSum;
			maxSum = resMaxSum;
			return true;
		}
		return false;
	}

	// Performs 'Unconditional Greedy' partition
	//	@ssCnt: count of subsets
	void UGreedy(ss_id ssCnt)
	{
		int i = 0, shift = 1;

		for (const auto& n : numbers) {
			finalResult.Bins[i].AddNumb(n);
			i += shift;
			if (i / ssCnt > 0) { i--; shift = -1; }	// last subset, flip to reverse round order
			else if (i < 0) { i++; shift = 1; }		// first subset, flip to direct  round order
		}
		SetRange(finalResult);		// set minSum, maxSum, sumDiff
		finalResult.SetIDs();		// set subsets ID
	}

	// Performs 'Greedy' partition
	//	@ssCnt: count of subsets
	void Greedy(ss_id ssCnt)
	{
		ss_id i, k;		// subsets cyclic index, index of subset with minimum sum

		for(effPartition::numb_it it=numbers.begin(); it!=numbers.end(); it++) {
			if(it->BinIInd)		continue;		// needs for SGreedy, for Greedy itself is redundant
			minSum = currResult.Bins[k=0].sumVal;
			for(i=1; i<ssCnt; i++)				// loop through the bins from the second one
				if(currResult.Bins[i].sumVal <= minSum)
					minSum = currResult.Bins[k=i].sumVal;
			currResult.Bins[k].TakeVal(it, ssCnt - k);
		}

		if(SetRange(currResult))			// is the current result better than after previous UGreedy() call?
			finalResult.Fill(numbers, 0, currResult.SumDiff);		// outfit final result
	}

    // Wrapper to 'Greedy' partition 
    //	@ssCnt: count of subsets
	void WrapGreedy(ss_id ssCnt)
    {
        numbers.Reset();
		currResult.Clear();
        Greedy(ssCnt);
    }

	// Performs 'Sequential stuffing Greedy' partition
	//	@ssCnt: count of subsets
	void SGreedy(ss_id ssCnt)
	{
		ss_id i;			// reversed bin index started with 1
		ss_it sit;			// bins iterator
		size_t freeCnt;		// count of unallocated numbers
		int k = 1;			// delta multiplicator
		float avrUp;		// raised average sum
		const float up = numbers.GetMinUp();	// delta above average sum
		effPartition::numb_it nit;				// numbers iterator

		// loop through the numbers until count of unallocated numbers becomes less then half count of bins
		do {
			freeCnt = numbers.size();
			numbers.Reset();
			currResult.Clear();
			avrUp = avrSum + up * k++;
			i = ssCnt;
			for(sit=currResult.Bins.begin(); sit!=currResult.Bins.end(); i--, sit++)
				for(nit=numbers.begin(); nit!=numbers.end(); nit++)
					if(nit->IsFitted(avrUp - sit->sumVal)) {
						sit->TakeVal(nit, i);
						freeCnt--;
					}
		}
		// this heuristic contition provided satisfactory inaccuracy in a single pass in most of cases
		while(freeCnt > size_t(ssCnt/2));	
		
		// distribute remaining unallocated numbers by Greed protocol
		// Checking for freeCnt==0 can be omitted, since as nit happens very rarely
		Greedy(ssCnt);
	}

	// Performs 'Dynamic Search Tree' ('perfect') partition
	//	@currnit: number's iterator, from which the cycle continues
	//	@invInd: current inverse index of subset
	void DSTree(effPartition::numb_it currnit, ss_id invInd)
	{
		if(complete)	return;
		effPartition::numb_it nit;
		const ss_it sit = currResult.Bins.end() - 1 - invInd;	// curent bins iterator
		
		if(++callCnt == callLimit)	RaiseComplFlag(limFlag);

		if(invInd) {		// not last bin
			// loop through the numbers starting from the current one
			for(nit=currnit; nit!=numbers.end(); nit++)		
				if(!nit->BinIInd && nit->Val + sit->sumVal < maxSum) {	// accuracy selection
					sit->TakeVal(nit, invInd);				// charge number's value
					if(nit+1 != numbers.end())				// checkup just to avoid blank recursive invoke
						DSTree(nit+1, invInd);				// try to fit next number to the same bin
					if(sit->sumVal > minSum)				// bin is full?
						DSTree(numbers.begin(), invInd-1);	// try to fit unallocated numbers to the next bin
					sit->ElimVal(nit);						// discharge number's value
				}
		}
		else {				// last bin
			// accumulate sum for the last bin
			for(nit=numbers.begin(); nit!=numbers.end(); nit++)
				// zero invIndex means that number belongs to the last bin
				if(!nit->BinIInd)
					sit->sumVal += nit->Val;
			if(SetRange(currResult)) {			// is inaccuracy better than the previous one?
				// keep current number's BinIInds as the standby one
				standbyNumbers.Copy(numbers);
				// keep standbySumDiff for the next standby selection
				lastSumDiff = standbySumDiff = sumDiff;	
				if(IsResultPerfect())		RaiseComplFlag(perfFlag);
			}
			else if(currResult.SumDiff < standbySumDiff) {	// should we keep current result as standby?
				// keep current numbers as the standby one
				standbyNumbers.Copy(numbers);
				standbySumDiff = currResult.SumDiff;
			}
			sit->sumVal = 0;				// clear last bin sum
		}
	}

	// Performs iterative 'Dynamic Search Tree' partition
	//	@ssCnt: count of subsets
	void ISTree(ss_id ssCnt)
	{
		// initial range expansion around average 
		numb_val up = avrSum < numbers[0].Val ? (numbers[0].Val - int(avrSum) + 2) : 1;
		if(up > avrSum)	up = int(avrSum) - 1;
		standbySumDiff = SUM_MAX;		// undefined standby inaccuracy
		lastSumDiff = finalResult.SumDiff;
		do {
			minSum = ss_sum(avrSum - up);
			maxSum = ss_sum(avrSum + up);
			sumDiff = maxSum - minSum;
			complete = 0;
			callCnt = 0;
			numbers.Reset();
			currResult.Clear();
			DSTree(numbers.begin(), ssCnt - 1);
			if(IsCompleteByPerfect() || (up *= 2) >= minSum)	// increase and checkup range expansion
				break;
			if(currResult.SumDiff > standbySumDiff)	// is current inaccuracy worse than standby one?
				break;					// use last fitted result
		}
		while(lastSumDiff != currResult.SumDiff); {	// until previous and current inaccuracy are different
			SetRange(finalResult);
			finalResult.Fill(standbyNumbers, 1, standbySumDiff);
		}

	}

	// Performes partitioning
	//	@ind: index of partition method to call
	//	@ssCnt: count od subsets
	//	return: true if result is 'perfect'
	bool DoPartition(int ind, const ss_id ssCnt)
	{
		(this->*methods[ind])(ssCnt);
		return IsResultPerfect();
	}

public:
	// Creates an instance and performs a partition
	//	@numbers: identified values to be distributed
	//	@res: final result
	//	@avr: average value sum
	//	@limMult: DSTree() call's limit multiplier; if 0 then omit DSTree method invoking
	Partition(effPartition::IdNumbers& numbers, effPartition::Result& res, float avr, UINT limMult)
		: numbers(numbers), finalResult(res), avrSum(avr), sumDiff(SUM_MAX),
		isAvrFract(int(avrSum) != avrSum), callLimit(ULLONG(minCallLimit * limMult))
	{
		const ss_id ssCnt = finalResult.SubsetCount();
		currResult.Init(numbers.size(), ssCnt);		// doesn't need if finished by UGreedy,
													// but this is extremely unlikely
		if(limMult)	standbyNumbers.Init(numbers);	// doesn't need if finished by Greedy or SGreedy,
													// but this is unlikely
        int i = 0;
        methods[i++] = &Partition::UGreedy;
		methods[i++] = &Partition::WrapGreedy;
		methods[i++] = &Partition::SGreedy;
		methods[i++] = &Partition::ISTree;
        int cnt = limMult ? 4 : 3;

        // for the degenerate case numbers.size()<=ssCnt method UGreedy() is comepletely enough
		if(numbers.size() <= ssCnt)	cnt = 1;
        for (i = 0; i < cnt; i++)
            if (DoPartition(i, ssCnt))	return;
	}
};

pMethod Partition::methods[4];

// ********** effPartition *********

// Initializes subsets by partitioned values
//	@vals: numbers to be distributed
//	@ssCnt: count of subsets
//	@perfect: true if partition should be 'perfect'
//	@limMult: DSTree call's limit multiplier; if 0 then omit DSTree method invoking
void effPartition::Init(IdNumbers& numbers, ss_id ssCnt, UINT limMult)
{
	result.Init(numbers.size(), ssCnt);
	if(ssCnt == 0) return;
	numbIDWidth = DigitsCnt(numbers.MaxID());
	if(ssCnt > 1)	numbers.Sort();		// sorted by chrom length, but disrupts chromosome order
	Partition p(numbers, result, avr = numbers.AvrSum(ssCnt), limMult);
}

// Creates numbers partition by identified values, with sums sorted in descending order
//	@chrFiles:chrom sizes container to be distributed
//	@ssCnt: count of subsets;
//	if 0 then creates an empty partition with undefined (maximum type's value) inaccuracy
//	@limMult: STree() call's limit multiplier -
//	it increases the limit of 1 million recursive invokes by limMult times;
//	if 0 then omit DSTree method invoking (fast, but not 'perfect')
effPartition::effPartition(const ChromSizesExt& cSizes, ss_id ssCnt, UINT limMult)
{
	IdNumbers numbs(cSizes);
	Init(numbs, ssCnt, limMult);
}

// Prints subsets
//	@prSumDiff: if true then prints sum diference (in absolute and relative)
//	@prNumbsCnt: maximum count of printed number IDs or 0 if print all (by default)
void effPartition::Print(bool prSumDiff, size_t prNumbsCnt) const
{
	result.Print(numbIDWidth, prNumbsCnt);
	if(prSumDiff)
		cout << SignDbg << "partition inaccuracy: " 
			 //<< Inacc() << " (" << setprecision(3) << RelInacc() << "%)\n";
			 << setprecision(3) << RelInacc() << "%\n";
}
