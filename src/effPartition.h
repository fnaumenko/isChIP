/**********************************************************
effPartition.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
Last modified: 21.03.2019

Effective Partition Problem (optimization version) solution
Effectively distributes a set of input positive integers into k subsets,
such that the difference between the subset sums is minimized.
***********************************************************/

#pragma once
#include <vector>
#include <algorithm>    // sort
#include "Data.h"		// ChromFiles

static const char* sThread = "thread ";
static const char* SignDbg = "## ";	// Marker of output debug info

typedef chrid	numb_id;	// type of number's ID
typedef chrlen	numb_val;	// type of number's value
typedef thrid ss_id;		// type of subset's ID
typedef unsigned long	ss_sum;		// type of subset's sum

class Partition;

// Creates numbers partition among subsets, possibly according equally number values.
class effPartition
{
	friend class Partition;
	struct	Result;
public:
	struct	Subset;
	class	IdNumbers;

	// Represent the identified number to be partitioned
	struct IdNumber
	{
		friend class	Partition;
		friend struct	Subset;
		friend class	IdNumbers;
		friend struct	Result;
	private:
		ss_id	BinIInd;	// bin's inverse index: count_of_bins-1 for the first and 0 for the last one
		numb_id	 Id;		// number's ID
		numb_val Val;		// number's value

		// Returns true if number is not allocated and its value is fitted to the upper limit
		//	@val: permitted sum value
		inline bool IsFitted(float val) const { return !BinIInd && Val <= val; }

	public:
		// Constructor by number's ID and value
		inline IdNumber(numb_id id, numb_val val) : BinIInd(0), Id(id), Val(val) {}

		// Returns number's value
		inline numb_val Value() const { return Val; }
	
		// Returns number's ID
		inline numb_val ID() const { return Id; }

		// for sorting by descent
		inline bool operator < (const IdNumber& numb) const { return numb.Val < Val; }
	};

	typedef std::vector<IdNumber>::iterator numb_it;
	typedef std::vector<IdNumber>::const_iterator numb_cit;
	typedef std::vector<numb_id> numb_ids;
	typedef std::vector<numb_id>::const_iterator numb_id_cit;

	// Represents specialized container of identified numbers
	class IdNumbers : public std::vector<IdNumber>
	{
		friend class Partition;
		friend class effPartition;
		
		// minimum tolerance: part of minimal number value defines the start accuracy
		// with which the free space of the average sum is measured.
		// It is a heuristic value provides satisfactory inaccuracy in a single pass in most of cases
		static const int minTol = 20;

		inline IdNumbers() {}

		// Constructor by chrom sizes
		IdNumbers(const ChromSizesExt& cSizes);

		// Initializes instance by numbers
		void Init(const IdNumbers& numbers);

		// Copies numbers
		void Copy(const IdNumbers& numbers)	{
			memcpy(data(), numbers.data(), numbers.size()*sizeof(IdNumber));
		}

		// Sorts numbers in value descending order
		inline void Sort() { sort(begin(), end()); }

		// Returns minimum accuracy with which the available space of the average sum should be measured
		float GetMinUp() const;

		// Returns average values sum
		//	@ssCnt: number of subsets
		float AvrSum(ss_id ssCnt) const;

		// returns maximum number's ID
		numb_id MaxID() const;

		// Marks all numbers as unallocated
		void Reset();
	};

	// Represents partition
	struct Subset
	{
		friend class Partition;
		friend class effPartition;
	private:
		ss_id	 id;		// subset's ID
		ss_sum	 sumVal;	// sum of number values
		numb_ids numbIDs;	// subset number IDs container

		inline static bool SortByDescend(Subset s1, Subset s2) { return s1.sumVal > s2.sumVal; }
		inline static bool SortByAscend (Subset s1, Subset s2) { return s1.sumVal < s2.sumVal; }

		// Creates an empty Subset with reserved capacity
		//	@iCnt: expected count of numbers 
		inline Subset(numb_id iCnt) : id(0), sumVal(0) { numbIDs.reserve(iCnt); }

		// Clears instance: removes all number IDs, sets the sum of number values to zero
		void Clear() { sumVal = 0; numbIDs.clear(); }

		// Adds number
		//	@it: number's iterator
		void AddNumb(const numb_it& it);

		// Takes into account (charges) number's value
		//	@it: number's iterator
		//	@binInd: bin's inverse index
		void TakeVal(const numb_it& it, ss_id binInd);

		// Eliminates (discharges) number's value
		//	@it: number's iterator
		void ElimVal(const numb_it& it);

		// Prints subset
		//	@sumW: maximum sum width (count of digits) to align sums left
		//	@valW: maximum number value width (count of digits) to align vals right
		//	@prNumbsCnt: maximum number of printed number IDs or 0 if all (by default)
		void Print(BYTE sumW, BYTE valW, size_t prNumbsCnt = 0) const;

	public:
		// Gets subset's ID
		inline ss_id ID() const { return id; }

		// Gets sum of subset numbers
        inline ss_sum Sum() const { return sumVal; }

		// Gets number IDs container
		inline const numb_ids& NumbIDs() const { return numbIDs; }
	};
	
	typedef std::vector<Subset>::const_iterator ss_cit;

private:
	// Represents obtained subsets and their difference between the subset sums
	struct Result
	{
		ss_sum				SumDiff;	// maximum difference between the subset sums
		std::vector<Subset>	Bins;		// subsets

		// Returns count of subsets.
		inline ss_id SubsetCount() const { return ss_id(Bins.size()); }

		// Initializes instance: creates empty bins with reserved capacity
		//	@nCnt: count of numbers
		//	@ssCnt: count of subsets
		void Init(size_t nCnt, ss_id ssCnt);

		// Sorts subsets in descending order and sets subsets ID starting with 1
		void SetIDs();

		// Clears instance: removes all subsets number IDs and resets the sum of the number values
		void Clear();

		// Outfits this result
		//	@numbers: numbers to be distributed
		//	@ind: index of the first subset: 1 for IGreedy() or 0 for STree()
		//	@diff: difference between the subset sums
		void Fill(IdNumbers& numbers, char ind, numb_val diff);

		// Sets the minimum/maximum subset sums and their difference for this instance
		//	@minSum: minimum subset sum to be retrieved
		//	@maxSum: maximum subset sum to be retrieved
		//	@return: true if sum difference is updated
		bool SetSumDiff(ss_sum& minSum, ss_sum& maxSum);

		// Sorts subsets by their sum
		//	@ascend: if true then in in ascending order, otherwise in descending order
		void Sort(bool ascend);

		// Prints sorted in descending order instance
		//	@valW: maximum number value width (count of digits) to align vals right
		//	@prNumbsCnt: maximum count of printed number IDs or 0 if all (by default)
		void Print(BYTE valW, size_t prNumbsCnt = 0) const;
	};

	Result	result;	// final partition
	float	avr;	// average sum value among subsets
	BYTE	numbIDWidth;	// maximum width (count of digits) of number ID

	// Initializes subsets by partitioned values
	//	@vals: numbers to be distributed
	//	@ssCnt: count of subsets
	//	@perfect: true if partition should be 'perfect'
	//	@limMult: DSTree call's limit multiplier; if 0 then omit DSTree method invoking
	void Init(IdNumbers& numbers, ss_id ssCnt, UINT limMult);

public:
	// Creates numbers partition by identified values, with sums sorted in descending order
	//	@cSizes: chrom sizes container to be distributed
	//	@ssCnt: count of subsets;
    //	if 0 then creates an empty partition with undefined (maximum type's value) inaccuracy
	//	@limMult: STree() call's limit multiplier -
	//	it increases the limit of 1 million recursive invokes by limMult times;
	//	if 0 then omit DSTree method invoking (fast, but not 'perfect')
	effPartition(const ChromSizesExt& cSizes, ss_id ssCnt, UINT limMult = 1);

	// Returns the count of subsets.
	inline ss_id SubsetCount() const { return result.SubsetCount(); }

	// Returns a reference to the subset at position n
	inline Subset& operator[](ss_id n) { return result.Bins[n]; }
	//inline const Subset& operator[](ss_id n) const { return result.Bins[n]; }

	// Returns a reference to the subsets container
	inline const std::vector<Subset>& Subsets() const { return result.Bins; }

	// Returns average summary value among subsets.
	inline float AvrSum() const { return avr; }

	// Returns inaccuracy: the difference between maximum and minimum summary value among subsets.
	inline ss_sum Inacc() const { return result.SumDiff; }

	// Returns relative inaccuracy: the difference between maximum and minimum subset summary value
	// in percentage to average summary value among subsets.
	inline float RelInacc() const { return 100.f * Inacc() / avr; }

	// Sorts subsets by their sum
	//	@ascend: if true then in in ascending order, otherwise in descending order
	//void Sort(bool ascend = true) { if(Inacc()) result.Sort(ascend); }

	// Prints subsets
	//	@prSumDiff: if true tnen prints sum diference (in absolute and relative)
	//	@prNumbsCnt: maximum count of printed number IDs or 0 if print all (by default)
	void Print(bool prSumDiff = true, size_t prNumbsCnt = 0) const;
};
