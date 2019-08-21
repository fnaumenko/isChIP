/**********************************************************
Data.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 16.06.2019
-------------------------
Provides common data functionality
***********************************************************/
#pragma once

#include "TxtFile.h"
#include <algorithm>    // std::sort

#ifndef _NO_MAP
	#include <map>
#endif	// _NO_MAP

#define	CID(it)	(it)->first

typedef pair<ULONG, ULONG> p_ulong;	// double chrlen

class ChromSizes;

// Basic class for objects keeping in Tab File
class Obj
{
public:
	enum eInfo {	// defines types of outputted info
		iNONE,	// nothing printed: it is never pointed in command line
		iLAC,	// laconic:		print file name if needs
		iNM,	// name:		print file name
		iEXT,	// standard:	print file name and items number
		iSTAT,	// statistics:	print file name, items number and statistics
	};

protected:
		// 'Spotter' handles items ambiguities and represents ambiguities statistics
	class Spotter
	{
	public:
		// Enum 'eCase' defines all possible ambiguous cases.
		// Check _CasesCnt value
		enum eCase	{
			DUPL,		// duplicated features
			CROSS,		// crossed features
			ADJAC,		// adjacent features
			COVER,		// coverage feature by another
			SHORT,		// too short features
			DIFFSZ,		// different size of reads
			//PE_ONLY,	// except PE reads only
			SCORE,		// filtered by score features
			EXCEED,		// start or stop position exceeded chrom length
			NEGL_CHR	// belong to negligible chromosome
		};

		// Enum 'eAction' defines all possible reactions for ambiguous.
		enum eAction /*:BYTE*/ {	// commented type doesn't compiled under Linux GNU
			ACCEPT,			// accept ambiguous feature/read "as is", without treatment
			HANDLE,			// handle ambiguous feature/read
			OMIT,			// omit ambiguous feature/read with alarm warning
			OMIT_SILENT,	// omit ambiguous feature/read without alarm warning
			ABORTING		// abort execution via exception
		};

	private:
		// pointer to reaction
		typedef int	(Obj::Spotter::*ReportCase)(eCase spotter);

		struct Msg {
			const char*	TotalAlarm;	// supplementary message, added to case message in statistics
			const char*	StatInfo;	// short text of ambiguous used in statistics
			const string LineAlarm;	// text of ambiguous used in line alarm
		};
		static Msg	_Msgs[];				// texts of all ambiguous
		static const char*	_ActionMsgs[];	// messages followed by reactions
		static const ReportCase	_Actions[];	// reactions
		static const BYTE	_CasesCnt = 9;	// count of cases of feature/read ambiguities
	
		struct Case {
			eAction Action;
			chrlen	Count;
			
			inline eAction TickAction() { Count++; return Action; }
		};
		Case	_cases[_CasesCnt];
		DataFile*	 _file;			// current reading file
		const	FT::fType _fType;	// type of data
		const	eInfo _info;		// input type of info
		mutable chrlen _count;		// count of discovered ambiguous
		const	bool _alarm;		// true if message should be printed
		mutable bool _alarmPrinted;	// true if warning was printed
#ifdef _BIOCC
		short	_treatcID;			// treated chrom ID: -1 initial, cID if only one chrom was treated.
									// UnID if more then one chrom was treated
#endif
		// ***** actions
		inline int Accept(eCase ambg)	{ return 1; }
		inline int Handle(eCase ambg)	{ PrintLineAlarm(ambg); return 0; }
		inline int Omit  (eCase ambg)	{ PrintLineAlarm(ambg); return -1; }
		inline int OmitQuiet(eCase ambg){ return -1; }
		inline int Abort (eCase ambg)	{ ThrowExcept(_Msgs[ambg].LineAlarm); return -1; }

		// Get action message 
		const inline char* Message(eCase spotter) const { return _ActionMsgs[_cases[spotter].Action]; }
		
		inline const string& ItemTitle(chrlen cnt = 1) const { return FT::ItemTitle(_fType, cnt!=1); }

		// Throws exception with given code, contained file name (if needed)
		inline const void ThrowExcept (Err::eCode code)	const {	
			Err(code, _file->LineNumbToStr(code).c_str()).Throw();
		}

		// Throws exception with given message, contained file name (if needed)
		inline const void ThrowExcept (const string& msg) const { 
			Err(msg, _file->LineNumbToStr().c_str()).Throw();
		}

		// Gets count of ambiguities
		chrlen Count() const;

		// Print given spotter as alarm
		//	@spotter: given spotter
		void PrintLineAlarm(eCase spotter) const;

		// Prints case statistics
		//	@spotter: spotter's case
		//	@allCnt: total count of ambiguities
		//	@total: if true then prints total warning case
		void PrintCaseStat(eCase spotter, chrlen allCnt, bool total=false) const;

		// Prints items with specifying chrom
		//	@cID: readed chromosome's ID or Chrom::UnID if all
		//	@prAcceptItems: if true then prints number of accepted items
		//	@itemCnt: count of accepted items after treatment
		void PrintItems(chrid cID, bool prAcceptItems, long itemCnt) const;

	public:
		bool unsortedItems;		// true if items are unsorted
		bool noCheck;			// true if neither collect statistics nor apply any ambig filter 
		bool wasPrinted;		// true if something (except warnings) was printed
		chrlen chrLen;			// length of readed chromosome

		// Creates an instance with omitted COVER, SHORT, SCORE and NEGL cases;
		// Features cases by default:
		// omitted DUPL cases and handled CROSS and ADJAC 
		Spotter (eInfo info, bool alarm, FT::fType format,
			eAction dupl = OMIT,
			eAction crossANDadjac = HANDLE,
			eAction diffsz = ACCEPT	// in fact for Features it even doesn't check
		);

		// Sets current Tab File
		//inline void SetFile (TabFile& file) { _file = &file; }
		inline void SetFile (DataFile& file) { _file = &file; }
		
		// Gets current Tab File
		//inline TabFile& File () const	{ return *_file; }
		inline DataFile& File () const	{ return *_file; }

		inline Obj::eInfo Info() const	{ return _info; }

		inline FT::fType FileType() const	{ return _fType; }

		//inline bool IsAlarmPrinted() const	{ return _alarmPrinted;	}
#ifdef _BIOCC
		// Remember treated chrom.
		void	SetTreatedChrom(chrid cid);

		// Gets treated chrom: set only one treated chrom, otherwise UnID.
		void KeepTreatedChrom() const { if(_treatcID != vUNDEF)	Chrom::SetCustomID(chrid(_treatcID)); }
#endif
		// Initializes given Region by second and third current reading line positions, with validating
		//	@rgn: Region that should be initialized
		//	@prevStart: previous start position
		//	return: true if Region was initialized successfully
		bool InitRegn(Region& rgn, chrlen prevStart);

		// Adds statistics and print given spotter as alarm (if permitted)
		//	@spotter: given spotter
		//	return: treatment code: 1 - accept, 0 - handle, -1 - omit
		inline int TreatCase(eCase spotter) {
			return (this->*_Actions[_cases[spotter].TickAction()])(spotter);
		}

		// Prints statistics.
		//	@cID: readed chromosome's ID or Chrom::UnID if all
		//	@title: string at the beginning; if NULL then this instance is used while initialization
		//	and don't feeds line
		//	@itemCnts: pair of count of all items AND count of accepted items 
		void Print(chrid cID, const char* title, const p_ulong& itemCnts);

		// Sets supplementary message, added to case message in statistics
		//	@spotter: given case
		//	@msg: supplementary message
		static inline void SetSupplAlarm(eCase spotter, const char* msg) {
			_Msgs[spotter].TotalAlarm = msg;
		}
	};	//***** end of Spotter

	bool _isBad;		// sign of invalidity
	bool _EOLneeded;	// true if EOL 

	inline Obj() : _isBad(false), _EOLneeded(false) {}

	// Initializes new instance by tab file name.
	//	@title: title printed before file name
	//	@fName: name of file
	//	@spotter: ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	@isInfo: true if file info hpuld be printed
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@scoreInd: index of 'score' filed; is set for FBED only
	void Init	(const char* title, const string& fName, Spotter& spotter, 
		ChromSizes& cSizes, bool isInfo, bool abortInvalid, BYTE scoreInd = 5);

	// Initializes child instance from tab file
	//	return: numbers of all and initialied items for given chrom
	virtual p_ulong InitDerived(Spotter&, const ChromSizes&) = 0;

	// Gets item's title.
	//	@pl: true if plural form
	virtual const string & ItemTitle(bool pl=false) const = 0;

	// Prints EOL if needs.
	//	@printEOL: true if EOL should be printed
	void PrintEOL(bool printEOL);

public:
	// Returns true if instance is invalid.
	inline bool IsBad()	const { return _isBad; }

	// Returns true if something was printed during initialization without EOL
	inline bool EOLNeeded() const { return _EOLneeded; }
};

static const string range_out_msg = "myMap[]: invalid key ";

// 'myMap' is implemented in 2 different ways: based on map or on vector
template <typename K, typename T> class myMap
{
protected:
#ifdef _NO_MAP
	typedef pair<K,T> Item;
	typedef vector<Item> Items;

	static inline bool Compare (const Item& i1, const Item& i2) {
		return i1.first < i2.first;
	}
#else
	typedef map<K,T> Items;
#endif	// _NO_MAP

private:
	Items _items;	// storage of key-value pairs

public:
	typedef typename Items::iterator Iter;			// iterator
	typedef typename Items::const_iterator cIter;	// constant iterator

	// Returns a random-access constant iterator to the first element in the container
	inline cIter cBegin() const { return _items.begin(); }
	// Returns the past-the-end constant iterator.
	inline cIter cEnd()	  const { return _items.end(); }

	// Returns a random-access iterator to the first element in the container
	inline Iter Begin()	{ return _items.begin(); }
	// Returns the past-the-end iterator.
	inline Iter End()	{ return _items.end(); }

protected:
	// Returns count of elements.
	inline size_t Count() const { return _items.size(); }

	// Adds class type to the collection without checking key.
	// Avoids unnecessery copy constructor call
	//	return: class type collection reference
	T& AddElem(K key, const T& val) {
	#ifdef _NO_MAP
		_items.push_back(Item(key, val));
		return (_items.end()-1)->second;
	#else
		return _items[key] = val;
	#endif	// _NO_MAP
	}

	// Adds empty class type to the collection without checking key
	//	return: class type collection reference
	inline T& AddEmptyElem(K key) { return AddElem(key, T()); }

public:
	// Copies entry. There is no option for _NO_MAP !!!
	void Assign(const myMap& map) {	
		_items = map._items;
	}
	
	// Reserves container's capacity
	inline void Reserve (K cnt) {
	#ifdef _NO_MAP
		if(cnt > 1)		_items.reserve(cnt);
	#endif
	}

	// Sorts container if possible
	inline void Sort() {
	#ifdef _NO_MAP
		sort(Begin(), End(), Compare);
	#endif
	}

	// Returns constant reference to the item at its key
	const T& At(const K& key) const {
	#ifdef _NO_MAP
		cIter it = GetIter(key);
		if( it == cEnd() )
			throw std::out_of_range (range_out_msg + BSTR(key));
		return it->second;
	#else
		return _items.at(key);
		//return _items[key];			// inserts a new element if no any element matched key
	#endif	// _NO_MAP
	}

	// Returns reference to the item at its key
	T& At(const K key) {
	#ifdef _NO_MAP
		Iter it = GetIter(key);
		if( it == End() )
			throw std::out_of_range (range_out_msg + BSTR(key));
		return it->second;
	#else
		return _items.at(key);
		//return _items[key];			// inserts a new element if no any element matched key
	#endif	// _NO_MAP
	}

	T& operator[] (K key) {
	#ifdef _NO_MAP
		Iter it = GetIter(key);
		return it == cEnd() ? AddEmptyElem(key) : it->second;
	#else
		return _items[key];
	#endif	// _NO_MAP
	}

	inline const T& operator[] (K key) const { return At(key); }

	// Searches the container for a key and returns an iterator to the element if found,
	// otherwise it returns an iterator to end (the element past the end of the container)
	Iter GetIter(K key) {
	#ifdef _NO_MAP
		for(Iter it = Begin(); it != End(); it++)
			if( it->first == key )	return it;
		return End();
	#else
		return _items.find(key);
	#endif	// _NO_MAP
	}

	// Searches the container for a key and returns a constant iterator to the element if found,
	// otherwise it returns an iterator to cEnd (the element past the end of the container)
	const cIter GetIter(K key) const {
	#ifdef _NO_MAP
		for(cIter it = cBegin(); it != cEnd(); it++)
			if( it->first == key )	return it;
		return cEnd();
	#else
		return _items.find(key);
	#endif	// _NO_MAP
	}

	// Adds value type to the collection without checking cID
	inline void AddVal(K key, const T & val) {
	#ifdef _NO_MAP
		_items.push_back(Item(key, val));
	#else
		_items[key] = val;
	#endif	// _NO_MAP
	}

	// Removes from the container an element key
	inline void Erase(K key) { 
	#ifdef _NO_MAP
		_items.erase(GetIter(key));
	#else
		_items.erase(key);
	#endif	// _NO_MAP
	}

	// Clear content
	inline void Clear() { _items.clear(); }

	// Insert value type to the collection: adds new value or replaces existed
	//void InsertVal(K key, const T & val) {
	//#ifdef _NO_MAP
	//	for(Iter it = Begin(); it != End(); it++)
	//		if( it->first == key ) {
	//			it->second = val;
	//			return;
	//		}
	//	_items.push_back( Item(key, val) );
	//#else
	//	_items[key] = val;
	//#endif	// _NO_MAP
	//}

	// Returns true if element with key exists in the container, and false otherwise.
	bool FindItem (K key) const {
	#ifdef _NO_MAP
		for(cIter it = cBegin(); it != cEnd(); it++)
			if( it->first == key )	
				return true;
		return false;
	#else
		return _items.count(key) > 0;
	#endif	// _NO_MAP
	}
};

// 'ChromData' implements sign whether chrom is involved in processing, and chrom's data itself
template <typename T> struct ChromData
{

	bool Treated;	// true if chrom is involved in processing
	T	 Data;		// chrom's data

	inline ChromData() : Treated(true), Data(T()) {}

	inline ChromData(const T& data) : Treated(true), Data(data) {}
};

// Basic class for chromosomes collection; keyword 'abstract' doesn't compiled in gcc
template <typename T> class Chroms : public myMap<chrid, ChromData<T> >
{
public:
	// Returns true if chromosome should be treated
	inline bool IsTreated(typename myMap<chrid,ChromData<T> >::cIter it) const { 
		return it->second.Treated;
	}
	
	inline const T& Data(typename myMap<chrid,ChromData<T> >::cIter it) const { return it->second.Data; }
	
	inline T& Data(typename myMap<chrid,ChromData<T> >::Iter it) { return it->second.Data; }

	// Returns count of chromosomes.
	inline chrid ChromCount() const { return chrid(this->Count()); }	// 'this' required in gcc

	// Gets count of treated chroms
	//chrid TreatedCount() const 
	//{
	//	typename myMap<chrid,ChromData<T> >::cIter it;
	//	chrid res = 0;
	//
	//	for(cIter it=this->cBegin(); it!=this->cEnd(); res += IsTreated(it++));
	//	return res;
	//}

	// Returns true if chromosome cID exists in the container, and false otherwise.
	inline bool FindChrom(chrid cID) const { return this->FindItem(cID); }

	// Adds value type to the collection without checking cID
	inline void AddValue(chrid cid, const T& val) { this->AddVal(cid, ChromData<T>(val)); }

	// Sets common chromosomes as 'Treated' in both of this instance and given Chroms.
	// Objects in both collections have to have second field as Treated.
	//	@obj: compared Chroms object
	//	@printWarn: if true print warning - uncommon chromosomes
	//	@throwExcept: if true then throw exception if no common chroms finded, otherwise print warning
	//	return: count of common chromosomes
	chrid	SetCommonChroms(Chroms<T>& obj, bool printWarn, bool throwExcept)
	{
		typename myMap<chrid,ChromData<T> >::Iter it;
		chrid commCnt = 0;

		// set treated chroms in this instance
		for(it = this->Begin(); it != this->End(); it++)	// 'this' required in gcc
			if( it->second.Treated = obj.FindChrom(CID(it)) )
				commCnt++;
			else if( printWarn )
				Err(Chrom::Absent(CID(it), "second file")).Warning();
		// set false treated chroms in a compared object
		for(it = obj.Begin(); it != obj.End(); it++)
			if( !FindChrom(CID(it)) ) {
				it->second.Treated = false;
				if( printWarn )
					Err(Chrom::Absent(CID(it), "first file")).Warning();
			}
		if( !commCnt )	Err("no common chromosomes").Throw(throwExcept);
		return commCnt;
	}
};

// 'ItemIndexes' representes a range of chromosome's features/reads indexes,
struct ItemIndexes
{
	chrlen	FirstInd;	// first index in items container
	chrlen	LastInd;	// last index in items container

	inline ItemIndexes(chrlen firstInd=0, chrlen lastInd=1)
		: FirstInd(firstInd), LastInd(lastInd-1) {}
		
	// Returns count of items
	inline size_t ItemsCount() const { return LastInd - FirstInd + 1; }
};

// Abctract class 'BaseItems' implements methods for creating items collection from file
class BaseItems : public Obj, public Chroms<ItemIndexes>
/*
 * Container of features/reads is complemented in derived class.
 * Skips comments in bed-file while reading.
 * strongly needs keyword 'abstract' but it doesn't compiled by GNU g++
 */
{
private:
	// Checks if chromosome is uniq and adds it to the container
	//	@cID: chroms id
	//	@firstInd: first item index
	//	@lastInd: last item index
	//	@fname: file name for exception message
	void AddChrom(chrid cID, chrlen firstInd, chrlen lastInd, const string& fname);

protected:
	static const BYTE _FieldsCnt = 6;	// count of fields readed from file

	// Initializes instance from tab file
	//	@spotter: spotter to control ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	return: numbers of all and initialied items for given chrom
	p_ulong InitDerived	(Spotter& spotter, const ChromSizes& cSizes);
	 
	// Initializes size of positions container.
	virtual void ReserveItemContainer(ULONG initSize) = 0;

	 // Shrink size of positions container.
	virtual void ShrinkItemContainer() = 0;

	// Checks the last element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@spotter: possible ambiguities
	//  return: false if some ambiguous has found; true if alright
	virtual bool CheckLastPos(const Region& rgn, Spotter& spotter) = 0;

	// Adds item to the container.
	//	@rgn: Region with mandatory fields
	//	@file: file to access to additionally fields
	//	return: true if Read was added successfully
	//virtual bool AddItem(const Region& rgn, TabFile& file) = 0;
	virtual bool AddItem(const Region& rgn, DataFile& file) = 0;

	// Sorts and rechecks items for the ambiguities.
	//	@spotter:  ambiguous filters
	virtual void SortItems(Spotter& spotter) = 0;

	// Sets and returns count of all Features/Reads.
	virtual void SetAllItemsCount() = 0;

	// Gets count of Features/Reads for chromosome
	virtual size_t ItemsCount(chrid cID) const = 0;

	// Gets count of all Features/Reads.
	virtual size_t AllItemsCount() const = 0;

	// Decreases Feature/Read position without checkup indexes.
	// Used in Shrink() method only.
	//	@cID: chromosome's ID
	//	@eInd: index of read
	//	@shift: decrease entry's start position on this value
	//	@regnEnd: region's end position to control Read location
	//	@return: true if Feature/Read is insinde this region
	//virtual bool DecreasePos(chrid cID, chrlen eInd, chrlen shift, chrlen regnEnd) = 0;

public:
	// Prints items name and count, adding chrom name if the instance holds only one chrom
	void PrintItemCount() const;

#ifdef DEBUG
	virtual void PrintItem(chrlen itemInd) const = 0;
#endif

public:
	//inline const ItemIndexes & operator[] (chrid cID) const { return At(cID); }

//#ifdef _READDENS
//	// Shifts elements's positions to collaps the 'holes' between external regions.
//	//	@cID: chromosome's ID
//	//	@regns: external (defined) regions
//	void ShrinkByID(chrid cID, const DefRegions &regns);
//#endif	// _READDENS

#ifdef DEBUG
	void PrintChrom() const;

	// Prints BaseItems with limited or all items
	//	@itemCnt: first number of items for each chromosome; all by default
	void Print(chrlen itemCnt=0) const;
#endif
};

// Abstract class 'Items' implements BaseItems's virtual container of features/reads
template <typename I> class Items : public BaseItems
{
protected:
	typedef typename vector<I>::iterator ItemsIter;

public:
	typedef typename vector<I>::const_iterator cItemsIter;

private:
	// Initializes size of positions container.
	inline void ReserveItemContainer(ULONG size) { _items.reserve(size); }

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: possible ambiguities
	//  return: true if item should be accepted; otherwise false
	virtual bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter) = 0;

	// Checks the last element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@spotter: possible ambiguities
	//  return: true if item should be accepted; otherwise false
	inline bool CheckLastPos(const Region& rgn, Spotter& spotter) {
		return CheckPrevPos(rgn, _items.end() - 1, spotter);
	}

	// Shrink size of positions container.
	inline void ShrinkItemContainer() {
	// actually shrink_to_fit() defined in C++11, 
	// but for any case replace shrink_to_fit() by swap()
#ifdef OS_Windows	
		_items.shrink_to_fit();
#else
		vector<I>(_items).swap(_items);
#endif
	}

	// Gets region by container iterator.
	virtual const Region Regn(cItemsIter it) const = 0;

	// Sorts and rechecks items for the ambiguities.
	//	@spotter:  ambiguous filters
	void SortItems(Spotter& spotter)
	{
		ULONG	rmvCnt = 0;		// counter of removed items in current chrom
		ItemsIter it, itLast;

		for(Iter cit=Begin(); cit!=End(); cit++) {
			// reduce indexes after previous chrom recheck
			it		= _items.begin() + (Data(cit).FirstInd -= rmvCnt);
			itLast	= _items.begin() + (Data(cit).LastInd  -= rmvCnt);
			_itemsCnt -= rmvCnt;
			// sort items for current chrom
			sort(it, itLast+1, I::CompareByStartPos);
			// recheck current chrom
			for(rmvCnt=0, it++; it<=itLast; it++)
				if( !CheckPrevPos(Regn(it), it-1, spotter) ) {
					// ambiguous situation: remove current item
					_items.erase(it--);
					Data(cit).LastInd--;
					itLast--;
					rmvCnt++;
				}
		}
	}

protected:
	vector<I> _items;	// vector of features/reads included in bed-file
	ULONG	_itemsCnt;	// count of all items. Initially it is equal to featrs.size(),
						// but may be reduced during extention

	//inline vector<I>* ItemsPoint() const { return &_items; };
	//inline BYTE	SizeOfItem() const { return sizeof(I); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const cItemsIter cItemsBegin(cIter cit) const { return _items.begin() + Data(cit).FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome'sconstant  iterator
	inline const cItemsIter cItemsEnd(cIter cit) const { return _items.begin() + Data(cit).LastInd + 1;	}

	// Returns an iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const ItemsIter ItemsBegin(cIter cit) { return _items.begin() + Data(cit).FirstInd; }

	// Returns an iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const ItemsIter ItemsEnd(cIter cit) { return _items.begin() + Data(cit).LastInd + 1;	}

	// Sets count of all features/reads.
	virtual void SetAllItemsCount()	{ _itemsCnt = _items.size(); }

	// Gets count of all features/reads.
	inline size_t AllItemsCount() const { return _itemsCnt; }

	// Gets count of items for chromosome or all by default.
	inline size_t ItemsCount(chrid cID) const {
		return cID==Chrom::UnID ? AllItemsCount() : At(cID).Data.ItemsCount();
	}

	// Gets item.
	//	@it: chromosome's iterator
	//	@iInd: index of item
	inline const I& Item(cIter it, chrlen iInd) const {
		return _items[Data(it).FirstInd + iInd];
	}

	// Returns item.
	//	@cInd: index of chromosome
	//	@iInd: index of item
	inline const I& Item(chrid cID, chrlen iInd) const { return Item(GetIter(cID), iInd); }

#ifdef DEBUG
	inline void PrintItem(chrlen itemInd) const { _items[itemInd].Print(); }
#endif

};

#if !defined _ISCHIP && !defined _WIGREG

// 'Reads' represents a collection of reads.
class Reads : public Items<Read>
{
private:
	readlen	_readLen;	// length of Read + 1
	bool	_strand;	// last readed Read strand
#if defined _VALIGN || defined _FRAGDIST
	bool	_isPEonly;	// true if only paired-end reads are acceptable
	bool	_paired;	// true if Reads are paired-end
#endif
#ifdef _VALIGN
	float	_minScore;	// score threshold: Reads with score <= _minScore are skipping
	float	_maxScore;	// maximum score along Reads
#endif

	// Gets a copy of region by container iterator.
	inline Region const Regn(cItemsIter it) const { return Region(it->Pos, it->Pos + ReadLen()); }

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: possible ambiguities
	//  return: true if item should be accepted; otherwise false
	bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter);

	// Adds Read to the container.
	//	@rgn: Region with mandatory fields
	//	@file: file to access to additionally fields
	//	return: true if Read was added successfully
	bool AddItem(const Region& rgn, DataFile& file);

	// Decreases Read's start position without checkup indexes.
	//	@cID: chromosome's ID
	//	@rInd: index of read
	//	@shift: decrease read's start position on this value
	//	@rgEnd: region's end position to control Read location
	//	@return: true if Read is insinde this region
	//bool DecreasePos(chrid cID, chrlen rInd, chrlen shift, chrlen rgEnd);

public:
	// Creates new instance from bed-file with output info.
	// Invalid instance will be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: name of bed-file
	//	@cSizes: chrom sizes to control the chrom length exceedeng; if NULL, no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@abortInval: true if invalid instance should abort excecution
	//	@alarm: true if warning messages should be printed 
	//	@PEonly: true if only paired-end reads are acceptable
	//	@acceptDupl: true if duplicates are acceptable 
	//	@minScore: score threshold (Reads with score <= minScore are skipping)
	Reads(const char* title, const string& fName, ChromSizes& cSizes, eInfo info,
		bool printfName, bool abortInval, bool alarm, bool PEonly, bool acceptDupl=true, int minScore=vUNDEF )
		: _readLen(0)
#if defined _VALIGN || defined _FRAGDIST
		, _isPEonly(PEonly), _paired(false)
#endif
#ifdef _VALIGN
		, _minScore(minScore), _maxScore(0)
#endif
	{ 
		FT::fType type = FT::GetType(fName.c_str());
		if(type == FT::BED)	type = FT::ABED;
		else if(type != FT::BAM)
			Err("wrong extension", printfName ? fName : NULL).Throw(abortInval);
		Spotter spotter(info, alarm, type,
			acceptDupl ? Spotter::ACCEPT : Spotter::OMIT_SILENT,	// duplicated reads
			Spotter::ACCEPT,		// crossed & adjacent reads: typical
			//ignoreDiffSize ? Spotter::OMIT : Spotter::ABORTING	// different Read size
			//Spotter::OMIT			// omit different Read size
			Spotter::ACCEPT			// accept different Read size
		);
		Init(title, fName, spotter, cSizes, info > iLAC || printfName, abortInval);
	}

#if defined _VALIGN || defined _FRAGDIST
	// Return true if Reads are paired-end
	inline bool	IsPE()	const { return _paired; }
#endif
#ifdef _VALIGN
	// Gets maximum score
	//inline Read::rNameType ReadNameType()	const { return _rNameType; }

	inline bool IsPosInName() const { return true; }

	// Gets maximum score
	inline float MaxScore()		const { return _maxScore; }
#endif
	// Gets an item's title
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const { return FT::ItemTitle(FT::ABED, pl); };

	// Gets length of Read.
	inline readlen ReadLen()		const { return _readLen; }

	// Returns Read's start position.
	//	@cID: chromosome's ID
	//	@rInd: index of read
	inline chrlen ReadPos(chrid cID, chrlen rInd=0) const { return Item(cID, rInd).Pos;	}

	// Returns iterator to first Read.
	//	@cit: chrom's const iterator
	inline cItemsIter ReadsBegin(cIter cit) const { return cItemsBegin(cit); }

	// Returns iterator to first Read.
	//	@cID: chrom's ID
	inline cItemsIter ReadsBegin(chrid cID) const { return ReadsBegin(GetIter(cID)); }

	// Returns iterator to last Read.
	//	@cit: chrom's const iterator
	inline cItemsIter ReadsEnd(cIter cit)	const { return cItemsEnd(cit); }

	inline cItemsIter ReadsEnd(chrid cID)	const { return ReadsEnd(GetIter(cID)); }

	// Gets count of Reads for chromosome or all by default.
	inline size_t Count(chrid cID=Chrom::UnID) const	{ return ItemsCount(cID); }
};

#endif	// !_ISCHIP && !_WIGREG

#ifdef _FEATURES
#ifdef _ISCHIP
struct Featr : public Region
{
	//const string	Name;	// features's name
	float	Score;			// features's score

	inline Featr(const Region& rgn, float score=1) : Score(score), Region(rgn) {}

#ifdef DEBUG
	void Print() const {
		cout << Start << TAB << End << TAB;
		//if( Name )	cout << Name << TAB;
		cout << Score << EOL;
	}
#endif	// DEBUG
};
#else	// NO _ISCHIP
typedef Region	Featr;
#endif	// _ISCHIP

// 'Features' represents a collection of crhoms features
class Features : public Items<Featr>
{
private:
	static const BYTE _MinFieldsCnt = 3;	// minimum count of fields readed from file

#ifdef _ISCHIP
	readlen	_minFtrLen;		// minimal length of feature
	float	_maxScore;		// maximal feature score after reading
#elif defined _BIOCC
	// these vars needed to get warning if user call Reads instead of Features (without -a option)
	long	_fLen;		// feature's length. 'long' since we compare the difference with Read length
	bool	_unifLen;	// true if all features have the same length
#endif	// _ISCHIP

	// Sets new end position on the feature if necessary.
	//	@rgn: current feature
	//	@end: potential new end position
	//	@treatCaseRes: result of treatment this spotter
	//	return: true if spotter is permitted (feature is valid)
	bool CorrectItemsEnd(Region& rgn, chrlen end, int treatCaseRes);

	// Gets item's title
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl=false) const	{ return FT::ItemTitle(FT::BED, pl); }
	
	// Gets a copy of Region by container iterator.
	inline Region const Regn(cItemsIter it) const { return *it; }

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: possible ambiguities
	//  return: true if item should be accepted; otherwise false
	bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter);

	// Adds feature to the container
	//	@rgn: Region with mandatory fields
	//	@file: file file to access to additionally fields
	//	return: true if Read was added successfully
	bool AddItem(const Region& rgn, DataFile& file);

	// Decreases Feature's positions without checkup indexes.
	//	@cID: chromosome's ID
	//	@fInd: index of Feature
	//	@shift: decrease Feature's positions on this value
	//	@rgEnd: region's end position to control feature location
	//	@return: true if Feature is insinde this region
	//bool DecreasePos(chrid cID, chrlen fInd, chrlen shift, chrlen rgEnd);

#ifdef _ISCHIP
	// Scales defined score through all features to the part of 1.
	void ScaleScores ();
#endif

public:
#ifdef _ISCHIP
	// Creates new instance by bed-file name
	// Invalid instance wil be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@scoreInd: index of 'score' field
	//	@bsLen: length of binding site: shorter features would be omitted
	//	@alarm: true if info about ambiguous lines is printed during initialization
	Features(const char* title, const string& fName, ChromSizes& cSizes, eInfo info,
		bool printfName, BYTE scoreInd, readlen bsLen, bool alarm)
		: _minFtrLen(bsLen), _maxScore(vUNDEF)
	{
		Spotter spotter(info, alarm, FT::BED);
		Init(title, fName, spotter, cSizes, info > iLAC || printfName, true, scoreInd);
		ScaleScores();
	}
#else
	// Creates new instance by bed-file name
	// Invalid instance will be completed by throwing exception.
	//	@title: title printed before file name or NULL
	//	@fName: name of bed-file
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@info: type of feature ambiguties that should be printed
	//	@printfName: true if file name should be printed unconditionally, otherwise deneds on info
	//	@abortInvalid: true if invalid instance should abort excecution
	//	@alarm: true if warning messages should be printed 
	Features(const char* title, const string& fName, ChromSizes& cSizes, eInfo info,
		bool printfName, bool abortInvalid, bool alarm)
#ifdef _BIOCC
		: _fLen(0), _unifLen(true)
#endif
	{
		Spotter spotter(info, alarm, FT::BED);
		Init(title, fName, spotter, cSizes, info > iLAC || printfName, abortInvalid);
	}
#endif	//  _ISCHIP
	
	// Gets chromosome's feature by ID
	//	@cID: chromosome's ID
	//	@fInd: feature's index, or first feature by default
	//inline const Featr& Feature(chrid cID, chrlen fInd=0) const { return Item(cID, fInd); }

	// Gets chromosome's feature by iterator
	//	@it: chromosome's iterator
	//	@fInd: feature's index, or first feature by default
	inline const Featr& Feature(cIter it, chrlen fInd=0) const { return Item(it, fInd); }

	// Returns count of chrom's features by its iter
	//	@it: chromosome's iterator
	inline chrlen Count(cIter it) const { return Data(it).ItemsCount(); }

	// Gets chromosome's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@it: chromosome's iterator
	//	@multiplier: 1 for numerics, 0 for letters
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	chrlen EnrRegLength(cIter it, BYTE multiplier, float fLen) const;

	// Gets chrom's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@cID: chromosome's ID
	//	@multiplier: 1 for numerics, 0 for nameds
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	//	return: chrom's total enriched regions length, or 0 if chrom is absent
	chrlen EnrRegLength(chrid cID, BYTE multiplier, float fLen) const;

	// Return min feature length
	chrlen GetMinFeatureLength() const;

	// Return min distance between features boundaries
	chrlen GetMinDistance() const;

	// Expands all features positions on the fixed length in both directions.
	// If extended feature starts from negative, or ends after chrom length, it is fitted.
	//	@extLen: distance on which Start should be decreased, End should be increased,
	//	or inside out if it os negative
	//	@cSizes: chrom sizes
	//	@info: displayed info
	//	return: true if positions have been changed
	bool Extend(chrlen extLen, const ChromSizes& cSizes, eInfo info);

	// Checks whether all features length exceed given length, throws exception otherwise.
	//	@len: given control length
	//	@lenDefinition: control length definition to print in exception message
	//	@sender: exception sender to print in exception message
	void CheckFeaturesLength(chrlen len, const string& lenDefinition, const char* sender);

#ifndef _ISCHIP

	// Returns iterator to first feature.
	//	@cit: chrom's const iterator
	inline cItemsIter FeaturesBegin(cIter cit) const { return cItemsBegin(cit); }

	// Returns iterator to last feature.
	//	@cit: chrom's const iterator
	inline cItemsIter FeaturesEnd(cIter cit) const	{ return cItemsEnd(cit); }

	// Gets the ordinary total length of all chromosome's features
	//	@it: chromosome's iterator
	inline chrlen FeaturesLength(cIter it) const {
		return chrlen(EnrRegLength(it, 0, 0));
	}

	// Gets the ordinary total length of all chromosome's features
	//	@cID: chromosome's ID
	//inline chrlen FeaturesLength(chrid cID) const {
	//	return chrlen(EnrRegLength(cID, 0, 0));
	//}

	// Copies features coordinates to external DefRegions.
	void FillRegions(chrid cID, Regions& regn) const;

#endif	// _ISCHIP

#ifdef _BIOCC
	// Gets count of features for chromosome or all by default.
	//	@cID: chromosome's ID
	//	return: count of features for existed chrom or 0, otherwise count of all features
	chrlen Count(chrid cID = Chrom::UnID) const;
	
	// Returns true if all features have the same length
	inline bool	SameFeaturesLength() const { return _unifLen; }

	friend class JointedBeds;	// to access GetIter(chrid)
#endif
};
#endif	// _FEATURES

// 'ChromSize' represents real and defined effective chrom lengths
struct ChromSize
{
	chrlen Real;			// real (actual) chrom length
#ifdef _ISCHIP
	mutable chrlen Defined;	// defined effective chrom length;
							// 'effective' means double length for autosomes, single one for somatic
	
	// Sets chrom's effective (treated) real length as defined
	inline chrlen SetEffDefined(bool autosome) const { return bool(Defined = (Real << int(autosome))); }
#endif

	inline ChromSize(chrlen size = 0) : Real(size)
#ifdef _ISCHIP
		, Defined(0)
#endif
	{}
};

// 'ChromSizes' represented chrom sizes with additional file system binding attributes
// Holds path to reference genome and to service files
class ChromSizes : public Chroms<ChromSize>
{
	string	_ext;		// FA files real extention; if empty then instance is initialized by service dir
	string	_gPath;		// ref genome path
	string	_sPath;		// service path
#ifdef _BIOCC
	mutable genlen	_gsize;		// size of whole genome
#endif

	// Returns length of common prefix before abbr chrom name of all file names
	//	@fName: full file name
	//	@extLen: length of file name's extention
	//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
	static int	CommonPrefixLength(const string& fName, BYTE extLen);

	// Initializes chrom sizes from file
	void Read(const string& fName);

	// Saves chrom sizes to file
	//	@fName: full file name
	void Write(const string& fName) const;

	// Fills external vector by chrom IDs relevant to file's names found in given directory.
	//	@cIDs: filling vector of chrom's IDs
	//	@gName: path to reference genome
	//	return: count of filled chrom's IDs
	BYTE GetChromIDs(vector<chrid>& cIDs, const string& gName);

	// Initializes the paths
	//	@gPath: reference genome directory
	//	@sPath: service directory
	//	@prMsg: true if print message about service fodler and chrom.sizes generation
	void SetPath(const string& gPath, const char* sPath, bool prMsg);

	// returns true if service path is defined
	inline bool IsServAvail() const { return _sPath.size(); }

protected:
	inline chrlen Length(cIter it) const { return Data(it).Real; }

public:
	inline const string& RefExt() const { return _ext; }
	
	// Creates and initializes an instance.
	//	@gName: reference genome directory or chrom.sizes file
	//	@sPath: service directory
	//	@prMsg: true if print message about service fodler and chrom.sizes generation
	ChromSizes(const char* gName, const char* sPath, bool prMsg);

	inline ChromSizes() { _ext = _gPath = _sPath = strEmpty; }

#if defined _READDENS || defined _BIOCC || defined _VALIGN
	// Initializes empty instance by SAM header data
	//	@cCnt: chroms count
	void Init(const string& samHeader, chrid cCnt);
#endif
	inline bool IsFilled() const { return Count(); }

	// Return true if chrom.sizes are defined explicitly, by user
	inline bool IsExplicit() const { return !_gPath.length(); }		// false if path == strEmpty
	
	// Returns reference directory
	inline const string& RefPath() const { return _gPath; }

	// Returns service directory
	inline const string& ServPath() const { return _sPath; }

	// Returns service directory
	inline const bool IsServAsRef() const { return _gPath == _sPath; }

	// Returns full ref chrom name by chrom ID 
	inline const string RefName(chrid cid) const { return _gPath + Chrom::AbbrName(cid); }

	// Returns full service chrom name by chrom ID 
	inline const string ServName(chrid cid) const { return _sPath + Chrom::AbbrName(cid); }

	inline chrlen operator[] (chrid cID) const { return At(cID).Data.Real; 
	}

#ifdef _BIOCC
	// Gets total size of genome.
	genlen GenSize() const;
#endif
#ifdef DEBUG
	void Print() const;
#endif
};

#ifdef _ISCHIP

// 'ChromSizesExt' provides additional functionality for the ChromSizes
class ChromSizesExt : public ChromSizes
{
	mutable chrid	_treatedCnt;	// number of treated chromosomes

	// Gets chrom's effective (treated) real length: a double length for autosomes, a single somatic
	//	@it: ChromSizes iterator
	inline chrlen SetEffLength(cIter it) const { return Data(it).SetEffDefined(IsAutosome(CID(it))); }

public:
	// Creates and initializes an instance.
	//	@gName: reference genome directory
	//	@sPath: service directory
	//	@printMsg: true if print message about chrom.sizes generation (in case of reference genome)
	inline ChromSizesExt(const char* gName, const char* sPath, bool printMsg)
		: ChromSizes(gName, sPath, printMsg) {}
	
	// Returns true if chrom is autosome, false for heterosome
	inline static bool IsAutosome(chrid cID) { return Chrom::IsAutosome(cID); }

	// Gets chrom's defined effective (treated) length
	//	@it: ChromSizes iterator
	chrlen DefEffLength(cIter it) const;

	// Gets count of treated chromosomes.
	inline chrid TreatedCount() const	{ return _treatedCnt; }

	// Sets actually treated chromosomes according template and custom chrom
	//	@templ: template bed or NULL
	//	return: number of treated chromosomes
	chrid	SetTreated	(bool statedAll, const BaseItems* const templ);
	
	// Prints threated chroms short names
	void	PrintTreatedChroms() const;
};

#endif	// _ISCHIP

// 'RefSeq' represented reference chromosome as array of nucleotides
class RefSeq
{
private:
	char*	_seq;		// the nucleotides buffer
	chrlen	_len,		// length of chromosome shortened by Read::Len
			_gapLen;	// total length of gaps
	Region	_effDefRgn;	// effective defined region (except 'N' at the begining and at the end)

	// Initializes instance and/or chrom's defined regions
	//	@fName: file name
	//	@rgns: chrom's defined regions: ripe or new
	//	@fill: if true fill sequence and def regions, otherwise def regions only
	//	return: true if chrom def regions are stated
	bool Init(const string& fName, ChromDefRegions& rgns, bool fill);

public:
	static bool	LetGaps;	// if true then include gaps at the edges of the ref chrom while reading
	static bool	StatGaps;	// if true count sum gaps for statistic output

	inline ~RefSeq()	{ if(_seq) delete [] _seq; }

	// Gets chrom legth
	inline chrlen Length()	const { return _len + Read::Len; }

	// Gets Read on position or NULL if exceeding the chrom length
	const char* Read(chrlen pos) const { return pos > _len ? NULL : _seq + pos;	}

#if defined _ISCHIP || defined _VALIGN

	// Creates a stub instance (for sampling cutting)
	//	@len: chrom length
	inline RefSeq(chrlen len) : _seq(NULL) { _effDefRgn.Set(0, _len = len - Read::Len); }

	// Creates and fills new instance
	RefSeq(chrid cID, const ChromSizes& cSizes);

#endif
#ifdef _ISCHIP
	// Gets count of nucleotides outside of defined region
	inline chrlen UndefLength() const { return Length() - _effDefRgn.Length(); }

	//inline float UndefLengthInPerc() const { return 100.f * (Length() - _effDefRgn.Length()) / Length(); }

	// Gets feature with defined region (without N at the begining and at the end) and score of 1
	inline const Featr DefRegion() const { return Featr(_effDefRgn); }
	
	// Gets start position of first defined nucleotide
	inline chrlen Start()	const { return _effDefRgn.Start; }
	
	// Gets end position of last defined nucleotide
	inline chrlen End()		const { return _effDefRgn.End; }
	
	// Gets total length of gaps
	inline chrlen GapLen()	const { return _gapLen; }

#elif defined _READDENS || defined _BIOCC

	// Creates an empty instance and fills chrom's defined regions
	//	@fName: FA file name with extension
	//	@rgns: new chrom's defined regions
	//	@minGapLen: minimal length which defines gap as a real gap
	inline RefSeq(const string& fName, ChromDefRegions& rgns, short minGapLen);

#endif	// _ISCHIP

//#if defined _FILE_WRITE && defined DEBUG 
//	// Saves instance to file by fname
//	void Write(const string & fname, const char *chrName) const;
//#endif
};

#if defined _READDENS || defined _BIOCC

// 'DefRegions' represents chrom's defined regions
// initialized from ChromSizes (statically, at once)
// or from .fa files (dynamically, by request).
class DefRegions : public Chroms<Regions>
{
private:
	ChromSizes&	_cSizes;
	const chrlen		_minGapLen;	// minimal allowed length of gap
#ifdef _BIOCC
	const bool	_singleRgn;	// true if this instance has single Region for each chromosome
#endif

public:
	// Creates an instance by genome name, from chrom sizes file or genome.
	//	@cSizes: chrom sizes
	//	@minGapLen: minimal length which defines gap as a real gap
	DefRegions(ChromSizes& cSizes, chrlen minGapLen);

	inline ChromSizes& ChrSizes() { return _cSizes; }

	// Gets chrom's size by chrom's iterator
	inline chrlen Size(cIter it)	const { return Data(it).LastEnd(); }

	// Gets chrom's size by chrom's ID
	inline chrlen Size(chrid cID)	const { return At(cID).Data.LastEnd(); }

	// Gets chrom regions by chrom ID; lazy for real chrom regions
	const Regions& operator[] (chrid cID);

#ifdef _BIOCC
	// Copying constructor: creates empty copy!
	inline DefRegions(const DefRegions& gRgns) :
		_cSizes(gRgns._cSizes), _minGapLen(gRgns._minGapLen), _singleRgn(gRgns._singleRgn) {}

	// Gets miminal size of chromosome: for represented chromosomes only
	chrlen MinSize() const;

	// Gets total genome's size: for represented chromosomes only
	genlen GenSize() const;

	// Gets true if this instance has single Region for each chromosome
	// (i.e. initialized from chrom.sizes file)
	inline bool	SingleRegions() const { return _cSizes.IsExplicit(); }

	// Adds chromosomes and regions without check up
	inline void AddChrom (chrid cID, const Regions& rgns) {	AddVal(cID, rgns); }
#endif	// _BIOCC

#ifdef DEBUG
	void Print() const;
#endif
};
#endif	// _READDENS || _BIOCC

#if defined _ISCHIP || defined _FRAGDIST
#include <queue>


// 'Freq' represents a fragment's length frequency statistics
class FragFreq : map<chrlen,ULONG>
{
	typedef map<chrlen,ULONG>::const_iterator citer;

	// Simple Moving Average splicer
	class SMA
	{
		short	_size;		// sliding subset
		size_t	_count;		// total count of adding elements
		ULONG	_sum;		// sum of adding elements
		queue<ULONG> _q;

	public:
		// Set splicer base (subset length) and clear it
		void SetSize(BYTE halfBase);

		// Add element and return average
		float Push(ULONG x);
	};

	// Simple Moving Median splicer
	class SMM
	{
		short	_size, _middle;	// sliding subset, the middle of subset
		citer	_end;			// end of external collection
		vector<ULONG> _v;

	public:
		// Set splicer base (subset length)
		//	@end: 'end' iterator of external collection
		void SetSize(BYTE halfBase, citer end);

		// Add element and return median
		ULONG Push(citer it);
	};

	// Calculate and print called lognormal distribution parameters
	void CalcDistrParams(dostream& s) const;

public:
	inline FragFreq() {}

	FragFreq(const char* fname);

	// Adds fragment to statistics
	//	@len: frag's length
	inline void AddFrag(fraglen len) { (*this)[len]++; }

	//inline citer Begin() const	{ return begin(); }
	//inline citer End()	 const	{ return end(); }

	// Calculate and print dist params
	//	@s: print stream
	void Print(dostream& s, bool prDistr = true) const;
};

#endif	// _ISCHIP
