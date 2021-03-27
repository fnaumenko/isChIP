/**********************************************************
Data.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 27.03.2021
-------------------------
Provides common data functionality
***********************************************************/
#pragma once

#include "TxtFile.h"
#include <algorithm>    // std::sort
#include <map>

#define	CID(it)	(it)->first

typedef pair<ULONG, ULONG> p_ulong;	// double chrlen

class ChromSizes;

// Basic class for objects keeping in Tab File
class Obj
{
public:
	enum class eInfo {	// defines types of outputted info
		NONE,	// nothing printed: it is never pointed in command line
		LAC,	// laconic:		print file name if needs
		NM,		// name:		print file name
		STD,	// standard:	print file name and items number
		STAT,	// statistics:	print file name, items number and statistics
	};

//protected:
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
		enum class eAction /*:BYTE*/ {	// commented type doesn't compiled under Linux GNU
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
		DataInFile*	 _file;				// current reading file
		const	FT::eType _fType;		// type of data; can't use _file!
		const	eInfo	_info;			// input type of info
		mutable chrlen	_count;			// count of discovered ambiguous
		const	bool	_alarm;			// true if message should be printed
		mutable bool	_alarmPrinted;	// true if warning was printed
#ifdef _BIOCC
		short	_treatcID;			// treated chrom ID: -1 initial, cID if only one chrom was treated.
#endif								// UnID if more then one chrom was treated

		// ***** actions
		inline int Accept(eCase ambg)	{ return 1; }
		inline int Handle(eCase ambg)	{ PrintLineAlarm(ambg); return 0; }
		inline int Omit  (eCase ambg)	{ PrintLineAlarm(ambg); return -1; }
		inline int OmitQuiet(eCase ambg){ return -1; }
		inline int Abort (eCase ambg)	{ ThrowExcept(_Msgs[ambg].LineAlarm); return -1; }

		// Get action message 
		const inline char* Message(eCase spotter) const { return _ActionMsgs[int(_cases[spotter].Action)]; }
		
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
		//	@ecase: spotter's case
		void PrintLineAlarm(eCase ecase) const;

		// Prints case statistics
		//	@a: spotter's case
		//	@allCnt: total count of ambiguities
		//	@total: if true then prints total warning case
		void PrintCaseStat(eCase ecase, chrlen allCnt, bool total=false) const;

		// Prints items with specifying chrom
		//	@cID: readed chromosome's ID or Chrom::UnID if all
		//	@prAcceptItems: if true then prints number of accepted items
		//	@itemCnt: count of accepted items after treatment
		void PrintItems(chrid cID, bool prAcceptItems, long itemCnt) const;

	public:
		bool unsortedItems;		// true if items are unsorted
		bool noCheck;			// true if neither collect statistics nor apply any ambig filter 
		bool hasPrinted;		// true if something (except warnings) has being printed
		chrlen chrLen;			// length of readed chromosome
#ifdef _WIG
		chrlen lastEnd;			// last end position: used for BEDGRAPH initialization only
#endif

		// Creates an instance with omitted COVER, SHORT, SCORE and NEGL cases;
		// Features cases by default:
		// omitted DUPL cases and handled CROSS and ADJAC 
		Spotter (FT::eType format, eInfo info, bool alarm = false, 
			eAction dupl = eAction::OMIT,
			eAction crossANDadjac = eAction::HANDLE,
			eAction diffsz = eAction::ACCEPT	// in fact for Features it even doesn't check
		);

		// Sets current Tab File
		inline void SetFile (DataInFile& file) { _file = &file; }
		
		// Gets current Tab File
		inline DataInFile& File () const	{ return *_file; }

		inline Obj::eInfo Info() const	{ return _info; }

		inline FT::eType FileType() const	{ return _fType; }

		//inline bool IsAlarmPrinted() const	{ return _alarmPrinted;	}
#ifdef _BIOCC
		// Remember treated chrom.
		void	SetTreatedChrom(chrid cid);

		// Gets single treated chrom if it's defined, otherwise UnID.
		void KeepTreatedChrom() const { if(_treatcID != vUNDEF)	Chrom::SetCustomID(chrid(_treatcID)); }
#endif
		// Initializes given Region by second and third current reading line positions, with validating
		//	@rgn: Region that should be initialized
		//	return: true if Region was initialized successfully
		bool InitRegn(Region& rgn);

		// Adds statistics and print given spotter as alarm (if permitted)
		//	@spotter: given spotter
		//	return: treatment code: 1 - accept, 0 - handle, -1 - omit
		inline int TreatCase(eCase spotter)  {
			return (this->*_Actions[int(_cases[spotter].TickAction())])(spotter);
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
protected:
	bool _isBad;		// sign of invalidity
	bool _EOLneeded;	// true if LF 

	inline Obj() : _isBad(false), _EOLneeded(false) {}

	// Initializes new instance by tab file name.
	//	@title: title printed before file name
	//	@fName: name of file
	//	@spotter: temporary values & ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	@isInfo: true if file info should be printed
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@scoreNumb: number of 'score' filed (for FBED and BedGraph)
	void Init	(const char* title, const string& fName, Spotter& spotter, 
		ChromSizes& cSizes, bool isInfo, bool abortInvalid, BYTE scoreNumb = 5);

	// Initializes child instance from tab file
	//	return: numbers of all and initialied items for given chrom
	virtual p_ulong InitDerived(Spotter&, const ChromSizes&) = 0;

	// Gets item's title.
	//	@pl: true if plural form
	virtual const string & ItemTitle(bool pl = false) const = 0;

	// Prints LF if needs and flash stdout
	//	@printEOL: true if LF should be printed explicitly
	void PrintEOL(bool printEOL);

public:
	// Returns true if instance is invalid.
	inline bool IsBad()	const { return _isBad; }

	// Returns true if something was printed during initialization without LF
	inline bool EOLNeeded() const { return _EOLneeded; }
};

static const string range_out_msg = "myMap[]: invalid key ";

// 'myMap' is implemented in 2 different ways: based on map or on vector
template <typename K, typename T> class myMap
{
protected:
	typedef map<K,T> Items;

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
	inline T& AddElem(K key, const T& val) { return _items[key] = val;	}

	// Adds empty class type to the collection without checking key
	//	return: class type collection reference
	inline T& AddEmptyElem(K key) { return AddElem(key, T()); }

public:
	// Copies entry. There is no option for _NO_MAP !!!
	inline void Assign(const myMap& map) {	_items = map._items; }
	
	// Returns constant reference to the item at its key
	inline const T& At(const K& key) const { return _items.at(key); }

	// Returns reference to the item at its key
	inline T& At(const K key) { return _items.at(key); }

	inline T& operator[] (K key) { return _items[key]; }

	inline const T& operator[] (K key) const { return At(key); }

	// Searches the container for a key and returns an iterator to the element if found,
	// otherwise it returns an iterator to end (the element past the end of the container)
	inline Iter GetIter(K key) { return _items.find(key); }

	// Searches the container for a key and returns a constant iterator to the element if found,
	// otherwise it returns an iterator to cEnd (the element past the end of the container)
	inline const cIter GetIter(K key) const { return _items.find(key);	}

	// Adds value type to the collection without checking cID
	inline void AddVal(K key, const T & val) { _items[key] = val; }

	// Removes from the container an element key
	inline void Erase(K key) { _items.erase(key); }

	// Clear content
	inline void Clear() { _items.clear(); }

	// Insert value type to the collection: adds new value or replaces existed
	//void InsertVal(K key, const T & val) { _items[key] = val; }

	// Returns true if element with key exists in the container, and false otherwise.
	inline bool FindItem (K key) const { return _items.count(key) > 0;	}
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
	// Returns true if chromosome by iterator should be treated
	inline bool IsTreated(typename myMap<chrid,ChromData<T> >::cIter it) const { 
		return it->second.Treated;
	}
	
	// Returns true if chromosome by ID should be treated
	inline bool IsTreated(chrid cID) const { return IsTreated(this->GetIter(cID));	}

	inline const T& Data(typename myMap<chrid,ChromData<T> >::cIter it) const { return it->second.Data; }
	
	inline T& Data(typename myMap<chrid,ChromData<T> >::Iter it) { return it->second.Data; }

	// Returns count of chromosomes.
	inline chrid ChromCount() const { return chrid(this->Count()); }	// 'this' required in gcc

	// Gets count of treated chroms
	//chrid TreatedCount() const 
	//{
	//	chrid res = 0;
	//
	//	for(auto it=this->cBegin(); it!=this->cEnd(); res += IsTreated(it++));
	//	return res;
	//}

	// Returns true if chromosome cID exists in the container, and false otherwise.
	inline bool FindChrom(chrid cID) const { return this->FindItem(cID); }

	// Adds value type to the collection without checking cID
	inline void AddValue(chrid cid, const T& val) { this->AddVal(cid, ChromData<T>(val)); }

#ifdef	_READDENS
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
		for (it = this->Begin(); it != this->End(); it++)	// 'this' required in gcc
			if (it->second.Treated = obj.FindChrom(CID(it)))
				commCnt++;
			else if (printWarn)
				Err(Chrom::Absent(CID(it), "second file")).Warning();
		// set false treated chroms in a compared object
		for (it = obj.Begin(); it != obj.End(); it++)
			if (!FindChrom(CID(it))) {
				it->second.Treated = false;
				if (printWarn)
					Err(Chrom::Absent(CID(it), "first file")).Warning();
			}
		if( !commCnt )	Err("no common " + Chrom::Title(true)).Throw(throwExcept);
		return commCnt;
	}
#endif	// _READDENS
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
//private:
protected:
	// Checks if items are initialized, chromosome is uniq, and adds it to the instance
	//	@cID: chroms id
	//	@firstInd: first item index
	//	@spotter: spotter to close items
	void AddChrom(chrid cID, chrlen firstInd, const Spotter& spotter);

protected:
	static const BYTE _FieldsCnt = 6;	// count of fields readed from file

	// Initializes generalized BED instance from tab file.
	//	@spotter: spotter to control ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	return: numbers of all and initialied items for given chrom
	// It's separated because it's invoked not only in BaseItems::InitDerived, but also in Cover::InitDerived
	p_ulong InitBed(Spotter& spotter, const ChromSizes& cSizes);

	// Initializes instance from tab file: proxy
	//	@spotter: spotter to control ambiguities
	//	@cSizes: chrom sizes to control chrom length exceedeing
	//	return: numbers of all and initialied items for given chrom
	inline p_ulong InitDerived(Spotter& spotter, const ChromSizes& cSizes) { return InitBed(spotter, cSizes); }

	// Initializes size of positions container.
	virtual void ReserveItems(ULONG initSize) = 0;

	// Checks the last element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@spotter: possible ambiguities
	//  return: false if some ambiguous has found; true if alright
	virtual bool CheckLastPos(const Region& rgn, Spotter& spotter) = 0;

	// Adds item to the container.
	//	@rgn: Region with mandatory fields
	//	@spotter: temporary values & ambiguities
	//	return: true if item was added successfully
	virtual bool AddItem(const Region& rgn, Spotter& spotter) = 0;

	// Adds last items for current chrom to the container. For WIG only.
	//	@spotter: used do get last item end
	//	return: count of added items
	virtual UINT FinishItems(const Spotter& spotter) = 0;

	// Fills items from intermediate container
	virtual void FillChromItems(const Spotter& spotter) = 0;

	// Gets total count of items
	virtual size_t ItemsCount() const = 0;

	// Gets count of items for chrom
	virtual size_t ItemsCount(chrid cID) const = 0;

	// Shrinks items, sorts chroms and items if there are not sorted and prints message
	virtual void SortIfNecessary(Spotter& spotter, ULONG estItemCnt) = 0;

public:
	// Prints items name and count, adding chrom name if the instance holds only one chrom
	//	@prLF: if true then print line feed
	void PrintItemCount(bool prLF = true) const;

#ifdef _DEBUG
	virtual void PrintItem(chrlen itemInd) const = 0;

	void PrintChrom() const;

	// Prints BaseItems with limited or all items
	//	@itemCnt: first number of items for each chromosome; all by default
	void Print(const char* title, chrlen itemCnt=0) const;
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

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: temporary values & ambiguities
	//  return: true if item should be accepted; otherwise false
	virtual bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter) = 0;

	// Checks the last element for the new potential start/end positions for all possible ambiguous.
	//	@rgn: checked start/stop positions
	//	@spotter: temporary values & ambiguities
	//  return: true if item should be accepted; otherwise false
	inline bool CheckLastPos(const Region& rgn, Spotter& spotter) {
		return CheckPrevPos(rgn, _items.end() - 1, spotter);
	}

	// Gets a copy region by container iterator.
	virtual const Region Regn(cItemsIter it) const = 0;

	// Shrinks items, sorts chroms and items if there are not sorted and prints message
	void SortIfNecessary(Spotter& spotter, ULONG estItemCnt)
	{
		if (estItemCnt / ItemsCount() > 2)
			// shrink items.
			// actually shrink_to_fit() is defined in C++11, 
			// but for any case replace shrink_to_fit() by swap()
#ifdef OS_Windows	
			_items.shrink_to_fit();
#else
			vector<I>(_items).swap(_items);
#endif
		if (spotter.unsortedItems) {
			const bool prInfo = spotter.Info() > Obj::eInfo::NM;
			if (prInfo) {
				Err(ItemTitle(false) + " sorting...", spotter.File().CondFileName()).Throw(false, false);
				spotter.hasPrinted = true;
			}
			SortItems(spotter);
			if (prInfo)	dout << " done";
		}
	}

	// Sorts and rechecks items for the ambiguities.
	//	@spotter: temporary values & ambiguities
	void SortItems(Spotter& spotter)
	{
		ULONG	rmvCnt = 0;		// counter of removed items in current chrom
		ItemsIter it, itLast;

		for(Iter cit=Begin(); cit!=End(); cit++) {
			// reduce indexes after previous chrom recheck
			it		= _items.begin() + (Data(cit).FirstInd -= rmvCnt);
			itLast	= _items.begin() + (Data(cit).LastInd  -= rmvCnt);
			// sort items for current chrom
			sort(it++, itLast+1, I::CompareByStartPos);
			// recheck current chrom
			while(it <= itLast)
				if( CheckPrevPos(Regn(it), it-1, spotter) )	// is ambiguous?
					it++;
				else {		// ambiguous processed, then remove current item
					it = _items.erase(it);
					rmvCnt++;
					itLast = _items.begin() + --Data(cit).LastInd;
				}
		}
	}

	// Empty implementation of BaseItems<> method
	void FillChromItems(const Spotter& spotter) {}

protected:
	vector<I> _items;	// vector of bed-file items

	//inline vector<I>* ItemsPoint() const { return &_items; };
	//inline BYTE	SizeOfItem() const { return sizeof(I); }

	// Initializes size of positions container.
	inline void ReserveItems(ULONG size) { _items.reserve(size); }

	// Gets item.
	//	@it: chromosome's iterator
	//	@iInd: index of item
	inline const I& Item(cIter it, chrlen iInd) const {	return _items[Data(it).FirstInd + iInd]; }

	// Returns item.
	//	@cInd: index of chromosome
	//	@iInd: index of item
	inline const I& Item(chrid cID, chrlen iInd) const { return Item(GetIter(cID), iInd); }

#ifdef _DEBUG
	inline void PrintItem(chrlen itemInd) const { _items[itemInd].Print(); }
#endif

public:
	// Gets total count of items
	inline size_t ItemsCount() const { return _items.size(); }

	// Gets count of items for chrom
	//	@cit: chrom's iterator
	inline size_t ItemsCount(cIter cit) const { return cit->second.Data.ItemsCount(); }

	// Gets count of items for chrom
	//	@cit: chrom's ID
	inline size_t ItemsCount(chrid cID) const { return At(cID).Data.ItemsCount(); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const cItemsIter ItemsBegin(cIter cit) const { return _items.begin() + Data(cit).FirstInd; }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cID: chromosome's ID
	inline const cItemsIter ItemsBegin(chrid cID) const { return ItemsBegin(GetIter(cID)); }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome'sconstant  iterator
	inline const cItemsIter ItemsEnd(cIter cit) const { return _items.begin() + Data(cit).LastInd + 1; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cID: chromosome's ID
	inline const cItemsIter ItemsEnd(chrid cID) const { return ItemsEnd(GetIter(cID)); }

	// Returns an iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const ItemsIter ItemsBegin(cIter cit) { return _items.begin() + Data(cit).FirstInd; }

	// Returns an iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome's constant iterator
	inline const ItemsIter ItemsEnd(cIter cit) { return _items.begin() + Data(cit).LastInd + 1; }
};

#if !defined _ISCHIP && !defined _WIGREG

class BasicReads
{
private:
	//readlen	_readLen;	// constant length of Read + 1
	//bool	_strand;	// last readed Read strand
	char	_isStrandPres = -1;	// 1 if strand mark is present on the right position in the input data

protected:
	void CheckStrand(const Obj::Spotter& spotter)
	{
		if (_isStrandPres < 0) {
			_isStrandPres = spotter.File().IsItemHoldStrand();	// initialize once
			if (!_isStrandPres)
				Err("missing strand character", spotter.File().CondFileName()).Throw();
		}
	}
};

// 'Reads' represents a collection of reads.
class Reads : BasicReads, public Items<Read>
{
private:
	readlen	_readLen;	// constant length of Read + 1
	bool	_strand;	// last readed Read strand
	//char	_isStrandPres = -1;	// 1 if strand mark is present on the right position in the input data
#if defined _VALIGN || defined _CALLDIST
	bool	_isPEonly;	// true if only paired-end reads are acceptable
	bool	_paired;	// true if Reads are paired-end
#endif
#ifdef _VALIGN
	float	_minScore;	// score threshold: Reads with score <= _minScore are skipping
	float	_maxScore;	// maximum score along Reads
#endif

	// Gets a copy of region by container iterator.
	// Items<> abstract method implementation.
	inline Region const Regn(cItemsIter it) const { return Region(it->Pos, it->Pos + ReadLen()); }

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	// Items<> abstract method implementation.
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: temporary values & ambiguities
	//  return: true if item should be accepted; otherwise false
	bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter);

	// Adds Read to the container.
	// Abstract BaseItems<> method implementation.
	//	@rgn: Region with mandatory fields
	//	@spotter: temporary values & ambiguities
	//	return: true if Read was added successfully
	bool AddItem(const Region& rgn, Spotter& spotter);

	// Abstract BaseItems<> method empty implementation
	inline UINT FinishItems(const Spotter&) { return 0;  }

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
	//	@printfName: true if file name should be printed unconditionally, otherwise depends on info
	//	@abortInval: true if invalid instance should abort excecution
	//	@alarm: true if warning messages should be printed 
	//	@PEonly: true if only paired-end reads are acceptable
	//	@acceptDupl: true if duplicates are acceptable 
	//	@minScore: score threshold (Reads with score <= minScore are skipping)
	Reads(const char* title, const string& fName, ChromSizes& cSizes, eInfo info,
		bool printfName, bool abortInval, bool alarm, bool PEonly, bool acceptDupl = true, 
		int minScore = vUNDEF);

#if defined _VALIGN || defined _CALLDIST
	// Return true if Reads are paired-end
	inline bool	IsPE()	const { return _paired; }
#endif
#ifdef _VALIGN
	// Gets maximum score
	//inline Read::rNameType ReadNameType()	const { return _rNameType; }

	inline bool IsPosInName() const { return true; }

	// Gets maximum score
	inline float MaxScore()	const { return _maxScore; }
#endif
	// Gets an item's title
	// Obj abstract method implementation.
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const { return FT::ItemTitle(FT::eType::ABED, pl); };

	// Gets length of Read.
	inline readlen ReadLen() const { return _readLen; }

	// Returns Read's start position.
	//	@cID: chromosome's ID
	//	@rInd: index of read
	inline chrlen ReadPos(chrid cID, chrlen rInd=0) const { return Item(cID, rInd).Pos;	}
};

#endif	// !_ISCHIP && !_WIGREG

#ifdef _FEATURES
struct ValRegion : public Region
{
	float	Value;			// features's score

	inline ValRegion(const Region& rgn, float val = 1) : Value(val), Region(rgn) {}

#ifdef _DEBUG
	inline void Print() const { cout << Start << TAB << End << TAB << Value << LF;	}
#endif	// _DEBUG
};
#ifdef _ISCHIP
typedef ValRegion	Featr;
#else	// NO _ISCHIP
typedef Region	Featr;
#endif	// _ISCHIP

// 'Features' represents a collection of crhoms features
class Features : public Items<Featr>
{
#ifdef _ISCHIP
	readlen	_minFtrLen;		// minimal length of feature
	float	_maxScore;		// maximal feature score after reading
#elif defined _BIOCC
	// these vars needed to get warning if user call Reads instead of Features (without -a option)
	long	_fLen = 0;			// feature's length. 'long' since we compare the difference with Read length
	char	_isStrandPres = -1;	// 1 if strand mark is present on the right position in the input data
#endif	// _ISCHIP

	// Sets new end position on the feature if necessary.
	//	@rgn: current feature
	//	@end: potential new end position
	//	@treatCaseRes: result of treatment this spotter
	//	return: true if spotter is permitted (feature is valid)
	bool CorrectItemsEnd(Region& rgn, chrlen end, int treatCaseRes);

	// Gets item's title
	// Obj abstract method implementation.
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl=false) const	{ return FT::ItemTitle(FT::eType::BED, pl); }
	
	// Gets a copy of Region by item's iterator.
	// Abstract BaseItems<> method implementation.
	inline Region const Regn(cItemsIter it) const { return *it; }

	// Checks the element for the new potential start/end positions for all possible ambiguous.
	// Items<> abstract method implementation
	//	@rgn: checked start/stop positions
	//	@it: iterator reffering to the compared element
	//	@spotter: possible ambiguities
	//  return: true if item should be accepted; otherwise false
	bool CheckPrevPos(const Region& rgn, ItemsIter it, Spotter& spotter);

	// Adds feature to the container
	// Abstract BaseItems<> method implementation.
	//	@rgn: Region with mandatory fields
	//	@spotter: temporary values & ambiguities
	//	return: true if Read was added successfully
	bool AddItem(const Region& rgn, Spotter& spotter);

	// Abstract BaseItems<> method empty implementation
	inline UINT FinishItems(const Spotter&) { return 0; }

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
		Spotter spotter(FT::eType::BED, info, alarm);
		Init(title, fName, spotter, cSizes, info > eInfo::LAC || printfName, true, scoreInd);
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
	{
		Spotter spotter(FT::eType::BED, info, alarm);
		Init(title, fName, spotter, cSizes, info > eInfo::LAC || printfName, abortInvalid);
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

	// Gets the ordinary total length of all chromosome's features
	//	@it: chromosome's iterator
	inline chrlen FeaturesLength(cIter it) const { return chrlen(EnrRegLength(it, 0, 0)); }

	// Gets the ordinary total length of all chromosome's features
	//	@cID: chromosome's ID
	//inline chrlen FeaturesLength(chrid cID) const { return chrlen(EnrRegLength(cID, 0, 0)); }

	// Copies features coordinates to external DefRegions.
	void FillRegions(chrid cID, Regions& regn) const;

#endif	// _ISCHIP

#ifdef _BIOCC
	// Returns true if strand mark is present in the file data
	inline bool	IsStrandPres() const { return _isStrandPres; }

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
	chrid GetChromIDs(vector<chrid>& cIDs, const string& gName);

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
	//	checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
	ChromSizes(const char* gName, const char* sPath, bool prMsg, bool checkGRef = false);

	inline ChromSizes() { _ext = _gPath = _sPath = strEmpty; }

#if defined _READDENS || defined _BIOCC || defined _VALIGN
	// Initializes empty instance by SAM header data
	void Init(const string& samHeade);
#endif
	inline bool IsFilled() const { return Count(); }

	// Return true if chrom.sizes are defined explicitly, by user
	//inline bool IsExplicit() const { return !_gPath.length(); }		// false if path == strEmpty
	
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
		: _treatedCnt(0), ChromSizes(gName, sPath, printMsg, true) {}
	
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
	chrid	_ID;		// chrom ID
	char*	_seq;		// the nucleotides buffer
	chrlen	_len,		// length of chromosome shortened by Read::FixedLen
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
	inline chrlen Length()	const { return _len + Read::FixedLen; }

	// Gets Read on position or NULL if exceeding the chrom length
	const char* Read(chrlen pos) const { return pos > _len ? NULL : _seq + pos;	}

	const char* Read(chrlen pos, readlen len) const { return pos + len > Length() ? NULL : _seq + pos; }

#if defined _ISCHIP || defined _VALIGN

	// Creates a stub instance (for sampling cutting)
	//	@len: chrom length
	inline RefSeq(chrlen len) : _ID(Chrom::UnID), _seq(NULL), _gapLen(0)
	{ _effDefRgn.Set(0, _len = len - Read::FixedLen); }

	// Creates and fills new instance
	RefSeq(chrid cID, const ChromSizes& cSizes);

#endif
#ifdef _ISCHIP
	inline chrid ID() const { return _ID; }

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
	ChromSizes&		_cSizes;
	const chrlen	_minGapLen;	// minimal allowed length of gap
	//bool			_isEmpty = true;	// true if regions are not initialized by real values
#ifdef _BIOCC
	const bool		_singleRgn = true;	// true if this instance has single Region for each chromosome
#endif


public:
	// Creates an instance by genome name, from chrom sizes file or genome.
	//	@cSizes: chrom sizes
	//	@minGapLen: minimal length which defines gap as a real gap
	inline DefRegions(ChromSizes& cSizes, chrlen minGapLen)
		: _cSizes(cSizes), _minGapLen(minGapLen)
	{ Init(); }

	void Init();

	// Returns true if regions are not initialized
	inline bool IsEmpty() const {	return !Count(); }

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

	// Adds chromosomes and regions without check up
	inline void AddChrom (chrid cID, const Regions& rgns) {	AddVal(cID, rgns); }
#endif	// _BIOCC

#ifdef DEBUG
	void Print() const;
#endif
};
#endif	// _READDENS || _BIOCC

#if defined _ISCHIP || defined _CALLDIST
#include <array>

typedef pair<float, float> fpair;

// 'LenFreq' represents a fragment's/read's length frequency statistics
class LenFreq : map<fraglen,ULONG>
{
public:
	// combined type of distribution
	enum /*class*/ eCType {		// not class to have a cast to integer by default
		NORM = 1 << 0,
		LNORM = 1 << 1,
		GAMMA = 1 << 2,
		CNT = 3,
	};

private:
	using dtype = int;	// consecutive distribution type: just to designate dist type, used as an index
	using spoint = pair<unsigned, ULONG>;	// initial sequence point
	using point = pair<unsigned, float>;	// distribution point 

	// Returns combined distribution type by consecutive distribution type
	inline static eCType GetCType(dtype type) { return eCType(1 << type); }

	// Returns consecutive distribution type by combined distribution type
	const static dtype GetDType(eCType ctype) { return RightOnePos(int(ctype)); }

	enum class eSpec {	// distribution specification
		CLEAR,		// normal quality;	exclusive
		SMOOTH,		// complementary
		MODUL,		// modulated;	complementary
		EVEN,		// exclusive
		CROP,		// cropped;	exclusive
		HCROP,		// heavily cropped;		exclusive
		SDEFECT,	// slightly defective;	exclusive
		DEFECT		// defective; exclusive
	};

	static const char* sTitle[];
	static const string sSpec[];
	static const string sParams;
	static const string sInaccurate;
	const fraglen smoothBase = 1;	// splining base for the smooth distribution
	static const int hRatio = 2;	// ratio of the summit height to height of the measuring point
	static const float lghRatio;	// log of ratio of the summit height to height of the measuring point

	// Moving window (Sliding subset)
	class MW : protected vector<ULONG>
	{
	protected:
		// Constructor by half of moving window length ("base")
		MW(unsigned base) { insert(begin(), Size(base), 0); }

		// Add last value and pop the first one (QUEUE functionality)
		void PushVal(ULONG x);

	public:
		// Returns length of moving window
		//	@base: moving window half-length
		static fraglen Size(fraglen base) { return 2 * base + 1; }
	};

	// Simple Moving Average splicer
	// https://en.wikipedia.org/wiki/Moving_average
	class MA : public MW
	{
		ULONG	_sum = 0;		// sum of adding values

	public:
		// Constructor by half of moving window length ("base")
		inline MA(fraglen base) : MW(base) {}

		// Add value and return average
		float Push(ULONG x);
	};

	// Simple Moving Median splicer
	class MM : public MW
	{
		vector<ULONG> _ss;		// sorted moving window (sliding subset)
		ULONG(MM::* _push)(ULONG) = &MM::GetMedian;	// pointer to GetMedian function: real or empty (stub)

		//	Return stub median
		inline ULONG GetMedianStub(ULONG x) { return x; }

		//	Add value and return median
		ULONG GetMedian(ULONG x);

	public:
		// Constructor
		//	@base: moving window half-length;
		//	if 0 then empty instance (initialized by empty Push() method)
		MM(fraglen base);

		// Add value and return median
		inline ULONG Push(ULONG x) { return (*this.*_push)(x); }
	};

	// Keeps distribution params: PCC, mean(alpha), sigma(beta)
	struct DParams
	{
	private:	
		static const float UndefPCC;
	public:
		float	PCC = 0;			// Pearson correlation coefficient
		fpair	Params{ 0,0 };		// mean(alpha), sigma(beta)

		bool operator >(const DParams& dp) const { return PCC > dp.PCC; }

		inline bool IsUndefPcc() const { return PCC == UndefPCC; };

		inline void SetUndefPcc() { PCC = UndefPCC; };
	};

	// 'AllDParams' represents a collection of restored distribution params for all type of distribution
	class AllDParams
	{
		// 'QualDParams' keeps restored distribution parameters
		struct QualDParams
		{
			eCType	Type;			//  combined distribution type
			const char* Title;		// distribution type title
			DParams dParams;		// PCC, mean(alpha), sigma(beta)

			// Sets combined type and title by consecutive type
			void SetTitle(dtype type) { Type = GetCType(type); Title = sTitle[type]; }

			// Returns true if distrib parameters set
			inline bool IsSet() const { return dParams.PCC != 0; }

			// Prints restored distr parameters
			//	@s: print stream
			//	@maxPCC: masimum PCC to print relative PCC percentage
			void Print(dostream& s, float maxPCC) const;
		};

		array<QualDParams, eCType::CNT>	_allParams;
		bool _sorted = false;

		// Returns true if distribution parameters set in sorted instance
		bool IsSetInSorted(eCType ctype) const;

		// Returns number of distribution parameters set in sorted instance
		int SetCntInSorted() const;

		// Returns DParams by combined distribution type
		inline DParams& Params(eCType ctype) { return _allParams[GetDType(ctype)].dParams; }

		// Sorts in PCC descending order
		void Sort();

	public:
		// Default constructor
		AllDParams();

		// Set distribution parameters by type
		//	@type: consecutive distribution type
		//	@dp: PCC, mean(alpha) & sigma(beta)
		inline void SetParams(dtype type, const DParams& dp) { _allParams[type].dParams = dp; }

		// Clear normal distribution if its PCC is less then lognorm PCC by the threshold
		void ClearNormDistBelowThreshold(float thresh) {
			if (Params(eCType::LNORM).PCC / Params(eCType::NORM).PCC > thresh)
				Params(eCType::NORM).PCC = 0;
		}

		// Returns parameters of distribution with the highest (best) PCC
		//	@dParams: returned PCC, mean(alpha) & sigma(beta)
		//	return: consecutive distribution type with the highest (best) PCC
		dtype GetBestParams(DParams& dParams);

		// Prints sorted distibutions params
		//	@s: print stream
		void Print(dostream& s);
	};
	
	// Returns specification string by specification type
	inline static const string Spec(eSpec s) { return "Distribution " + sSpec[int(s)]; }
	
	// Returns true if consecutive type is represented in combo cType
	inline static bool IsType(eCType cType, dtype type) { return cType & (1 << type); }

	// Returns true if combo type is represented in combo cType
	inline static bool IsType(eCType cType, eCType type) { return cType & type; }

	// Calls distribution parameters by consecutive distribution type
	//	@keypts: key pointers: X-coord of highest point, X-coord of right middle hight point
	//	@params: returned mean(alpha) & sigma(beta)
	//	Defined as a member of the class only to use the private short name lghRatio
	static void (*SetParams[])(const fpair& keypts, fpair& params);

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	AllDParams	_allParams;			// distributions parameters
#ifdef _DEBUG
	mutable vector<point> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Returns estimated moving window half-length ("base")
	//	return: estimated base, or 0 in case of degenerate distribution
	fraglen GetBase();

	// Defines key pointers
	//	@base: moving window half-length
	//	@summit: returned X,Y coordinates of spliced (smoothed) summit
	//	return: key pointers: X-coord of highest point, X-coord of right middle hight point
	fpair GetKeyPoints(fraglen base, point& summit) const;

	// Compares this sequence with calculated one by given mean & sigma, and returns PCC
	//	@type: consecutive distribution type
	//	@dParams: returned PCC, input mean(alpha) & sigma(beta)
	//	@Mode: X-coordinate of summit
	//	@full: if true then correlate from the beginning, otherwiase from summit
	//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
	void CalcPCC(dtype type, DParams& dParams, fraglen Mode, bool full = true) const;

	// Calculates distribution parameters
	//	@type: consecutive distribution type
	//	@keyPts: key pointers: X-coord of highest point, X-coord of right middle hight point
	//	@dParams: returned PCC, mean(alpha) & sigma(beta)
	//	@Mode: X-coordinate of summit
	void SetPCC(dtype type, const fpair& keypts, DParams& dParams, fraglen Mode) const;

	// Calculates and print called distribution parameters
	//	@type: consecutive distribution type
	//	@base: moving window half-length
	//	@summit: returned X,Y coordinates of best spliced (smoothed) summit
	void CallParams(dtype type, fraglen base, point& summit);

	// Prints original distribution features
	//	@s: print stream
	//	@base: moving window half-length
	//	@summit: X,Y coordinates of spliced (smoothed) summit
	void PrintTraits(dostream& s, fraglen base, const point& summit);

	// Prints original distribution as a set of <frequency>-<size> pairs
	//	@s: print stream
	void PrintSeq(dostream& s) const;

public:
	// Default constructor
	LenFreq() {}

	// Constructor by pre-prepared frequency distribution file
	//	@fname: name of pre-prepared frequency distribution file
	LenFreq(const char* fname);

	// Returns true if distribution has not enough size
	inline bool IsDegenerate() const { return size() < 5; }

	// Adds fragment/read to statistics
	//	@len: frag's length
	inline void AddLen(fraglen len) { (*this)[len]++; }

	// Calculate and print dist fpair
	//	@s: print stream
	//	@type: combined type of distribution
	//	@prDistr: if true then print distribution additionally
	void Print(dostream& s, eCType type, bool prDistr = true);
};

#endif	// _ISCHIP