#ifndef HT_H
#define HT_H
#include <vector>
#include <iostream>
#include <cmath>

typedef size_t HASH_INDEX_T;


// Complete - Base Prober class
template <typename KeyType>
struct Prober {
    // Data members
    HASH_INDEX_T start_;    // initial hash location, h(k)
    HASH_INDEX_T m_;        // table size
    size_t numProbes_;      // probe attempts for statistic tracking
    static const HASH_INDEX_T npos = (HASH_INDEX_T)-1; // used to indicate probing failed
    void init(HASH_INDEX_T start, HASH_INDEX_T m, const KeyType& key) 
    {
        (void) key;  // avoid unused argument warnings since base class doesn't use key
        start_ = start;
        m_ = m;
        numProbes_ = 0;
    }
    HASH_INDEX_T next() {
        throw std::logic_error("Not implemented...should use derived class");
    }
};

// Almost Complete - Fill in the if statement below.
template <typename KeyType> // can be any type of key
struct LinearProber : public Prober<KeyType> { // inherits from prober .. linear type

    HASH_INDEX_T next() // int next
    {
        // Complete the condition below that indicates failure
        // to find the key or an empty slot
        if(this->numProbes_ == this->m_ ) { // if reaches table size amount of tries
            return this->npos; 
        }
        HASH_INDEX_T loc = (this->start_ + this->numProbes_) % this->m_; // h(k) + 1 ,  h(k) + 2 ... mod m
        this->numProbes_++; // update # of tries
        return loc;
    }
};

// To be completed
template <typename KeyType, typename Hash2>
struct DoubleHashProber : public Prober<KeyType> // inherits from prober .. doublehasher type
{
    Hash2 h2_;              /// h2(k)
    HASH_INDEX_T dhstep_;   /// Stepsize to use for double hash probing
    HASH_INDEX_T modulus;

    /// Moduli to use for double hashing as table increases (resizes)
    static const HASH_INDEX_T DOUBLE_HASH_MOD_VALUES[]; // mod values to use for double hashing given below
    /// The number of elements in the array above
    static const int DOUBLE_HASH_MOD_SIZE;

    //==================================
    // Add data members, as desired
    //==================================

private:
    // Complete
    HASH_INDEX_T findModulusToUseFromTableSize(HASH_INDEX_T currTableSize) // size of hashtable is passed
    { // pass function with current table size and will return next biggest mod to use that is smaller than table size 
        HASH_INDEX_T modulus = DOUBLE_HASH_MOD_VALUES[0]; 
        // find the modulus that is just smaller than the table size
        for(int i=0; i < DOUBLE_HASH_MOD_SIZE && DOUBLE_HASH_MOD_VALUES[i] < currTableSize; i++)
        { // mod value should be less than current table size
            modulus = DOUBLE_HASH_MOD_VALUES[i];
        }
        return modulus;
    }
public:
    /**
     * @brief Construct a new Double Hash Prober
     *     Accepts information that must be provided outside the hash table 
     *      (i.e. apart from internal hash table info/parameters)
     *        
     * 
     * @param h2  h2(k) - Object with an operator()(const KeyType&) defined for it
     */
    DoubleHashProber(const Hash2& h2 = Hash2()) : h2_(h2) {} // 
    
    /**
     * @brief Supplies info the hash table must provide
     * 
     * @param start Starting location for probing (i.e. h1(k))
     * @param m     Table size
     * @param key   Key (in case further hashing is necessary)
     */
    // Complete
    void init(HASH_INDEX_T start, HASH_INDEX_T m, const KeyType& key) 
    { // starting idx for probe sequence (where linear collided), hash table size, the key to be inserted/deleted
        Prober<KeyType>::init(start, m, key);
        HASH_INDEX_T modulus = findModulusToUseFromTableSize(m); // gives mod
        // Compute probe stepsize given modulus and h2(k) 
        dhstep_ = modulus - h2_(key) % modulus; // gives step size
    }

    // To be completed
    HASH_INDEX_T next() // should return either nullptr or 0 if reached upper limit
    { // if reached end size m, resize using new mod
     // else: use m2 - h2(k)%m2 to hash (the dhstep_)

      if (this->numProbes_ == this->m_) {
        return this->npos;
      }
      // h(k,i) = [ h1(k) + i*h2(k) ] mod m
      HASH_INDEX_T idx = (this->start_ + this->numProbes_* this->dhstep_) % this-> m_;
      this->numProbes_ +=1;
      return idx;
      // the next() of doublehasher will have a step size  oh dhstep which is calculated in init()

    }
};

// Initialization of static array (do not alter)
template <typename KeyType, typename Hash2>
const HASH_INDEX_T DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_VALUES[] =
{
    7, 19, 43, 89, 193, 389, 787, 1583, 3191, 6397, 12841, 25703, 51431, 102871,
    205721, 411503, 823051, 1646221, 3292463, 6584957, 13169963, 26339921, 52679927,
    105359939, 210719881, 421439749, 842879563, 1685759113
};

// Initialization of static array size (do not alter)
template <typename KeyType, typename Hash2>
const int DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_SIZE = 
    sizeof(DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_VALUES)/sizeof(HASH_INDEX_T);

// Hash Table Interface
template<
    typename K, 
    typename V, 
    typename Prober = LinearProber<K>,
    typename Hash = std::hash<K>, 
    typename KEqual = std::equal_to<K> >
class HashTable
{
public:
    typedef K KeyType;
    typedef V ValueType;
    typedef std::pair<KeyType, ValueType> ItemType;
    typedef Hash Hasher;
    struct HashItem {
        ItemType item;
        bool deleted;
        HashItem(const ItemType& newItem){
            item = newItem;
            deleted = false;
        }
    };

    /**
     * @brief Construct a new Hash Table object
     * 
     * @param resizeAlpha Loading factor threshold at which the table should resize
     * @param prober Probing object of type Prober
     * @param hash Hash functor that supports hash(key) and returns a HASH_INDEX_T
     * @param kequal Functor that checks equality of two KeyType objects
     */
    HashTable(
        double resizeAlpha = 0.4, 
        const Prober& prober = Prober(),
        const Hasher& hash = Hasher(), 
        const KEqual& kequal = KEqual());

    /**
     * @brief Destroy the Hash Table object and delete all remaining
     *         key,value pairs
     * 
     */
    ~HashTable();

    /**
     * @brief Returns true if the table has no non-deleted key,value pairs,
     *        and false otherwise
     * 
     */
    bool empty() const;

    /**
     * @brief Returns number of (non-deleted) key,value pairs in the table
     * 
     * @return size_t 
     */
    size_t size() const;

    /**
     * @brief Inserts a new item into the map, or, if an item with the
     *        given key already exists, it updates the Value of that item
     *        with the second value of the pair, p
     * 
     * @param p Pair to insert  
     * @throw std::logic_error If no free location can be found
     */
    void insert(const ItemType& p);

    /**
     * @brief Removes (marks as deleted) the item with the given key.  
     * Does nothing if an item with the given key does not exist.
     * 
     * @param key 
     */
    void remove(const KeyType& key);

    /**
     * @brief Finds an item with the given key and returns a pointer 
     * to the key,value pair
     * 
     * @param key 
     * @return ItemType const* nullptr is returned if the key does not exist
     */
    ItemType const * find(const KeyType& key) const;
    ItemType * find(const KeyType& key);

    /**
     * @brief Returns the value corresponding to the given key
     * 
     * @param key 
     * throw std::out_of_range if the key does not exist
     * @return ValueType Value associated with key
     */
    const ValueType& at(const KeyType& key) const;
    ValueType& at(const KeyType& key);
    const ValueType& operator[](const KeyType& key) const;
    ValueType& operator[](const KeyType& key);

    // Debug / Performance functions
    void reportAll(std::ostream& out) const;
    void clearTotalProbes() { totalProbes_ = 0; }
    size_t totalProbes() const { return totalProbes_; }
private:
    /**
     * @brief Helper routine to find a given key
     * 
     * @param key 
     * @return HashItem* returns nullptr if key does not exist
     */
    HashItem * internalFind(const KeyType& key) const;
    /**
     * @brief Performs the probing sequence and returns the index
     * of the table location with the given key or the location where
     * key can be inserted (i.e. the index now contains nullptr) but is
     * available.
     * 
     * @param key 
     * @return returns npos is the key does not exist and
     * no free location is available
     */
    HASH_INDEX_T probe(const KeyType& key) const;

    // Constant to signify an invalid hash location is being returned
    static const HASH_INDEX_T npos = Prober::npos;

    /**
     * @brief Resizes the hash table replacing the old with a new
     * table of the next prime size given in CAPACITIES.  Must rehash
     * all non-deleted items while freeing all deleted items.
     * 
     * Must run in O(m) where m is the new table size
     * 
     * @throws std::logic_error if no more CAPACITIES exist
     */
    void resize();

    // Data members
    std::vector<HashItem*> table_; // actual hash table
    Hasher hash_;   
    KEqual kequal_;
    mutable Prober prober_;  // mutable allows const member functions to modify this member
    // debug/performance counters
    mutable size_t totalProbes_; // mutable allows const member functions to modify this member
    // prime capacities to be used when resizing/rehashing is needed
    static const HASH_INDEX_T CAPACITIES[];
    HASH_INDEX_T mIndex_;  // index to CAPACITIES
    double resizeAlpha;

    // ADD MORE DATA MEMBERS HERE, AS NECESSARY

};

// ----------------------------------------------------------------------------
//                           Hash Table Implementation
// ----------------------------------------------------------------------------

// Static array of prime table sizes
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
const HASH_INDEX_T HashTable<K,V,Prober,Hash,KEqual>::CAPACITIES[] =
    {
        11, 23, 47, 97, 197, 397, 797, 1597, 3203, 6421, 12853, 25717, 51437, 102877,
        205759, 411527, 823117, 1646237, 3292489, 6584983, 13169977, 26339969, 52679969,
        105359969, 210719881, 421439783, 842879579, 1685759167
    };

// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HashTable<K,V,Prober,Hash,KEqual>::HashTable(
    double resizeAlpha, const Prober& prober, const Hasher& hash, const KEqual& kequal)
       :  resizeAlpha(resizeAlpha), hash_(hash), kequal_(kequal), prober_(prober)
{ // init hashtable
    // Initialize any other data members as necessary
    mIndex_ = 0; // idx of table size (CAPACITIES) vector
    totalProbes_ = 0; // init
  
    table_.resize(CAPACITIES[mIndex_]); // first table size in list of sizes

    

}

// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HashTable<K,V,Prober,Hash,KEqual>::~HashTable()
{ // delete/deallocate hashtable vector
  // delete table_;
  for (auto it = table_.begin(); it != table_.end(); it++) { // delete the pointers in table_
    delete *it;
    *it = nullptr;

  }
  table_.clear();
}

// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
bool HashTable<K,V,Prober,Hash,KEqual>::empty() const
{ // return true if table_ is empty
  // ask runtime requirement

  for (auto it = table_.begin(); it != table_.end(); it++) {
    if (*it != NULL && (*it)->deleted ==false) {
      return false;
    } // if there is a value, leave loop and return false right away
  }

  return true;
}

// To be completed 
// update/store size as a data member 
//update the size in insert/remove

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
size_t HashTable<K,V,Prober,Hash,KEqual>::size() const
{ // return # of items STORED
  size_t count = 0;
  for (auto it = table_.begin(); it != table_.end(); it++) {
    if (*it != NULL && !(*it)->deleted) {
      count +=1;
    }
  }
  return count;

}

// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::insert(const ItemType& p)
{ // allocate HashItem - DONE
  // find free loc ( != nullptr) ( numprobes != m) do NOT stop at a bool deleted item (no reusing) -> deleted = false only - DONE
  // if impossible, throw logic_error (usually not the case)
  // if loading factor >= alpha BEFORE the insert, resize capacity - DONE
  //  call prober::init WHICH calls next() in its own function ?

  //BEFORE INSERT**
  // check if loading factor is greater than or equal to alpha
  if (size()/ static_cast<double>(CAPACITIES[mIndex_]) >= resizeAlpha) {
    std::cout << "resize from " << CAPACITIES[mIndex_] << std:: endl;
    resize();
    std::cout <<"resized to " << CAPACITIES[mIndex_] << std::endl;
  }

  // HASH_INDEX_T idx = probe(p.first);
  


  //allocate prober and find empty spot
  //void init(HASH_INDEX_T start, HASH_INDEX_T m, const KeyType& key) 

    HASH_INDEX_T modulus = CAPACITIES[mIndex_]; // does this->m_ method work? or table_.size() m_ is Prober's member
    HASH_INDEX_T start = hash_(p.first) % modulus;
    prober_.init(start, modulus, p.first); //p.first is the key of item p

    HASH_INDEX_T idx = prober_.next(); // should return loc (might be empty or might be collision) (linear OR double hash)
    while (idx != this->npos) { //  while given a valid idx of hashtable
      // if (table_[idx] == NULL || table_[idx]->deleted == true) { // make sure there hasn't been a deleted item here before
      //   if (table_[idx]==NULL) {
      //     table_[idx] = new HashItem(p); // allocate new hashitem
      //   } else {
      //     table_[idx]->item = p;
      //     table_[idx]->deleted = false;
      //   }
      //   return;
      if (table_[idx] != NULL) {
        if (table_[idx]->deleted == false && kequal_(table_[idx]->item.first, p.first)) {
          table_[idx]->item.second = p.second; // update value
          return;
        }
      } else {
        table_[idx] = new HashItem(p); // if null, allocate new item
        return;
      }

      // if it isn't empty space, call next() to probe again
      idx = prober_.next();
    }


    // throw std::logic_error("failed insert");


}

// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::remove(const KeyType& key) 
{ // mark an item as "deleted"
  // cannot find() or operate[] these values
  // counts toward loading factor table_.size()
  // 


  // calculate mod
  HASH_INDEX_T modulus = CAPACITIES[mIndex_]; // does this->m_ work? or table_.size() m_ is Prober's member
  // start index is hashed key
  HASH_INDEX_T start = hash_(key) % modulus;
  //prober initialize
  prober_.init(start, modulus, key); //p.first is the key of item p

  HASH_INDEX_T idx = prober_.next();
  while (idx != Prober::npos) { // loop until invalid is returned
    // check if spot is full, never been deleted, and is the item we are searching
    if (table_[idx] != NULL && table_[idx]->deleted == false && kequal_(table_[idx]->item.first,key)) {
        table_[idx]->deleted = true;
        return;
    }
    idx = prober_.next(); // if doesn't remove in this try, probe again 
  }

  // throw std::logic_error("failed to remove");

}


// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ItemType const * HashTable<K,V,Prober,Hash,KEqual>::find(const KeyType& key) const
{
    HASH_INDEX_T h = this->probe(key);
    if((npos == h) || nullptr == table_[h] ){
        return nullptr;
    }
    return &table_[h]->item;
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ItemType * HashTable<K,V,Prober,Hash,KEqual>::find(const KeyType& key)
{
    HASH_INDEX_T h = this->probe(key);
    if((npos == h) || nullptr == table_[h] ){
        return nullptr;
    }
    return &table_[h]->item;
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
const typename HashTable<K,V,Prober,Hash,KEqual>::ValueType& HashTable<K,V,Prober,Hash,KEqual>::at(const KeyType& key) const
{
    HashItem const * item = this->internalFind(key);
    if(item == nullptr) { throw std::out_of_range("Bad key"); }
    return item->item.second;
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ValueType& HashTable<K,V,Prober,Hash,KEqual>::at(const KeyType& key)
{
    HashItem * item = this->internalFind(key);
    if(item == nullptr) { throw std::out_of_range("Bad key"); }
    return item->item.second;
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
const typename HashTable<K,V,Prober,Hash,KEqual>::ValueType& HashTable<K,V,Prober,Hash,KEqual>::operator[](const KeyType& key) const
{
    return this->at(key);
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ValueType& HashTable<K,V,Prober,Hash,KEqual>::operator[](const KeyType& key)
{
    return this->at(key);
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::HashItem* HashTable<K,V,Prober,Hash,KEqual>::internalFind(const KeyType& key) const
{
    HASH_INDEX_T h = this->probe(key);
    if((npos == h) || nullptr == table_[h] ){
        return nullptr;
    }
    return table_[h];
}

//************
// To be completed
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::resize() // after deleting item in null set it to nullptr
{ 
  // check if loading factor is greater than or equal to alpha
  // if (size()/CAPACITIES[mIndex_] < 0.4) {
  //   throw std::logic_error("no more space");
  // }    
  // resizing operation
  // to next prime in CAPACITIES

  if (mIndex_ + 1 < sizeof(CAPACITIES)/sizeof(CAPACITIES[0])) {
    mIndex_ +=1;
  }
  else {
    throw std::logic_error("no more capacities");
  }
  // size_t newcapacity = findModulusToUseFromTableSize(this->m_);

  // size_t newsize = CAPACITIES[newcapacity]; // given new table size from CAPACITIES[]

  // create new table with larger size, set all to NULL
  std::vector<HashItem*> resizedtable(CAPACITIES[mIndex_ ], NULL);
  // std::vector<HashItem*> oldtable = std::move(table_);
  // table_ = std::move(resizedtable);
  // // size_t prev = CAPACITIES[mIndex_];
  // mIndex_ +=1;

  // iterate thru each item in table_ 
  for (auto item : table_) {
    if (item && !item->deleted) { // check if it exists and that it hasn't been deleted before
      HASH_INDEX_T newidx = hash_(item->item.first) % CAPACITIES[mIndex_]; // find new idx
      while (resizedtable[newidx] != NULL) { // if idx within new table is already taken
        newidx = (newidx + 1) % CAPACITIES[mIndex_]; // find new spot
      }
      resizedtable[newidx] = item; // update with item
    }
    else {
      delete item; // if it has been deleted, don't move to new table
      item = nullptr;
    }
  }
  table_ = std::move(resizedtable);


}

//******************
// Almost complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HASH_INDEX_T HashTable<K,V,Prober,Hash,KEqual>::probe(const KeyType& key) const
{
    HASH_INDEX_T h = hash_(key) % CAPACITIES[mIndex_];
    prober_.init(h, CAPACITIES[mIndex_], key);

    HASH_INDEX_T loc = prober_.next(); 
    HASH_INDEX_T firstloc = loc;
    totalProbes_++;

    while(Prober::npos != loc)
    {
        if(table_[loc] == NULL ) {
          
            return loc;
        }
        // fill in the condition for this else if statement which should 
        // return 'loc' if the given key exists at this location
        else if(table_[loc] != NULL && table_[loc]->deleted == false && kequal_(table_[loc]->item.first, key)) {
            return loc;
        }
        loc = prober_.next();
        totalProbes_++;

      if (loc == firstloc) {
        break;
      }
    }

    return npos;
}

// Complete
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K, V, Prober, Hash, KEqual>::reportAll(std::ostream& out) const
{
	for(HASH_INDEX_T i = 0; i < CAPACITIES[mIndex_]; ++i)
	{
		if(table_[i] != nullptr)
		{
			out << "Bucket " << i << ": " << table_[i]->item.first << " " << table_[i]->item.second << std::endl;
		}
	}
}

#endif