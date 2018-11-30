//Start with a random point
//#define STARTPOINT1
//Start with the point that is furthest from all
#define STARTPOINT2

//Every CHECKRATE iterations, it can verify the streamlinging calculations
//It might be useful when rounding errors accumulate
#define CHECKRATE 100000
#define PRINTRATE 1000000

//This is useful if ones wants to keep all selected vertices in an ordered list
//This way, one can go through all selected vertices in O(k) instead of O(n)
//Adding a vertex is done in O(n/k) --- in our implementation (in average)
//Dropping a vertex is done in O(1).
//The wasted O(n/k) of the add operation is compensated by several gains (each iteration) of order (O(n)-O(k))
//We could have use std::set, but it looks quite inefficient to go through these point using (motr complex) iterators
//Let it uncommented if you want the ordered list to be used
//#define DISABLE_FAST_ORDERED_LIST

//Only to make the source code more easy to follow
//This is like #define NULL 0
#define NO_VTX_SELECTED 0

/* Useful Tweaks */
//You can remove all asserts (and other checks) with line bellow
#define NDEBUG

