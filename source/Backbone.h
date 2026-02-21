#ifndef BACKBONE_H
#define BACKBONE_H 1

#include <vector>

#include "BinaryVolume.h"
#include "Geometry.h"

using namespace std;

/* Backbone class for skeletonized segments
 *
 * Backbone stores the vertebra(e) associated with a segment of meat, as well as the length of the backbone to save a bit of computation.
 */
template <class T = voxelType>
class Backbone{
private:
	double len;
	vector<T> vertebrae;
	
public:
	
	/* Default constructor
	 */
	Backbone() : len(-1) {}
	
	/* Copy constructor
	 */
	template <class S>
	Backbone(const vector<S> &backbone) : len(-1), vertebrae(backbone) {}
	
	/* Single vertebra constructor.
	 * @i the <BinaryVolume> index of the single vertebra
	 *
	 * Constructs a <Backbone> from with the single vertebra `i`.
	 */
	template <class S>
	Backbone(S i) : len(-1), vertebrae(1, i) {}
	
	
	///////////////////////////////// Just the same as std::vector:
	
	/* Returns true if there are no vertebrae.
	 *
	 * Returns true if there are no vertebrae.
	 * @return true if there are no vertebrae
	 */
	bool empty() const { return vertebrae.empty(); }
	
	/* Clears all vertebrae.
	 */
	void clear(){ vertebrae.clear(); }
	
	/* Returns the number of vertebrae.
	 *
	 * Returns the number of vertebrae.
	 * @return the number of vertebrae
	 */
	unsigned long size() const { return vertebrae.size(); }
	
	/* Direct vertebrae access.
	 * @i the index of the vertebra in `vertebrae`
	 *
	 * Accesses the vertebra at `i` directly for reading and writing.
	 * @return reference to the vertebra at `i`
	 */
	template <class S>
	T& operator[](S i){ return vertebrae[i]; }
	
	/* Vertebrae access.
	 * @i the index of the vertebra in `vertebrae`
	 *
	 * Accesses the vertebra at `i` for reading.
	 * @return a copy of the vertebra at `i`
	 */
	template <class S>
	T operator[](S i) const { return vertebrae[i]; }
	
	/* Direct access to the back of vertebrae.
	 *
	 * Accesses the last vertebra directly for reading and writing.
	 * @return reference to the last vertebra
	 */
	T& back(){ return vertebrae.back(); }
	
	/* Access to the back of vertebrae.
	 *
	 * Accesses the last vertebra directly for reading.
	 * @return a copy of the last vertebra
	 */
	T back() const { return vertebrae.back(); }
	
	/* Direct access to the list of all vertebrae.
	 *
	 * Accesses the list of all vertebrae directly for reading and writing.
	 * @return reference to the list of all vertebrae
	 */
	vector<T>& getVertebrae() { return vertebrae; }
	
	/* Access to the list of all vertebrae.
	 *
	 * Accesses the list of all vertebrae for reading.
	 * @return a copy of the list of all vertebrae
	 */
	vector<T> getVertebrae() const { return vertebrae; }
	
	
	///////////////////////////////// Similar to std::vector:
	
	/* Adds `i` to the list of vertebrae.
	 * @i <BinaryVolume> index to add to vertebrae
	 *
	 * Adds `i` to the back of the list of vertebrae.
	 */
	template <class S>
	void add(S i){ vertebrae.push_back(i); }
	
	/* Adds `v` to the list of vertebrae.
	 * @v list of <BinaryVolume> indices to add to vertebrae
	 *
	 * Adds `v` to the back of the list of vertebrae.
	 */
	template <class S>
	void add(const vector<S> &v){ vertebrae.insert(vertebrae.end(), v.begin(), v.end()); }
	
	/* Adds <Backbone> `b` to the list of vertebrae.
	 * @b <Backbone> to add to vertebrae
	 *
	 * Adds <Backbone> `b` to the back of the list of vertebrae.
	 */
	template <class S>
	void add(const Backbone<S> &b){ vertebrae.insert(vertebrae.end(), b.vertebrae.begin(), b.vertebrae.end()); }
	
	/* Adds `i` to the list of vertebrae at the index `insert_index`.
	 * @i <BinaryVolume> index to add to vertebrae
	 * @insert_index index in `vertebrae` where `i` should be inserted
	 *
	 * Adds `i` to the list of vertebrae at index `insert_index`.
	 */
	template <class S, class R>
	void addAt(S i, R insert_index){
		if(insert_index < vertebrae.size())
			vertebrae.insert(vertebrae.begin() + insert_index, i);
		else
			vertebrae.push_back(i);
	}
	
	/* Removes the first instance of the <BinaryVolume> index `v` from the list of vertebrae.
	 * @v <BinaryVolume> index to remove from the list of vertebrae
	 *
	 * Removes the first instance of the <BinaryVolume> index `v` from the list of vertebrae (returns after the first instance is found).
	 */
	// assumes i appears only once in vertebrae
	template <class S>
	void remove(S v){
		for(typename vector<T>::iterator it(vertebrae.begin()); it != vertebrae.end(); it++){
			if(*it == v){
				vertebrae.erase(it);
				break;
			}
		}
	}
	
	/* Removes the vertebra at index `i`
	 * @i the index in `vertebrae` to remove
	 *
	 * Removes the vertebra at index `i` (in `vertebrae`) if the index exists.
	 */
	void removeIndex(unsigned long i){ if(i < vertebrae.size()) vertebrae.erase(vertebrae.begin() + i); }
	
	
	///////////////////////////////// Not the same as std::vector:
	
	/* Returns a <Backbone> with reversed vertebrae.
	 *
	 * Returns a <Backbone> that is identical in every way to this <Backbone> reversed.
	 * @return a <Backbone> with reversed vertebrae
	 */
	Backbone<> reversed() const {
		Backbone<> rev(*this); // in case there are other quantities or states that need to be copied (some day in the future...)
		
		unsigned long offset(vertebrae.size() - 1);
		for(unsigned long i(0); i < vertebrae.size(); i++)
			rev[offset - i] = vertebrae[i];
		
		return rev;
	}
	
	/* Returns the length of the <Backbone> (avoids recalculation if possible).
	 * @voxdim the physical dimensions of a voxel
	 * @voldim the extent of (number of voxels along) each axis
	 *
	 * Returns the length of the <Backbone>.  Does not recalculate the length if `len` is valid.
	 * @return the length of the <Backbone>
	 */
	// assumes vertebrae are ordered and voxdim are constant (to avoid recomputing)
	template <class S, class R>
	S length(const S (&voxdim)[3], const R (&voldim)[3]){
		if(len > 0)
			return len;
		
		if(vertebrae.empty())
			return 0;
		
		if(vertebrae.size() == 1)
			return (voxdim[0] + voxdim[1] + voxdim[2])/6.0; // half of the average voxel dimension
		
		len = 0;
		for(unsigned int i(1); i < vertebrae.size(); i++)
			len += separation3D(vertebrae[i - 1], vertebrae[i], voxdim, voldim);
		return len; // basically distance if we were to walk along the backbone (straightened out)
	}
	
	/* Returns the length of the <Backbone> (always recalculates).
	 * @voxdim the physical dimensions of a voxel
	 * @voldim the extent of (number of voxels along) each axis
	 *
	 * Returns the length of the <Backbone>.  Always recalculates the length.
	 * @return the length of the <Backbone>
	 */
	template <class S, class R>
	S newLength(const S (&voxdim)[3], const R (&voldim)[3]){
		len = -1;
		return length(voxdim, voldim);
	}
	
	/* Organizes the backbone by ordering the vertebrae in a sequence from tip to tip.
	 * @B <BinaryVolume> with volume extent information
	 *
	 * Organizes the backbone by ordering the vertebrae in a sequence from tip to tip.  Assumes that there are no more than two vertebrae with less than two neighbors, and that there are no vertebrae with more than two neighbors, and that the backbone is a single connected component.
	 */
	void organize(const BinaryVolume &B){
		
		if(vertebrae.size() < 3) // nothing to organize
			return;
		
		vector<T> ub(vertebrae), // the unorganized backbone
				neighbors(2*vertebrae.size(), 0); // list of neighbor indices in ub
				// critically, assumes at most 2 neighbors, linear chain
		vector<char> numNeighbors(vertebrae.size(), 0); // number of neighbors for each vertebra in ub (always less than 27)
		for(unsigned long v(0); v < vertebrae.size(); v++){
			if(numNeighbors[v] > 1)
				continue;
			for(unsigned long w(v + 1); w < vertebrae.size(); w++){
				if(inNeighborhood(vertebrae[v], vertebrae[w], B)){
					neighbors[2*v + numNeighbors[v]] = w;
					neighbors[2*w + numNeighbors[w]] = v;
					numNeighbors[v]++;
					numNeighbors[w]++;
					if(numNeighbors[v] > 1)
						w = vertebrae.size();
				}
			}
		}
		
		// find a tip
		unsigned long tipIndex(0);
		for(unsigned long v(0); v < vertebrae.size(); v++){
			if(numNeighbors[v] < 2){
				tipIndex = v;
				v = vertebrae.size();
			}
		}
		if(numNeighbors[tipIndex] > 1){
			cerr << "\n Warning Backbone::organize(): no tip found in Backbone with vertebra(e) " << makeString(vertebrae) << endl;
			return;
		}
		
		if(numNeighbors[tipIndex] < 1) // the tip has no neighbor
			vertebrae = vector<T>(1, ub[tipIndex]);
		
		// walk along neighbors until another tip
		unsigned long numVertebraeOrganized(2),
				prevIndex(tipIndex),
				curIndex(neighbors[2*tipIndex]),
				nextIndex;
		vertebrae[0] = ub[prevIndex];
		vertebrae[1] = ub[curIndex];
		while(numNeighbors[curIndex] > 1){
			nextIndex = (prevIndex == neighbors[2*curIndex]) ? neighbors[2*curIndex + 1] : neighbors[2*curIndex];
			vertebrae[numVertebraeOrganized] = ub[nextIndex];
			prevIndex = curIndex;
			curIndex = nextIndex;
			numVertebraeOrganized++;
		}
		
		if(numVertebraeOrganized < vertebrae.size()){
			cerr << "\n Warning Backbone:organize(): Did not include all vertebrae!" << endl;
			cout << "\n Warning Backbone:organize(): Did not include all vertebrae!" << endl;
			
			len = -1;
			if(numVertebraeOrganized < 1)
				vertebrae.clear();
			else
				vertebrae.erase(vertebrae.begin() + numVertebraeOrganized - 1, vertebrae.end());
		}
	}
	
};

#endif
