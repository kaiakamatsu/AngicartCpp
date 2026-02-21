#ifndef BINARYVOLUME_H
#define BINARYVOLUME_H 1

#include <vector>
#include <sstream>
#include <string>

using namespace std;

// these are type aliases
using voxelType = unsigned long long; // single index to identify each voxel (8 bytes (64 bits integer))
using branchType = vector<vector<voxelType> >::size_type; // get's the size type and indexing type for vector
using vecType = vector<voxelType>::size_type; // get's the size type and indexing type for vector

/* The <BinaryVolume> class holds true/false information in a 3D voxel space
 *
 * Implements a binary volume to hold true/false in a 3D voxel space, utilizing the bitset nature of std::vector<bool>.
 */
struct BinaryVolume{
	
private:
	vector<bool> bits; // flat 1D storage of all voxels
	voxelType totSize, size[3];
	
public:
	
	/* Default constructor
	 */
	BinaryVolume() : totSize(0) {
		size[0] = size[1] = size[2] = 0;
	}
	
	/* Copy constructor
	 */
	BinaryVolume(const BinaryVolume &B) : bits(B.bits), totSize(B.totSize) { // vector knows how to copy itself and make sure there are no memory leaks.
		size[0] = B.size[0];
		size[1] = B.size[1];
		size[2] = B.size[2];
	}
	
	/* Constructor for a volume with specific extent and optional truthiness
	 * @x extent of the volume in the *x* dimension
	 * @y extent of the volume in the *y* dimension
	 * @z extent of the volume in the *z* dimension
	 * @tf optional truthiness setting (default false)
	 */
	template <class T>
	BinaryVolume(T x, T y, T z, bool tf = false) : bits(x*y*z, tf), totSize(x*y*z) {
		size[0] = x;
		size[1] = y;
		size[2] = z;
	}
	
	/* Constructor for a volume with specific extent and optional truthiness
	 * @givenSize array with three elements to define the extent of the volume in each dimensions
	 * @tf optional truthiness setting (default false)
	 */
	BinaryVolume(const unsigned int givenSize[3], bool tf = false) : bits(givenSize[0]*givenSize[1]*givenSize[2], tf), totSize(givenSize[0]*givenSize[1]*givenSize[2]) {
		size[0] = givenSize[0];
		size[1] = givenSize[1];
		size[2] = givenSize[2];
	}
	
	/* Checks whether the bit corresponding to voxel (`x`, `y`, `z`) is (true)
	 * @x the *x* coordinate of the voxel of interest
	 * @y the *y* coordinate of the voxel of interest
	 * @z the *z* coordinate of the voxel of interest
	 *
	 * Returns the state of the corresponding voxel (`x`, `y`, `z`).  Voxels outside the volume are assumed to be false.
	 *
	 * @return the state of the bit corresponding to voxel (`x`, `y`, `z`)
	 */
	template <class T>
	bool is(T x, T y, T z) const {
		if(x >= size[0])
			return false;
		if(y >= size[1])
			return false;
		if(z >= size[2])
			return false;
		return bits[indexOf(x, y, z)];
	}
	
	 /* Checks whether the bit corresponding to voxel `i` is (true)
	 * @i the index in this <BinaryVolume> of the voxel of interest
	 *
	 * Returns the state of the corresponding voxel `i`.  Voxels outside the volume are assumed to be false.
	 *
	 * @return the state of the bit corresponding to voxel `i`
	 */
	template <class T>
	bool is(T i) const {
		if(i < totSize)
			return bits[i];
		return false;
	}
	
	/* Sets the bit corresponding to voxel (`x`, `y`, `z`) to the desired state (default true)
	 * @x the *x* coordinate of the voxel of interest
	 * @y the *y* coordinate of the voxel of interest
	 * @z the *z* coordinate of the voxel of interest
	 * @tf optional state to set (default true)
	 *
	 * Sets the state of the corresponding voxel (`x`, `y`, `z`) to state `tf`.
	 */
	template <class T>
	void set(T x, T y, T z, bool tf = true){
		if(x < size[0] && y < size[1] && z < size[2])
			bits[indexOf(x, y, z)] = tf;
	}
	
	 /* Sets the bit corresponding to voxel `i` to the desired state (default true)
	 * @i the index in this <BinaryVolume> of the voxel of interest
	 * @tf optional state to set (default true)
	 *
	 * Sets the state of the corresponding voxel `i` to state `tf`.
	 */
	template <class T>
	void set(T i, bool tf = true){
		if(i < totSize)
			bits[i] = tf;
	}
	
	/* Sets the bit corresponding to voxel (`x`, `y`, `z`) to true
	 * @x the *x* coordinate of the voxel of interest
	 * @y the *y* coordinate of the voxel of interest
	 * @z the *z* coordinate of the voxel of interest
	 *
	 * Sets the state of the corresponding voxel (`x`, `y`, `z`) to true.
	 */
	template <class T>
	void t(T x, T y, T z){
		set(x, y, z, true);
	}
	
	/* Sets the bit corresponding to voxel `i` to true
	 * @i the index in this <BinaryVolume> of the voxel of interest
	 *
	 * Sets the state of the corresponding voxel `i` to true.
	 */
	template <class T>
	void t(T i){
		set(i, true);
	}
	
	/* Sets the bits corresponding to voxels in `v` to true
	 * @v a list of the indices in this <BinaryVolume> to set to true
	 *
	 * Sets the state of the corresponding voxels in `v` to true.
	 */
	template <class T>
	void t(const vector<T> &v){
		for(unsigned long long i(0); i < v.size(); i++)
			t(v[i]);
	}
	
	 /* Sets the bit corresponding to voxel (`x`, `y`, `z`) to false
	 * @x the *x* coordinate of the voxel of interest
	 * @y the *y* coordinate of the voxel of interest
	 * @z the *z* coordinate of the voxel of interest
	 *
	 * Sets the state of the corresponding voxel (`x`, `y`, `z`) to false.
	 */
	template <class T>
	void f(T x, T y, T z){
		set(x, y, z, false);
	}
	
	/* Sets the bit corresponding to voxel `i` to false
	 * @i the index in this <BinaryVolume> of the voxel of interest
	 *
	 * Sets the state of the corresponding voxel `i` to false.
	 */
	template <class T>
	void f(T i){
		set(i, false);
	}
	
	/* Sets the bits corresponding to voxels in `v` to false
	 * @v a list of the indices in this <BinaryVolume> to set to false
	 *
	 * Sets the state of the corresponding voxels in `v` to false.
	 */
	template <class T>
	void f(const vector<T> &v){
		for(unsigned long long i(0); i < v.size(); i++)
			f(v[i]);
	}
	
	/* Sets all bits to false.
	 */
	void clear(){
		bits = vector<bool>(totSize, false);
	}
	
	/* Sets all of bits to true.
	 */
	void fill(){
		bits = vector<bool>(totSize, true);
	}
	
	/* Finds the first bit at or after index `i` that is true
	 * @i the index in the <BinaryVolume> to start the search
	 *
	 * Sequentially searches through the bits in the <BinaryVolume> starting at `i` until the end is reached or a bit in the true state is found.  If the end is reached without finding a true bit, the total number of voxels in the volume is returned.
	 *
	 * @return the index of the first bit in the true state at or after index `i`; if none found, then returns the total <BinaryVolume> volume
	 */
	template <class T>
	voxelType findFirstAtOrAfter(T iT) const {
		voxelType i(iT), iMax(totSize);
		while(i < iMax){
			if(bits[i])
				return i;
			i++;
		}
		return totSize;
	}
	
	 /* Finds the first bit at or after index `i` that is false
	 * @i the index in the <BinaryVolume> to start the search
	 *
	 * Sequentially searches through the bits in the <BinaryVolume> starting at `i` until the end is reached or a bit in the false state is found.  If the end is reached without finding a false bit, the total number of voxels in the volume is returned.
	 *
	 * @return the index of the first bit in the false state at or after index `i`; if none found, then returns the total <BinaryVolume> volume
	 */
	template<class T>
	voxelType findFirstFalseAtOrAfter(T i) const {
		voxelType iMax(totalSize());
		while(i < iMax){
			if(bits[i])
				i++;
			else
				return i;
			
		}
		return totSize;
	}
	
	/* Returns the total number of voxels in the <BinaryVolume>.
	 *
	 * Returns the total number of voxels in the <BinaryVolume>, corresponding to size[0]*size[1]*size[3].
	 * @return the total number of voxels in the <BinaryVolume>
	 */
	unsigned long long totalSize() const {
		return totSize;
	}
	
	/* Returns the total number of voxels in the true state in the <BinaryVolume>.
	 *
	 * Returns the total number of voxels in the true state in the <BinaryVolume> by inspecting each bit.
	 * @return the total number of voxels in the true state in the <BinaryVolume>
	 */
	voxelType totalTrue() const {
		voxelType totTrue(0);
		for(voxelType i(0); i < totSize; i++){
			if(bits[i])
				totTrue++;
		}
		return totTrue;
	}
	
	/* Returns the index of the point (`x`, `y`, `z`).
	 * @x the *x* coordinate of the voxel of interest
	 * @y the *y* coordinate of the voxel of interest
	 * @z the *z* coordinate of the voxel of interest
	 *
	 * Returns the index of the point (`x`, `y`, `z`).
	 * @return the index of the point (`x`, `y`, `z`)
	 */
	template <class T>
	voxelType indexOf(T x, T y, T z) const { // it is critical that we understand this is a 1D array to make sense of this indexing. 
		return x*size[1]*size[2] + y*size[2] + z;
	}
	
	/* Flips every bit in `bits`.
	 */
	void flip(){
		bits.flip();
	}
	
	/* Returns the *x* coordinate of the index `i`.
	 * @i the index of interest
	 *
	 * Returns the *x* coordinate of the index `i`.
	 * @return the *x* coordinate of the index `i`
	 */
	template <class T>
	voxelType x(T i) const {
		return i/(size[1]*size[2]);
	}
	
	/* Returns the *y* coordinate of the index `i`.
	 * @i the index of interest
	 *
	 * Returns the *y* coordinate of the index `i`.
	 * @return the *y* coordinate of the index `i`
	 */
	template <class T>
	voxelType y(T i) const {
		return (i%(size[1]*size[2]))/size[2];
	}
	
	/* Returns the *z* coordinate of the index `i`.
	 * @i the index of interest
	 *
	 * Returns the *z* coordinate of the index `i`.
	 * @return the *z* coordinate of the index `i`
	 */
	template <class T>
	voxelType z(T i) const {
		return i%size[2];
	}
	
	/* Returns the size of the <BinaryVolume> in the `i`-th direction.
	 *
	 * Returns the size of the <BinaryVolume> in the `i`-th direction.
	 * @return the size of the <BinaryVolume> in the `i`-th direction.
	 */
	template <class T>
	voxelType getSize(T i) const {
		return size[i];
	}
	
};

#endif
