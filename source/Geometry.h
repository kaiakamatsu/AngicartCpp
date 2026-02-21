#ifndef GEOMETRY_H
#define GEOMETRY_H 1

#include <cmath>

#include "BinaryVolume.h"
#include "Lumens.h"
#include "utilMinSurfTests.h"

/* Returns the absolute difference between `x` and `y`.
 * @x the value to subtract `y` from
 * @y the value to subtract from `x`
 *
 * Returns the absolute difference between `x` and `y` (useful for unsigned types).
 *
 * @return the absolute difference between `x` and `y`.
 */
template <class T, class S>
T absDiff(T x, S y){
	if(x > y)
		return x - y;
	return y - x;
}

/* Converts the <BinaryVolume> index `i` into coordinate triplet (`x`, `y`, `z`).
 * @i <BinaryVolume> index of interest
 * @x the *x*-coordinate of the voxel of interest
 * @y the *y*-coordinate of the voxel of interest
 * @z the *z*-coordinate of the voxel of interest
 * @volDim the extent of the volume in all three directions
 *
 * Converts the <BinaryVolume> index `i` into coordinate triplet (`x`, `y`, `z`) given volume extent from `volDim`.
 *
 * @return coordinate triplet (`x`, `y`, `z`) of <BinaryVolume> index `i`
 */
template <class T, class S, class R>
bool xyz(const T i, S &x, S &y, S &z, const R (&volDim)[3]){
	if(volDim[0] == 0)
		return false;
	if(volDim[1] == 0)
		return false;
	if(volDim[2] == 0)
		return false;
	x = i/(volDim[1]*volDim[2]);
	y = (i%(volDim[1]*volDim[2]))/volDim[2];
	z = i%volDim[2];
	return true;
}

/* Returns the square of the length between (`x`, `y`, `z`) and (`X`, `Y`, `Z`).
 * @x the *x* coordinate of the first point
 * @y the *y* coordinate of the first point
 * @z the *z* coordinate of the first point
 * @X the *x* coordinate of the second point
 * @Y the *y* coordinate of the second point
 * @Z the *z* coordinate of the second point
 *
 * Returns the square of the length between (`x`, `y`, `z`) and (`X`, `Y`, `Z`).
 *
 * @return the square of the length between (`x`, `y`, `z`) and (`X`, `Y`, `Z`)
 */
double separationSq3D(double x, double y, double z, double X, double Y, double Z){
	double dx(x - X),
			dy(y - Y),
			dz(z - Z);
	return dx*dx + dy*dy + dz*dz;
}

/* Returns the square of the physical length between voxels `i` and `j`.
 * @i an index in the <BinaryVolume> `B`
 * @j another index in the <BinaryVolume> `B`
 * @voxDim the dimensions of a single voxel
 * @volDim the extent of the volume in all three directions
 *
 * Returns the square of the physical length between voxels `i` and `j` given voxel size `voxDim` and volume extent from `volDim`.
 *
 * @return the square of the physical length between voxels `i` and `j`.
 */
template <class T, class S, class R>
double separationSq3D(T i, T j, const S (&voxDim)[3], const R (&volDim)[3]){
	voxelType xi, yi, zi, xj, yj, zj;
	xyz(i, xi, yi, zi, volDim);
	xyz(j, xj, yj, zj, volDim);
	double dx(voxDim[0]*absDiff(xi, xj)),
			dy(voxDim[1]*absDiff(yi, yj)),
			dz(voxDim[2]*absDiff(zi, zj));
	return dx*dx + dy*dy + dz*dz;
}

/* Returns the square of the physical length between voxels `i` and `j`.
 * @i an index in the <BinaryVolume> `B`
 * @j another index in the <BinaryVolume> `B`
 * @voxDim the dimensions of a single voxel
 * @B <BinaryVolume> that encodes volume extent
 *
 * Returns square of the physical length between voxels `i` and `j` given voxel size `voxDim` and volume extent from <BinaryVolume> `B`.
 *
 * @return the square of the physical length between voxels `i` and `j`.
 */
template <class T, class S>
double separationSq3D(T i, T j, const S (&voxDim)[3], const BinaryVolume &B){
	voxelType volDim[3] = {B.getSize(0), B.getSize(1), B.getSize(2)};
	return separationSq3D(i, j, voxDim, volDim);
}

/* Returns the physical length between voxels `i` and `j`.
 * @i an index in the <BinaryVolume> `B`
 * @j another index in the <BinaryVolume> `B`
 * @voxDim the dimensions of a single voxel
 * @volDim the extent of the volume in all three directions
 *
 * Returns the physical length between voxels `i` and `j` given voxel size `voxDim` and volume extent from `volDim`.
 *
 * @return the physical length between voxels `i` and `j`.
 */
template <class T, class R>
double separation3D(T i, T j, const double (&voxDim)[3], const R (&volDim)[3]){
	return sqrt(separationSq3D(i, j, voxDim, volDim));
}

/* Returns the physical length between voxels `i` and `j`.
 * @i an index in the <BinaryVolume> `B`
 * @j another index in the <BinaryVolume> `B`
 * @voxDim the dimensions of a single voxel
 * @B <BinaryVolume> that encodes volume extent
 *
 * Returns the physical length between voxels `i` and `j` given voxel size `voxDim` and volume extent from <BinaryVolume> `B`.
 *
 * @return the physical length between voxels `i` and `j`.
 */
template <class T>
double separation3D(T i, T j, const double (&voxDim)[3], const BinaryVolume &B){
	return sqrt(separationSq3D(i, j, voxDim, B));
}

/* Determines if indices `i` and `j` are within each other's 3x3x3 neighborhood.
 * @i an index in the <BinaryVolume> `B`
 * @j another index in the <BinaryVolume> `B`
 * @B the <BinaryVolume> that corresponds to the indices `i` and `j`
 *
 * Deternmines the neighborliness of the voxels with indices `i` and `j` within the <BinaryVolume> `B` by checking the absolute difference in each (x, y, z) coordinate.
 *
 * @return true if `i` and `j` are neighbors
 */
template <class T>
bool inNeighborhood(T i, T j, const BinaryVolume &B){
	return absDiff(B.x(i), B.x(j)) < 2 && absDiff(B.y(i), B.y(j)) < 2 && absDiff(B.z(i), B.z(j)) < 2;
}

/* Constructs a list of the blocks in the 3x3x3 neighborhood of `i` -- excluding the central point itself -- in <BinaryVolume> `B`.
 * @B the <BinaryVolume> that indicates the volume extent
 * @i the index in <BinaryVolume> `B` of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at `i` and adds all neighbors to the list for output.  Note that the central point `i` is intentionally excluded from this list.
 *
 * @return a list of voxel indices in <BinaryVolume> `B` that are neighbors of the voxel `i`
 */
vector<voxelType> blocksInNeighborhood(const BinaryVolume &B, voxelType index){
	voxelType x(B.x(index)), y(B.y(index)), z(B.z(index)), xStart(0), yStart(0), zStart(0);
	if(x > 1)
		xStart = x - 1;
	if(y > 1)
		yStart = y - 1;
	if(z > 1)
		zStart = z - 1;
	vector<voxelType> nh;
	nh.reserve(26);
	for(voxelType i(xStart); i < x + 2 && i < B.getSize(0); i++){
		for(voxelType j(yStart); j < y + 2 && j < B.getSize(1); j++){
			for(voxelType k(zStart); k < z + 2 && k < B.getSize(2); k++){
				if(i == x && j == y && k == z) // might be unnecessary when optimized...
					continue;
				nh.push_back(B.indexOf(i, j, k));
			}
		}
	}
	return nh;
}

/* Counts the number of vessel blocks in the 3x3x3 neighborhood of (`x`, `y`, `z`) -- excluding the central point itself -- that are true in <BinaryVolume> `B`.
 * @B <BinaryVolume> that indicates which voxels qualify as neighbors
 * @i the index of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at (`x`, `y`, `z`) and count qualifying neighbors (i.e., those that are true in <BinaryVolume> `B`).  Note that the central point (`x`, `y`, `z`) is intentionally excluded from this list.
 *
 * @return the number of voxels in <BinaryVolume> `B` that qualify as neighbors of the voxel (`x`, `y`, `z`) but are not the voxel itself
 */
unsigned int numVesselBlocksInNeighborhood(const BinaryVolume &B, voxelType index){
	voxelType x(B.x(index)), y(B.y(index)), z(B.z(index)), xStart(0), yStart(0), zStart(0);
	if(x > 1)
		xStart = x - 1;
	if(y > 1)
		yStart = y - 1;
	if(z > 1)
		zStart = z - 1;
	unsigned int numVesselBlocks(0);
	for(voxelType i(xStart); i < x + 2 && i < B.getSize(0); i++){
		for(voxelType j(yStart); j < y + 2 && j < B.getSize(1); j++){
			for(voxelType k(zStart); k < z + 2 && k < B.getSize(2); k++){
				if(i == x && j == y && k == z) // might be unnecessary when optimized...
					continue;
				if(B.is(i, j, k))
					numVesselBlocks++;
			}
		}
	}
	return numVesselBlocks;
}

/* Constructs a list of the blocks in the 3x3x3 neighborhood of (`x`, `y`, `z`) -- excluding the central point itself -- that are true in <BinaryVolume> `B`.
 * @B the <BinaryVolume> that indicates which voxels qualify as neighbors
 * @x the *x*-coordinate of the voxel of interest
 * @y the *y*-coordinate of the voxel of interest
 * @z the *z*-coordinate of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at (`x`, `y`, `z`) and adds qualifying neighbors (i.e., those that are true in <BinaryVolume> `B`) to the list for output.  Note that the central point (`x`, `y`, `z`) is intentionally excluded from this list.
 *
 * @return a list of voxel indices in <BinaryVolume> `B` that qualify as neighbors of the voxel (`x`, `y`, `z`) but are not the voxel itself
 */
template <class T>
vector<voxelType> vesselBlocksInNeighborhood(const BinaryVolume &B, T x, T y, T z){
	// defines indices of the 3 x 3 x 3 neighborhood while respecting the boundaries of the volume
	voxelType xStart((x > 0) ? x - 1 : 0),
			yStart((y > 0) ? y - 1 : 0),
			zStart((z > 0) ? z - 1 : 0),
			xEnd((x + 2 > B.getSize(0)) ? B.getSize(0) : x + 2),
			yEnd((y + 2 > B.getSize(1)) ? B.getSize(1) : y + 2),
			zEnd((z + 2 > B.getSize(2)) ? B.getSize(2) : z + 2);
	vector<voxelType> vesselBlockIndices;
	vesselBlockIndices.reserve(26); // pre-allocates memory for up to 26 neighbors, excludes itself. 
	for(voxelType i(xStart); i < xEnd; i++){
		for(voxelType j(yStart); j < yEnd; j++){
			for(voxelType k(zStart); k < zEnd; k++){
				if(B.is(i, j, k)){ // checks if the voxel is true in the binary volume
					if(i != x) // these three checks ensure the point (x, y, z) is not included (itself.)
						vesselBlockIndices.push_back(B.indexOf(i, j, k)); // indexOf gives you the index in the 1D binary array that corresponds to the voxel (i, j, k)
					else if(j != y)
						vesselBlockIndices.push_back(B.indexOf(i, j, k));
					else if(k != z)
						vesselBlockIndices.push_back(B.indexOf(i, j, k));
				}
			}
		}
	}
	return vesselBlockIndices;
}

/* Constructs a list of the blocks in the 3x3x3 neighborhood of (`x`, `y`, `z`) -- excluding the central point itself -- that are false in <BinaryVolume> `B`.
 * @B the <BinaryVolume> that indicates which voxels qualify as neighbors
 * @x the *x*-coordinate of the voxel of interest
 * @y the *y*-coordinate of the voxel of interest
 * @z the *z*-coordinate of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at (`x`, `y`, `z`) and adds qualifying neighbors (i.e., those that are false in <BinaryVolume> `B`) to the list for output.  Note that the central point (`x`, `y`, `z`) is intentionally excluded from this list.
 *
 * @return a list of nonvessel indices in <BinaryVolume> `B` that qualify as neighbors of the voxel (`x`, `y`, `z`) but are not the voxel itself
 */
template <class T>
vector<voxelType> nonvesselBlocksInNeighborhood(const BinaryVolume &B, T x, T y, T z){
	voxelType xStart((x > 0) ? x - 1 : 0),
			yStart((y > 0) ? y - 1 : 0),
			zStart((z > 0) ? z - 1 : 0),
			xEnd((x + 2 > B.getSize(0)) ? B.getSize(0) : x + 2),
			yEnd((y + 2 > B.getSize(1)) ? B.getSize(1) : y + 2),
			zEnd((z + 2 > B.getSize(2)) ? B.getSize(2) : z + 2);
	vector<voxelType> vesselBlockIndices;
	vesselBlockIndices.reserve(26);
	for(voxelType i(xStart); i < xEnd; i++){
		for(voxelType j(yStart); j < yEnd; j++){
			for(voxelType k(zStart); k < zEnd; k++){
				if(!B.is(i, j, k)){
					if(i != x) // these three checks ensure the point (x, y, z) is not included
						vesselBlockIndices.push_back(B.indexOf(i, j, k));
					else if(j != y)
						vesselBlockIndices.push_back(B.indexOf(i, j, k));
					else if(k != z)
						vesselBlockIndices.push_back(B.indexOf(i, j, k));
				}
			}
		}
	}
	return vesselBlockIndices;
}

/* Constructs a list of the blocks in the 3x3x3 neighborhood of `i` -- excluding the central point itself -- that are true in <BinaryVolume> `B`.
 * @B the <BinaryVolume> that indicates which voxels qualify as neighbors
 * @i the index in <BinaryVolume> `B` of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at `i` and adds qualifying neighbors (i.e., those that are true in <BinaryVolume> `B`) to the list for output.  Note that the central point `i` is intentionally excluded from this list.
 *
 * @return a list of voxel indices in <BinaryVolume> `B` that qualify as neighbors of the voxel `i`
 */
vector<voxelType> vesselBlocksInNeighborhood(const BinaryVolume &B, voxelType index){
	return vesselBlocksInNeighborhood(B, B.x(index), B.y(index), B.z(index));
}

/* Constructs a list of the blocks in the 3x3x3 neighborhood of `i` -- excluding the central point itself -- that are false in <BinaryVolume> `B`.
 * @B <BinaryVolume> that indicates which voxels qualify as neighbors
 * @i index in <BinaryVolume> `B` of the voxel of interest
 *
 * Checks every voxel in the 3x3x3 neighborhood centered at `i` and adds qualifying neighbors (i.e., those that are false in <BinaryVolume> `B`) to the list for output.  Note that the central point `i` is intentionally excluded from this list.
 *
 * @return a list of nonvessel indices in <BinaryVolume> `B` that qualify as neighbors of the voxel `i`
 */
vector<voxelType> nonvesselBlocksInNeighborhood(const BinaryVolume &B, voxelType index){
	return nonvesselBlocksInNeighborhood(B, B.x(index), B.y(index), B.z(index));
}

/* Finds the indices in `v` that are neighbors of `i` (includes `i` if `i` is in `v`).
 * @B <BinaryVolume> that indicates which voxels qualify as neighbors
 * @i index in <BinaryVolume> `B` of the voxel of interest
 * @v list of possible neighbor voxels
 *
 * Finds the indices in `v` that are neighbors of `i` (includes `i` if `i` is in `v`).
 *
 * @return a list of the voxels in `v` that are neighbors of `i` (includes `i` if `i` is in `v`)
 */
vector<voxelType> findNeighbors(const BinaryVolume &B, voxelType i, const vector<voxelType> &v){
	vector<voxelType> nh;
	for(unsigned int j(0); j < v.size(); j++){
		if(inNeighborhood(i, v[j], B))
			nh.push_back(v[j]);
	}
	return nh;
}

/* Finds the radius for a cylinder with the given volume and length.
 * @vol cylinder volume
 * @len cylinder length
 *
 * Computes the radius of a cylinder with volume `vol` and length `len`.
 *
 * @return cylinder radius
 */
template <class S, class T>
double radFromVolLen(S vol, T len){
	return sqrt((vol/len)/acos(-1.0));
}

#endif
