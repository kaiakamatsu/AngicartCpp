#define TRY_WX 0	//	compiles with wxWidgets window if TRY_WX is set to 1; this can be partially avoided by setting TRY_WX to 0

#include <algorithm>
#include <chrono>
#include <future>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <atomic>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "BinaryVolume.h"
#include "Backbone.h"
#include "Lumens.h"
#include "utilMinSurfTests.h"
#include "Geometry.h"

#include "lodepng.h"
#include "pngMinSurf.h"

unsigned int ccSelected(0), bbSelected(0), vertSelected(0), labelStatus(0);
bool segIsSelected(false), updateSelectedStatus(false);

chrono::steady_clock::time_point startTime;
unsigned int numLocThreads(1u); // set from command line or default (hardware_concurrency - 1) in main 

vector<vector<vector<branchType> > > branchpointsGlobal;
vector<vector<Backbone<> > > backbonesGlobal;
vector<vector<double> > volumesGlobal, aveRadSolidGlobal, aveRadSurfGlobal;
vector<vector<vector<double> > > aveRadSolidByVertebraGlobal, aveRadSurfByVertebraGlobal;
vector<vector<vector<bool> > > branchAtFrontGlobal;
double voxdimsGlobal[3] = {0.0, 0.0, 0.0};
voxelType voldimsGlobal[3] = {0, 0, 0};
double defaultPointColor[] = {0.0, 0.0, 0.0, 0.2};
string lengthUnitGlobal("units");

#if TRY_WX == 1
#include "wxMinSurfTests.h"
#endif

/* Updates the displayed status message.
 * @newStatus the new status message
 *
 * Updates the displayed status message to cout and cerr, as well as to the wx window, if applicable.
 */
void generalStatusUpdate(string newStatus){
	cout << endl << newStatus << endl;
	cerr << endl << newStatus << endl;
#if TRY_WX == 1
	wxGetApp().myFrame->updateStatus(newStatus);
#endif
}

/* Produces a <BinaryVolume> from the intensities in <Lumens> `L` based on the threshold `thresh`.
 * @L <Lumens> intensities
 * @thresh intensity threshold
 *
 * Simplifies the intensity information from <Lumens> `L` by setting voxels with intensity greater than `thresh` to true and less than `thresh` to false in the returned <BinaryVolume> `B`.
 *
 * @return A <BinaryVolume> `B`
 */
BinaryVolume threshThrash(const Lumens &L, double thresh){
	double minL(L.minLumen()), maxL(L.maxLumen());
	BinaryVolume B(L.size[0], L.size[1], L.size[2]);
	if(minL == maxL)
		return B;
	for(unsigned int x(0); x < L.size[0]; x++){
		for(unsigned int y(0); y < L.size[1]; y++){
			for(unsigned int z(0); z < L.size[2]; z++){
				if((L.lumens[x][y][z] - minL)/(maxL - minL) > thresh) // looks at normalized intensity of each voxel
					B.t(x, y, z); // sets to true, see BinaryVolume.h for more details.
			}
		}
	}
	return B;
}

/* Determines the connectivity of the voxels `nh` in <BinaryVolume> `B`.
 * @nh the set of voxels of interest for determining connectivity
 * @B the <BinaryVolume> that contains `nh`
 *
 * Searches a connected component of the voxels of `nh` in <BinaryVolume> `B` to see if it spans `nh`.
 *
 * @return true if the specified neighborhood is indeed a single connected component
 */
bool isSingleConnectedComponent(const vector<voxelType> &nh, const BinaryVolume &B){
	vector<voxelType> f;
	f.reserve(nh.size());
	f.push_back(nh[0]);
	vector<bool> visited(nh.size(), false);
	visited[0] = true;
	for(unsigned int fi(0); fi < f.size(); fi++){
		for(unsigned int j(1); j < nh.size(); j++){ // can start with 1 because nh[0] is the seed
			if(visited[j])
				continue;
			
			if(inNeighborhood(f[fi], nh[j], B)){
				visited[j] = true;
				f.push_back(nh[j]);
			}
		}
	}
	return f.size() == nh.size();
}

/* Finds the center of mass of the given set of indices (with identical weights)
 * @B <BinaryVolume> that encodes volume extent
 * @vi list of indices in <BinaryVolume> `B` of the object of interest
 *
 * Finds the center of mass in <BinaryVolume> `B` of the object specified by the voxel indices in the list `vi`.  Note that each voxel has the same weight (mass).
 *
 * @return the index in <BinaryVolume> `B` of the voxel that contains the center of mass of the voxels in `vi`.
 */
template <class T>
voxelType centerOfMass(const BinaryVolume &B, const vector<T> &vi){
	if(vi.size() < 1)
		return B.totalSize();
	unsigned int xSum(0), ySum(0), zSum(0);
	for(unsigned int i(0); i < vi.size(); i++){
		xSum += B.x(vi[i]);
		ySum += B.y(vi[i]);
		zSum += B.z(vi[i]);
	}
	return B.indexOf(xSum/vi.size(), ySum/vi.size(), zSum/vi.size());
}

/* Finds the voxen in the given set of indices that is nearest the center of mass of the set of indices (with identical weights)
 * @B the <BinaryVolume> that encodes volume extent
 * @vi the list of indices in <BinaryVolume> `B` of the object of interest
 *
 * Finds the voxel in the given list of voxels `vi` that is nearest to the center of mass in <BinaryVolume> `B` of the object specified by `vi`.  For empty sets, the total size of `B` is returned.  Note that each voxel has the same weight (mass).
 *
 * @return the index in <BinaryVolume> `B` of the voxel in `vi` that is nearest the center of mass of the voxels in `vi`.
 */
template <class T>
voxelType centerOfMassInSet(const BinaryVolume &B, const vector<T> &vi){
	if(vi.empty())
		return B.totalSize();
	
	voxelType com(centerOfMass(B, vi));
	if(com == B.totalSize())
		return com;
	
	voxelType bestVoxel(vi[0]);
	double bestDistanceSq(separationSq3D(com, vi[0], voxdimsGlobal, B)); // EuclideanDistanceSq(B, com, vi[0]));
	for(unsigned int i(1); i < vi.size(); i++){
		double thisDistanceSq(separationSq3D(com, vi[i], voxdimsGlobal, B));//EuclideanDistanceSq(B, com, vi[i]));
		if(bestDistanceSq > thisDistanceSq){
			bestDistanceSq = thisDistanceSq;
			bestVoxel = vi[i];
		}
	}
	return bestVoxel;
}

/* Finds the minimum separation between voxel `i` and the closest voxel in `v`.
 * @i <BinaryVolume> index of the voxel of interest
 * @v set of voxels to compare with `i`
 *
 * Finds the minimum separation between voxel `i` and the closest voxel in `v` by an exhaustive search.
 * @return the minimum separation between voxel `i` and the closest voxel in `v`
 */
double minimumSeparation3D(voxelType i, const vector<voxelType> &v){
	if(v.empty())
		return -1.0;
	
	double minSepSq(separationSq3D(i, v[0], voxdimsGlobal, voldimsGlobal));
	for(unsigned long j(0); j < v.size(); j++){
		double sepSq(separationSq3D(i, v[j], voxdimsGlobal, voldimsGlobal));
		if(minSepSq > sepSq)
			minSepSq = sepSq;
	}
	
	return sqrt(minSepSq);
}

/* Finds the average solid and surface radius of the meat associated with <Backbone> `organizedBackbone`.
 * @B <BinaryVolume> that indicates vasculature
 * @sMap the meat associated with <Backbone> `organizedBackbone`
 * @organizedBackbone the <Backbone> of interest
 * @aveRadSolid (output) average radius of the solid of the meat associated with <Backbone> `organizedBackbone`.
 * @aveRadSurf (output) average radius of the surface of the meat associated with <Backbone> `organizedBackbone`.
 *
 * Finds the average solid and surface radius of the meat associated with <Backbone> `organizedBackbone`.  Any voxel that does not have a half-filled neighborhood (with 14 other voxels) is considered a surface voxel.
 */
void averageRadius(const BinaryVolume &B, const map<voxelType, double> &sMap, const Backbone<> &organizedBackbone, double &aveRadSolid, double &aveRadSurf){
	aveRadSolid = aveRadSurf = 0;
	if(organizedBackbone.empty())
		return;
	
	double solidVol(0), surfVol(0);
	for(map<voxelType, double>::const_iterator it(sMap.begin()); it != sMap.end(); it++){  // for every voxel, it->second has the fraction of that voxel that belongs to the backbone of interest. 
		double minSep(minimumSeparation3D(it->first, organizedBackbone.getVertebrae())); // distance from the voxel to the closest vertebra voxel
		aveRadSolid += minSep*it->second; 
		solidVol += it->second;
		if(numVesselBlocksInNeighborhood(B, it->first) < 14){ // valid discriminator because internal gaps have been filled
			aveRadSurf += minSep*it->second;
			surfVol += it->second;
		}
	}
	
	if(solidVol > 0)
		aveRadSolid /= solidVol;
	if(surfVol > 0)
		aveRadSurf /= surfVol;
}

/* Finds the average solid and surface radius profile along <Backbone> `organizedBackbone`.
 * @B <BinaryVolume> that indicates vasculature
 * @sMap the meat associated with <Backbone> `organizedBackbone`
 * @organizedBackbone the <Backbone> of interest
 * @aveRadSolidByVertebra (output) average solid radius at each vertebra in `organizedBackbone`
 * @aveRadSurfByVertebra (output) average surface radius at each vertebra in `organizedBackbone`
 *
 * A version of averageRadius that computes and stores the average solid radius and average surface radius at every point along the backbone, instead of averaging over all voxels in the backbone's territory.
 * Each voxel in `sMap` is assigned to its nearest vertebra in `organizedBackbone`, and weighted averages are computed per vertebra using the same voxel-fraction weighting and surface criterion as the scalar version.
 */
void averageRadiusByVertebra(const BinaryVolume &B, const map<voxelType, double> &sMap, const Backbone<> &organizedBackbone, vector<double> &aveRadSolidByVertebra, vector<double> &aveRadSurfByVertebra){
	aveRadSolidByVertebra.clear(); // this is a reference to a specific vector of doubles for this backbone so clearing this is safe and will not affect the other backbones being processed by other threads. 
	aveRadSurfByVertebra.clear();
	if(organizedBackbone.empty())
		return;
	
	vector<voxelType> vertebrae(organizedBackbone.getVertebrae());
	aveRadSolidByVertebra = vector<double>(vertebrae.size(), 0.0);
	aveRadSurfByVertebra = vector<double>(vertebrae.size(), 0.0);
	vector<double> solidVol(vertebrae.size(), 0.0), surfVol(vertebrae.size(), 0.0);
	
	for(map<voxelType, double>::const_iterator it(sMap.begin()); it != sMap.end(); it++){
		unsigned int nearestIndex(0);
		double minSepSq(separationSq3D(it->first, vertebrae[0], voxdimsGlobal, voldimsGlobal));
		for(unsigned int i(1); i < vertebrae.size(); i++){
			double sepSq(separationSq3D(it->first, vertebrae[i], voxdimsGlobal, voldimsGlobal));
			if(minSepSq > sepSq){
				minSepSq = sepSq;
				nearestIndex = i;
			}
		} // find the nearest vertebra voxel to the current voxel
		
		double minSep(sqrt(minSepSq));
		aveRadSolidByVertebra[nearestIndex] += minSep*it->second;
		solidVol[nearestIndex] += it->second;
		if(numVesselBlocksInNeighborhood(B, it->first) < 14){
			aveRadSurfByVertebra[nearestIndex] += minSep*it->second;
			surfVol[nearestIndex] += it->second;
		}
	}
	
	for(unsigned int i(0); i < vertebrae.size(); i++){
		if(solidVol[i] > 0.0)
			aveRadSolidByVertebra[i] /= solidVol[i];
		if(surfVol[i] > 0.0)
			aveRadSurfByVertebra[i] /= surfVol[i];
	} // aveRadSolidByVertebra[i] corresponds to the average solid radius at the i-th voxel along the backbone. 
}

/* Grows the sphere specified by `seed` to determine if the appropriate sphere is disjoint from existing spheres
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that specifies unvisited voxels
 * @fB <BinaryVolume> that specifies which voxels have not been visited by the frontier
 * @seed the index in <BinaryVolume> `B` that seeds the sphere
 * @critFrac the critical fraction of new voxels that are not vasculature that terminates growth
 * @r the radius of the sphere
 * @falsifiedPoints the <BinaryVolume> indices for the set of voxels that belong to the sphere and were falsified in the frontier search
 *
 * Grows a sphere centered at `seed` while at least `critFrac` fraction of voxels in the newly-grown shell (e.g., between `r=3` and `r=4`) are vascular.  Returns a boolean to indicate whether or not the sphere is disjoint, and modifies r to be the final radius of the sphere and modifies falsifiedPoints to be the points flipped to false.
 * @return true if the sphere is disjoint
 */
bool disjointCriticalSphere(const BinaryVolume &B, BinaryVolume &uB, BinaryVolume &fB, voxelType seed, double critFrac, double &r, vector<voxelType> &falsifiedPoints
#if TRY_WX == 1
		, const map<voxelType, unsigned int> &pm
#endif
		){

	// uB and fB are and can be modified in this function. 
	// r is also updated
	// uB is necessary so we can find non-overlapping spheres. 
	
	if(!uB.is(seed)) // seed already visited
		return false;
	
	
	falsifiedPoints.clear();
	vector<voxelType> f(1, seed); // frontier list, starting with seed
	vector<double> d(1, 0); // distances of each voxel in frontier from seed
	
	fB.fill(); // mark all voxels in fB as true (available)
	fB.f(f[0]); // mark seed as visited in fB
	
#if TRY_WX == 1
	map<voxelType, unsigned int>::const_iterator pmit(pm.find(f[0]));
	if(pmit != pm.end()){
		map<unsigned int, vector<double> > uc;
		unsigned int pmIndex(pmit->second);
		uc[pmIndex] = vector<double>(4, 1.0); uc[pmIndex][0] = uc[pmIndex][1] = 0.0; uc[pmIndex][3] = 0.1;
		wxGetApp().myFrame->m_canvas->updateCols(uc);
	}
#endif
	
	double frac(1.0), rInc(sqrt(voxdimsGlobal[0]*voxdimsGlobal[0] + voxdimsGlobal[1]*voxdimsGlobal[1] + voxdimsGlobal[2]*voxdimsGlobal[2]));
	// rInc is the increment in radius based on the voxel dimensions.
	voxelType trueCount(1), count(1);
	r = -rInc/2;
	while(frac > critFrac && count > 0){ // Loop continues while the fraction of vascular voxels in the shell is above critFrac
		
		r = r + rInc;
		trueCount = count = 0;
		
		while(d[0] - r < 1.0e-9){ // while the distance from the seed to the frontier is less than the radius, we can grow the sphere.
			
			count++;
			if(B.is(f[0]))
				trueCount++;
			vector<voxelType> nh(vesselBlocksInNeighborhood(fB, f[0]));
			
			for(unsigned int i(0); i < nh.size(); i++){
				
				if(B.is(nh[i]) && !uB.is(nh[i])){ // has been visited before (not unvisited)
					vector<voxelType> trueFrontiers;
					for(unsigned int i(0); i < f.size(); i++){
						if(B.is(f[i]))
							trueFrontiers.push_back(f[i]);
					}
#if TRY_WX == 1
					setDefaultColor(trueFrontiers, pm);
#endif
					// fB.t(trueFrontiers);
					return false; // If a neighboring voxel is vascular but already visited, the sphere would overlap another sphere, so this is a failed sphere
				}// in what follows, nh[i] has not been visited before in any sphere
				
				// insert nh[i] and c into the sorted lists f and d, respectively
				double c(separation3D(seed, nh[i], voxdimsGlobal, B)); // distance from seed to neighbor
				unsigned int hi((unsigned int)d.size()), lo(0), mid(hi/2);
				while(hi - lo > 1){ // binary search to find the correct position to insert the neighbor
					if(d[mid] < c)
						lo = mid;
					else if(d[mid] == c)
						hi = lo = mid;
					else // d[mid] > c
						hi = mid;
					mid = lo + (hi - lo)/2;
				}
				
				unsigned int j(lo);
				if(d[lo] < c)
					j = hi;
				f.insert(f.begin() + j, nh[i]); // insert neighbor into frontier list
				d.insert(d.begin() + j, c); // insert distance to seed
				fB.f(nh[i]); // visited by frontier
				
#if TRY_WX == 1
				if(uB.is(nh[i])){
					
						map<voxelType, unsigned int>::const_iterator pmit(pm.find(nh[i]));
					if(pmit != pm.end()){
						map<unsigned int, vector<double> > uc;
						unsigned int pmIndex(pmit->second);
						
						uc[pmIndex] = vector<double>(4, 1.0); uc[pmIndex][0] = uc[pmIndex][1] = 0.0; uc[pmIndex][3] = 0.7;
						
						wxGetApp().myFrame->m_canvas->updateCols(uc);
					
					}
				}
#endif
				
			} // end of iterating over neighbors
			
			if(uB.is(f[0])){
				uB.f(f[0]); // set to visited
				falsifiedPoints.push_back(f[0]);
#if TRY_WX == 1
				map<voxelType, unsigned int>::const_iterator pmit(pm.find(f[0]));
				if(pmit != pm.end()){
					map<unsigned int, vector<double> > uc;
					unsigned int pmIndex(pmit->second);
					uc[pmIndex] = vector<double>(4, 1.0); uc[pmIndex][0] = uc[pmIndex][2] = 0.0; uc[pmIndex][3] = 0.7;
					wxGetApp().myFrame->m_canvas->updateCols(uc);
				}
#endif
			}
			f.erase(f.begin()); // remove the processed voxel from the frontier list
			d.erase(d.begin()); // remove the processed voxel from the distances list
			if(d.empty()){ // no frontier left
				vector<voxelType> trueFrontiers;
				for(unsigned int i(0); i < f.size(); i++){
					if(B.is(f[i]))
						trueFrontiers.push_back(f[i]);
				}
#if TRY_WX == 1
				setDefaultColor(trueFrontiers, pm);
#endif
				// fB.t(trueFrontiers);
				return false;
			}
		}
		frac = double(trueCount)/count; // update the fraction of vascular voxels in the shell
	} // end of while loop for sphere expansion
	vector<voxelType> trueFrontiers;
	for(unsigned int i(0); i < f.size(); i++){
		if(B.is(f[i]))
			trueFrontiers.push_back(f[i]);
	}
	
#if TRY_WX == 1
	setDefaultColor(trueFrontiers, pm);
#endif
	// fB.t(trueFrontiers);
	return true;
}

/* Roughly vectorizes the vasculature in <BinaryVolume> `B` into a map of seed to radius of spheres.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> of vasculature that has not yet been visited by a sphere (i.e., intersphere space)
 * @critFrac the critical fraction of new voxels that are not vasculature that terminates growth
 * @rSeed the red value of the color for the seed region
 * @gSeed the green value of the color for the seed region
 * @bSeed the blue value of the color for the seed region
 * @aSeed the alpha value of the color for the seed region
 *
 * Roughly vectorizes the vasculature in <BinaryVolume> `B` into a map of index of the seed to radius of spheres by choosing random seed points until too many consecutive seed points fail.  Modifies `uB` to be the unvisited (i.e., intersphere) vasculature.
 * @return a map of the <BinaryVolume> index of the seed to the radius of the disjoint spheres
 */
map<voxelType, double> sphereCoarsen(const BinaryVolume &B, BinaryVolume &uB, double critFrac
#if TRY_WX == 1
		, const map<voxelType, unsigned int> &pm
#endif
		, double rSeed = 0.0, double gSeed = 1.0, double bSeed = 1.0, double aSeed = 1.0){
	chrono::steady_clock::time_point sphereTime(getNow());
	map<voxelType, double> rv; // map seed voxel index to the radius of the sphere 
	unsigned int consecutiveFailures(0), maxSeedFailures(100000), maxConsecutiveFailures(200), reachedConsecutiveFailures(0); // how many times to try before returning the map.
	double minSeedMovement(1.0); // minimum distance a seed must move during centering
	BinaryVolume sB(uB), // for untried seeds
		fB(B.getSize(0), B.getSize(1), B.getSize(2), true); // for possible sphere volumes
	voxelType numUntriedSeeds(sB.totalTrue());

	// Portable RNG: full range on all platforms (avoids problems with Windows RAND_MAX=32767)
	// previous implementation had rand()%sB.totalSize()
	//std::uniform_int_distribution<voxelType> dist(0, sB.totalSize() - 1); // works but different across systems even with the same random seed number.
	std::mt19937 gen(112358u); // 112358u is another srandSet. mt19937 produces a random 32-bit integer when gen() is called
	uint64_t range = sB.totalSize();   // assumes range <= 2^32. promotes to uint64 to avoid overflow

	while(consecutiveFailures < maxConsecutiveFailures && numUntriedSeeds > 0){
		// Loop continues until too many consecutive seed attempts fail or there are no untried seeds left
		// WARNING: runtime is sensitive to this number. but I think it is worth it for the quality of the coarsening
		if(reachedConsecutiveFailures < consecutiveFailures)
			reachedConsecutiveFailures = consecutiveFailures;
		
#if TRY_WX == 1
		int paddedDigits(1 + (int)floor(log10(maxConsecutiveFailures)));
		wxGetApp().myFrame->updateStatus("Spheroidal segmentation... "
			+ paddedInt(consecutiveFailures, paddedDigits) + " up to "
			+ paddedInt(reachedConsecutiveFailures, paddedDigits) + " of "
			+ makeString(maxConsecutiveFailures) + "."
			+ "  Found " + makeString(rv.size()) + " points so far..."
			);
#endif		
		// sample a random seed
		uint32_t raw = gen();
		uint64_t value = (uint64_t(raw) * range) >> 32; // right shift by 32 bits is effectively dividing by 2^32
		voxelType seed(sB.findFirstAtOrAfter(value)), consecutiveSeedFailures(0);
		//voxelType seed(sB.findFirstAtOrAfter(dist(gen))), consecutiveSeedFailures(0); // (PREVIOUS IMPLEMENTATION, not wrong just inconsistent acros machine)
		//voxelType seed(sB.findFirstAtOrAfter(rand() % sB.totalSize())), consecutiveSeedFailures(0); // (PREVIOUS IMPLEMENTATION, WRONG FOR WINDOWS.)
		while(seed == sB.totalSize() && consecutiveSeedFailures < maxSeedFailures){
			// Loop continues until a valid seed is found or too many consecutive seed attempts fail
			// runtime is also likely sensitive to this number. 
			consecutiveSeedFailures++;
			uint32_t raw = gen();
			uint64_t value = (uint64_t(raw) * range) >> 32;
			seed = sB.findFirstAtOrAfter(value);
			//seed = sB.findFirstAtOrAfter(dist(gen)); // (PREVIOUS IMPLEMENTATION, not wrong just inconsistent acros machine)
			//seed = sB.findFirstAtOrAfter(rand() % sB.totalSize()); // (PREVIOUS IMPLEMENTATION, WRONG FOR WINDOWS.)
		}
		if(seed == sB.totalSize()){ // if no valid seed found, then this iteration was a failure
			seed = sB.findFirstAtOrAfter(0); // try to find a single seed anywhere
			if(seed == sB.totalSize()){
				numUntriedSeeds = 0;
				consecutiveFailures = maxConsecutiveFailures; // no seed exists anywhere; stop the loop
				continue;
			}
		}
		
		// move seed to central location within vasculature
		double seedMovement(minSeedMovement);
		bool seedFailed(false); // seed fails from disjointCriticalSphere()
		double r(-1.0), rPrev(-1.0);
		vector<voxelType> falsifiedPoints;
		bool centering(true); // flag indicating that the center of the vessel has not yet been found
		while(centering && !seedFailed){
			sB.f(seed); // mark seed as used
			r = -1.0;
			// try placing a sphere around the seed, r is radius of the sphere and falsifiedPoints are the voxels included in this sphere
			seedFailed = !disjointCriticalSphere(B, uB, fB, seed, critFrac, r, falsifiedPoints
#if TRY_WX == 1
				, pm
#endif
				);
			if(seedFailed){
				uB.t(falsifiedPoints); // if failed put the voxels in the sphere back into the unvisited binary volume
#if TRY_WX == 1
				setDefaultColor(falsifiedPoints, pm);
#endif
			}else{
				voxelType newSeed(centerOfMassInSet(B, falsifiedPoints)); // find the center of mass of the voxels in the sphere
				seedMovement = separation3D(seed, newSeed, voxdimsGlobal, B); // EuclideanDistance(B, seed, newSeed);
				centering = seedMovement >= minSeedMovement && r > 1.0*rPrev; // if the seed moved enough and the radius is larger than the previous radius, then we can center the sphere on the new seed
				if(centering){
					uB.t(falsifiedPoints); // put the voxels in the sphere back into the unvisited binary volume
#if TRY_WX == 1
					setDefaultColor(falsifiedPoints, pm);
#endif
					seed = newSeed; // update the seed to the new seed
				}
			}
			rPrev = r; // note that r is updated in disjointCriticalSphere
		} // while loop ends when the sphere is centered or the seed fails to place a sphere.
		if(seedFailed){
			uB.t(falsifiedPoints); // just in case make
#if TRY_WX == 1
			setDefaultColor(falsifiedPoints, pm);
#endif
			consecutiveFailures++;
		}else{
			rv[seed] = r; // add the seed and radius to the map (centered and successfully placed sphere)
			sB.f(falsifiedPoints); // these are not in a sphere so they have to be false in the untried seeds
			consecutiveFailures = 0;
#if TRY_WX == 1
			setPointSetColor(falsifiedPoints, pm, intrasphereColor[0], intrasphereColor[1], intrasphereColor[2], intrasphereColor[3]);
			setPointSetColor(vector<voxelType>(1, seed), pm, rSeed, gSeed, bSeed, aSeed);
			double x(0.0), y(0.0), z(0.0);
			normalizedCoord(seed, x, y, z);
			wxGetApp().myFrame->m_canvas->addNumber(seed, x, y, z);
#endif
		}
	}
	cerr << "\n\t sphereCoarsen() took " << niceTimeSince(sphereTime) << endl;
	return rv;
}

/* Creates a map with a key of intersphere seed and value of a list of the adjacent intrasphere seeds
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that indicates intersphere vasculature
 * @rv rough vectorization of the vasculature (keys being the intrasphere seeds)
 *
 * Explores all of the connected vascular space.  First, an intersphere seed is chosen.  This is expanded to find all bordering intrasphere voxels.  The intrasphere space is then explored in order to find the corresponding intrasphere seeds.  In the end, a map with keys being intersphere seeds and values being a list of adjacent intrasphere seeds is returned.
 *
 * @return a map with keys being intersphere seeds and values being a list of adjacent intrasphere seeds
 */
map<voxelType, vector<voxelType> > adjacentSpheres(const BinaryVolume &B, const BinaryVolume &uB, map<voxelType, double> &rv
#if TRY_WX == 1
		, const map<voxelType, unsigned int> &pm
#endif
		){
	
	map<voxelType, vector<voxelType> > adj;
	BinaryVolume cB(uB), // unvisited during exploration of intersphere space
		sB(B); // unvisited during exploration of intrasphere space
	voxelType seed(cB.findFirstAtOrAfter(0)); // start at a voxel not in a sphere
	while(seed < cB.totalSize()){
		vector<voxelType> F, f(1, seed), centers; // F is for intraspheres, f for unexplored structure
		for(unsigned int fi(0); fi < f.size(); fi++){
			
			vector<voxelType> nh(vesselBlocksInNeighborhood(cB, f[fi]));
			f.insert(f.end(), nh.begin(), nh.end()); // flowing into the intersphere space
			cB.f(nh); // mark as visited
#if TRY_WX == 1
			setPointColor(f[fi], pm, 0.0, 1.0, 0.0, 0.2);
#endif
			nh = vesselBlocksInNeighborhood(B, f[fi]); // get the neighbors of the current voxel that are true in the binary volume
			for(unsigned int i(0); i < nh.size(); i++){
				if(!uB.is(nh[i])){ // if they are in a sphere
					F.push_back(nh[i]); // add to the intrasphere list
					sB.f(nh[i]); // mark as visited in the intrasphere space
					if(rv.find(nh[i]) != rv.end()) // is center
						pushUnique(nh[i], centers);//centers.push_back(nh[i]);
#if TRY_WX == 1
					else // is not center
						setPointColor(nh[i], pm, 1.0, 0.0, 0.0, 0.2);
#endif
				}
			}
		} // explored intersphere space from current seed. 
		
		vector<voxelType> intrasphereVisited(F);
		for(unsigned int Fi(0); Fi < F.size(); Fi++){
			
			if(rv.find(F[Fi]) != rv.end()) // replacing rv.find a few lines below
				pushUnique(F[Fi], centers); // add to the centers list
			
			vector<voxelType> nh(vesselBlocksInNeighborhood(sB, F[Fi])); // not visited in intrasphere (starts as all of B)
			for(unsigned int i(0); i < nh.size(); i++){
				if(!uB.is(nh[i])){ // is intrasphere (not unvisited -> visited)
					F.push_back(nh[i]); // add to the intrasphere list
					sB.f(nh[i]); // mark as visited in the intrasphere space
					intrasphereVisited.push_back(nh[i]);
				}
			}
		} // explored intrasphere space
		
		sB.t(intrasphereVisited); // mark as unvisited again for next iteration which can identify the same sphere adjacent to a different seed from another intersphere space. 
		if(!centers.empty()) // don't want to create empty values for a key in adj
			adj[seed] = centers;

		seed = cB.findFirstAtOrAfter(seed + 1);
	}
	
	return adj;
}

/* Creates a new <BinaryVolume> with the smallest connected components removed.
 * @B <BinaryVolume> that indicates vasculature
 * @keepAtLeastThisManyLargest the (minimum) number of largest connected components to allow
 *
 * Creates a new <BinaryVolume> with the smallest connected components removed.  The new <BinaryVolume> can have more than `keepAtLeastThisManyLargest` connected components if there are multiple connected components that are the same size.
 * @return a <BinaryVolume> with the smallest connected components removed.
 */
BinaryVolume removeSmallestConnectedComponents(const BinaryVolume &B, unsigned int keepAtLeastThisManyLargest
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm
#endif
		){
	
	if(keepAtLeastThisManyLargest == 0)
		keepAtLeastThisManyLargest = 1; // pointless to be empty
	BinaryVolume uB(B); // the unvisited BinaryVolume
	vector<voxelType> seeds;
	vector<unsigned long> cc_sizes;
	voxelType seed(uB.findFirstAtOrAfter(0)); // Finds index of the first bit at or after index `i` that is true
	while(seed < uB.totalSize()){
		uB.f(seed); // set the seed index to false
		vector<voxelType> f(1, seed); // vector of size 1 with the seed index
#if TRY_WX == 1
		setPointSetColor(f, pm, 1.0, 0.0, 0.0, 0.1);
#endif
		for(unsigned int fi(0); fi < f.size(); fi++){ // interate over indices of vector f
			vector<unsigned int long long> nh(vesselBlocksInNeighborhood(uB, f[fi])); // get the neighbors of the current voxel that are true in the binary volume
			uB.f(nh); // set the neighbors to false
			f.insert(f.end(), nh.begin(), nh.end()); //vector.insert(position, first, last): insert the neighbors at the end of the vector f
#if TRY_WX == 1
			setPointSetColor(nh, pm, 1.0, 0.0, 0.0, 0.1);
#endif
		} // this for loop ends when a connected component ends (f grows until there is a disconnectivity, i.e., no more true neighbors.)
		seeds.push_back(seed); // add the seed index to the seeds vector
		cc_sizes.push_back(f.size()); // add the size of the current connected component to the cc_sizes vector
		seed = uB.findFirstAtOrAfter(seed + 1); // find the next seed index that is true
	}
	
	vector<unsigned long> sorted_sizes(cc_sizes);
	sort(sorted_sizes.rbegin(), sorted_sizes.rend()); // sort the sizes in descending order
	unsigned long minimumSize(sorted_sizes.back());
	if(keepAtLeastThisManyLargest < sorted_sizes.size())
		minimumSize = sorted_sizes[keepAtLeastThisManyLargest];
	
	uB = B;
	BinaryVolume lccsB(B); // initialize the binary volume to be returned, it will contain the largest connected components.
	for(unsigned int si(0); si < seeds.size(); si++){
		if(cc_sizes[si] >= minimumSize)
			continue; // continue if the current connected component is larger than the minimum size
		uB.f(seeds[si]);
		lccsB.f(seeds[si]);
		vector<voxelType> f(1, seeds[si]);
#if TRY_WX == 1
		setPointSetColor(f, pm, 0.0, 0.0, 0.5, 0.2); // setDefaultColor(f, pm);
#endif
		for(unsigned int fi(0); fi < f.size(); fi++){
			vector<unsigned int long long> nh(vesselBlocksInNeighborhood(uB, f[fi]));
			uB.f(nh);
			lccsB.f(nh); // sets the neighbors to false 
			f.insert(f.end(), nh.begin(), nh.end());
#if TRY_WX == 1
			setPointSetColor(nh, pm, 0.0, 0.0, 0.5, 0.2); // setDefaultColor(nh, pm);
#endif
		} // this loop ends when all indices of this connected component are set to false. 
	}
	
	return lccsB;
}

/* Finds the frontier of the intersphere region that contains `seed`.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that indicates intersphere vasculature
 * @seed seed voxel for the intersphere region
 *
 * Finds the frontier of the intersphere region that contains `seed`.  The frontier is in some intrasphere region.
 * @return the sets of voxels that are on the frontier of the intersphere region that contains `seed`.
 */
vector<voxelType> visitedFrontier(const BinaryVolume &B, const BinaryVolume &uB, BinaryVolume &cB, voxelType seed){
	if(!uB.is(seed))
		return vector<voxelType>();
	
	vector<voxelType> f(1, seed), F;
	cB.f(f);
	for(unsigned long fi(0); fi < f.size(); fi++){
		vector<voxelType> nh(vesselBlocksInNeighborhood(cB, f[fi]));
		for(unsigned int i(0); i < nh.size(); i++){
			if(uB.is(nh[i]))
				f.push_back(nh[i]);
			else
				F.push_back(nh[i]);
			cB.f(nh[i]);
		}
	}
	return F;
}

/* Finds the intrasphere seed in `rv` that identifies the intrasphere region connected with `seed`.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that indicates intersphere vasculature
 * @rv rough vectorization map with keys being the seeds to intrasphere regions
 * @seed <BinaryVolume> index that the intrasphere region contains
 *
 * Finds the intrasphere seed in `rv` that identifies the intrasphere region connected with `seed` by expanding the region from `seed` until a key in `rv` is found.
 * @return the <BinaryVolume> index of the intrasphere seed in `rv` that identifies the intrasphere region connected with `seed`
 */
voxelType findVisitedCenter(const BinaryVolume &B, const BinaryVolume &uB, const map<voxelType, double> &rv, voxelType seed){
	if(!B.is(seed) || uB.is(seed)) // not vasculature or seed itself is not visited
		return B.totalSize();
	
	vector<voxelType> f(1, seed);
	set<voxelType> v;
	v.insert(seed);
	for(unsigned long i(0); i < f.size(); i++){
		vector<voxelType> nh(vesselBlocksInNeighborhood(B, f[i]));
		for(unsigned int n(0); n < nh.size(); n++){
			if(uB.is(nh[n]))
				continue;
			
			if(rv.find(nh[n]) != rv.end())
				return nh[n];
			
			if(v.insert(nh[n]).second)
				f.push_back(nh[n]);
		}
	}
	
	return B.totalSize();
}

/* Finds the largest false component in <BinaryVolume> `B` with the assumption that the component is at least half of the volume of `B`
 * @B the <BinaryVolume> where the inside is set to true
 *
 * Finds the voxels that are outside of voxels that are inside (true) in <BinaryVolume> `B`, assuming the outside of `B` is at least half of `B`.  Note that the output <BinaryVolume> `O` is no necessarily `B.`<flip>() if `B` has internal holes.
 *
 * @return the <BinaryVolume> `O` that  corresponds to the outisde of <BinaryVolume> `B`
 */
BinaryVolume findOutside(const BinaryVolume &B){
	
	if(B.getSize(0) < 2 || B.getSize(1) < 2 || B.getSize(2) < 2){ // 2D slice: all nonvascular voxels are the outside
		BinaryVolume O(B);
		O.flip();
		return O;
	}
	
	
	BinaryVolume uO(B); // unvisited outside
	uO.flip();
	voxelType seed(uO.findFirstAtOrAfter(0));
	while(seed < uO.totalSize()){
		BinaryVolume O(B.getSize(0), B.getSize(1), B.getSize(2), false);
		vector<voxelType> f(1, seed);
		
		
		for(unsigned int fi(0); fi < f.size(); fi++){
			vector<voxelType> nh(vesselBlocksInNeighborhood(uO, f[fi]));
			f.insert(f.end(), nh.begin(), nh.end());
			uO.f(nh);
			O.t(nh);
		}
		if(f.size() > B.totalSize()/2) // assumption about half of voxels outside of vasculature
			return O;
		
		O.clear();
		seed = uO.findFirstAtOrAfter(seed + 1); // B.findFirstFalseAtOrAfter(seed + 1);
	}
	
	// did not find obvious contiguous outside: return blank
	return BinaryVolume(B.getSize(0), B.getSize(1), B.getSize(2), false);
}

/* Constructs the map with key of intrasphere seed to store the list of intersphere seeds (inverse of `adj`).
 * @adj map that lists the intrasphere seeds that are adjacent to (the map's key) intersphere seeds
 *
 * Constructs a map of the adjacent intersphere seeds from the map `adj` of seeds adjacent to a intersphere seed.
 * @return a map with keys being intrasphere seeds that store a list of intersphere seeds
 */
map<voxelType, vector<voxelType> > sphereConnectionSeeds(map<voxelType, vector<voxelType> > &adj){
	map<voxelType, vector<voxelType> > scs;
	for(map<voxelType, vector<voxelType> >::iterator it(adj.begin()); it != adj.end(); it++){
		for(unsigned int i(0); i < it->second.size(); i++)
			scs[it->second[i]].push_back(it->first); // pushUnique(it->first, scs[it->second[i]]);
	}
	
	return scs;
}

/* Removes false tips from `adj`.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that indicates intersphere vasculature
 * @rv rough vectorization map with keys being the seeds to intrasphere regions
 * @adj map with keys being intersphere seeds and values being the list of adjacent intrasphere seeds
 * @scs map with keys being intrasphere seeds and values being the list of adjacent intersphere seeds
 *
 * Removes false tips from `adj` by considering the connectivity of each intersphere seed.  If there is only a single adjacent intrasphere region, then the intersphere region is likely inappropriate.
 */
void removeFalseTipsByContext(const BinaryVolume &B, BinaryVolume &uB, map<voxelType, double> &rv, map<voxelType, vector<voxelType> > &adj,
		map<voxelType, vector<voxelType> > &scs
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm
#endif
		){
	// absorb singly-connected intersphere seeds into adjacent intrasphere
	vector<voxelType> badAdj;
	for(map<voxelType, vector<voxelType> >::iterator it(adj.begin()); it != adj.end(); it++){
		if(it->second.size() == 1){ // only a single component
			vector<voxelType> f(1, it->first);
			uB.f(f[0]);
			for(unsigned int fi(0); fi < f.size(); fi++){
				vector<voxelType> nh(vesselBlocksInNeighborhood(uB, f[fi]));
				uB.f(nh); // note the invalid intersphere space are marked as visited so we don't use it later. 
				f.insert(f.end(), nh.begin(), nh.end());
			}
			removeFrom(it->first, scs[it->second[0]]);
			badAdj.push_back(it->first);
#if TRY_WX == 1
			setPointSetColor(f, pm, intrasphereColor[0], intrasphereColor[1], intrasphereColor[2], intrasphereColor[3]);
#endif
		}
	}
	
	// remove bad intersphere seeds from adj
	for(unsigned int i(0); i < badAdj.size(); i++)
		adj.erase(badAdj[i]);
}

/* Fills the internal voxels that do not contact the outside of <BinaryVolume> `B` and returns the number of holes filled.
 * @B the <BinaryVolume> to be filled
 * @O the <BinaryVolume> outisde of `B`
 *
 * Fills the internal holes (i.e., the false voxels that do not contact the outside of `B`) in <BinaryVolume> `B` by taking the <flip>()'d outside <BinaryVolume> `O`.
 *
 * @return the number of holes filled
 */
voxelType fillInternalFalses(BinaryVolume &B, const BinaryVolume &O){
	voxelType prevTrue(B.totalTrue());
	B = O;
	B.flip(); // this flips all bits in the binary volume
	return B.totalTrue() - prevTrue;
}

/* Finds the connected components of `frontier`.
 * @frontier the set of voxels of interest
 * @B <BinaryVolume> that contains `frontier`
 *
 * Finds the connected components of `frontier` and returns a vector with the list of voxels in each component.
 * @return a vector with lists of the voxels in each of the connected components of `frontier`
 */
vector<vector<voxelType> > segregateConnectedComponents(const vector<voxelType> &frontier, const BinaryVolume &B){
	vector<voxelType> fRemains(frontier);
	vector<vector<voxelType> > ccs;
	while(fRemains.size() > 0){ // to get all connected components
		ccs.push_back(vector<voxelType>(1, fRemains[0]));
		vector<voxelType> fAdd(1, fRemains[0]);
		fRemains.erase(fRemains.begin());
		while(fAdd.size() > 0 && fRemains.size() > 0){ // to get all points in connected component
			for(int i(0); i < (int)fRemains.size(); i++){
				if(inNeighborhood(fAdd[0], fRemains[i], B)){
					fAdd.push_back(fRemains[i]);
					ccs.back().push_back(fAdd.back());
					fRemains.erase(fRemains.begin() + i);
					i--;
				}
			}
			fAdd.erase(fAdd.begin());
		}
	}
	return ccs;
}

/* Finds the voxel in `f` that is farthest from any nonvascular voxel.
 * @B <BinaryVolume> that indicates vasculature
 * @f set of voxels of interest
 *
 * Finds the voxel in `f` that is farthest from any nonvascualr voxel.
 * @return the voxel in `f` that is farthest from any false point in `B`
 */
voxelType findFarthestInside(const BinaryVolume &B, const vector<voxelType> &f){
	
	// initialize distances
	vector<double> d(f.size(), -1.0);
	for(unsigned int i(0); i < f.size(); i++){
		vector<voxelType> nh(nonvesselBlocksInNeighborhood(B, f[i]));
		if(nh.empty())
			continue;
		d[i] = separation3D(nh[0], f[i], voxdimsGlobal, voldimsGlobal);
		for(unsigned int j(1); j < nh.size(); j++){
			double sep(separation3D(nh[j], f[i], voxdimsGlobal, voldimsGlobal));
			if(d[i] > sep)
				d[i] = sep;
		}
	}

	// adjacencies and/or distances in f
	vector<vector<unsigned int> > adjf(f.size(), vector<unsigned int>());
	vector<vector<double> > dadjf(f.size(), vector<double>());
	for(unsigned int i(0); i < f.size(); i++){
		for(unsigned int j(i + 1); j < f.size(); j++){
			if(inNeighborhood(f[i], f[j], B)){
				double sep(separation3D(f[i], f[j], voxdimsGlobal, voldimsGlobal));
				adjf[i].push_back(j);
				dadjf[i].push_back(sep);
				adjf[j].push_back(i);
				dadjf[j].push_back(sep);
			}
		}
	}
	
	// find shortest distances (somewhat inefficiently)
	bool changed = true;
	while(changed){
		changed = false;
		for(unsigned int i(0); i < f.size(); i++){
			for(unsigned int j(0); j < adjf[i].size(); j++){
				unsigned int k(adjf[i][j]);
				if(d[k] >= 0.0){ // some path to outside is known
					double dist(d[k] + dadjf[i][j]);
					if(d[i] < 0.0 || d[i] > dist){
						d[i] = dist;
						changed = true;
					}
				}
			}
		}
	}

	// choose smallest (favoring center of tied regions)
	double farthest(d[0]);
	for(unsigned int i(0); i < f.size(); i++){
		if(farthest < d[i])
			farthest = d[i];
	}
	vector<voxelType> farthests;
	for(unsigned int i(0); i < f.size(); i++){
		if(d[i] == farthest)
			farthests.push_back(f[i]);
	}

	return centerOfMassInSet(B, farthests);
}

/* Finds the critical vertebrae for each intrasphere region.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that identifies intersphere vasculature
 * @rv rough vectorization map with keys being the seeds to intrasphere regions
 * @adj map with keys being intersphere seeds and values being the list of adjacent intrasphere seeds
 * @scs map with keys being intrasphere seeds and values being the list of adjacent intersphere seeds
 * 
 * Finds the critical vertebrae for each intrasphere region, with the order of vertebrae for each region corresponding to the order of regions in `adj`.
 * @return map of critical vertebrae as values for each intrasphere region seed key
 */
map<voxelType, vector<voxelType> > findCriticalVertebrae(const BinaryVolume &B, const BinaryVolume &uB, map<voxelType, double> &rv, map<voxelType, vector<voxelType> > &adj, map<voxelType, vector<voxelType> > &scs){
	
	map<voxelType, vector<voxelType> > critVert(adj); // organized such that it maps intersphere seeds to intrasphere seeds
	
	BinaryVolume cB(B); // allocate once here so it doesn't need to allocated again and again for visitedFrontier
	
	for(map<voxelType, vector<voxelType> >::iterator it(adj.begin()); it != adj.end(); it++){ // interate over intersphere seeds -> adjacent intrasphere centers
		vector<vector<voxelType> > vfcc(segregateConnectedComponents(visitedFrontier(B, uB, cB, it->first), B)); // visited frontier connected components for the intersphere seed it->first
		// visitedFrontier contains voxels at the interface of intersphere and intrasphere voxels for this seed and vfcc contains the connected components of these voxels as a list. 
		// pretty sure this is important for finding junctions. 
		map<voxelType, vector<voxelType> > vfcc_matched; // vgcc_matched[i] is the list of vfcc elements that belong to intrasphere seed i
		for(unsigned int i(0); i < vfcc.size(); i++){
			voxelType centerLabel(findVisitedCenter(B, uB, rv, vfcc[i][0])); // intrasphere center connected to this interface voxel
			if(centerLabel < B.totalSize())
				vfcc_matched[centerLabel].insert(vfcc_matched[centerLabel].end(), vfcc[i].begin(), vfcc[i].end()); // maps center to the list of voxels at the interface of this center and the intersphere space. 
		}
		for(map<voxelType, vector<voxelType> >::iterator itf(vfcc_matched.begin()); itf != vfcc_matched.end(); itf++){
			unsigned int i(0);
			bool foundIndex(false);
			while(i < it->second.size() && !foundIndex){ // finding index in adj, so looking at it (NOT itf)
				foundIndex = it->second[i] == itf->first; // it->second[i] is the intrasphere center, itf->first is the center in vfcc_matched
				if(!foundIndex)
					i++;
			}
			if(foundIndex)
				critVert[it->first][i] = findFarthestInside(B, itf->second); // find the farthest voxel inside vasculature within itf->second which is the list of voxels at the interface
			else{
				adj[it->first].push_back(itf->first);
				scs[itf->first].push_back(it->first);
				critVert[it->first].push_back(findFarthestInside(B, itf->second));
			}
				
		}
	}
	
	return critVert; // maps intersphere seeds to the voxel that interfaces a intrasphere space, notably the interfaces are disconnected so we can find junctions. 
}

/* Finds the meat of the **intrasphere** connected component defined by the given voxel `seed` in <BinaryVolume> `B`.
 * @B the <BinaryVolume> where meat is indicated by true
 * @uB the <BinaryVolume> that indicates what true voxels in <BinaryVolume> `B` are still **unvisited**
 * @sB the <BinaryVolume> that can be modified for exploration of the meat from `seed` (allocated once outside this method)
 * @seed the index in <BinaryVolume> `B` to seed the connected meat region
 *
 * Finds the meat in <BinaryVolume> `B` of the **visited** connected component that is defined by <BinaryVolume> `uB` and contains the voxel `seed`.
 *
 * @return a list of indices that are the meat defined in `B` that are **visited**
 */
vector<voxelType> intraspherePoints(const BinaryVolume &B, const BinaryVolume &uB, BinaryVolume &sB, voxelType seed){
	if(uB.is(seed)) // unvisited: not a sphere
		return vector<voxelType>();
	
	vector<voxelType> s(1, seed); // , f(1, seed)
	//BinaryVolume sB(B);
	sB.f(seed);
	for(unsigned int si(0); si < s.size(); si++){
		vector<voxelType> nh(vesselBlocksInNeighborhood(sB, s[si]));
		sB.f(nh);
		for(unsigned int i(0); i < nh.size(); i++){
			if(!uB.is(nh[i])) // not unvisited = visited: nh[i] is in the sphere
				s.push_back(nh[i]);
		}
	}
	return s;
}

/* Finds the meat of the **intersphere** connected component defined by the given voxel `seed` in <BinaryVolume> `B`.
 * @B the <BinaryVolume> where meat is indicated by true
 * @uB the <BinaryVolume> that indicates what true voxels in <BinaryVolume> `B` are visited
 * @sB the <BinaryVolume> that can be modified for exploration of the meat from `seed` (allocated once outside this method)
 * @seed the index in <BinaryVolume> `B` to seed the connected meat region
 *
 * Finds the meat in <BinaryVolume> `B` of the **un**visited connected component that is defined by <BinaryVolume> `uB` and contains the voxel `seed`.
 *
 * @return a list of indices that are the meat defined in `B` that are **un**visited
 */
vector<voxelType> interspherePoints(const BinaryVolume &B, const BinaryVolume &uB, BinaryVolume &sB, voxelType seed){
	if(!uB.is(seed)) // visited inside a sphere
		return vector<voxelType>();
	
	vector<voxelType> s(1, seed);
	
	sB.f(seed);
	for(unsigned int si(0); si < s.size(); si++){
		vector<voxelType> nh(vesselBlocksInNeighborhood(sB, s[si]));
		sB.f(nh);
		for(unsigned int i(0); i < nh.size(); i++){
			if(uB.is(nh[i])) // unvisited: is not in a sphere
				s.push_back(nh[i]);
		}
	}
	return s;
}

/* Erodes the given small set of voxels to yield a backbone that connects the vertebrae in `cv`.
 * @B <BinaryVolume> that indicates vasculature
 * @s the set of voxels to erode
 * @cv the set of critical vertebrae that must be included in the final backbone
 *
 * Erodes the given small set of voxels `s` to yield a backbone that connects the critical vertebrae in `cv`.  Voxels nearer to nonvascular voxels are eroded first.
 * @return an eroded set of voxels that make up a backbone that connects the vertebrae in `cv`.
 */
Backbone<> erodeForBackbone(const BinaryVolume &B, vector<voxelType> s, const vector<voxelType> &cv){ // , const BinaryVolume &O
	// add missing critical voxels
	for(unsigned int i(0); i < cv.size(); i++)
		pushUnique(cv[i], s);
	
	// initialize distances that are adjacent to the outside of B
	vector<double> d(s.size(), -1.0);
	for(unsigned int i(0); i < s.size(); i++){ // iterating over voxels to erode
		vector<voxelType> nh(nonvesselBlocksInNeighborhood(B, s[i]));
		if(nh.empty())
			continue;
		
		// search for nearest outisde voxel
		d[i] = separation3D(nh[0], s[i], voxdimsGlobal, voldimsGlobal);
		for(unsigned int j(1); j < nh.size(); j++){
			double sep(separation3D(nh[j], s[i], voxdimsGlobal, voldimsGlobal));
			if(d[i] > sep) // finding minimum distance to an outside voxel
				d[i] = sep;
		}
	}
	
	// build graph with adjacencies and distances in s
	vector<vector<voxelType> > adjf(s.size(), vector<voxelType>()); // adjacent index in s
	vector<vector<double> > dadjf(s.size(), vector<double>()); // adjacent distance corresponding to adjf
	map<voxelType, vector<voxelType> > nhs; // for each voxel in s stores neighbors also in s
	
	// iterate over upper triangular part of s x s matrix
	for(unsigned long i(0); i < s.size(); i++){
		nhs[s[i]];
		map<voxelType, vector<voxelType> >::iterator nhsIt(nhs.find(s[i]));
		for(unsigned long j(i + 1); j < s.size(); j++){ // note i and j are the indices in s/d, and the loop structure ensures that i does not equal j
			if(inNeighborhood(s[i], s[j], B)){
				adjf[i].push_back(j);
				double sep(separation3D(s[i], s[j], voxdimsGlobal, voldimsGlobal));
				dadjf[i].push_back(sep);
				nhsIt->second.push_back(s[j]); // nhs[s[i]]
				adjf[j].push_back(i);
				dadjf[j].push_back(sep);
				nhs[s[j]].push_back(s[i]);
			}
		}
	}
	
	// find shortest distances to outside voxels (somewhat redundantly, but this doesn't take too much time)
	// this is sort of like Bellman-Ford dynamic programming, utilizing the shortest path of the neighboring voxels to compute its own shortest path. 
	// note that at this point d already contains the distances to outside for the voxels lining the boundary of the vasculatue. 
	bool changed = true;
	while(changed){
		changed = false;
		for(unsigned int i(0); i < s.size(); i++){ // index in s/d
			if(d[i] == 0.0) // the shortest path for boundary voxels is trivial
				continue;
			for(unsigned int j(0); j < adjf[i].size(); j++){ // index in adjacency
				unsigned long k(adjf[i][j]); // index in s/d
				if(d[k] >= 0.0){ 
					double dist(d[k] + dadjf[i][j]); // total distance from outside voxel to voxel i if shortest path goes through voxel j
					if(d[i] < 0.0 || d[i] > dist){ // if we found a new minimum distance to an outside voxel change 
						d[i] = dist;
						changed = true;
					}
				}
			}
		}
	}
	// note that these distances will likely come in handy when defining radius of the backbone. 
	
	// organize the voxels into a distance-sorted structure
	map<double, vector<voxelType> > toErode; // maps distance to list of voxels at that distance
	for(unsigned int i(0); i < s.size(); i++)
		toErode[d[i]].push_back(s[i]);
	
	// erode locally
	bool eroded(true);
	while(eroded){
		eroded = false;
		for(map<double, vector<voxelType> >::iterator it(toErode.begin()); it != toErode.end(); it++){ // it->first is the distance, it->second is the list of voxels at that distance
			for(int i(0); i < (int)it->second.size(); i++){ // i is an int so that it can still be decremented when i is zero
				voxelType potRem(it->second[i]); // voxel to potentially remove
				if(isIn(potRem, cv)) // do not remove critical voxels
					continue;
				
				map<voxelType, vector<voxelType> >::iterator nhsIt(nhs.find(potRem)); // nhs[potRem] is the list of neighbors of potRem in s
				bool shouldRemove(nhsIt->second.size() < 2); // if the voxel has less than 2 neighbors in s, it can be removed
				if(!shouldRemove)
					shouldRemove = isSingleConnectedComponent(nhsIt->second, B);
				if(shouldRemove){
					eroded = true;
					for(unsigned int j(0); j < nhsIt->second.size(); j++)
						removeFrom(potRem, nhs[nhsIt->second[j]]); // remove from current neighbors neighborhoods
					// don't care if nhs get updated, so don't bother erasing: nhs.erase(potRem);
					removeFrom(potRem, it->second); // remove potRem from toErode
					i--;
				}
			}
			
		}
	}
	
	// translate into vector and erode globally
	Backbone<> backbone;
	for(map<double, vector<voxelType> >::iterator it(toErode.begin()); it != toErode.end(); it++)
		backbone.add(it->second); // add list of voxels to backbone (this was distance sorted but we add all of them here in a loop)
	
	eroded = true;
	while(eroded){
		eroded = false;
		for(unsigned long i(0); i < backbone.size(); i++){
			voxelType potRem(backbone[i]);
			if(isIn(potRem, cv)) // do not remove critical voxels
				continue;
			
			
			map<unsigned long long, vector<unsigned long long> >::iterator nhsIt(nhs.find(potRem));
			bool shouldRemove(nhsIt->second.size() < 2);
			if(shouldRemove)
				backbone.removeIndex(i);
			else{
				backbone.removeIndex(i);
				shouldRemove = isSingleConnectedComponent(backbone.getVertebrae(), B);
				if(!shouldRemove)
					backbone.addAt(potRem, i); // if removing the voxel would break the connectivity, add it back
			}
			if(shouldRemove){
				eroded = true;
				for(unsigned int j(0); j < nhsIt->second.size(); j++)
					removeFrom(potRem, nhs[nhsIt->second[j]]); // remove from current neighbors neighborhoods
				// don't care if nhs get updated, so don't bother erasing: nhs.erase(potRem);
				i--;
			}
		}
	}
	
	return backbone;
}

/* Erodes each intra- and inter-sphere region to (unorganized) backbones for the subset of regions assigned to this thread.
 * @adj map with keys being intersphere seeds and values being the list of adjacent intrasphere seeds
 * @scs map with keys being intrasphere seeds and values being the list of adjacent intersphere seeds
 * @critVert map of critical vertebrae for each region 
	// critVert maps intersphere seeds to the voxel that interfaces a intrasphere space, notably the interfaces are disconnected so we can find junctions.
 * @backbones list of all backbones
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that determines unvisited (intersphere) vasculature
 * @sBintra <BinaryVolume> for exploring intrasphere segments
 * @sBinter <BinaryVolume> for exploring intersphere segments
 *
 * Erodes each intra- and inter-sphere region to yield (unorganized) backbones for the subset of regions assigned to this thread using <erodeForBackbone>.
 * @return a map with key being the label for the intra-/inter-sphere region and value being the eroded (unorganized) backbone segments for each region
 */
bool erodeIntraIntersphereBackbones(map<voxelType, vector<voxelType> > &adj, map<voxelType, vector<voxelType> > &scs, map<voxelType, vector<voxelType> > &critVert, map<voxelType, Backbone<> > &backbones, const BinaryVolume &B, const BinaryVolume &uB, BinaryVolume &sBintra, BinaryVolume &sBinter, unsigned int th, unsigned int numThreads
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm
#endif
		){

	// intrasphere
	map<voxelType, vector<voxelType> >::iterator itscs(scs.begin());
	for(unsigned int i(0); i < th && itscs != scs.end(); i++) // each thread starts at it's offset (th = 1 works on items 1, numThreads + 1, 2numThreads + 1, etc.)
		itscs++;
	while(itscs != scs.end()){
		vector<voxelType> s(intraspherePoints(B, uB, sBintra, itscs->first)), // intrasphere meat
				cv; // critical voxels
#if TRY_WX == 1
		setPointSetColor(s, pm, intrasphereColor[0], intrasphereColor[1], intrasphereColor[2], 0.03);
#endif
		
		for(unsigned int i(0); i < itscs->second.size(); i++){ // iterate over adjacent intersphere seeds
			voxelType intersphereSeed(itscs->second[i]);
			map<voxelType, vector<voxelType> >::iterator adjit(adj.find(intersphereSeed)); // to a list of adjacent intrasphere seeds
			unsigned long critIndex(adjit->second.size()); // size of the list of adjacent intrasphere seeds
			
			// looking for the index of the current intrasphere seed in the adjacency vector
			for(unsigned int j(0); j < critIndex; j++){
				if(adjit->second[j] == itscs->first)
					critIndex = j; // index of the current intrasphere seed in the adjacency vector
			}
			
			if(critIndex >= adjit->second.size())
				cout << "\n Warning erodeIntraIntersphereBackbones(): could not find " << itscs->first << " in adj[intersphereSeed] = " << makeString(adjit->second) << endl;
			else
				cv.push_back(critVert[intersphereSeed][critIndex]); // add the critical voxel that interfaces the intersphere and intrasphere spaces of interest
				// roughly imagine starting at an intersphere seed and finding all paths to a neighboring intrasphere seed
				// critical vertebrae values lie at the interface of the intra and inter sphere spaces
				// cv contains all of the critical voxels that must remain after erosion to ensure that the connectivity map does not change. 
		}
		
		if(cv.empty())
			cout << "\n Warning erodeIntraIntersphereBackbones(): intrasphere cv for scs[itrv->first=" << makeString(itscs->first) << "] = " << makeString(itscs->second) << " is empty!" << endl;
		else{
			backbones[itscs->first] = erodeForBackbone(B, s, cv);
#if TRY_WX == 1
			setPointSetColor(backbones[itscs->first].getVertebrae(), pm, intrasphereColor[0], intrasphereColor[1], intrasphereColor[2], 0.5);
#endif
		}

		for(unsigned int i(0); i < numThreads && itscs != scs.end(); i++)
			itscs++;
	}
	
	// intersphere
	map<voxelType, vector<voxelType> >::iterator itcrit(critVert.begin());
	for(unsigned int i(0); i < th && itcrit != critVert.end(); i++)
		itcrit++;
	while(itcrit != critVert.end()){
		if(itcrit->second.empty())
			cout << "\n Warning erodeIntraIntersphereBackbones(): intersphere cv for intersphere is empty!" << endl;
		else{
			vector<voxelType> s(interspherePoints(B, uB, sBinter, itcrit->first));
#if TRY_WX == 1
			setPointSetColor(s, pm, intersphereColor[0], intersphereColor[1], intersphereColor[2], 0.03);
#endif
			backbones[itcrit->first] = erodeForBackbone(B, s, itcrit->second);
#if TRY_WX == 1
			setPointSetColor(backbones[itcrit->first].getVertebrae(), pm, intersphereColor[0], intersphereColor[1], intersphereColor[2], 0.5);
#endif
		}

		for(unsigned int i(0); i < numThreads && itcrit != critVert.end(); i++)
			itcrit++;
	}
	
	return true;
}

/* Erodes each intra- and inter-sphere region to (unorganized) backbones for each region.
 * @B <BinaryVolume> that indicates vasculature
 * @uB <BinaryVolume> that determines unvisited (intersphere) vasculature
 * @adj map with keys being intersphere seeds and values being the list of adjacent intrasphere seeds
 * @scs map with keys being intrasphere seeds and values being the list of adjacent intersphere seeds
 * @critVert map of critical vertebrae for each region
 *
 * Erodes each intra- and inter-sphere region to yield (unorganized) backbones for each region using <erodeIntraIntersphereBackbones>.  Any empty backbones are removed.
 * @return a map with key being the label for the intra-/inter-sphere region and value being the eroded (unorganized) backbone segments for each region
 */
map<voxelType, Backbone<> > backboneErosion(const BinaryVolume &B, const BinaryVolume &uB, // const BinaryVolume &O,
		map<voxelType, vector<voxelType> > &adj,
		map<voxelType, vector<voxelType> > &scs,
		map<voxelType, vector<voxelType> > &critVert
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm
#endif
		){
	
	map<voxelType, Backbone<> > backbones;
	BinaryVolume sB(B), sBa(B); // allocate once here instead of to each call of intraspherePoints and interspherePoints; intrasphere and intersphere regions are disjoint
	for(map<voxelType, vector<voxelType> >::iterator it(critVert.begin()); it != critVert.end(); it++) // it->first is start node it->second is the end node of the critical vertebra
		backbones[it->first]; // initialize the backbones corresponding to each critical vertebra

	std::cout << "[backboneErosion] launching " << numLocThreads << " thread(s)." << std::endl;
	vector<future<bool> > locThreads;
	for(unsigned int th(0); th < numLocThreads; th++)
		locThreads.push_back(async(launch::async, erodeIntraIntersphereBackbones, ref(adj), ref(scs), ref(critVert), ref(backbones), ref(B), ref(uB), ref(sB), ref(sBa), th, numLocThreads
#if TRY_WX == 1
		, ref(pm)
#endif
		)); // launches parallel processes of erodeIntraIntersphereBackbones
	
	for(unsigned int th(0); th < locThreads.size(); th++)
		locThreads[th].get();

	// at this point the backbones contains a backbone segment for each intrasphere and intersphere region. 
	
	// erase empty backbones
	vector<voxelType> emptyBackbones;
	for(map<voxelType, Backbone<> >::iterator it(backbones.begin()); it != backbones.end(); it++){
		if(backbones[it->first].empty())
			emptyBackbones.push_back(it->first);
	}
	for(unsigned long i(0); i < emptyBackbones.size(); i++)
		backbones.erase(emptyBackbones[i]);
	
	return backbones;
}

/* Segments, organizes, and produces connectivity information from raw sets of voxels that make up backbones.
 * @B <BinaryVolume> that indicates vasculature
 * @rawBackbones the list of raw (unorganized and not segmented) backbones
 * @critVert list of critical vertebrae in each <Backbone> that should not be eroded
 * @branchpoints sets of adjacent backbones after segmentation (i.e., the  connectivity information)
 *
 * Segments, organizes, and produces connectivity information from raw sets of voxels that make up backbones.  Segmentation involves identifying branching points and the connections between branching points.
 * @return the list of organized, segmented backbones that correspond with the other output `branchpoints` (indices in `backbones` that are adjacent).
 */
vector<Backbone<> > segmentBackbones(const BinaryVolume &B, map<voxelType, Backbone<> > &rawBackbones,
		map<voxelType, vector<voxelType> > &critVert, vector<vector<branchType> > &branchpoints){
	
	if(rawBackbones.size() == 0)
		return vector<Backbone<> >();
	
	if(rawBackbones.size() == 1){
		Backbone<> onlyBackbone(rawBackbones.begin()->second);
		onlyBackbone.organize(B); // orderSegmentedBackbone(B, onlyBackbone);
		return vector<Backbone<> >(1, onlyBackbone.getVertebrae()); // organized backbone by assuming a linear chain and walking from tip to tip. 
	}

	/* FROM KAI, OVERKILL COMMENTING: but this part is kind of tricky so it will be worth the space to write this out.
	
	BIG PICTURE: For each seed (key), you are reorganizing one raw backbone into multiple linear segments and recording how those segments connect at branchpoints.

	understanding the data structure below is key to understanding the code. In the future, all code with unconventional data structures should be described in detail in comments or documentation. 

	for each seed, we build two things: 
		individualOrganizedBackbones[seed] of type map<voxelType, vector<Backbone<> > >
			seed → vector of Backbone segments (linear chains of voxels)
			indexing occurs like this: individualOrganizedBackbones[seed][branchID]
		individualBranchPoints[seed] of type map<voxelType, vector<vector<branchType>>>
			seed → list of branchpoints
			Each inner vector represents one branchpoint event
			individualBranchPoints[seed][i] contains segment or branch IDs that MEET at the branchpoint. (IMPORTANT, this was conceptually difficult for me to understand when I was reading the code without comments.)
		
	when the algorithm reaches an inline branchpoint (the frontier voxel we are exploring also happens to be a critical vertebra voxel, I think this can happen as a function of the stochastic nature of sphere formation)
		for individualOrganizedBackbones[seed]:
			at an inline branchpoint, the frontier voxel that happens to also be a critical voxel is contained in its own backbone.
			in a regular branchpoint the neighbors become the starting tip voxel for the next branch
		for individualBranchPoints[seed]:
			at an inline branchpoint, [ oldSegmentID, f0SegmentID, neighborSegmentID ]
			at a regular branchpoint, any number of segments or branchIDs that meet at the branchpoint
	*/
	
	// organize each seed individually
	map<voxelType, vector<vector<branchType> > > individualBranchPoints;
	map<voxelType, vector<Backbone<> > > individualOrganizedBackbones;
	map<voxelType, vector<voxelType> > tipSeeds; // keeps track of tips so that they can easily be found for stitching
	for(map<voxelType, Backbone<> >::iterator it(rawBackbones.begin()); it != rawBackbones.end(); it++){ // each unsorted backbone segment should be strictly arborescent (i.e., no reticulations)
		// rawBackbones is a map of (intersphere or intrasphere) seed to the corresponding backbone segment of that region. 
		Backbone<> unorganizedBackbone(it->second); // get the backbone segment for this seed
		
		unsigned long tipIndex(0); // need to find a tip to start; could find from critVert, but critVert is not conveniently organized for spheres
		vector<voxelType> nh(findNeighbors(B, unorganizedBackbone[tipIndex], unorganizedBackbone.getVertebrae()));
		while(nh.size() > 2 && tipIndex < unorganizedBackbone.size()){
			tipIndex++;
			nh = findNeighbors(B, unorganizedBackbone[tipIndex], unorganizedBackbone.getVertebrae());
		}
		if(tipIndex >= unorganizedBackbone.size()){
			cout << "\n Warning segmentBackbones(): could not find a tip for seed " << it->first << " in backbone " << makeString(unorganizedBackbone.getVertebrae()) << endl;
			continue;
		}
		
		// explore segment backbone starting from a tip (unorganizedBackbone[tipIndex])
		vector<voxelType> f(1, unorganizedBackbone[tipIndex]), bi(1, 0); // f is the frontier, bi is the segment index
		individualOrganizedBackbones[it->first]; // initialize the organized backbone for this seed
		map<voxelType, vector<Backbone<> > >::iterator indOrgBackboneIt(individualOrganizedBackbones.find(it->first)); // get the iterator for (seed, backbone) pair for this seed
		indOrgBackboneIt->second.push_back(Backbone<>(f[0])); // append backbone with just the tip to the vector of backbones for this seed
		tipSeeds[f[0]].push_back(it->first); // link tip to the corresponding seed
		unorganizedBackbone.removeIndex(tipIndex); // mark tip as visited by removing from unorganizedBackbone which is the backbone segment for this seed
		individualBranchPoints[it->first]; // initialize the branchpoints for this seed
		map<voxelType, vector<vector<branchType> > >::iterator indBPIt(individualBranchPoints.find(it->first)); // iterator for (seed, branchpoints) pair for this seed
		while(!f.empty()){
			nh = findNeighbors(B, f[0], unorganizedBackbone.getVertebrae()); // explore from tip along this backbone, only among unvisited skeleton voxels
			
			if(nh.size() > 1){ // branchpoint
				indBPIt->second.push_back(vector<branchType>(1, bi[0])); // add to the branchpoints in the (seed, (segments or branches, branchpoints) ) for this seed, bi[0] is the label for the current branchpoint
				f.erase(f.begin());
				bi.erase(bi.begin());
				for(unsigned int i(0); i < nh.size(); i++){
					f.push_back(nh[i]); // add neighbor to frontier
					bi.push_back((unsigned int)indOrgBackboneIt->second.size()); // add a new segment index because each neighbor at a branchpoint will become a new segment, named based on the number of backbones in indOrgBackboneIt
					indBPIt->second.back().push_back(bi.back()); // add new label to the current branch's branchpoint record
					indOrgBackboneIt->second.push_back(Backbone<>(nh[i])); // vector<voxelType>(1, nh[i]) // create and add new backbone segment starting at that neighbor.
					unorganizedBackbone.remove(nh[i]); // removeFrom(nh[i], unorganizedBackbone); // mark neighbor as visited
				} // indBPIt->second.back() now contains the segment or branch IDs that MEET at the branchpoint (if the branchpoint has 3 neighbors, it will have 3 entries)
			}else if(nh.size() == 1){ // skeleton segment might continue
				if(f[0] != indOrgBackboneIt->second[0][0] && isIn(f[0], critVert)){ // inline branchpoint
					indBPIt->second.push_back(vector<branchType>(1, bi[0])); // add new (segment, branchpoint) for this seed, containing the current branchpoint index 
					indOrgBackboneIt->second[bi[0]].remove(f[0]); // removeFrom(f[0], individualOrganizedBackbones[it->first][bi[0]]); // remove f[0] from previous branch because it becomes a new branch here
					// remember that indOrgBackboneIt->second has the backbone for each branch 
					indBPIt->second.back().push_back((unsigned int)indOrgBackboneIt->second.size()); // add new branch index to the recently added (segment, branchpoint) 
					indOrgBackboneIt->second.push_back(Backbone<>(f[0])); // vector<voxelType>(1, f[0]) // create new branch with f[0] as the tip
					bi[0] = (unsigned int)indOrgBackboneIt->second.size(); // neighbor will be in new segment (note that (unsigned int)indOrgBackboneIt->second.size() just incremented by 1)
					indBPIt->second.back().push_back(bi[0]); // add neighbor's segment to branch point
					indOrgBackboneIt->second.push_back(Backbone<>(nh[0])); // vector<voxelType>(1, nh[0])
					pushUnique(it->first, tipSeeds[f[0]]);
				}else // skeleton segment continues
					indOrgBackboneIt->second[bi[0]].add(nh[0]); // individualOrganizedBackbones[it->first][bi[0]].push_back(nh[0]); // we just keep traversing the neighbors here
				f[0] = nh[0]; // now next frontier is the neighbor
				unorganizedBackbone.remove(nh[0]);//removeFrom(nh[0], unorganizedBackbone); // mark neighbor as visited
			}else{ // tip
				pushUnique(it->first, tipSeeds[f[0]]);//tipSeeds[f[0]].push_back(it->first);
				f.erase(f.begin());
				bi.erase(bi.begin());
			}
		}
	}
	
	// create vector of backbones into one big vector of individual backbones
	vector<Backbone<> > backbones;
	branchpoints.clear(); // I think should be empty at this point already. 
	map<unsigned long, unsigned long> offsets; // the number of backbones already entered ahead of this segment: for seed x, individualOrganizedBackbones[x][i] is backbones[offsets[x] + i]
	for(map<voxelType, vector<Backbone<> > >::iterator it(individualOrganizedBackbones.begin()); it != individualOrganizedBackbones.end(); it++){ // for each seed, there are multiple backbone segments
		offsets[it->first] = (unsigned long)backbones.size();
		backbones.insert(backbones.end(), it->second.begin(), it->second.end()); // insert all the backbone segments
		for(unsigned int i(0); i < individualBranchPoints[it->first].size(); i++){ // for each of the brachpoints corresponding to this seed
			branchpoints.push_back(vector<branchType>());
			for(unsigned int j(0); j < individualBranchPoints[it->first][i].size(); j++) // for each of the branches or segment IDs that meet at this branchpoint
				branchpoints.back().push_back(offsets[it->first] + individualBranchPoints[it->first][i][j]); // offset + segment ID = global index of the branch
		}
	}
	
	// stitch tips that have at least two segments by grouping into a branchpoint
	for(map<voxelType, vector<voxelType> >::iterator it(tipSeeds.begin()); it != tipSeeds.end(); it++){ // iterating over tips (tip, seeds that have this tip)
		if(it->second.size() < 2) // actual tip or stranded component
			continue;
		
		branchpoints.push_back(vector<branchType>()); // init a branchpoint vector
		for(unsigned int i(0); i < it->second.size(); i++){
			voxelType seed(it->second[i]); // for each of the seeds that have this tip
			map<voxelType, vector<Backbone<> > >::iterator indOrgBackboneIt(individualOrganizedBackbones.find(seed)); // backbones corresponding to this seed
			unsigned long tipIndividualBackboneIndex((unsigned long)indOrgBackboneIt->second.size()); // number of backbones for this seed
			for(unsigned int j(0); j < indOrgBackboneIt->second.size(); j++){ // over each of the individual seed's backbones
				unsigned int kmax((unsigned int)indOrgBackboneIt->second[j].size()); // number of voxels in this backbone
				for(unsigned int k(0); k < kmax; k++){ // over each voxel in the backbone
					if(indOrgBackboneIt->second[j][k] == it->first){ // if this voxel is the tip
						tipIndividualBackboneIndex = j;
						k = kmax; // skip the rest of the loop
						j = (unsigned int)indOrgBackboneIt->second.size(); // skip the rest of the loop
					}
				}
			}
			
			if(tipIndividualBackboneIndex >= indOrgBackboneIt->second.size())
				cout << "\n Warning segmentBackbones(): tip point " << it->first << " not found in backbones from seed " << seed << endl;
			else
				branchpoints.back().push_back(offsets[seed] + tipIndividualBackboneIndex); // add this to the branchpoint vector so that the segments that share this tip are connected via this branchpoint
		} // after this for loop, we are done looking for the branches connected at this tip. 
	}
	
	// combine backbones in branch points with only two backbone segments for final version
	// these backbones should always share the critical vertex that defined them
	for(long i(0); i < (long)branchpoints.size(); i++){ // for each element in branchpoints (which is a vector of the segments that meet at a particular branchpoint)
		if(branchpoints[i].size() == 2){ // really these only come from tip stitching when exactly two segments meet at a shared tip
			branchType h(branchpoints[i][0]), t(branchpoints[i][1]);
			Backbone<> head, tail;
			// each of these statements make sure that end of head equals start of tail
			if(backbones[h][0] == backbones[t][0]){ 
				head = backbones[h].reversed();
				tail = backbones[t];
			}else if(backbones[h].back() == backbones[t][0]){ 
				head = backbones[h];
				tail = backbones[t];
			}else if(backbones[h][0] == backbones[t].back()){
				head = backbones[t];
				tail = backbones[h];
			}else if(backbones[h].back() == backbones[t].back()){
				head = backbones[h];
				tail = backbones[t].reversed();
			}else{
				cout << "\n Warning segmentBackbones(): backbones at branchpoint " << i
					<< " share no critical point" << "\n\t backbones are:"
					<< "\n\t\t" << makeString(backbones[h].getVertebrae())
					<< "\n\t\t" << makeString(backbones[t].getVertebrae())
					<< endl;
				continue;
			}
			
			// perform stitching
			if(tail.size() == 1){
				// keep head unchanged and simply clear tail without adding to head, this single tail voxel is already at the end of head
			}else{
				bool overlapping(true);
				while(overlapping && tail.size() > 1){
					tail.removeIndex(0);
					if(tail.size() == 1) // if tail happens to be length 2, tail[1] indexing will cause problems. we break here to avoid this. 
						break;
					vector<voxelType> nh(findNeighbors(B, tail[1], head.getVertebrae()));
					overlapping = nh.size() > 0;
				}
				// there is one highly unlikely edge case where actually the whole tail is overlapping with head, which could mean that this remaining voxel is already in head
				if(tail.size() > 0 && isIn(tail[0], head.getVertebrae()))
					tail.removeIndex(0);
				if(tail.size() > 0)
					head.add(tail);
			}
			backbones[h] = head;
			backbones[t].clear();

			// rewrite instances of t in branchpoints to h
			for(unsigned int j(0); j < branchpoints.size(); j++){
				for(unsigned int k(0); k < branchpoints[j].size(); k++){
					if(branchpoints[j][k] == t)
						branchpoints[j][k] = h;
				}
			}

			branchpoints.erase(branchpoints.begin() + i);
			i--;
		}
	}

	
	// condense backbones so that it is compact (no empty backbones) and adjust branchPoints accordingly
	for(int i(0); i < (int)backbones.size(); i++){
		if(backbones[i].size() == 0){
			backbones.erase(backbones.begin() + i);
			for(unsigned int j(0); j < branchpoints.size(); j++){
				for(unsigned int k(0); k < branchpoints[j].size(); k++){
					if((int)branchpoints[j][k] > i)
						branchpoints[j][k]--;
				}
			}
			i--;
		}
	}

	return backbones;
}

/* Unsnarls <Backbone> `backbone` locally by removing unnecessary single-voxel diversions
 * @backbone <Backbone> to be unsnarled
 * @B <BinaryVolume> that indicates vasculature
 * @voxdims the extent of a single voxel
 * @voldims the number of voxels in each dimension
 *
 * Unsnarls (i.e., straightens out) each <Backbone> (in place) by removing unnecessary single-voxel detours (if possible).  Ends are held constant.
 */
void localUnsnarl(Backbone<> &backbone, const BinaryVolume &B, const double (&voxdims)[3], const voxelType (&voldims)[3]){
	bool changed(false);
	for(unsigned int w(2); w < backbone.size(); w++){
		unsigned int v(w - 1), u(w - 2);
		vector<voxelType> nh(vesselBlocksInNeighborhood(B, backbone[v]));
		double lenBest(separation3D(backbone[u], backbone[v], voxdims, B) + separation3D(backbone[v], backbone[w], voxdims, B));
		long nBest(-1);
		for(unsigned int n(0); n < nh.size(); n++){
			if(!inNeighborhood(nh[n], backbone[u], B))
				continue;
			if(!inNeighborhood(nh[n], backbone[w], B))
				continue;
			
			double len(separation3D(backbone[u], nh[n], voxdims, B) + separation3D(nh[n], backbone[w], voxdims, B));
			if(lenBest < 0 || lenBest > len){ // if the neighbor minimizes the length of the backbone (striaghten out backbone locally)
				lenBest = len;
				nBest = n;
			}
		}
		if(nBest >= 0){
			backbone[v] = nh[nBest];
			changed = true;
		}
	}
	
	if(changed)
		backbone.newLength(voxdims, voldims); // recompute length with straightened out backbone
}

/* Unsnarls backbones locally by calling <localUnsnarl>
 * @backbones list (possibly of list(s)) of <Backbone> objects
 * @B <BinaryVolume> that indicates vasculature
 * @voxdims the extent of a single voxel
 * @voldims the number of voxels in each dimension
 *
 * Calls <localUnsnarl> for each <Backbone> in `backbones`.
 */
template <class T>
void localUnsnarl(vector<T> &backbones, const BinaryVolume &B, const double (&voxdims)[3], const voxelType (&voldims)[3]){
	for(unsigned int i(0); i < backbones.size(); i++)
		localUnsnarl(backbones[i], B, voxdims, voldims);
}

/* Determines the fraction of the voxel `v` that belongs to the <Backbone> identified by `lab`.
 * @B <BinaryVolume> that indicated vasculature
 * @vertebraMap map of all voxels (keys) in any <Backbone> with a mapping to the label of the <Backbone> (values)
 * @v index in `B` of interest
 * @lab label of the <Backbone> of interest
 *
 * Searches for the nearest <Backbone> to voxel `v` and evaluates how much (if any) of `v` belongs to the <Backbone> identified by `lab`.
 * @return the fraction of the voxel `v` that belongs to the <Backbone> identified by `lab`
 */
template <class L>
double belongingnessTo(const BinaryVolume &B, const map<voxelType, L > &vertebraMap, voxelType v, const L &lab){
	typename map<voxelType, L>::const_iterator it(vertebraMap.find(v));
	if(it != vertebraMap.end()){
		if(it->second == lab)
			return 1;
		return 0;
	}
	
	map<voxelType, double> exploredBulk;
	//frontier propagates by distance
	map<double, vector<voxelType> > frontier;
	frontier[0].push_back(v);
	while(frontier.size() > 0){
		for(unsigned long f(0); f < frontier.begin()->second.size(); f++){
			voxelType vox(frontier.begin()->second[f]);
			if(exploredBulk.find(vox) == exploredBulk.end()){ // then this is the shortest path; otherwise there has already been another (as short or shorter) path
				it = vertebraMap.find(vox);
				if(it == vertebraMap.end()){ // not a vertebra
					vector<voxelType> nh(vesselBlocksInNeighborhood(B, vox));
					for(unsigned long n(0); n < nh.size(); n++){
						if(exploredBulk.find(nh[n]) == exploredBulk.end()) // not in exploredBulk; otherwise, don't add nh[n] to anything
							frontier[frontier.begin()->first + separation3D(vox, nh[n], voxdimsGlobal, B)].push_back(nh[n]); // duplicates will be handled later
					}
				}else{ // the first vertebra is found
					unsigned long labCount((it->second == lab) ? 1 : 0), anyVertCount(1); 
					for(unsigned long g(f + 1); g < frontier.begin()->second.size(); g++){
						it = vertebraMap.find(frontier.begin()->second[g]);
						if(it != vertebraMap.end()){
							anyVertCount++;
							if(it->second == lab)
								labCount++;
						}
					}
					return double(labCount)/anyVertCount;
				}
				
				exploredBulk[vox] = frontier.begin()->first;
			}
		}
		frontier.erase(frontier.begin());
	}
	/*// frontier propagates by adjacency
	vector<voxelType> frontier(1, v);
	for(unsigned long f(0); f < frontier.size(); f++){
		voxelType vox(frontier[f]);
		if(exploredBulk.find(vox) == exploredBulk.end()){ // no duplicates
			it = vertebraMap.find(vox);
			if(it == vertebraMap.end()){ // not a vertebra
				vector<voxelType> nh(vesselBlocksInNeighborhood(B, vox));
				for(unsigned long n(0); n < nh.size(); n++){
					if(exploredBulk.find(nh[n]) == exploredBulk.end()) // not in exploredBulk; otherwise, don't add nh[n] to anything
						frontier.push_back(nh[n]); // duplicates will be handled later
				}
			}else{ // the first vertebra is found
				unsigned long labCount((it->second == lab) ? 1 : 0), anyVertCount(1); // all or nothing (for now...)
				return double(labCount)/anyVertCount;
			}
		
			exploredBulk[vox]; // just needs to exist
		}
	}*/
	
	return 0; // shouldn't ever be reached, but it'll be here just in case
}

/* Updates distance map `dist` with first list in `frontier`.
 * @B <BinaryVolume> that indicates vasculature
 * @frontier voxels in the frontier, sorted by distance: key is distance from source, value is (unsorted, not necessarily unique) list of voxels
 * @dist list that gets updated of shortest distance to points: key is voxel, value is shortest distance found to voxel so far
 * @expandDist list that does not get updated of shortest distance to points from a possibly different source: key is voxel, value is shortest distance found to voxel so far
 * @tol tolerance for qualification of sources as equidistant
 *
 * Updates distance map `dist` with first list in `frontier`, taking a step in a simple implementation of Dijkstra's algorithm.  Returns the list of newly identified voxels with shortest distances (i.e., the updates to `dist` from `frontier`).
 * @return list of newly-established voxels with shortest distances
 */
vector<voxelType> expandMeatFrontier(const BinaryVolume &B, map<double, vector<voxelType> > &frontier, map<voxelType, double> &dist, map<voxelType, double> &expandDist, double tol){
	vector<voxelType> acceptedFrontier;
	if(frontier.empty())
		return acceptedFrontier;
	
	for(unsigned int f(0); f < frontier.begin()->second.size(); f++){
		voxelType vox(frontier.begin()->second[f]);
		map<voxelType, double>::iterator it(dist.find(vox));
		bool shouldAccept(it == dist.end());
		if(!shouldAccept)
			shouldAccept = frontier.begin()->first < it->second;
		if(shouldAccept && &dist != &expandDist){ // also check expandDist
			it = expandDist.find(vox);
			if(it != expandDist.end())
				shouldAccept = frontier.begin()->first < it->second + tol;
		}
		if(shouldAccept){
			dist[vox] = frontier.begin()->first;
			acceptedFrontier.push_back(vox);
			vector<voxelType> nh(vesselBlocksInNeighborhood(B, vox));
			for(unsigned long n(0); n < nh.size(); n++){
				double nDist(frontier.begin()->first + separation3D(vox, nh[n], voxdimsGlobal, B));
				it = dist.find(nh[n]);
				bool shouldAdd(it == dist.end());
				if(!shouldAdd)
					shouldAdd = it->second > nDist;
				if(shouldAdd && &dist != &expandDist){ // also check expandDist
					it = expandDist.find(nh[n]);
					if(it != expandDist.end())
						shouldAdd = it->second + tol > nDist;
				}
				if(shouldAdd)
					frontier[nDist].push_back(nh[n]);
			}
		}
	}
	frontier.erase(frontier.begin());
	return acceptedFrontier;
}

/*Finds the meat associated with <Backbone> `backbones[b]`.
 * @B <BinaryVolume> that indicates vasculature
 * @vertebraMap map of all voxels (keys) in any <Backbone> with a mapping to the label of the <Backbone> (values)
 * @backbones list of all <Backbone> objects
 * @b index in `backbones` of the <Backbone> of interest
 * @tol tolerance for qualification of sources as equidistant
 *
 * Finds the meat associated with <Backbone> `backbones[b]` by starting with the vertebra(e) of `backbones[b]` and growing the list of voxels that are nearest it.  When another <Backbone> vertebra is encountered, its influence is added to the competing expansion.  Once the frontier for `backbones[b]` vanishes, the meat is grown to double the farthest association in order to verify that there are not other vertebrae to complete with `backbones[b]`.
 * @return a map with the fraction of voxels that are associated with <Backbone> `backbone`
 */
map<voxelType, double> segmentMeatFromAdjacencies(const BinaryVolume &B, const map<voxelType, unsigned long> &vertebraMap, const vector<Backbone<> > &backbones, unsigned long b, double tol = 1.0e-3){
	// initialize frontiers
	map<voxelType, double> dist, otherDist;
	map<double, vector<voxelType> > frontier, otherFrontier;
	// note that key for frontier is the distance from the source vertebra. 
	for(unsigned long v(0); v < backbones[b].size(); v++) // initialize the frontier with the vertebrae of the backbone of interest
		frontier[0].push_back(backbones[b][v]);
	
	// expand frontiers
	// think of this like a multisource Dijkstra's algorithm. 
	double lastFrontierDist(0);
	while(!frontier.empty()){
		bool expandOtherFrontier(otherFrontier.empty() ? false : otherFrontier.begin()->first <= frontier.begin()->first); // expand whichever frontier has the closest voxel from the source vertebra
		while(expandOtherFrontier){
			expandMeatFrontier(B, otherFrontier, otherDist, dist, tol);
			expandOtherFrontier = otherFrontier.empty() ? false : otherFrontier.begin()->first <= frontier.begin()->first;
		}
		
		lastFrontierDist = frontier.begin()->first;
		vector<voxelType> expansion(expandMeatFrontier(B, frontier, dist, otherDist, tol));
		for(unsigned long i(0); i < expansion.size(); i++){
			map<voxelType, unsigned long>::const_iterator it(vertebraMap.find(expansion[i]));
			if(it != vertebraMap.end()){
				if(it->second != b) // if expansion touches a vertebra from another backbone, this triggers competition. 
					otherFrontier[0].push_back(expansion[i]); // add to other frontier at distance 0 because it it at the vertebra of a different backbone.
			}
		}
	}
	
	// What if there is another backbone nearby, but separated by a slightly longer path, so it was never directly touched? 
	// Then we would incorrectly assign its territory to backbone b. 
	// This block of code addresses this problem. 
	map<voxelType, double> checkDist, emptyDist;
	map<double, vector<voxelType> > checkFrontier;
	for(map<voxelType, double>::iterator it(dist.begin()); it != dist.end(); it++)
		checkFrontier[it->second].push_back(it->first); // rebuilds a frontier from all voxels claimed by b
	bool expandCheckFrontier(checkFrontier.empty() ? false : checkFrontier.begin()->first <= 2*lastFrontierDist); // expand out to 2 times the maximum distance from the source vertebra. 
	while(expandCheckFrontier){
		expandMeatFrontier(B, checkFrontier, checkDist, emptyDist, tol);
		expandCheckFrontier = checkFrontier.empty() ? false : checkFrontier.begin()->first <= 2*lastFrontierDist;
	}
	for(map<voxelType, double>::iterator it(checkDist.begin()); it != checkDist.end(); it++){ // now look for other vertebrae that are within 2 times the maximum distance from the source vertebra. 
		map<voxelType, unsigned long>::const_iterator vIt(vertebraMap.find(it->first));
		if(vIt != vertebraMap.end()) if(vIt->second != b)
			otherFrontier[0].push_back(it->first);
	}
	bool expandOtherFrontier(otherFrontier.empty() ? false : otherFrontier.begin()->first <= lastFrontierDist); // expand the newly added competing vertebrae
	while(expandOtherFrontier){
		expandMeatFrontier(B, otherFrontier, otherDist, dist, tol);
		expandOtherFrontier = otherFrontier.empty() ? false : otherFrontier.begin()->first <= lastFrontierDist;
	}
	
	map<voxelType, double> sMap; // final ownership decisions 
	for(map<voxelType, double>::iterator it(dist.begin()); it != dist.end(); it++){
		map<voxelType, double>::iterator oIt(otherDist.find(it->first));
		if(oIt == otherDist.end())
			sMap[it->first] = 1; // no other vertebrae reached this voxel, fully belongs to b
		else if(oIt->second + tol >= it->second)
			sMap[it->first] = belongingnessTo(B, vertebraMap, it->first, b); // fractional ownership 
	}
	
	return sMap;
}

/* Finds the <BinaryVolume> index contained in `s` that is farthest from any other segment in `B`.
 * @s the meat of the segment of interest
 * @B <BinaryVolume> that indicates vasculature
 *
 * Finds the <BinaryVolume> index that is contained in `s` and is farthest from any other segment in `B`.  This is used to find the effective tip for the segment meat defined by `s`.
 * @return the <BinaryVolume> index that is contained in `s` that is farthest from any other segment in `B`.
 */
voxelType farthestFromOtherSegments(vector<voxelType> &s, const BinaryVolume &B){

	// adjacencies and/or distances in s
	vector<vector<voxelType> > adjf(s.size(), vector<voxelType>());
	map<voxelType, vector<voxelType> > nhs; // stored indices in B of the neighborhood of the key, which itself is an index in B
	vector<vector<double> > dadjf(s.size(), vector<double>());
	for(unsigned long i(0); i < s.size(); i++){
		nhs[s[i]];
		map<voxelType, vector<voxelType> >::iterator nhsIt(nhs.find(s[i]));
		for(unsigned long j(i + 1); j < s.size(); j++){
			if(inNeighborhood(s[i], s[j], B)){
				double sep(separation3D(s[i], s[j], voxdimsGlobal, voldimsGlobal));
				adjf[i].push_back(j);
				dadjf[i].push_back(sep);
				nhsIt->second.push_back(s[j]);
				adjf[j].push_back(i);
				dadjf[j].push_back(sep);
				nhs[s[j]].push_back(s[i]);
			}
		}
	}
	
	// initialize distances
	vector<double> d(s.size(), -1);
	bool foundBoundary(false);
	for(unsigned long i(0); i < s.size(); i++){
		vector<voxelType> nh(vesselBlocksInNeighborhood(B, s[i]));
		for(unsigned int j(0); j < nh.size(); j++){
			if(nhs.find(nh[j]) == nhs.end()){// if(!isIn(nh[j], s)) // is on boundary with other segments
				d[i] = 0;
				foundBoundary = true;
			}
		}
	}
	if(!foundBoundary){ // no boundary found; find farthest from center of mass
		voxelType cm(centerOfMassInSet(B, s));
		for(unsigned long i(0); i < s.size(); i++){
			if(s[i] == cm){
				d[i] = 0;
				i = s.size();
			}
		}
	}
	
	// find shortest distances (somewhat inefficiently)
	bool changed = true;
	while(changed){
		changed = false;
		for(unsigned long i(0); i < s.size(); i++){
			for(unsigned long j(0); j < adjf[i].size(); j++){
				unsigned long  k(adjf[i][j]);
				if(d[k] >= 0.0){ // some path to other segments is known
					double dist(d[k] + dadjf[i][j]);
					if(d[i] < 0.0 || d[i] > dist){
						d[i] = dist;
						changed = true;
					}
				}
			}
		}
	}

	double maxDist(0.0);
	for(unsigned int i(0); i < s.size(); i++){
		if(maxDist < d[i])
			maxDist = d[i];
	}
	for(unsigned int i(0); i < s.size(); i++){
		if(maxDist == d[i])
			return s[i];
	}
	return B.totalSize();
}

/* Converts the connected component and backbone index into a double for a unique label.
 * @backboneSetIndex the connected component to which the backbone belongs
 * @backboneIndex the index of the backbone within its connected component
 *
 * Converts the connected component and backbone index into a double for a unique label.  The index of the backbone within its connected component appears before the decimal point, and the connected component label appears after the  decimal point.
 * @return a double that can be used to uniquely label a backbone
 */
template <class T, class S, class R>
inline double backboneNumberName(T backboneSetIndex, S backboneIndex, R numBackboneSets){
	return backboneIndex + (backboneSetIndex + 1)/(double)pow(10.0, ceil(log10(numBackboneSets + 1)));
}

/* Creates a string to describe the subtree of segments.
 * @voxelVolume physical volume of a single voxel
 * @parent index of the parent segment of the <Backbone> of interest
 * @i index of the connected component in segregated backbones of the <Backbone> of interest (for naming purposes)
 * @numComponents number of connected components
 * @j the index within the connected component (`backbones`) of the <Backbone> of interest
 * @exclude list of the indices in `backbones` that should be excluded
 * @backbones list of <Backbone> objects in connected component `i`
 * @branchpoints sets of adjacent backbones after segmentation (i.e., the  connectivity information)
 * @volumes volume of each <Backbone> in `backbones`
 * @aveRad average radius of each <Backbone> in `backbones`
 *
 * Creates a string to describe the subtree of segments.  Ends with new line.
 */
string toStringSubtreeWithParent(double voxelVolume, unsigned long parent, unsigned long i, unsigned int numComponents, unsigned long j, vector<unsigned long long> &exclude,
                                 vector<Backbone<> > &backbones, // note that backbones is not const because of length() call
                                 const vector<vector<branchType> > &branchpoints,
                                 const vector<double> &volumes, const vector<double> &aveRad){
    vector<branchType> adj;
    for(unsigned long k(0); k < branchpoints.size(); k++){ // for each segment, get all adjacent segments
        if(isIn(j, branchpoints[k]))
            pushUniqueSet(branchpoints[k], adj);
    }
    for(unsigned long k(0); k < exclude.size(); k++)
        removeFrom(exclude[k], adj); // remove the parent segment (visited segment) from the list of adjacent segments
    stringstream strstr;
    double vol(voxelVolume*volumes[j]),
    len(backbones[j].length(voxdimsGlobal, voldimsGlobal));
    strstr << backboneNumberName(i, j, numComponents)
    << "\t" << vol << "\t" << len
    << "\t" << radFromVolLen(vol, len) << "\t" << aveRad[j];
    if(parent < backbones.size())
        strstr << "\t" << backboneNumberName(i, parent, numComponents);
    else
        strstr << "\tN/A"; // if there is no parent (roots will have no parent. parent argument is set to backbones.size() when this function is called with j as the index of the root in the current connected component.)
    strstr << "\t" << adj.size();
    for(unsigned long k(0); k < adj.size(); k++){
        strstr << "\t" << backboneNumberName(i, adj[k], numComponents);
        exclude.push_back(adj[k]);
    }
    strstr << endl;
    for(unsigned long k(0); k < adj.size(); k++) // add child strings
        strstr << toStringSubtreeWithParent(voxelVolume, j, i, numComponents, adj[k], exclude, backbones, branchpoints, volumes, aveRad);
		// this is like a depth-first search, but it avoids cycles by removing the visited segments
    return strstr.str();
}

/* Writes the vascular structural analysis with parent-child hierarchy to the specified file.
 * @outputTSVwithRootsFn the name of the file to write to
 * @lengthUnits name of length unit in images
 * @roots the index in `backbones` of the chosen root
 * @backbones the list of <BinaryVolume> indices of voxels for each backbone in each connected component
 * @branchpoints the list of branching points for each connected component
 * @volumes volume of meat for each segment associated with a backbone
 * @aveRad the average radius for each segment associated with a backbone
 *
 * Writes the vascular structural analysis to the specified file.  Uses `roots` to specify the parentless segment of each connected component.  Note that some segments may have only one child because of reticulations.
 */
void writeTSVwithRoots(string outputTSVwithRootsFn, string lengthUnits, const vector<unsigned long long> &roots,
                       vector<vector<Backbone<> > > &backbones, // note backbones is not const because of length() call
                       const vector<vector<vector<branchType> > > &branchpoints,
                       const vector<vector<double> > &volumes, const vector<vector<double> > &aveRad){
    
    double voxelVolume(voxdimsGlobal[0]*voxdimsGlobal[1]*voxdimsGlobal[2]);
    ofstream outFile(outputTSVwithRootsFn.c_str());
    outFile << " name"
    << "\t vol(cu." << lengthUnits << ")\t len(" << lengthUnits << ")"
    << "\t <r>_vl(" << lengthUnits << ")\t <r>_obs(" << lengthUnits <<")"
    << "\t par\t num_child\t children..."
    << endl;
    for(unsigned long i(0); i < backbones.size(); i++){ // for each connected component
        vector<voxelType> exclude(1, roots[i]);
        outFile << toStringSubtreeWithParent(voxelVolume, (unsigned int)backbones[i].size(), i, backbones.size(), roots[i], exclude, backbones[i], branchpoints[i], volumes[i], aveRad[i]);
    }
    outFile.close();
}

/* Writes the vascular structural analysis to the specified file.
 * @outputTSVfn the name of the file to write to
 * @lengthUnits name of length unit in images
 * @backbones the list of <BinaryVolume> indices of voxels for each backbone in each connected component
 * @branchpoints the list of branching points for each connected component
 * @volumes volume of meat for each segment associated with a backbone
 * @aveRad the average radius for each segment associated with a backbone
 *
 * Writes the vascular structural analysis to the specified file with adjacent segments.
 */
// NOTE: backbones input here is really the segregated backbones, segregated by connected components. 
void writeTSV(string outputTSVfn, string lengthUnits,
		vector<vector<Backbone<> > > &backbones, // note backbones is not const because of length() call
		const vector<vector<vector<branchType> > > &branchpoints,
		const vector<vector<double> > &volumes, const vector<vector<double> > &aveRad){
	double voxelVolume(voxdimsGlobal[0]*voxdimsGlobal[1]*voxdimsGlobal[2]);
	ofstream outFile(outputTSVfn.c_str());
	outFile << " name"
		<< "\t vol(cu." << lengthUnits << ")\t len(" << lengthUnits << ")"
		<< "\t <r>_vl(" << lengthUnits << ")\t <r>_obs(" << lengthUnits <<")"
		<< "\t num_adj\t adj..."
		<< endl;
	for(unsigned long i(0); i < backbones.size(); i++){ // for each connected component
		for(unsigned long j(0); j < backbones[i].size(); j++){ // for each backbone in the connected component
			vector<branchType> adj;
			for(unsigned long k(0); k < branchpoints[i].size(); k++){
				if(isIn(j, branchpoints[i][k]))
					pushUniqueSet(branchpoints[i][k], adj);
			}
			removeFrom(j, adj);
			double vol(voxelVolume*volumes[i][j]),
				len(backbones[i][j].length(voxdimsGlobal, voldimsGlobal));
			outFile << backboneNumberName(i, j, (unsigned int)backbones.size())
				<< "\t" << vol << "\t" << len
				<< "\t" << radFromVolLen(vol, len) << "\t" << aveRad[i][j]
				<< "\t" << adj.size();
			for(unsigned int k(0); k < adj.size(); k++)
				outFile << "\t" << backboneNumberName(i, adj[k], (unsigned int)backbones.size());
			outFile << endl;
		}
	}
	outFile.close();
}

/* Chooses a root for the connected component by selecting a backbone with the largest average radius.
 * @branchpoints list of branching points
 * @aveRad average radius of each segment in the connected component
 *
 * Searches for the largest average radius attained by any tip in the connected component.  Returns the index of the first segment that matches the largest average radius.  Note that there could be alternative means for choosing a root.  In case of an error, the size of aveRad is returned.
 * @return the index in `aveRad` of the root (size() of `aveRad` if error)
 */
voxelType assignRoot(const vector<vector<branchType> > &branchpoints, const vector<double> &aveRad){
	double maxRad(0.0);
	for(unsigned long i(0); i < aveRad.size(); i++){ // for each segment in the connected component
		unsigned int branchpointCount(0); // to keep track of the number of branchpoints that this segment is included in
		for(unsigned int j(0); j < branchpoints.size(); j++){ // for each branchpoint, each branchpoint is a vector of indices in `backbones` that are adjacent
			if(isIn(i, branchpoints[j])) // if the segment is included in the branchpoint
				branchpointCount++;
		}
		if(maxRad < aveRad[i] && branchpointCount < 2) // a root must be a tip with no more than one associated brach point
			maxRad = aveRad[i];
	}
	for(unsigned int i(0); i < aveRad.size(); i++){
		if(aveRad[i] == maxRad)
			return i; // return the segment with the largest average radius that is a tip. 
	}
	return (voxelType)aveRad.size();
}

/* Chooses a root for all supplied connected component by selecting a backbone with the largest average radius.
 * @branchpoints list of branching points for all connected components
 * @aveRad list of average radius of each segment in the connected components
 *
 * Searches for the largest average radius attained by any tip in the connected component.  Returns the index of the first segment that matches the largest average radius.  Note that there could be alternative means for choosing a root.  In case of an error, the size of aveRad is returned.
 * @return list of indices of the correesponding component in `aveRad` of the root (size() of component's `aveRad` if error)
 */
vector<voxelType> assignRoots(const vector<vector<vector<branchType> > > &branchpoints, const vector<vector<double> > &aveRad){
	vector<voxelType> roots;
	for(unsigned int i(0); i < aveRad.size(); i++) // for each connected component
		roots.push_back(assignRoot(branchpoints[i], aveRad[i]));
	return roots;
}

/* Attempts to extend the tip of <Backbone> `backbone` within the segment's meat `sMap`.
 * @B <BinaryVolume> that indicates vasculature
 * @backbone the <Backbone> of interest
 * @sMap the meat associated with `backbone`
 *
 * Attempts to extend the tip of <Backbone> `backbone` within the segment's meat `sMap` out to the voxel in `sMap` that is farthest from the meat of any adjacent segment.
 */
void extendTipInMeat(const BinaryVolume &B, Backbone<> &backbone, const map<voxelType, double> &sMap){
	vector<voxelType> s;
	s.reserve(sMap.size());
	for(map<voxelType, double>::const_iterator it(sMap.begin()); it != sMap.end(); it++)
		s.push_back(it->first);
	
	// erode using two critical vertebrae
	voxelType critTip(farthestFromOtherSegments(s, B)); // find critical vertebra farthest from branchpoint -- this is the target for the extension
	voxelType tipTip(backbone[0]);
	if(separationSq3D(critTip, backbone.back(), voxdimsGlobal, B) < separationSq3D(critTip, backbone[0], voxdimsGlobal, B))
		tipTip = backbone.back();
	vector<voxelType> cv(1, tipTip);
	cv.push_back(critTip);
	Backbone<> missingTip(erodeForBackbone(B, s, cv).getVertebrae());
	missingTip.organize(B);
	
	if(missingTip.size() > 1){ // includes more than just tipTip; stitch with existing backbone
		Backbone<> head, tail;
		if(missingTip[0] == tipTip && backbone[0] == tipTip){ // branchpoint at back; missingTip must be reversed
			head = missingTip.reversed();
			tail = backbone;
		}else if(missingTip[0] == tipTip && backbone.back() == tipTip){ // branchpoint at front
			head = backbone;
			tail = missingTip;
		}else if(missingTip.back() == tipTip && backbone[0] == tipTip){ // branchpoint at back
			head = missingTip;
			tail = backbone;
		}else if(missingTip.back() == tipTip && backbone.back() == tipTip){ // branchpoint at front; missingTip must be reversed
			head = backbone;
			tail = missingTip.reversed();
		}else{
			cout << "\n Error extendTips(): Could not properly match endpoints for backbone "
				<< makeString(backbone.getVertebrae()) << " to missingTip = " << makeString(missingTip.getVertebrae())
				<< "\n\t Leaving the Backbone as it is.\n"
				<< endl;
		}
		tail.removeIndex(0); head.add(tail);
		backbone = head;
		backbone.newLength(voxdimsGlobal, voldimsGlobal); // must recalculate length of backbones[i] since it has changed
	}
}

/* Analyzes the meat associated with each <Backbone> in `backbones` from the subset assigned to the thread, extending tips if possible.
 * @B <BinaryVolume> that indicates vasculature
 * @vertebraMap map of all voxels (keys) in any <Backbone> with a mapping to the label of the <Backbone> (values)
 * @backbones list of all <Backbone> objects
 * @isTip boolean list indicating if the associated <Backbone> in `backbones` is a tip
 * @branchpoints branching points with lists of the corresponding index in `backbones`
 * @vol list of the volume of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSolid list of the average radius (of the complete solid) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSurf list of the average radius (of the surface) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSolidByVertebra (output) per-vertebra average solid radius for each backbone
 * @aveRadSurfByVertebra (output) per-vertebra average surface radius for each backbone
 * @offset thread offset
 * @numThreads total number of threads
 *
 * Analyzes the meat associated with each <Backbone> in `backbones` from the subset assigned by `offset` and `numThreads`, extending tips if possible.  Updates the corresponding elements of `vol` and `aveRadSolid` and `aveRadSurf`.
 * @return true to indicate thread completion
 */
bool analyzeMeatAfterExtendingTipsParallel(const BinaryVolume &B, const map<voxelType, unsigned long> &vertebraMap, vector<Backbone<> > &backbones, const vector<bool> &isTip, const vector<vector<branchType> > &branchpoints,
		vector<double> &vol, vector<double> &aveRadSolid, vector<double> &aveRadSurf, vector<vector<double> > &aveRadSolidByVertebra, vector<vector<double> > &aveRadSurfByVertebra, unsigned int offset, unsigned int numThreads
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm, const vector<vector<unsigned char> > &backboneColors, BinaryVolume &Btemp
#endif
		){
	for(unsigned long i(offset); i < backbones.size(); i += numThreads){ // for each backbone in the subset assigned to this thread
		// get the meat associated with the backbone
		map<voxelType, double> sMap(segmentMeatFromAdjacencies(B, vertebraMap, backbones, i)); // key is voxel, value is fraction of each voxel that belongs to the backbone of interest. 
		
		// calculate volume of the segment
		vol[i] = 0;
		for(map<voxelType, double>::iterator sMapIt(sMap.begin()); sMapIt != sMap.end(); sMapIt++)
			vol[i] += sMapIt->second; // sum of the fractions of the voxels that belong to the backbone of interest
		
#if TRY_WX == 1
		vector<voxelType> sVec;
		sVec.reserve(sMap.size());
		for(map<voxelType, double>::iterator sMapIt(sMap.begin()); sMapIt != sMap.end(); sMapIt++)
			sVec.push_back(sMapIt->first);
		Btemp.f(sVec);
		setPointSetColor(sVec, pm, 1.0 - backboneColors[i][0], 1.0 - backboneColors[i][1], 1.0 - backboneColors[i][2], 0.03);
#endif
		// extend tips, if necessary
		// So after segmentation of the meat, each backbone owns a region described by sMap. Now we try to push the backbone’s endpoint further into that region just in case it was previously truncated.
		if(isTip[i])
			extendTipInMeat(B, backbones[i], sMap);
		
		// unsnarl the backbone
		double snarledLength(backbones[i].length(voxdimsGlobal, voldimsGlobal)); // basically distance if we were to walk along the backbone (straightened out)
		localUnsnarl(backbones[i], B, voxdimsGlobal, voldimsGlobal);
		cout << "\n" + makeString(i) + "\t" + makeString(snarledLength) + "\t" + makeString(backbones[i].length(voxdimsGlobal, voldimsGlobal));
		
#if TRY_WX == 1
		setPointSetColor(backbones[i].getVertebrae(), pm, backboneColors[i][0], backboneColors[i][1], backboneColors[i][2], 0.5);
#endif
		
		// find the average radius
		averageRadius(B, sMap, backbones[i], aveRadSolid[i], aveRadSurf[i]);
		averageRadiusByVertebra(B, sMap, backbones[i], aveRadSolidByVertebra[i], aveRadSurfByVertebra[i]);
	}
	
	return true;
}

/* Analyzes the meat associated with each <Backbone> in `backbones`, extending tips if possible.
 * @B <BinaryVolume> that indicates vasculature
 * @backbones list of all <Backbone> objects
 * @branchpoints branching points with lists of the corresponding index in `backbones`
 * @vol list of the volume of the meat associated with the corresponding <Backbone> in `backbones` in terms of number of voxels
 * @aveRadSolid list of the average radius (of the complete solid) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSurf list of the average radius (of the surface) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSolidByVertebra (output) per-vertebra average solid radius for each backbone
 * @aveRadSurfByVertebra (output) per-vertebra average surface radius for each backbone
 *
 * Analyzes the meat associated with each <Backbone> in `backbones`, extending tips if possible.  Updates the corresponding elements of `vol` and `aveRadSolid` and `aveRadSurf`.
 * @return true to indicate thread completion
 */
void analyzeMeatAfterExtendingTips(const BinaryVolume &B, vector<Backbone<> > &backbones, const vector<vector<branchType> > &branchpoints,
		vector<double> &vol, vector<double> &aveRadSolid, vector<double> &aveRadSurf, vector<vector<double> > &aveRadSolidByVertebra, vector<vector<double> > &aveRadSurfByVertebra
#if TRY_WX == 1
		, map<voxelType, unsigned int> &pm
#endif
		){
	aveRadSolidByVertebra.resize(backbones.size());
	aveRadSurfByVertebra.resize(backbones.size());
	// find tips
	vector<bool> isTip(backbones.size(), false);
	for(unsigned int i(0); i < branchpoints.size(); i++){
		for(unsigned int j(0); j < branchpoints[i].size(); j++)
			isTip[branchpoints[i][j]] = !isTip[branchpoints[i][j]]; // tips are listed exactly once; nontips are listed exactly twice
	} // tips would be true after this because it is flipped only once 
	
	// construct map of each vertebra
	map<voxelType, unsigned long> vertebraMap;
	for(unsigned int i(0); i < backbones.size(); i++){
		for(unsigned int v(0); v < backbones[i].size(); v++)
			vertebraMap[backbones[i][v]] = i; // highest i will win any tie
	} // maps voxel to the index of the backbone to which it belongs
	
#if TRY_WX == 1
	vector<vector<vector<unsigned char> > > segCol(segColors(vector<vector<Backbone<> > >(1, backbones), vector<vector<vector<branchType> > >(1, branchpoints)));
	vector<vector<unsigned char> > backboneColor(segCol[0]);
	//setAllPointsTo(0.5, 0.5, 0.5, 1.0);
#endif
	
	cout << "\n i\t snarL\t unsnarL" << endl;
	
	BinaryVolume Btemp(B); // !!!!
	vector<future<bool> > locThreads;
	for(unsigned int th(0); th < numLocThreads; th++)
		locThreads.push_back(async(launch::async, analyzeMeatAfterExtendingTipsParallel, ref(B), ref(vertebraMap), ref(backbones), ref(isTip), ref(branchpoints), ref(vol), ref(aveRadSolid), ref(aveRadSurf), ref(aveRadSolidByVertebra), ref(aveRadSurfByVertebra), th, numLocThreads
#if TRY_WX == 1
		, ref(pm), ref(backboneColor), ref(Btemp) // !!!!
#endif
		));
	for(unsigned int th(0); th < locThreads.size(); th++)
		locThreads[th].get();
	cerr << "\n Btemp.totalTrue()*voxDim = " << Btemp.totalTrue()*voxdimsGlobal[0]*voxdimsGlobal[1]*voxdimsGlobal[2] << endl;
	cout << endl;
}

/* Returns the segregated list of <Backbones> `backbones` into its connected components.
 * @backbones the list of <Backbone> objects
 * @branchpoints sets of adjacent backbones after segmentation (i.e., the  connectivity information)
 * @volumes list of the volume of the meat associated with the corresponding <Backbone> in `backbones` in terms of number of voxels
 * @aveRadSolid list of the average radius (of the complete solid) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSurf list of the average radius (of the surface) of the meat associated with the corresponding <Backbone> in `backbones`
 * @aveRadSolidByVertebra per-vertebra average solid radius for each backbone
 * @aveRadSurfByVertebra per-vertebra average surface radius for each backbone
 * @segregatedVolumes segregated `volumes`
 * @segregatedAveRadSolid segregated aveRadSolid
 * @segregatedAveRadSurf segregated aveRadSurf
 * @segregatedAveRadSolidByVertebra (output) segregated per-vertebra average solid radius
 * @segregatedAveRadSurfByVertebra (output) segregated per-vertebra average surface radius
 * @segregatedBranchpoints segregated branching points
 *
 * Returns the segregated list of <Backbones> `backbones` into its connected components.
 * @ return the segregated list of <Backbones> `backbones` into its connected components
 */
vector<vector<Backbone<> > > segregateBackbones(const vector<Backbone<>> &backbones, const vector<vector<branchType> > &branchpoints, const vector<double> &volumes, const vector<double> &aveRadSolid, const vector<double> &aveRadSurf, const vector<vector<double> > &aveRadSolidByVertebra, const vector<vector<double> > &aveRadSurfByVertebra, vector<vector<double> > &segregatedVolumes, vector<vector<double> > &segregatedAveRadSolid, vector<vector<double> > &segregatedAveRadSurf, vector<vector<vector<double> > > &segregatedAveRadSolidByVertebra, vector<vector<vector<double> > > &segregatedAveRadSurfByVertebra, vector<vector<vector<branchType> > > &segregatedBranchpoints){
	
	// construct backbone adjacency lists from branchpoints
	vector<vector<voxelType> > adj(backbones.size(), vector<voxelType>());
	for(unsigned int i(0); i < branchpoints.size(); i++){
		for(unsigned int j(0); j < branchpoints[i].size(); j++){
			for(unsigned int k(j + 1); k < branchpoints[i].size(); k++){
				adj[branchpoints[i][j]].push_back(branchpoints[i][k]); // duplicates are not a problem
				adj[branchpoints[i][k]].push_back(branchpoints[i][j]);
			}
		}
	}
	
	// construct connected components and mapping for branching points
	vector<bool> included(backbones.size(), false);
	vector<unsigned long> branchpointMapping(backbones.size(), backbones.size());
	vector<unsigned long> ccMapping(backbones.size(), backbones.size()); // maps each backbone to the connected component index it belongs to (same value means same connected component, each index represents a backbone.)
	vector<vector<unsigned long> > ccs; // list of connected components, each component is a list of backbone indices.
	vector<unsigned long> ccSizes;
	for(unsigned long i(0); i < backbones.size(); i++){
		if(included[i])
			continue;
		
		// classic BFS like algorithm
		included[i] = true;
		branchpointMapping[i] = 0;
		ccMapping[i] = ccs.size();
		ccs.push_back(vector<unsigned long>(1, i));
		for(unsigned long c(0); c < ccs.back().size(); c++){
			unsigned long j(ccs.back()[c]);
			for(unsigned long k(0); k < adj[j].size(); k++){
				unsigned long n(adj[j][k]);
				if(!included[n]){
					included[n] = true;
					branchpointMapping[n] = ccs.back().size();
					ccMapping[n] = ccMapping[i];
					ccs.back().push_back(n);
				}
			}
		}
		ccSizes.push_back(ccs.back().size());
	}
	
	// sort ccs based on ccSizes; bubble works well enough
	vector<unsigned long> ccInd(ccs.size());
	for(unsigned long cc(0); cc < ccInd.size(); cc++)
		ccInd[cc] = cc;
	bool changed(true);
	while(changed){
		changed = false;
		for(unsigned long cc(1); cc < ccSizes.size(); cc++){
			if(ccSizes[ccInd[cc - 1]] < ccSizes[ccInd[cc]]){
				unsigned long ccTemp(ccInd[cc - 1]);
				ccInd[cc - 1] = ccInd[cc];
				ccInd[cc] = ccTemp;
				changed = true;
			}
		}
	}
	
	vector<unsigned long> ccIndInv(ccs.size());
	for(unsigned long cc(0); cc < ccInd.size(); cc++)
		ccIndInv[ccInd[cc]] = cc;
	
	
	vector<vector<unsigned long> > ccsCopy(ccs);
	for(unsigned long cc(0); cc < ccs.size(); cc++)
		ccs[cc] = ccsCopy[ccInd[cc]];
	for(unsigned long i(0); i < ccMapping.size(); i++)
		ccMapping[i] = ccIndInv[ccMapping[i]];
	
	// create segregated data structures
	vector<vector<Backbone<> > > segregatedBackbones(ccs.size(), vector<Backbone<> >());
	segregatedVolumes = segregatedAveRadSolid = segregatedAveRadSurf = vector<vector<double> >(ccs.size(), vector<double>());
	segregatedAveRadSolidByVertebra = segregatedAveRadSurfByVertebra = vector<vector<vector<double> > >(ccs.size(), vector<vector<double> >());
	for(unsigned long cc(0); cc < ccs.size(); cc++){
		segregatedBackbones[cc] = vector<Backbone<> >(ccs[cc].size(), Backbone<>());
		segregatedVolumes[cc] = segregatedAveRadSolid[cc] = segregatedAveRadSurf[cc] = vector<double>(ccs[cc].size(), 0);
		segregatedAveRadSolidByVertebra[cc] = vector<vector<double> >(ccs[cc].size());
		segregatedAveRadSurfByVertebra[cc] = vector<vector<double> >(ccs[cc].size());
		for(unsigned long c(0); c < ccs[cc].size(); c++){
			unsigned long i(ccs[cc][c]);
			segregatedBackbones[cc][c] = backbones[i];
			segregatedVolumes[cc][c] = volumes[i];
			segregatedAveRadSolid[cc][c] = aveRadSolid[i];
			segregatedAveRadSurf[cc][c] = aveRadSurf[i];
			segregatedAveRadSolidByVertebra[cc][c] = aveRadSolidByVertebra[i];
			segregatedAveRadSurfByVertebra[cc][c] = aveRadSurfByVertebra[i];
		}
		
	}
	segregatedBranchpoints = vector<vector<vector<branchType> > >(ccs.size(), vector<vector<branchType> >());
	for(unsigned long i(0); i < branchpoints.size(); i++){
		unsigned long cc(ccMapping[branchpoints[i][0]]);
		segregatedBranchpoints[cc].push_back(branchpoints[i]);
		for(unsigned long j(0); j < branchpoints[i].size(); j++)
			segregatedBranchpoints[cc].back()[j] = branchpointMapping[branchpoints[i][j]];
	}
	
	/*// !!!! temp segregatedBranchpoints check
	for(unsigned long i(0); i < segregatedBranchpoints.size(); i++){
		for(unsigned long j(0); j < segregatedBranchpoints[i].size(); j++){
			cout << "\n i=" << i << " j=" << j << ":" << makeString(segregatedBranchpoints[i][j]) << ":";
			for(unsigned long x(0); x < segregatedBranchpoints[i][j].size(); x++){
				Backbone<> X(segregatedBackbones[i][segregatedBranchpoints[i][j][x]]);
				if(X.empty())
					continue;
				for(unsigned int y(x + 1); y < segregatedBranchpoints[i][j].size(); y++){
					Backbone<> Y(segregatedBackbones[i][segregatedBranchpoints[i][j][y]]);
					if(Y.empty())
						continue;
					unsigned int bestX(0), bestY(0);
					double bestSep(separation3D(X[bestX], Y[bestY], voxdimsGlobal, voldimsGlobal));
					for(unsigned long vx(0); vx < X.size(); vx++){
						for(unsigned long vy(0); vy < Y.size(); vy++){
							double sep(separation3D(X[vx], Y[vy], voxdimsGlobal, voldimsGlobal));
							if(bestSep > sep)
								bestSep = sep;
						}
					}
					cout << " " << bestSep;
				}
			}
		}
		cout << endl;
	}*/
	
	return segregatedBackbones;
}

/* Saves a sanity check visualization of the binary volume before backbone erosion.
 * @B the <BinaryVolume> to visualize
 * @outputBase base filename for output files (will append "_vessels.png" and "_vessels_coords.txt")
 *
 * Saves a PNG visualization and a text file with voxel coordinates for debugging/visualization purposes.
 * Only runs when TRY_WX == 0 (non-GUI mode).
 */
 void saveVesselSanityCheck(const BinaryVolume &B, const string &outputBase){
	#if TRY_WX == 0
		generalStatusUpdate("Saving sanity check visualization of vessels before backbone erosion...");
		writePNGBinaryVolume(B, outputBase + "_vessels.png");
		
		// COMMENTED OUT - but maybe useful in the future for figure visualizations
		// Also save 3D voxel coordinates for Python visualization
		//ofstream voxelCoordsOut((outputBase + "_vessels_coords.txt").c_str());
		//voxelCoordsOut << "# Vessel voxel coordinates (not physical coordinates) - one per line\n";
		//voxelCoordsOut << "# Total voxels: " << B.totalTrue() << "\n";
		//voxelCoordsOut << "# Volume dimensions: " << B.getSize(0) << " x " << B.getSize(1) << " x " << B.getSize(2) << "\n";
		
		//voxelType i(0);
		//unsigned long long voxelCount(0);
		//while(i < B.totalSize()){
			//i = B.findFirstAtOrAfter(i);
			//if(i >= B.totalSize())
				//break;
			//voxelType xCoord(B.x(i)), yCoord(B.y(i)), zCoord(B.z(i));
			//voxelCoordsOut << xCoord << "\t" << yCoord << "\t" << zCoord << "\n";
			//voxelCount++;
			//i++; // move to next index
		//}
		//voxelCoordsOut.close();
		
		generalStatusUpdate("Sanity check image saved to " + outputBase + "_vessels.png");
		//generalStatusUpdate("Vessel coordinates saved to " + outputBase + "_vessels_coords.txt (" + makeString(voxelCount) + " voxels)");
	#endif
}

/* Analyzes the vascular structure from images imStart to imEnd in imageDir (expecting the file naming convention XXXXX.png)
 * @imageDir the directory path that contains the images of interest
 * @imStart the start index of the image to be included
 * @imEnd the end index of the image to be included
 * @outputTSVfnBase output file base name
 * @voxdims the dimensions of the voxels
 * @lengthUnit the length unit of voxdims
 * @thresh the intensity threshold to segregate vascular and nonvascular voxels
 * @critFrac the critical fraction of vascular voxels within a given sphere surface frontier that is required to continue growing the sphere
 *
 * Analyzes the vasculature in the specified images; identifies vascular images from voxel intensities; removes small components; segments the vasculature into disjoint intrasphere and intersphere spaces; erodes vascular meat to find backbones; analyzes structure and connectivity.
*/

void analyzeVascularStructure(string imageDir, int imStart, int imEnd,
		string outputTSVfnBase,
		double voxdims[], string lengthUnit,
		double thresh, double critFrac){

	chrono::steady_clock::time_point analyzeTime(getNow()), tempTime;
	srandSet(112358);
	
#if TRY_WX == 1
	wxGetApp().myFrame->updateStatus("Loading data...");
#endif
	
	// first read the PNG images into a Lumens object which is a 3D array of doubles. 
	// see pngMinSurf.h and Lumens.h for more details. 
	// TODO: revisit the Lumens.h assignment operator to ensure no memory leaks
	Lumens L(readPNGImages(imageDir, imStart, imEnd));
	if(L.totalSize() == 0){
		generalStatusUpdate("\n No data to analyze after trying to load from " + imageDir);
		return;
	}
	BinaryVolume B(threshThrash(L, thresh));
	cout << "\n B.totalSize() = " << B.totalSize() << "\n B.totalTrue() = " << B.totalTrue() << endl;

	// sanity check added by Kai to visualize projection of vessel mask and save vessel mask coordinates for downstream visualization. 
	saveVesselSanityCheck(B, outputTSVfnBase); // TODO: for efficient checks of optimal threshold, we can add an option to stop here
	
	lengthUnitGlobal = lengthUnit;
	for(unsigned int i(0); i < 3; i++){
		voxdimsGlobal[i] = voxdims[i];
		voldimsGlobal[i] = B.getSize(i); // size of volume from binary volume class.
	}
	
	// finding largest connected components
#if TRY_WX == 1
	generalStatusUpdate("Initializing points...");
	map<voxelType, unsigned int > pm;
	vector<double> pc(pointsWithColors(B, pm));
	wxGetApp().myFrame->m_canvas->newPoints(pc);
	generalStatusUpdate("Showing all points.  Finding largest connected components...");
#else
	generalStatusUpdate("Finding largest connected components...");
#endif
	BinaryVolume lccsB = removeSmallestConnectedComponents(B, 9 // lccsB is another Binary volume that will contain the 9 largest connected components. 
#if TRY_WX == 1
		, pm
#endif
		);
#if TRY_WX == 1
	generalStatusUpdate("Showing all points. Initializing points for largest ccs... (expect no visualization interaction)");
	pc = pointsWithColors(lccsB, pm);
	wxGetApp().myFrame->m_canvas->newPoints(pc);
	generalStatusUpdate("Showing largest connected components.  Finding outside...");
#else
	generalStatusUpdate("Finding outside...");
#endif
	
	// filling gaps
	B = lccsB; // update the binary volume to only contain the largest connected components.
	generalStatusUpdate("Filling gaps...");
	voxelType numFillings(fillInternalFalses(B, findOutside(B))); 
	// findOutside(B) gives you another binary volume where the largest connected component among the non-vascular parts are set to true. 
	// then fillInternalFalses fills all holes in the vascular volume by taking the background and simply flipping all bits. 
	// critically, findOutside(B) should give us the true background - to make sure this is the case we make sure the largest connected component returned is > 50% of the entire volume. 
	cout << "\n Filled " << numFillings << " hole" << makePluralSuffix(numFillings) << "."
			<< "\n B.totalSize() = " << B.totalSize() << "\n B.totalTrue() = " << B.totalTrue() << endl;
	
#if TRY_WX == 1
	generalStatusUpdate("Initializing points with gaps filled... (expect no visualization interaction)");
	pc = pointsWithColors(B, pm);
	 wxGetApp().myFrame->m_canvas->newPoints(pc);
#endif

	//cout << "RAND_MAX = " << RAND_MAX << "  B.totalSize() = " << B.totalSize() << endl;
	//RAND_MAX = 32767  B.totalSize() = 16252928. This was causing issues for windows machines in sphereCoarsen function. 
	// rand()%sB.totalSize() would never return a majority of indices because B.totalSize() is larger than the max random value. 
	// what this means algorithmically is that a majority of the volume was left without spheres because the seed was never placed in the indices > 32767. 
	// downstream I think this means that there is potentially a very large intersphere space without any spheres and the algorithm got stuck processing that single intersphere space. 

	// coarsending vascular structure
	generalStatusUpdate("Coarsening...");
	BinaryVolume uB(B);
	// Roughly vectorizes the vasculature in <BinaryVolume> `B` into a map of seed to radius of spheres.
	map<voxelType, double> rv(sphereCoarsen(B, uB, critFrac // rv for rough vectorization: key is voxel, value is radius of sphere
#if TRY_WX == 1
		, pm
#endif
		));
	cout << endl << rv.size() << " spheres for critFrac = " << critFrac << endl;
	
	// Sanity check: sphere volume vs vessel volume (dispersion check; radius in physical units)
	// a properly functioning sphereCoarsen should result in a sphere volume that covers at least most of the vessel volume. Note that spheres are nonoverlapping
	// also print max, min, average, and standard deviation of the sphere radii to understand the distribution of sphere radii. Save the sphere radii to an output file too. 
	{
		double voxVol(voxdimsGlobal[0] * voxdimsGlobal[1] * voxdimsGlobal[2]);
		double vesselVol(B.totalTrue() * voxVol);
		double totalSphereVol(0.0);
		const double fourThirdsPi(4.0 / 3.0 * acos(-1.0));
		double sum_r(0.0), sq_sum_r(0.0), max_r(0.0), min_r(0.0);
		vector<double> sphereRadii;
		for(map<voxelType, double>::const_iterator it(rv.begin()); it != rv.end(); it++){
			double r(it->second);
			totalSphereVol += fourThirdsPi * r * r * r;
			sum_r += r;
			sq_sum_r += r * r;
			if(r > max_r) max_r = r;
			if(r < min_r) min_r = r;
			sphereRadii.push_back(r);
		}
		double ratio((vesselVol > 0.0) ? (totalSphereVol / vesselVol) : 0.0);
		cout << "[sanity check] sphere volume / vessel volume (should be > 0.8) = " << ratio << endl;
		cout << "Max sphere radius = " << max_r << endl;
		cout << "Min sphere radius = " << min_r << endl;
		double mean_r = sum_r / sphereRadii.size();
		cout << "Average sphere radius = " << mean_r << endl;
		// Var[X] = E[X^2] - E[X]^2
		cout << "Standard deviation of sphere radii = " << std::sqrt(sq_sum_r / sphereRadii.size() - mean_r * mean_r) << endl;
		ofstream sphereRadiiFile(outputTSVfnBase + "_sphere_radii.tsv");
		sphereRadiiFile << "sphere_radius" << endl;
		for(double r : sphereRadii)
			sphereRadiiFile << r << endl;
		sphereRadiiFile.close();
	}
	
	// finding adjacent sphere
	generalStatusUpdate("Coarsened (" + makeString(rv.size()) + " points).  Finding adjacent spheres...");
	// uB at this point contains true where voxels are unvisited during the sphere coarsening. i.e. they are not in a sphere. 
	map<voxelType, vector<voxelType> > adj(adjacentSpheres(B, uB, rv
#if TRY_WX == 1
		, pm
#endif
		)); // adj is the list of intersphere connections between spheres: key is seed for intersphere meat, value is list of intrasphere seeds (centers)
	
	// removing obvious false tips
	generalStatusUpdate("Found adjacent spheres.  "
		+ makeString(adj.size()) + " connections to spheres. Removing false tips...");
	map<voxelType, vector<voxelType> > scs(sphereConnectionSeeds(adj)); // scs is the list of intersphere connections in terms of spheres: key is the intrasphere seed, value is the list of adjacent intersphere seeds (reverse of adj)
	unsigned long prevAdjTips(adj.size());
	removeFalseTipsByContext(B, uB, rv, adj, scs
#if TRY_WX == 1
		, pm
#endif
		); // note that we pass uB and edit it to mark the invalid intersphere space as visited so we don't use it later. 

#if TRY_WX == 1
	for(map<voxelType, vector<voxelType> >::iterator it(adj.begin()); it != adj.end(); it++){
		double x(0.0), y(0.0), z(0.0);
		normalizedCoord(it->first, x, y, z);
		wxGetApp().myFrame->m_canvas->addNumber(it->first, x, y, z, 0.0, 0.5, 0.0, 1.0);
	}
#endif
	
	generalStatusUpdate("Removed " + makeString(prevAdjTips - adj.size()) + " suspicious tip" + makePluralSuffix(prevAdjTips - adj.size()) + ". Finding critical vertebrae...");
	
	map<voxelType, vector<voxelType> > critVert(findCriticalVertebrae(B, uB, rv, adj, scs));
	// maps intersphere seeds to the voxel that interfaces a intrasphere space, notably the interfaces are disconnected so we can find junctions.
#if TRY_WX == 1
	for(map<voxelType, vector<voxelType> >::iterator it(critVert.begin()); it != critVert.end(); it++)
		setPointSetColor(it->second, pm, 0.0, 0.0, 0.5, 1.0);
#endif
	
	// eroding for raw backbones
	generalStatusUpdate("Found " + makeString(critVert.size()) + " critical vertebra(e).  "
		+ makeString(adj.size()) + " connection(s) between spheres. Eroding for backbones...");
	map<voxelType, Backbone<> > rawBackbones(backboneErosion(B, uB, adj, scs, critVert
#if TRY_WX == 1
		, pm
#endif
		));
	// rawBackbones is a map of (intersphere or intrasphere) seed to the corresponding backbone segment of that region. 
	// at this point it is not ordered. 
	
#if TRY_WX == 1
	//for(map<voxelType, Backbone<> >::iterator it(rawBackbones.begin()); it != rawBackbones.end(); it++)
	//	setPointSetColor(it->second.getVertebrae(), pm, 0.5, 0.0, 0.0, 1.0);
	for(map<voxelType, vector<voxelType> >::iterator it(critVert.begin()); it != critVert.end(); it++)
		setPointSetColor(it->second, pm, 0.0, 0.0, 0.5, 1.0);
#endif

	generalStatusUpdate("Found " + makeString(rawBackbones.size()) + " raw (not segmented or segregated) backbones. Segmenting backbones...");
	vector<vector<branchType> > branchpoints;
	vector<Backbone<> > backbones(segmentBackbones(B, rawBackbones, critVert, branchpoints));
	// `branchpoints` (indices in `backbones` that are adjacent) edited in segmentBackbones function
	// Segments, organizes, and produces connectivity information from raw sets of voxels that make up backbones.  Segmentation involves identifying branching points and the connections between branching points.
	
	generalStatusUpdate("Segmented " + makeString(backbones.size()) + " backbone" + makePluralSuffix(backbones.size()) + " with "
		+ makeString(branchpoints.size()) + " branch point" + makePluralSuffix(branchpoints.size()) + " (before segregation).");
	
	generalStatusUpdate("Analyzing meat (and extending tips)...");
	
	vector<double> volumes(backbones.size(), 0),
		aveRadSolid(backbones.size(), 0),
		aveRadSurf(backbones.size(), 0);
	vector<vector<double> > aveRadSolidByVertebra, aveRadSurfByVertebra; // for each backbone, stores a vector of radius values corresponding to the radius at each voxel along the backbone. 
	vector<vector<double> > segregatedVolumes,
		segregatedAveRadSolid,
		segregatedAveRadSurf;
	vector<vector<vector<double> > > segregatedAveRadSolidByVertebra, segregatedAveRadSurfByVertebra; // segregates the per-vertebra radius vectors for each backbone into connected components. 
	vector<vector<vector<branchType> > > segregatedBranchpoints; // segregates branchpoints by connected components (branchpoints indicate for each segment which backbones indices are adjacent).

	analyzeMeatAfterExtendingTips(B, backbones, branchpoints, volumes, aveRadSolid, aveRadSurf, aveRadSolidByVertebra, aveRadSurfByVertebra
#if TRY_WX == 1
			, pm
#endif
			);

	vector<vector<Backbone<> > > segregatedBackbones(segregateBackbones(backbones, branchpoints, volumes, aveRadSolid, aveRadSurf, aveRadSolidByVertebra, aveRadSurfByVertebra, segregatedVolumes, segregatedAveRadSolid, segregatedAveRadSurf, segregatedAveRadSolidByVertebra, segregatedAveRadSurfByVertebra, segregatedBranchpoints));
	backbonesGlobal = segregatedBackbones;
	branchpointsGlobal = segregatedBranchpoints;
	volumesGlobal = segregatedVolumes;
	aveRadSolidGlobal = segregatedAveRadSolid;
	aveRadSurfGlobal = segregatedAveRadSurf;
	aveRadSolidByVertebraGlobal = segregatedAveRadSolidByVertebra;
	aveRadSurfByVertebraGlobal = segregatedAveRadSurfByVertebra;
	
	// !!!! volume check
	double voxVol(voxdimsGlobal[0]*voxdimsGlobal[1]*voxdimsGlobal[2]),
			volB(B.totalTrue()*voxVol),
			volSum(0),
			segVolSum(0);
	for(unsigned long i(0); i < volumes.size(); i++)
		volSum += volumes[i];
	for(unsigned long i(0); i < segregatedVolumes.size(); i++){
		for(unsigned long j(0); j < segregatedVolumes[i].size(); j++)
			segVolSum += segregatedVolumes[i][j];
	}
	volSum *= voxVol;
	segVolSum *= voxVol;
	string volUnits(lengthUnitGlobal + "^3");
	cerr << "\n\t volume check:"
			<< "\n\t\t voxVol = " << voxVol << " " << volUnits
			<< "\n\t\t volB = " << volB << " " << volUnits
			<< "\n\t\t volSum = " << volSum << " " << volUnits
			<< "\n\t\t segVolSum = " << segVolSum << " " << volUnits
			<< endl;
	
	cout << "\n cc\t #backbones\t #branch points" << endl;
	for(unsigned int i(0); i < segregatedBackbones.size(); i++)
		cout << i << "\t" << segregatedBackbones[i].size() << "\t" << segregatedBranchpoints[i].size() << endl;
	
	vector<vector<vector<bool> > > branchAtFront(vector<vector<vector<bool> > >(segregatedBranchpoints.size(), vector<vector<bool> >()));
	for(unsigned int i(0); i < segregatedBranchpoints.size(); i++){ // over each cc
		
		branchAtFront[i] = vector<vector<bool> >(segregatedBranchpoints[i].size(), vector<bool>());
		
		for(unsigned int j(0); j < segregatedBranchpoints[i].size(); j++){ // over each branch point
			
			branchAtFront[i][j] = vector<bool>(segregatedBranchpoints[i][j].size(), false);
			
			voxelType nearBranchpoint(segregatedBackbones[i][segregatedBranchpoints[i][j][0]][0]),
				nh(segregatedBackbones[i][segregatedBranchpoints[i][j][1]][0]),
				nt(segregatedBackbones[i][segregatedBranchpoints[i][j][1]].back());
			
			if(!inNeighborhood(nearBranchpoint, nh, B) && !inNeighborhood(nearBranchpoint, nt, B))
				nearBranchpoint = segregatedBackbones[i][segregatedBranchpoints[i][j][0]].back();
			
			for(unsigned int k(0); k < segregatedBranchpoints[i][j].size(); k++){
				voxelType f(segregatedBackbones[i][segregatedBranchpoints[i][j][k]][0]), b(segregatedBackbones[i][segregatedBranchpoints[i][j][k]].back());
				double distFront(separation3D(nearBranchpoint, f, voxdimsGlobal, voldimsGlobal)), distBack(separation3D(nearBranchpoint, b, voxdimsGlobal, voldimsGlobal));
				branchAtFront[i][j][k] = distFront <= distBack;
			}
			
		}
		
	}
	branchAtFrontGlobal = branchAtFront;
	
	writeTSV(outputTSVfnBase + ".tsv", lengthUnit, segregatedBackbones, segregatedBranchpoints, segregatedVolumes, segregatedAveRadSurf); // inputs are already segregated by connected components. 
	// this is an undirected graph representation (cycles can be identified via graph traversal)
	
	vector<voxelType> roots(assignRoots(segregatedBranchpoints, segregatedAveRadSurf));
	writeTSVwithRoots(outputTSVfnBase + "_withRoots.tsv", lengthUnit, roots, segregatedBackbones, segregatedBranchpoints,  segregatedVolumes, segregatedAveRadSurf);
	// this is a tree-like representation (no cycles are represented here)
	
#if TRY_WX == 1

	vector<double> newNumbers;
	for(unsigned int i(0); i < segregatedBackbones.size(); i++){
		for(unsigned int j(0); j < segregatedBackbones[i].size(); j++){
			voxelType com(centerOfMassInSet(B, segregatedBackbones[i][j].getVertebrae()));
			double x(0.0), y(0.0), z(0.0);
			normalizedCoord(com, x, y, z);
			newNumbers.push_back(x);
			newNumbers.push_back(y);
			newNumbers.push_back(z);
			newNumbers.push_back(0.0); // r
			newNumbers.push_back(0.0); // g
			newNumbers.push_back(0.0); // b
			newNumbers.push_back(1.0); // a
			newNumbers.push_back(backboneNumberName(i, j, (unsigned int)segregatedBackbones.size()));
		}
	}
	wxGetApp().myFrame->m_canvas->newNumbers(newNumbers);
#endif


	// Per-backbone .dat: name (same as TSV), vertebra index, x, y, z, r_solid, r_surf. One block per backbone, blank line between.
	// IMPORTANT: this is NOT the x y z coordinates in physical units. These are the x y z indices in voxel space if we imagine the volume as a 3D grid.
	ofstream voxelOut(outputTSVfnBase + ".dat");
	for(unsigned int connectedComponent(0); connectedComponent < segregatedBackbones.size(); connectedComponent++){
		for(unsigned int backboneIndex(0); backboneIndex < segregatedBackbones[connectedComponent].size(); backboneIndex++){
			voxelOut << "# " << backboneNumberName(connectedComponent, backboneIndex, (unsigned int)segregatedBackbones.size()) << endl;
			voxelOut << "idx\tx\ty\tz\tr_solid\tr_surf" << endl;
			const vector<double> &radSolid = aveRadSolidByVertebraGlobal[connectedComponent][backboneIndex];
			const vector<double> &radSurf = aveRadSurfByVertebraGlobal[connectedComponent][backboneIndex];
			for(unsigned int vertebra(0); vertebra < segregatedBackbones[connectedComponent][backboneIndex].size(); vertebra++){
				voxelType xVox, yVox, zVox;
				xyz(segregatedBackbones[connectedComponent][backboneIndex][vertebra], xVox, yVox, zVox, voldimsGlobal);
				double rS = (vertebra < radSolid.size()) ? radSolid[vertebra] : 0.0;
				double rU = (vertebra < radSurf.size()) ? radSurf[vertebra] : 0.0;
				voxelOut << vertebra << "\t" << xVox << "\t" << yVox << "\t" << zVox << "\t" << rS << "\t" << rU << endl;
			}
			voxelOut << endl;
		}
	}
	voxelOut.close();
	
	string niceRuntime(niceTimeSince(analyzeTime));
	generalStatusUpdate("Output written to outputTSVfnBase = " + outputTSVfnBase + "...  [Runtime: " + niceRuntime + "]");    
	
}

void sphereCoarsenTest(){
	srandSet(112358);

//#if defined(__WXMSW__) || defined(__WINDOWS__)
//	string pathToImages("."); // Windows
//#elif defined(__WXOSX__)
	string pathToImages("C:/Users/kaiak/OneDrive/Documents/UCLA/SAVAGELAB/AngicartC++/AngicartCpp/sample_data");
//#else
//	string pathToImages(".");
//#endif
	
	double thresh(0.2), critFrac(0.16); // critFrac of 0.16 \approx 1/6 for 1 out of the six faces having a vessel neighbor
	//string imageDir(pathToImages + "/FITC-MCA0_N12_NIH001_s1_4"); int imStart(0), imEnd(66); double voxdims[] = {1.648, 1.648, 0.492};  thresh = 0.22; string lengthUnit("\u03BCm");
	string imageDir(pathToImages + "/FITC-MCA0_N12_PI001_s1"); int imStart(0), imEnd(61); double voxdims[] = {0.412, 0.412, 0.492}; thresh = 0.2; string lengthUnit("\u03BCm");
    //string imageDir(pathToImages); int imStart(0), imEnd(95); double voxdims[] = {195.3125, 195.3125, 200}; thresh = 0.85; string lengthUnit("\u03BCm");
	//string imageDir(pathToImages + "/Lung Images (Processed) - Study 2310 L2"); int imStart(1), imEnd(337); double voxdims[] = {1.0, 1.0, 1.0}; string lengthUnit("\u03BCm"); // um units not to scale
	
	string outputTSVfnBase("C:/Users/kaiak/OneDrive/Documents/UCLA/SAVAGELAB/AngicartC++/AngicartCpp/sample_outputs/FITC-MCA0_N12_PI001_s1");
	
	analyzeVascularStructure(imageDir, imStart, imEnd, outputTSVfnBase, voxdims, lengthUnit, thresh, critFrac);
}

#if TRY_WX != 1

int main(int argc, char *argv[]){
	startTime = getNow();
	srandSet(time(NULL));
	// default threads: hardware_concurrency - 1 (leave one core free), at least 1
	{
		unsigned int h(thread::hardware_concurrency());
		numLocThreads = (h <= 1u) ? 1u : (h - 1u);
	} // default max - 1 threads
	if(argc == 2){ // if 2 args are provided, use the second as the thread count
		int n(atoi(argv[1]));
		if(n >= 1) numLocThreads = (unsigned int)n;
	}
	if(argc == 1 || argc == 2){ // if no args or 2 args are provided, run the test
		generalStatusUpdate("\n Running sphereCoarsenTest()");
		sphereCoarsenTest();
		generalStatusUpdate("\n sphereCoarsenTest() complete.");
	}else if(argc == 10){ // keep 10 args option for reverse compatibility
		string imageDir(argv[1]);
		int imStart(atoi(argv[2])), imEnd(atoi(argv[3]));
		string outputTSVfnBase(argv[4]);
		double voxdims[] = {atof(argv[5]), atof(argv[6]), atof(argv[7])};
		string lengthUnit(argv[8]);
		double thresh(atof(argv[9])), critFrac(0.16);
		analyzeVascularStructure(imageDir, imStart, imEnd, outputTSVfnBase, voxdims, lengthUnit, thresh, critFrac);
	}else if(argc == 11){ // if 11 args are provided, use the 11th as the thread count
		int n(atoi(argv[10]));
		if(n >= 1) numLocThreads = (unsigned int)n;
		string imageDir(argv[1]);
		int imStart(atoi(argv[2])), imEnd(atoi(argv[3]));
		string outputTSVfnBase(argv[4]);
		double voxdims[] = {atof(argv[5]), atof(argv[6]), atof(argv[7])};
		string lengthUnit(argv[8]);
		double thresh(atof(argv[9])), critFrac(0.16);
		analyzeVascularStructure(imageDir, imStart, imEnd, outputTSVfnBase, voxdims, lengthUnit, thresh, critFrac);
	}else{
		cout << "\n\n Expected: no args or [thread_count] for sphereCoarsenTest(); or 10 args (image_directory image_start image_end output_base voxdim_x voxdim_y voxdim_z length_unit threshold) with optional 11th arg = thread_count."
			<< "\n\t string image_directory,\n\t int image_start_inclusive"
			<< "\n\t int image_end_inclusive"
			<< "\n\t string output_filename_base"
			<< "\n\t double voxel_dimension_x"
			<< "\n\t double voxel_dimension_y"
			<< "\n\t double voxel_dimension_z"
			<< "\n\t string length_unit_name"
			<< "\n\t double normalized_intensity_threshold"
			<< "\n\t [optional] int thread_count (default: hardware_concurrency - 1)"
			<< endl;
	}
	
	cout << "\n\n Runtime: " << niceTimeSince(startTime);
	cout << "\nProgram complete.\n";
	
	return 0;
}

#endif
