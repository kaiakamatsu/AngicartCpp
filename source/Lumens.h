#ifndef LUMENS_H
#define LUMENS_H 1

// note from Kai: I think this might be easier implemented as a 1D array of doubles, but I inherited the code like this so let's use as is. No need to fix if it is not broken. 

struct Lumens{ // really this is an alternative to a vector of vector of vectors
	unsigned int size[3]; //stores three dimensions
	double ***lumens; // triple pointer, think of an array in an array in an array

    // default constructor and destructor, not used, just an example
	Lumens(){ size[0] = size[1] = size[2] = 0; lumens = new double**[1]; lumens[0] = new double*[1]; lumens[0][0] = new double[1]; }
	~Lumens(){
		for(unsigned int i(0); i < size[0]; i++){
			for(unsigned int j(0); j < size[1]; j++)
				delete[] lumens[i][j];
			delete[] lumens[i];
		}
        // outermost array I think should be deleted here but it is missing.
        //delete[] lumens; // WARNING: I added this, return here if we see any inconsistencies.
		size[0] = size[1] = size[2] = 0;
	}

    // main constructor, manually initializing the 3D array
	Lumens(unsigned int x, unsigned int y, unsigned int z){
		size[0] = x; size[1] = y, size[2] = z;
		lumens = new double**[x]; // array of x pointers
		for(unsigned int i(0); i < x; i++){
			lumens[i] = new double*[y]; // for each x pointer, we create an array of y pointers
			for(unsigned int j(0); j < y; j++){
				lumens[i][j] = new double[z]; // for each x pointer y pointer we allocate z doubles
				for(unsigned int k(0); k < z; k++) // for each x pointer y pointer z double we set to 0.0
					lumens[i][j][k] = 0.0;
			}
		}
	}

    // copy constructor, deep copy of the 3D array
	Lumens(const Lumens &L){
		unsigned int x(L.size[0]), y(L.size[1]), z(L.size[2]);
		size[0] = x; size[1] = y, size[2] = z;
		lumens = new double**[x];
		for(unsigned int i(0); i < x; i++){
			lumens[i] = new double*[y];
			for(unsigned int j(0); j < y; j++){
				lumens[i][j] = new double[z];
				for(unsigned int k(0); k < z; k++)
					lumens[i][j][k] = L.lumens[i][j][k];
			}
		}
	}

    // assignment operator, but we should revisit this to handle self-assignment and memory management
    // simply perform copy constructor of input and swap (which should handle self-assignment and memory management properly)
	Lumens& operator=(const Lumens &L){
		unsigned int x(L.size[0]), y(L.size[1]), z(L.size[2]);
		size[0] = x; size[1] = y, size[2] = z;
		lumens = new double**[x];
		for(unsigned int i(0); i < x; i++){
			lumens[i] = new double*[y];
			for(unsigned int j(0); j < y; j++){
				lumens[i][j] = new double[z];
				for(unsigned int k(0); k < z; k++)
					lumens[i][j][k] = L.lumens[i][j][k];
			}
		}
		return *this;
	}
	
	unsigned int totalSize() const ;
	unsigned int indexOf(unsigned int x, unsigned int y, unsigned int z) const ;
	unsigned int x(unsigned int i) const ;
	unsigned int y(unsigned int i) const ;
	unsigned int z(unsigned int i) const ;
	double minLumen() const ;
	double maxLumen() const ;
};

void normalizeLumens(Lumens &L);
Lumens simpleLumensCube(unsigned int cubeSide = 2);

unsigned int Lumens::totalSize() const {
    return size[0]*size[1]*size[2];
}

unsigned int Lumens::indexOf(unsigned int x, unsigned int y, unsigned int z) const {
    return x*size[1]*size[2] + y*size[2] + z;
}
unsigned int Lumens::x(unsigned int i) const {
    return i/(size[1]*size[2]);
}
unsigned int Lumens::y(unsigned int i) const {
    return (i%(size[1]*size[2]))/size[2];
}
unsigned int Lumens::z(unsigned int i) const {
    return i%size[2];
}
double Lumens::minLumen() const {
    double m(lumens[0][0][0]);
    for(unsigned int i(0); i < size[0]; i++){
        for(unsigned int j(0); j < size[1]; j++){
            for(unsigned int k(0); k < size[2]; k++){
                if(m > lumens[i][j][k])
                    m = lumens[i][j][k];
            }
        }
    }
    return m;
}
double Lumens::maxLumen() const {
    double m(lumens[0][0][0]);
    for(unsigned int i(0); i < size[0]; i++){
        for(unsigned int j(0); j < size[1]; j++){
            for(unsigned int k(0); k < size[2]; k++){
                if(m < lumens[i][j][k])
                    m = lumens[i][j][k];
            }
        }
    }
    return m;
}

void normalizeLumens(Lumens &L){
    if(L.size[0] < 1 || L.size[1] < 1 || L.size[2] < 1)
        return;
    double maxL(L.lumens[0][0][0]), minL(L.lumens[0][0][0]);
    for(unsigned int x(0); x < L.size[0]; x++){
        for(unsigned int y(0); y < L.size[1]; y++){
            for(unsigned int z(0); z < L.size[2]; z++){
                if(maxL < L.lumens[x][y][z])
                    maxL = L.lumens[x][y][z];
                if(minL > L.lumens[x][y][z])
                    minL = L.lumens[x][y][z];
            }
        }
    }
    if(maxL == minL){
        for(unsigned int x(0); x < L.size[0]; x++){
            for(unsigned int y(0); y < L.size[1]; y++){
                for(unsigned int z(0); z < L.size[2]; z++)
                    L.lumens[x][y][z] = 0.0;
            }
        }
        return;
    }
    for(unsigned int x(0); x < L.size[0]; x++){
        for(unsigned int y(0); y < L.size[1]; y++){
            for(unsigned int z(0); z < L.size[2]; z++)
                L.lumens[x][y][z] = (L.lumens[x][y][z] - minL)/(maxL - minL);
        }
    }
}

Lumens simpleLumensCube(unsigned int cubeSide){
    Lumens L(cubeSide, cubeSide, cubeSide);
    for(unsigned int x(0); x < L.size[0]; x++){
        for(unsigned int y(0); y < L.size[1]; y++){
            for(unsigned int z(0); z < L.size[2]; z++)
                L.lumens[x][y][z] = 1.0;
        }
    }
    L.lumens[0][0][0] = 0.0;
    return L;
}


#endif
