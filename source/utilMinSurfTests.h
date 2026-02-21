#ifndef _UTILMINSURFTESTS_H
#define _UTILMINSURFTESTS_H 1

#include <chrono>
#include <map>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

template <class T>
string makeString(const T x, streamsize prec = 8){
	stringstream strstr;
	strstr.precision(prec);
	strstr << x;
	return strstr.str();
}

template <class T>
string makeString(const vector<T> &x, streamsize prec = 8){
	if(x.size() == 0)
		return "()";
	stringstream strstr;
	strstr.precision(prec);
	strstr << "(" << x[0];
	for(unsigned int i(1); i < x.size(); i++)
		strstr << ", " << x[i];
	strstr << ")";
	return strstr.str();
}

string makeStringPos(const vector<unsigned int> &x, const vector<unsigned int> &y, const vector<unsigned int> &z){
	if(x.size() == 0)
		return "[]";
	stringstream strstr;
	strstr << "[(" << x[0] << "," << y[0] << "," << z[0] << ")";
	for(unsigned int i(1); i < x.size(); i++){
		strstr << "; (" << x[i] << "," << y[i] << "," << z[i] << ")";
	}
	strstr << "]";
	return strstr.str();
}

template<class T, class S>
string makeStringMap(const map<T, S> &m){
	string s("\n");
	for(typename map<T, S>::const_iterator it(m.begin()); it != m.end(); it++)
		s += makeString(it->first) + "\t" + makeString(it->second) + "\n";
	return s;
}

template <class T, class S>
string makeString(const pair<T, S> &p){
	return "(" + makeString(p.first) + ", " + makeString(p.second) + ")";
}

string makeStringOneLine(map<unsigned int, unsigned int> &m){
	stringstream strstr;
	for(map<unsigned int, unsigned int>::iterator it(m.begin()); it != m.end(); it++)
		strstr << makeString(it->first) << "(" << makeString(it->second) << ") ";
	return strstr.str();
}

template <class T>
string makePluralSuffix(T i){
	if(i == 1)
		return "";
	return "s";
}

string paddedInt(int x, int digits){
	string s(makeString(x));
	int thresh(10);
	for(int i(1); i < digits; i++){
		if(x < thresh)
			s = "0" + s;
		thresh *= 10;
	}
	return s;
}

double uniformRN(){
	return double(rand())/double(RAND_MAX);
}

template <class T>
void srandSet(T i){
	string srandSetMessage("\n Setting srand(" + makeString(i) + ").\n");
	cerr << srandSetMessage;
	cout << srandSetMessage; cout.flush();
	srand((unsigned int)i);
}

// returns index of first false; if not found, returns v.size()
unsigned int findFirstFalse(const vector<bool> &v){
	for(unsigned int i(0); i < v.size(); i++){
		if(!v[i])
			return i;
	}
	return (unsigned int)v.size();
}

// returns true if x is found in v
template <class T, class S>
bool isIn(const T &x, const vector<S> &v){
	for(unsigned int i(0); i < v.size(); i++){
		if(x == v[i])
			return true;
	}
	return false;
}

bool isIn(unsigned long long x, map<unsigned long long, vector<unsigned long long> > &m){
	for(map<unsigned long long, vector<unsigned long long> >::iterator it(m.begin()); it != m.end(); it++){
		for(unsigned int i(0); i < it->second.size(); i++){
			if(it->second[i] == x)
				return true;
		}
	}
	return false;
}

// returns true if v is all true
bool allTrue(const vector<bool> &v){
	return findFirstFalse(v) == v.size();
}

// returns true if v is all false
bool allFalse(const vector<bool> &v){
	for(unsigned int i(0); i < v.size(); i++){
		if(v[i])
			return false;
	}
	return true;
}

// returns true if added
template <class T>
bool pushUnique(T x, vector<T> &v){
	for(unsigned int i(0); i < v.size(); i++){
		if(x == v[i])
			return false;
	}
	v.push_back(x);
	return true;
}

// returns true if something from x added to v
template <class T>
bool pushUniqueSet(const vector<T> &x, vector<T> &v){
	bool pushed(false);
	for(unsigned int i(0); i < x.size(); i++)
		pushed = pushUnique(x[i], v) || pushed;
	return pushed;
}

// removes the first x from v that it finds, if any
template <class T, class S>
void removeFrom(T x, vector<S> &v){
	for(unsigned int i(0); i < v.size(); i++){
		if(x == v[i]){
			v.erase(v.begin() + i);
			return;
		}
	}
}

string niceTime(double t, bool roundSeconds){
	if(t >= 0.0){
		stringstream nT;
		int m(60);
		int h(60*m);
		int d(24*h);
		int w(7*d);
		int wt((int)floor(t/w));
		t -= w*wt;
		int dt((int)floor(t/d));
		t -= d*dt;
		int ht((int)floor(t/h));
		t -= h*ht;
		int mt((int)floor(t/m));
		t -= m*mt;
		if(wt > 0)
			nT << wt << "w ";
		if(wt > 0 || dt > 0)
			nT << dt << "d ";
		if(wt > 0 || dt > 0 || ht > 0)
			nT << ht << "h ";
		if(wt > 0 || dt > 0 || ht > 0 || mt > 0)
			nT << mt << "m ";
		if(roundSeconds)
			nT << ceil(t) << "s";
		else
			nT << t << "s";
		return nT.str();
	} else return "-(" + niceTime(-t, roundSeconds) + ")";
}

string niceTimeSince(time_t firstTime){
	return niceTime(difftime(time(NULL), firstTime), false);
}

inline chrono::steady_clock::time_point getNow(){ return chrono::steady_clock::now(); }

double getDiff(const chrono::steady_clock::time_point &t2, const chrono::steady_clock::time_point &t1){
	return (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1)).count();
}

double getDiffSince(const chrono::steady_clock::time_point &t1){
	return getDiff(getNow(), t1);
}

string niceTime(const chrono::steady_clock::time_point &endTime, const chrono::steady_clock::time_point &startTime, bool roundSeconds = false){
	return niceTime(getDiff(endTime, startTime), roundSeconds);
}

string niceTimeSince(const chrono::steady_clock::time_point &startTime, bool roundSeconds = false){
	return niceTime(getDiff(getNow(), startTime), roundSeconds);
}

// cycles x into the interval [0, 1]
template <class T>
double colorProfile(T x){
	x = 6*(x - floor(x));
	
	if(x < 1 || x >= 5)
		return 1;
	if(x < 2)
		return 2.0 - x;
	if(x < 4)
		return 0;
	// if(x < 5)
		return x - 4.0;
		
}

template <class T, class S>
void rainbowRGB(T x, S &r, S &g, S&b){
	r = colorProfile(x);
	g = colorProfile(x + 2.0/3.0);
	b = colorProfile(x + 4.0/3.0);
}

void rainbowColor(unsigned int i, unsigned int iMax, unsigned char &r, unsigned char &g, unsigned char &b){
	double R, G, B;
	rainbowRGB(double(i)/iMax, R, G, B);
	r = (unsigned char)floor(127*R);
	g = (unsigned char)floor(127*G);
	b = (unsigned char)floor(127*B);
}

#endif
