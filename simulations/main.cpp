//
//

#include <iostream>
#include <unordered_map>
#include <vector>
#include <map>
#include <unistd.h>
#include <sys/wait.h>

// settings for the experiment
#define BIT_LEN 2
#define HIST_SIZE 4 // make sure HIST_SIZE = 2^BIT_LEN

typedef unsigned long long ull;
typedef unsigned long ulong;
typedef unsigned short int ushort;


struct PMD_Histogram {
    // counts historgram of bitvector distributions
    ushort histogram[HIST_SIZE];
    
    PMD_Histogram() {memset(histogram, 0, sizeof(PMD_Histogram));}
    
    void incBitVectorCount(ulong binNumber) {
        assert(binNumber < HIST_SIZE);
        histogram[binNumber] ++;
    }

    inline const ushort& operator [](int index) const {return histogram[index];}
    
    void print() const {
        for (ulong i = 0; i < HIST_SIZE; i++) {
            std::cout << histogram[i] << ".";
        }
    }
    
    ull* getUll() const {return (ull*) &histogram;}
};

bool operator==(const PMD_Histogram& lhs, const PMD_Histogram& rhs) {
    switch (HIST_SIZE) {
        case 4:
            return *(lhs.getUll()) == *(rhs.getUll());
    }

    for (ulong i = 0; i < HIST_SIZE; i++) {
        if (lhs.histogram[i] != rhs.histogram[i]) return false;
    }
    return true;
}

bool operator< ( const PMD_Histogram& lhs, const PMD_Histogram& rhs ) {
    
    switch (HIST_SIZE) {
        case 4:
            return *(lhs.getUll()) < *(rhs.getUll());
    }
    
    for (ulong i = 0; i < HIST_SIZE; i++) {
        if (lhs.histogram[i] < rhs.histogram[i]) return true;
    }
    return false;
}


struct PMD_Histogram_Hash {
    
    std::hash<ushort> hasher;
    std::hash<ull> ullHasher;
    std::size_t operator()(PMD_Histogram const& histogram) const noexcept
    {
        switch (HIST_SIZE) {
            case 4:
                return ullHasher(*histogram.getUll());
        }
        std::size_t hashValue = 0;
        for (int i = 0; i < HIST_SIZE; i++)
        {
            hashValue ^= hasher(histogram[i]);
        }
        return hashValue;
    }
};


struct PMD_Probs {
    
    double qProb;
    double binProbs[HIST_SIZE][HIST_SIZE];
  
    const double* getBinProbs(ulong binNumber) const { return binProbs[binNumber];}
    
    static ushort countSetBits(ulong number) {
        ulong count = 0;
        while (number > 0) {
            // Brian Kernighan’s Algorithm to count set bits in a number
            number = number & (number - 1);
            count++;
        }
        return count;
    }
    
    // interface to probability vector for each hisgtrogram bin
    static void binProbabilities(ulong binNumber, double q, double* probs)
    {
        // clean the probs
        memset(probs, 0, sizeof(double) * HIST_SIZE);
        // compute p
        double p = 1 - q;
        // generate probs
        for (ulong i = 0; i < HIST_SIZE; i++) {
            // compute commong bits between i and binVector
            ulong diffeferentBits = countSetBits(binNumber ^ i);
            probs[i] = pow(q,diffeferentBits) * pow(p,BIT_LEN-diffeferentBits);
        }
    }
    
    PMD_Probs(double q) {
        qProb = q;
        for (ulong i = 0; i < HIST_SIZE; i++) {
            binProbabilities(i, q, binProbs[i]);
        }
    }

};

typedef std::unordered_map<PMD_Histogram, double, PMD_Histogram_Hash> PMD_Map;
typedef PMD_Map::iterator PMD_Iterator;
typedef std::map<PMD_Histogram, double> PMD_Order_MAP;
struct PMD_Deltas {
    double unitOverZero;
    double zeroOverUnit;
    PMD_Deltas():unitOverZero(0),zeroOverUnit(0) { }
    PMD_Deltas(double uoz, double zou):unitOverZero(uoz),zeroOverUnit(zou) { }
};

struct PMD_Distro {
    PMD_Map *theMap;
    ulong originalHistogram[HIST_SIZE];
    
    ~PMD_Distro() {delete theMap;}

    PMD_Distro(const PMD_Probs& probs, ulong* binCounts) {
        memcpy(originalHistogram, binCounts, sizeof(originalHistogram));
        theMap = NULL;
        buildMapFromScrantch(probs);
    }
    
    PMD_Distro(PMD_Distro& distro, const PMD_Probs& probs, ulong binNumber) {
        memcpy(originalHistogram, distro.originalHistogram, sizeof(originalHistogram));
        originalHistogram[binNumber] ++;
        theMap = distro.addProbabilityVector(probs.getBinProbs(binNumber));
    }

    void buildMapFromScrantch(const PMD_Probs& probs) {
        if (theMap) delete theMap;
        
        theMap = NULL;
        // add vectors in every bin one by one
        for (ulong i = 0; i < HIST_SIZE; i++) {
            // add each bin vector
            for (ulong j = 0; j < originalHistogram[i]; j++) {
                if (theMap == NULL) {
                    theMap = new PMD_Map();
                    initPMDMap(*theMap, probs.getBinProbs(i));
                }
                else {
                    PMD_Map * newMap = addProbabilityVector(probs.getBinProbs(i));
                    delete theMap;
                    theMap = newMap;
                }
            }
        }
    }
    
    // initliazes an empty PMD map with a single porbability vector
    void initPMDMap(PMD_Map &map, const double* porbabilites) {
        for (int i=0; i < HIST_SIZE; i++) {
            PMD_Histogram hist;
            hist.incBitVectorCount(i);
            map.insert({hist, porbabilites[i]});
        }
    }
    
    // add a probability vector to the distribution
    PMD_Map * addProbabilityVector(const double porbabilites[HIST_SIZE]) {
        PMD_Map * newMap = new PMD_Map();
        // walk over the old map and for every histogram in it compute convolution with the probability vector
        // and reinsert new histograms into a newMap
        for (PMD_Iterator itr = theMap->begin(); itr != theMap->end();  itr++) {
            double histProbability = itr->second;
            for (int i = 0; i < HIST_SIZE; i++)
            {
                // make a copy of the original histogram
                PMD_Histogram hist = itr->first;
                hist.incBitVectorCount(i);
                double newProb = histProbability * porbabilites[i];
                // we now need to add the newHistogram and its probability into the newMap
                PMD_Iterator itemItr = newMap->find(hist);
                if (itemItr == newMap->end()) {
                    // we did not find that histogram, so insert it
                    newMap->insert({hist, newProb});
                }
                else {
                    // otherwise add a new probability to the existing historam
                    itemItr->second += newProb;
                }
            } // ends bit vector probability walk
        } // ends original map walk
        // at this point newMap contains all the new probability for all the extended histograms
        // so swap the maps
        return newMap;
    }
    
    void computeDelta(const PMD_Probs& probs, double lambda, PMD_Deltas& deltas)
    {
        PMD_Distro unitDistro(*this, probs, HIST_SIZE-1);
        PMD_Distro zeroDistro(*this, probs, 0);
        
        // compute both probabilites
        ulong unitOverZeroBinCount = 0;
        double unitOverZeroProbs = 0;
        ulong zeroOverUnitBinCount = 0;
        double zeroOverUnitProbs = 0;
        bool printHistograms = false;
        
        PMD_Order_MAP orderedMap(unitDistro.theMap->begin(), unitDistro.theMap->end());
        for (PMD_Order_MAP::iterator unitItr =orderedMap.begin(); unitItr != orderedMap.end();  unitItr++) {
            PMD_Iterator zeroItr = zeroDistro.theMap->find(unitItr->first);
            assert(zeroItr != zeroDistro.theMap->end());
            if (printHistograms) unitItr->first.print();
            // now we have two iterators unitItr pinting to unitDistro historgram and zeroItr pointing to zeroDisto histogram
            // let's compare thier probabilities
            if (unitItr->second > lambda * zeroItr->second) {
                // we found a histogram which unit prob is above lambda * zero_prob
                unitOverZeroBinCount ++;
                unitOverZeroProbs += unitItr->second;
                if (printHistograms) std::cout << " UOZ " << unitItr->second / zeroItr->second;
            }
            if (zeroItr->second > lambda * unitItr->second) {
                // same thing but on the zero end
                zeroOverUnitBinCount ++;
                zeroOverUnitProbs += zeroItr->second;
                if (printHistograms) std::cout << " ZOU " << zeroItr->second / unitItr->second ;
            }
            if (printHistograms) std::cout << "\n";
        }
        if (printHistograms) {
            std::cout << "\n\n";
            std::cout << "Unit Over Zero " << unitOverZeroBinCount << " " << unitOverZeroProbs << "\n";
            std::cout << "Zero Over Unit " << zeroOverUnitBinCount << " " << zeroOverUnitProbs << "\n";
        }
        deltas.unitOverZero = unitOverZeroProbs;
        deltas.zeroOverUnit = zeroOverUnitProbs;
    }
    
    void print() {
        PMD_Order_MAP orderedMap(theMap->begin(), theMap->end());
        for (PMD_Order_MAP::iterator itr = orderedMap.begin(); itr != orderedMap.end();  itr++) {
            itr->first.print();
            std::cout << " " << *(itr->first.getUll()) << " " << itr->second << "\n";
        }
    }
};

void runTests(double q, double lambda, ulong multiple) {
    ulong N = 20;
    ulong binCounts[][HIST_SIZE] = {
      {20,0,0,0},
       {15,5,0,0},
        {15,0,5,0},
        {15,0,0,5},
        {10,5,5,0},
        {10,0,5,5},
        {10,0,0,10},
        {10,10,0,0},
        {0,0,10,10},
        {5,5,5,5},
        {5,0,5,10},
        {5,5,0,10},
        {0,5,5,10},
        {0,0,0,20},
    };

    
    PMD_Probs probs(q);
    ulong total = sizeof(binCounts) / (sizeof(ulong) * HIST_SIZE);
    ulong binTotalBits[HIST_SIZE];
    for (ulong i = 0; i < HIST_SIZE; i++) {
        binTotalBits[i] = PMD_Probs::countSetBits(i);
    }
    
    // multiple volume
    N *= multiple;
    for (ulong i = 0; i < total; i++) {
        for (ulong j = 0; j < HIST_SIZE; j++) {
            binCounts[i][j] *= multiple;
        }
    }
    
    for (ulong i = 0; i < total; i++) {
        std::cout << N << " " << q << " " << lambda << " ";
        PMD_Distro distro(probs, binCounts[i]);
        PMD_Deltas deltas;
        distro.computeDelta(probs,2, deltas);
        double expectations[HIST_SIZE];
        memset(expectations, 0, sizeof(expectations));
        for (ulong j = 0; j < HIST_SIZE; j++) {
            ulong binCount = binCounts[i][j];
            std::cout << binCount << "|";
            const double* j_probs =  probs.getBinProbs(j);
            // accumulate expectations
            for (ulong k=0; k < HIST_SIZE; k++) {
                expectations[k] += binCount * j_probs[k];
            }
        }
        std::cout << " ";
        
        double unitDist = 0;
        double zeroDist = 0;
        const double * unitProbs = probs.getBinProbs(HIST_SIZE-1);
        const double * zeroProbs = probs.getBinProbs(0);
        for (ulong k=0; k < HIST_SIZE; k++) {
            std::cout << expectations[k] << "|";
            unitDist += pow(expectations[k] / N - unitProbs[k],2);
            zeroDist += pow(expectations[k] / N - zeroProbs[k],2);
        }
        std::cout << " " << sqrt(zeroDist) << " " << sqrt(unitDist)  << " " << deltas.zeroOverUnit << " " <<  deltas.unitOverZero  << '\n' << std::flush;
;
    }
}

int main(int argc, const char * argv[]) {
    
    /*ulong binCounts[HIST_SIZE] = {0,0,0,10};
    double q = 0.2;
    PMD_Probs probs(q);
    PMD_Distro distro(probs, binCounts);
    PMD_Deltas deltas;
    distro.computeDelta(probs,2, deltas);
    std::cout << deltas.unitOverZero << " " <<  deltas.zeroOverUnit << '\n';
    */
    for (double qp = 0.1; qp <= 0.4; qp += 0.1 ) {
        for (double lambda = 2; lambda <= 4; lambda += 1 ) {
            pid_t childPid = fork();
            if (childPid == 0) {
                // i am a child run the thing
                std::cerr << "running child for " << " " << qp << " " << lambda << '\n' << std::flush;
                for (ulong multiple = 1; multiple < 9; multiple += 2) {
                    runTests(qp, lambda, multiple);
                    std::cerr << "done with " << 20*multiple << " " << qp << " " << lambda << '\n' << std::flush;
                }
                exit(0);
            }
        }
        // ok i am a parent keep making kids
    }
    // and now wait for all of them
    int stat_loc;
    pid_t deadChild;
    while ( (deadChild=wait(&stat_loc)) != -1) {
        std::cerr << deadChild << " is dead" << "\n";
    }
    return 0;
}
