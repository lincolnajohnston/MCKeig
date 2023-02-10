#ifndef MATERIAL
#define MATERIAL

#include <vector>

class Material {
    private:
        double totalXS;
        double scatterXS;
        double captureXS;
        double censusXS;
        double fissionXS;
        double nu;
        double beta;
        double decayConst;
    public:
        Material(double totalXS, double scatterXS, double captureXS, double censusXS, double fissionXS, double nu, double beta, double decayConst):totalXS(totalXS), 
            scatterXS(scatterXS), captureXS(captureXS), censusXS(censusXS), fissionXS(fissionXS), nu(nu), beta(beta), decayConst(decayConst){

        }

        double getTotalXS() {
            return totalXS;
        }
        double getScatterXS() {
            return scatterXS;
        }
        double getCaptureXS() {
            return captureXS;
        }
        double getCensusXS() {
            return censusXS;
        }
        double getFissionXS() {
            return fissionXS;
        }
        double getNu() {
            return nu;
        }
        double getBeta() {
            return beta;
        }
        double getDecayConst() {
            return decayConst;
        }

};




#endif