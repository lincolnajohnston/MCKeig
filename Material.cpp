#ifndef MATERIAL
#define MATERIAL

#include <vector>

class Material {
    private:
        int groups;
        std::vector<double> totalXS;
        std::vector<double> scatterXS;
        std::vector<std::vector<double>> scatterMatrixXS;
        std::vector<double> captureXS;
        std::vector<double> censusXS;
        std::vector<double> fissionXS;
        std::vector<double> nu;
        std::vector<double> beta;
        std::vector<double> decayConst;
    public:
        Material(int groups, std::vector<double> totalXS, std::vector<double> scatterXS, std::vector<std::vector<double>> scatterMatrixXS, std::vector<double> captureXS, std::vector<double> censusXS,
                 std::vector<double> fissionXS, std::vector<double> nu, std::vector<double> beta, std::vector<double> decayConst):groups(groups), totalXS(totalXS), 
            scatterXS(scatterXS), scatterMatrixXS(scatterMatrixXS), captureXS(captureXS), censusXS(censusXS), fissionXS(fissionXS), nu(nu), beta(beta), decayConst(decayConst){

        }

        double getTotalXS(int group) {
            return totalXS[group];
        }
        double getScatterXS(int group) {
            return scatterXS[group];
        }
        std::vector<double> getScatterMatrixXS(int group) {
            return scatterMatrixXS[group];
        }
        double getCaptureXS(int group) {
            return captureXS[group];
        }
        double getCensusXS(int group) {
            return censusXS[group];
        }
        double getFissionXS(int group) {
            return fissionXS[group];
        }
        double getNu(int group) {
            return nu[group];
        }
        double getBeta(int group) {
            return beta[group];
        }
        double getDecayConst(int group) {
            return decayConst[group];
        }

};




#endif