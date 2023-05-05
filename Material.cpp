#ifndef MATERIAL
#define MATERIAL

#include <vector>

class Material {
    private:
        int energy_groups;
        int delayed_groups;
        std::vector<double> totalXS;
        std::vector<double> scatterXS;
        std::vector<std::vector<double>> scatterMatrixXS;
        std::vector<double> captureXS;
        std::vector<double> censusXS;
        std::vector<double> fissionXS;
        std::vector<double> velo;
        std::vector<double> nu;
        std::vector<double> beta;
        double totalBeta;
        std::vector<double> decayConst;
    public:
        Material(int e_groups, int del_groups, std::vector<double> totalXS, std::vector<double> scatterXS, std::vector<std::vector<double>> scatterMatrixXS, std::vector<double> captureXS, std::vector<double> censusXS,
                 std::vector<double> fissionXS, std::vector<double> velo, std::vector<double> nu, std::vector<double> beta, std::vector<double> decayConst):energy_groups(e_groups), delayed_groups(del_groups), totalXS(totalXS), 
            scatterXS(scatterXS), scatterMatrixXS(scatterMatrixXS), captureXS(captureXS), censusXS(censusXS), fissionXS(fissionXS), velo(velo), nu(nu), beta(beta), decayConst(decayConst){
            totalBeta = 0;
            for (double group_beta : beta) {
                totalBeta += group_beta;
            }
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
        double getVelocity(int group) {
            return velo[group];
        }
        double getNu(int group) {
            return nu[group];
        }
        double getBeta() {
            return totalBeta;
        }
        double getBeta(int group) {
            return beta[group];
        }
        double getDecayConst(int group) {
            return decayConst[group];
        }
        double getNumDelayedGroups() {
            return delayed_groups;
        }
        int sampleDelayedGroup(Rand &rng) {
            double ksi = rng.getRand();
            double accumulator = 0;
            for(int g = 0; g < delayed_groups; g++) {
                accumulator += beta[g] / totalBeta;
                if (accumulator > ksi) {
                    return g;
                }
            }
        }

};




#endif