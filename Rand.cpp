#ifndef RAND
#define RAND

class Rand {
    public:

    static double getRand() {
        return ((double)rand())/(double)RAND_MAX;
    }
};

#endif