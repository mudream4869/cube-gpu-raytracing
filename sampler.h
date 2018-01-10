#include <random>

int getRand(int* prng){
    (*prng) = (*prng)*8763;
    (*prng) %= 1000000007;
    (*prng) ^= 0xdeadbeef;
    if((*prng) < 0) (*prng) = -(*prng);
    return *prng;
}

// Return Random in [0, 1)
double sample(int* prng){
    int a = getRand(prng)%12345678;
    return a/12345678.;
}

void samplesSquare(int* prng, double* x, double* y, int n){
    double delta = 1./n;
    for(int lx = 0;lx < n;lx++){
        for(int ly = 0;ly < n;ly++){
            x[lx*n + ly] = delta*(lx + sample(prng));
            y[lx*n + ly] = delta*(ly + sample(prng));
        }
    }
    return;
}
