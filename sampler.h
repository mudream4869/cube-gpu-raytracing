#include <random>

void sampleSquare(double& x, double& y){
    std::random_device r;
    std::uniform_int_distribution<int> uniform_dist(0, 1000);
    std::default_random_engine e1(r());

    int nx = uniform_dist(e1), ny = uniform_dist(e1);

    x = nx/1000., y = ny/1000.;
    return;
}


// |Sample| = n*n
void samplesSquare(double* x, double* y, int n){
    std::random_device r;
    std::uniform_int_distribution<int> uniform_dist(0, 1000);
    std::default_random_engine e1(r());
    
    double delta = 1./n;
     
    for(int lx = 0;lx < n;lx++){
        for(int ly = 0;ly < n;ly++){
            int nx = uniform_dist(e1), ny = uniform_dist(e1);
            x[lx*n + ly] = delta*lx + nx/1000./n;
            y[lx*n + ly] = delta*ly + ny/1000./n;
        }
    }

    return;
}

void sampleLine(double &x){
    std::random_device r;
    std::uniform_int_distribution<int> uniform_dist(0, 1000);
    std::default_random_engine e1(r());

    int nx = uniform_dist(e1);

    x = nx/1000.;
}
