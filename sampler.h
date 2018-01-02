#include <random>

void sampleSquare(double& x, double& y){
    std::random_device r;
    std::uniform_int_distribution<int> uniform_dist(1, 999);
    std::default_random_engine e1(r());

    int nx = uniform_dist(e1), ny = uniform_dist(e1);

    x = nx/1000., y = ny/1000.;
    return;
}

void sampleLine(double &x){
    std::random_device r;
    std::uniform_int_distribution<int> uniform_dist(1, 999);
    std::default_random_engine e1(r());

    int nx = uniform_dist(e1);

    x = nx/1000.;
}
