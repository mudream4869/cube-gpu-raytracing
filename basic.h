#include <iostream>
#include <limits>

double clamp_(double x, double l, double u){
    return (x < l) ? l : ((x > u) ? u : x);
}

// Return True if |x| < eps
bool isZero(double x){
    return std::abs(x) < std::numeric_limits<double>::min();
}

template<typename T>
struct Vec3 {
    T x,y,z;
    Vec3<T>(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}

    Vec3<T> operator + (const Vec3<T>& v) const { return Vec3<T>(x+v.x, y+v.y, z+v.z); }
    Vec3<T> operator - (const Vec3<T>& v) const { return Vec3<T>(x-v.x, y-v.y, z-v.z); }
    Vec3<T> operator * (const Vec3<T>& v) const { return Vec3<T>(x*v.x, y*v.y, z*v.z); }

    Vec3<T> operator * (T d) const { return Vec3<T>(x*d, y*d, z*d); }
    Vec3<T> operator / (T d) const { return Vec3<T>(x/d, y/d, z/d); }

    bool operator < (const Vec3<T>& v) const {
        return std::tie(x, y, z) < std::tie(v.x, v.y, v.z);
    }

    double len() const {
        return sqrt(x*x + y*y + z*z);
    }

    Vec3<T> normalize() const {
        double mg = sqrt(x*x + y*y + z*z);
        return Vec3<T>(x/mg,y/mg,z/mg);
    }

    void clamp(double lower_bound, double upper_bound) {
        x = clamp_(x, lower_bound, upper_bound); 
        y = clamp_(y, lower_bound, upper_bound); 
        z = clamp_(z, lower_bound, upper_bound); 
        return;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const Vec3<T>& v){
    o << "(" << v.x << "," << v.y << "," << v.z << ")";
    return o;
}

template<typename T>
T dot(const Vec3<T>& a, const Vec3<T>& b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;

struct Ray {
    Vec3d o,d;
    double t;
    Ray(const Vec3d& o, const Vec3d& d, double t = 0) : o(o), d(d.normalize()), t(t) {}

    Vec3d operator()(double _t) const{
        return o + d*_t;
    }
};
