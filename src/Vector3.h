#pragma once
#include <vector>
#include <fstream>
#include <cmath>
namespace anita{
    template <typename T>
    class Vector3{        
    public:
        const T x, y, z;
        Vector3(T x, T y, T z) : x(x), y(y), z(z){}        

        Vector3<T> operator+ (const Vector3<T>& B) const {
            return Vector3(x + B.x, y + B.y, z + B.z);
        }

        Vector3<T> operator+= (const Vector3<T>& B) const {
            return Vector3(x + B.x, y + B.y, z + B.z);
        }

        Vector3<T> operator- (const Vector3<T>& B) const {
            return Vector3<T>(x - B.x, y - B.y, z - B.z);
        }

        Vector3<T> operator-= (const Vector3<T>& B) const {
            return Vector3<T>(x - B.x, y - B.y, z - B.z);
        }

        T operator* (const Vector3<T>& B) const {
            return x*B.x + y*B.y + z*B.z;
        }
        
        Vector3<T> operator* (const T s) const {
            return Vector3<T>(s*x, s*y, s*z);
        }

        Vector3<T> operator*= (const T s) const {
            return Vector3<T>(s*x, s*y, s*z);
        }

        Vector3<T> operator/ (const T s) const {
            return Vector3<T>(x/s, y/s, z/s);
        }

        Vector3<T> operator/= (const T s) const {
            return Vector3<T>(x/s, y/s, z/s);
        }

        T magSqr() const {
            return (x*x + y*y + z*z);
        }

        T mag() const {
            return sqrt(x*x + y*y + z*z);
        }

        Vector3<T> norm() const {
            T magnitude = sqrt(x*x + y*y + z*z);
            if(magnitude == (T)0){
                return *this;
            }
            else{
                return Vector3(x/magnitude, y/magnitude, z/magnitude);
            }
        }

        Vector3<T> cross (const Vector3<T>& B) const {
            return Vector3<T>(
                y*B.z - z*B.y,
                z*B.x - x*B.z,
                x*B.y - y*B.x
            );
        }

        static Vector3<T> cross(const Vector3<T>& A, const Vector3<T>& B) {
            return Vector3<T>(
                A.y*B.z - A.z*B.y,
                A.z*B.x - A.x*B.z,
                A.x*B.y - A.y*B.x
            );
        }
        
    };
}
