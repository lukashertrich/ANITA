#pragma once
#include <vector>
#include <fstream>
namespace anita{
    template <typename T>
    class Vector2{        
    public:
        const T x, y;
        Vector2(T X, T Y) : x(X), y(Y){}        

        Vector2<T> operator+ (const Vector2<T>& B) const {
            return Vector2(x + B.x, y + B.y);
        }

        Vector2<T> operator+= (const Vector2<T>& B) const {
            return Vector2(x + B.x, y + B.y);
        }

        Vector2<T> operator- (const Vector2<T>& B) const {
            return Vector2<T>(x - B.x, y - B.y);
        }

        Vector2<T> operator-= (const Vector2<T>& B) const {
            return Vector2<T>(x - B.x, y - B.y);
        }

        T operator* (const Vector2<T>& B) const {
            return x*B.x + y*B.y;
        }
        
        Vector2<T> operator* (const T s) const {
            return Vector2<T>(s*x, s*y);
        }

        Vector2<T> operator*= (const T s) const {
            return Vector2<T>(s*x, s*y);
        }

        Vector2<T> operator/ (const T s) const {
            return Vector2<T>(x/s, y/s);
        }

        Vector2<T> operator/= (const T s) const {
            return Vector2<T>(x/s, y/s);
        }

        T magSqr() const {
            return (x*x + y*y);
        }

        T mag() const {
            return sqrt(x*x + y*y);
        }

        Vector2<T> norm() const {
            T magnitude = sqrt(x*x + y*y);
            if(magnitude == (T)0){
                return this;
            }
            else{
                return Vector2(x/magnitude, y/magnitude);
            }
        }
    };
}
