//
//  Vec.h
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH_Vec_h
#define SPH_Vec_h


#include<cmath>
#include<string>
#include<iostream>

template<int N>
class Vec{
    
private:
    double _elems[N];
    int _transpose;
    
public:
    Vec();
    Vec(double a);          //Initialize the vector with value 'a'
    Vec(double, double);
    Vec(double, double, double);
    Vec( const Vec&V);
    ~Vec(){};
    int length();
    double* elem_return();
    
//    Vec& operator=( Vec&);
    Vec& operator=( Vec);
    Vec& operator=( double&);
    void transpose(){_transpose *=-1;}
    int return_trans(){return _transpose;}
    double& operator[](int i) {return _elems[i];}
    double& operator()(int i) {return _elems[i];}
    void set(int i, double a) {_elems[i] = a;}
    Vec& operator+=(Vec&);
    Vec& operator-=(Vec&);
    Vec& operator*=(double& a);
    Vec& operator/=(double& a);
    double l2norm(){
        double sum = 0.0;
        for(int i =0 ; i < N; i++){
            sum += _elems[i];
        }
        return sqrt(sum);
    }
    double squaredNorm(){
        double sum = 0.0;
        for(int i =0 ; i < N; i++){
            sum += _elems[i];
        }
        return sum;
    }
    friend std::ostream& operator<<(std::ostream& output, Vec& A);
    
};

template<int N>
Vec<N>::Vec( ){
    _transpose = -1;
    for(int i = 0; i < N; i++){
        _elems[i]  = 0.0;
    }
}

template<int N>
Vec<N>::Vec(double a){
    _transpose = -1;
    for(int i = 0; i < N; i++){
        _elems[i]  = a;
    }
}

template<int N>
Vec<N>::Vec(double a, double b){
    _transpose = -1;
    
        _elems[0]  = a;
    _elems[1]  = b;
}

template<int N>
Vec<N>::Vec(double a, double b, double c){
    _transpose = -1;
    
    _elems[0]  = a;
    _elems[1]  = b;
    _elems[2] = c;
}

template<int N>
Vec<N>::Vec(const Vec<N> &V){
    _transpose = -1;
    for(int i = 0; i < N; i++){
        _elems[i]  = V._elems[i];
    }
}



//template<int N>
//Vec<N>& Vec<N>::operator =( Vec<N>& V){
//    if( this != &V)
//        for(int i = 0; i < N; i++){
//            _elems[i] = V._elems[i];
//        }
//    return *this;
//}

template<int N>
Vec<N>& Vec<N>::operator =( Vec<N> V){
    if( this != &V)
        for(int i = 0; i < N; i++){
            _elems[i] = V._elems[i];
        }
    return *this;
}

template<int N>
Vec<N>& Vec<N>::operator =( double& a){
    for(int i = 0; i < N; i++){
        _elems[i] = a;
    }
    return *this;
}

template<int N>
Vec<N>& Vec<N>::operator +=( Vec<N>& V){
    
    for(int i = 0; i < N; i++){
        _elems[i] += V[i];
    }
    return *this;
}

template<int N>
Vec<N>& Vec<N>::operator -=( Vec<N>& V){
    
    for(int i = 0; i < N; i++){
        _elems[i] -= V[i];
    }
    return *this;
}

template<int N>
Vec<N>& Vec<N>::operator *=( double& a){
    
    for(int i = 0; i < N; i++){
        _elems[i] *= a;
    }
    return *this;
}

template<int N>
Vec<N>& Vec<N>::operator /=( double& a){
    double b;
    if(a == 0){
        b = 1.0;
    }else{
        b = a;
    }
    for(int i = 0; i < N; i++){
        _elems[i] /= b;
    }
    return *this;
}

template<int N>
Vec<N> operator+( Vec<N>&u, Vec<N>&v){
    return Vec<N>(u) +=v;
}

template<int N>
Vec<N> operator-( Vec<N>&u, Vec<N>&v){
    return Vec<N>(u) -=v;
}

template<int N>
Vec<N> operator*( Vec<N>&u, double&a){
    return Vec<N>(u) *=a;
}

template<int N>
Vec<N> operator*(double&a ,Vec<N>&u){
    return Vec<N>(u) *=a;
}

template<int N>
Vec<N> operator/(double&a ,Vec<N>&u){
    return Vec<N>(u) /=a;
}

template<int N>
Vec<N> operator/( Vec<N>&u, double&a){
    return Vec<N>(u) /=a;
}

template<int N>
Vec<N> operator+(Vec<N>&u){
    return u;
}
template<int N>
Vec<N> operator-(Vec<N>&u){
    return Vec<N>(u) *=-1;
}

template<int N>
double operator*(Vec<N>&u, Vec<N>&v){
    double sum = 0.0;
    if(u.return_trans() == 1 && v.return_trans() == -1){
        for(int i = 0; i < N; i++){
            sum += u[i]*v[i];
        }
    }
    return sum;
}

template<int N>
std::ostream& operator << (std::ostream& output, Vec<N>& A){
    output << "[ ";
    for(int i = 0; i < N; i++){
        output << A[i] << std::endl;
    }
    output << "]" << std::endl;
    return output;
}

#endif
