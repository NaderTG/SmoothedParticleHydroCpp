//
//  Vec.h
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_Vec_h
#define SPH3_Vec_h


#include <cmath>
#include<string>
#include<iostream>

template< int N> class Vec{
    double _elems[N];
public:
    Vec();
    Vec(const double& );
    Vec(const double& a,const double& b){
        _elems[0] = a; _elems[1] = b;
    }  // constructor for 2-d Vecs
    Vec(const double& a,const double& b,const double& c){
        _elems[0] = a; _elems[1] = b; _elems[2] = c;
    }  // constructor for 3-d Vecs
    Vec(const Vec&);
    const Vec& operator=(const Vec&);
    const Vec& operator=(const double& );
    ~Vec(){}  //  destructor
    double&  operator[](int i) { return _elems[i]; }  //ith component
    const double&  operator()(int i) const{ return _elems[i]; }  //ith component
    void set(int i,const double&  a){ _elems[i] = a; }  //  change ith component
    const Vec& operator+=(const Vec&);
    const Vec& operator-=(const Vec&);
    const Vec& operator*=(const double& );
    const Vec& operator/=(const double& );
    
    double l2norm(){
        double sum = 0.0;
        for(int i =0 ; i < N; i++){
            sum += _elems[i]*_elems[i];
        }
        return sqrt(sum);
    }
    double squaredNorm(){
        double sum = 0.0;
        for(int i =0 ; i < N; i++){
            sum += _elems[i]*_elems[i];
        }
        return sum;
    }
    
    friend std::ostream& operator<<(std::ostream& output, Vec& A){
        output << "\n[\n";
        for(int i = 0; i < N; i++){
            output << A[i] <<", "<< std::endl;
        }
        output << "]\n" << std::endl;
        return output;
    }
    
};

template<int N>
Vec<N>::Vec(){
    for(int i = 0; i < N; i++)
        _elems[i] = 0.0;
}  //  constructor


template<int N>
Vec<N>::Vec(const double&  a ){
    for(int i = 0; i < N; i++)
        _elems[i] = a;
}  //  constructor

template<int N>
Vec<N>::Vec(const Vec<N>& v){
    for(int i = 0; i < N; i++)
        _elems[i] = v._elems[i];
}  //  copy constructor

template<int N>
const Vec<N>& Vec<N>::operator=(const Vec<N>& v){
    if(this != &v)
        for(int i = 0; i < N; i++)
            _elems[i] = v._elems[i];
    return *this;
}  //  assignment operator

template<int N>
const Vec<N>& Vec<N>::operator=(const double&  a){
    for(int i = 0; i < N; i++)
        _elems[i] = a;
    return *this;
}  //  assignment operator with a scalar argument

template<int N>
const Vec<N>& Vec<N>::operator+=(const Vec<N>&v){
    for(int i = 0; i < N; i++)
        _elems[i] += v(i);
    return *this;
}  //  adding a Vec to the current Vec

template<int N>
const Vec<N>& Vec<N>::operator-=(const Vec<N>&v){
    for(int i = 0; i < N; i++)
        _elems[i] -= v(i);
    return *this;
}  //  subtracting a Vec from the current Vec

template<int N>
const Vec<N>& Vec<N>::operator*=(const double&  a){
    for(int i = 0; i < N; i++)
        _elems[i] *= a;
    return *this;
}  //  multiplying the current Vec by a scalar

template<int N>
const Vec<N>& Vec<N>::operator/=(const double&  a){
    for(int i = 0; i < N; i++)
        _elems[i] /= a;
    return *this;
}  //  multiplying the current Vec by a scalar

template<int N>
const Vec<N> operator+(const Vec<N>&u, const Vec<N>&v){
    return Vec<N>(u) += v;
}  //  Vec plus Vec

template<int N>
const Vec<N> operator-(const Vec<N>&u, const Vec<N>&v){
    return Vec<N>(u) -= v;
}  //  Vec minus Vec

template<int N>
const Vec<N> operator*(const Vec<N>&u, const double&  a){
    return Vec<N>(u) *= a;
}  //  Vec times scalar

template<int N>
const Vec<N> operator*(const double&  a, const Vec<N>&u){
    return Vec<N>(u) *= a;
}  //  'T' times Vec

template<int N>
const Vec<N> operator/(const Vec<N>&u, const double&  a){
    return Vec<N>(u) /= a;
}  //  Vec times scalar

template<int N>
const Vec<N>& operator+(const Vec<N>&u){
    return u;
}  //  negative of a Vec

template<int N>
const Vec<N> operator-(const Vec<N>&u){
    return Vec<N>(u) *= -1;
}  //  negative of a Vec

template<int N>
const double operator*(const Vec<N>&u, const Vec<N>&v){
    double sum = 0;
    for(int i = 0; i < N; i++)
        sum += u(i) * +v(i);
    return sum;
}  //  Vec times Vec (inner product)


template<int N>
void print(const Vec<N>&v){
    printf("(");
    for(int i = 0;i < N; i++){
        printf("v[%d]=",i);
        print(v[i]);
    }
    printf(")\n");
}  //  printing a Vec



#endif
