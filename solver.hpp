#pragma once
#include <complex>
using namespace std;

namespace solver 
{


    class RealVariable 
    {
        double a;
        double b;
        double c;
    
    public:
        RealVariable(): a(0), b(1), c(0) {}
        RealVariable(double a, double b, double c)
        {
            this->a = a;
            this->b = b;
            this->c = c;
        }
        const double & getA() const 
        {
            return a; 
        }
        const double & getB() const 
        {
            return b;
        }
        const double & getC() const 
        {
            return c; 
        }
        friend RealVariable operator + (const RealVariable&,const double);
        friend RealVariable operator + (const double,const RealVariable&);
        
        friend RealVariable operator - (const RealVariable&,const double);
        friend RealVariable operator - (const double,const RealVariable&);
        RealVariable operator- () const;
        
        friend RealVariable operator * (const RealVariable&,const double);
        friend RealVariable operator * (const double,const RealVariable&);
        
        friend RealVariable operator / (const RealVariable&,const double);
        
        friend RealVariable operator ^ (const RealVariable&,const double);
        
        friend RealVariable operator + (const RealVariable&,const RealVariable&);
        friend RealVariable operator - (const RealVariable&,const RealVariable&);
        friend RealVariable operator * (const RealVariable&,const RealVariable&);
        friend RealVariable operator / (const RealVariable&,const RealVariable&);
        
        friend RealVariable operator== (const RealVariable& ,const RealVariable&);
        friend RealVariable operator== (const double n,const RealVariable&);
        friend RealVariable operator== (const RealVariable& , const double);

        bool operator! () const;

 
    };
    class ComplexVariable
    {
        complex<double> a;
        complex<double> b;
        complex<double> c;
        public:
        ComplexVariable(): a(0), b(1,0), c(0) {}

        ComplexVariable(complex<double> a, complex<double> b, complex<double> c) 
        {
            this->a = complex(a);
            this->b = complex(b);
            this->c = complex(c);
        }
        
        const complex<double> & get_a() const 
        {
            return a; 
        }
        const complex<double> & get_b() const
        {
            return b;
        }
        const complex<double> & get_c() const 
        {
            return c;
        }
        friend ComplexVariable operator+ (const ComplexVariable& x, const complex<double> n);
        friend ComplexVariable operator+ (const complex<double> n ,const ComplexVariable& x);
        friend ComplexVariable operator+ (const ComplexVariable& left,const ComplexVariable& right);
        friend ComplexVariable operator- (const ComplexVariable& x, const complex<double> n);
        friend ComplexVariable operator- (const complex<double> n ,const ComplexVariable& x);
        friend ComplexVariable operator- (const ComplexVariable& left,const ComplexVariable& right);
        ComplexVariable operator- () const;
        friend ComplexVariable operator* (const ComplexVariable& x, const complex<double> n);
        friend ComplexVariable operator* (const complex<double> n ,const ComplexVariable& x);
        friend ComplexVariable operator* (const ComplexVariable& left,const ComplexVariable& right);
        friend ComplexVariable operator^ (const ComplexVariable& x, const complex<double> n);
        friend ComplexVariable operator/ (const ComplexVariable& x, const complex<double> n);
        friend ComplexVariable operator/ (const ComplexVariable& left, const ComplexVariable& right);
        friend ComplexVariable operator== (const ComplexVariable& left,const ComplexVariable& right);
        friend ComplexVariable operator== (const complex<double> n,const ComplexVariable& right);
        friend ComplexVariable operator== (const ComplexVariable& left, const complex<double> n);
               bool operator! () const;
        friend string operator+ (const ComplexVariable& x, const string str);
        friend string operator+ (const string str ,const ComplexVariable& x);
    };
    double solve(const RealVariable& );
    double solve(const bool);
    

    complex<double> solve(const ComplexVariable&);
};