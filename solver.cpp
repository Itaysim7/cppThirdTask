#include "solver.hpp"
#include <iostream>
#include <string>

using namespace std;
namespace solver 
{
    RealVariable operator+ (const RealVariable& r, const double n) 
    {
        return RealVariable(r.a, r.b, r.c+n);
    }
    RealVariable operator+ (const double n,const RealVariable& r) 
    {
        return RealVariable(r.a, r.b, r.c+n);
    }
    RealVariable operator- (const RealVariable& r, const double n) 
    {
        return RealVariable(r.a, r.b, r.c-n);
    }
    RealVariable operator- (const double n,const RealVariable& r) 
    {
        return RealVariable(-r.a, -r.b, -r.c+n);
    }
    RealVariable RealVariable::operator- () const 
    {
        return (-1)* (*this);
    }
    RealVariable operator* (const RealVariable& r, const double n) 
    {
        return RealVariable(r.a*n, r.b*n, r.c*n);
    }
    RealVariable operator* (const double n,const RealVariable& r) 
    {
        return RealVariable(r.a*n, r.b*n, r.c*n);
    }
    RealVariable operator / (const RealVariable& r, const double n) 
    {
        if(n==0)
        {
            throw runtime_error("you can't divide by zero");
        }
        return RealVariable(r.a/n, r.b/n, r.c/n);
    }
    RealVariable operator + (const RealVariable& r1,const RealVariable& r2) 
    {
        return RealVariable(r1.a+r2.a, r1.b+r2.b, r1.c+r2.c);
    }
    RealVariable operator- (const RealVariable& r1,const RealVariable& r2) 
    {
        return RealVariable(r1.a-r2.a, r1.b-r2.b, r1.c-r2.c);
    }
    RealVariable operator* (const RealVariable& r1,const RealVariable& r2) 
    {
        if((r1.a>0&&r2.a>0)||(r1.a>0&&r2.b>0)||(r2.a>0&&r1.b>0))
        {
            throw runtime_error("power is bigger then 2");
        }
        return RealVariable(r1.a*r2.c+r2.a*r1.c+r1.b*r2.b, r1.b*r2.c+r2.b*r1.c, r1.c*r2.c);
    }
    RealVariable operator/ (const RealVariable& r1,const RealVariable& r2) 
    {
        if(!r2)
            throw runtime_error("divide by zero");
        if(!(r1-r2))
            return RealVariable(0,0,1);
        if(!(r2-r1))
            return RealVariable(0,0,-1);
        if(r1.a!=0&&r1.b==0&&r1.c==0&&r2.a!=0&&r2.b==0&&r2.c==0)//case r1 and r2 is type of d*(x)^2
            return RealVariable(0,0,r1.a/r2.a);
        if(r1.a==0&&r1.b!=0&&r1.c==0&&r2.a==0&&r2.b!=0&&r2.c==0)//case r1 and r2 is type of d*x
            return RealVariable(0,0,r1.b/r2.b);
        if(r1.c==0&&r2.a==0&&r2.b!=0&&r2.c==0)//case r1 is type of d*(x)^2+x and r2 is type of d*x
            return RealVariable(0,r1.a/r2.b,r1.b/r2.b);
        if(r2.b==0&&r2.c==0)// r2 is type of double
            return r1/r2.c;
        throw runtime_error("illegal divide");
    }
    RealVariable operator^ (const RealVariable& r, const double n)
    {
        if(n==1)
            return RealVariable(r.a,r.b,r.c);
        if(n==0)
            return RealVariable(0,0,1);
        if(n>2||n<0)
            throw runtime_error("illegal power");
        if(n==2&&r.a>0)
            throw runtime_error("illegal power");
        return RealVariable(r.b,0,r.c*r.c);
    }
    RealVariable operator== (const double n, const RealVariable& r)
    {
        return n-r;
    }
    RealVariable operator== (const RealVariable& r, const double n)
    {
         return r-n;
    }
    RealVariable operator== (const RealVariable& r1, const RealVariable& r2)
    {
        return r1-r2;
    }

    bool RealVariable::operator! () const
    {
        if(this->a==0&&this->b==0&&this->c==0)
        {
            return true;
        }
        return false;
    }
    
    // ///////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    
    
    ComplexVariable operator+ (const ComplexVariable& c, const complex<double> n) 
    {
        return ComplexVariable(c.a,c.b,c.c+n);
    }
    ComplexVariable operator+ (const complex<double> n , const ComplexVariable& c)
    {
        return ComplexVariable(c.a,c.b,c.c+n);
    }
    ComplexVariable operator- (const complex<double> n , const ComplexVariable& c)
    {
        return ComplexVariable(-c.a, -c.b , -c.c+n);
    }
    ComplexVariable operator- (const ComplexVariable& c, const complex<double> n)
    {
        return ComplexVariable(c.a,c.b,c.c-n);
    }
    ComplexVariable operator* (const ComplexVariable& c, const complex<double> n)
    {
        return ComplexVariable(c.a*n, c.b*n, c.c*n);
    }
    ComplexVariable operator* (const complex<double> n , const ComplexVariable& c)
    {
        return ComplexVariable(c.a*n, c.b*n, c.c*n);
    }
    ComplexVariable operator/ (const ComplexVariable& c, const complex<double> n)
    {
        if(n.real()==0&&n.imag()==0)
            throw invalid_argument("Division by zero is not valid");
        return ComplexVariable(c.a/n,c.b/n,c.c/n);
    }
    ComplexVariable operator+ (const ComplexVariable& c1, const ComplexVariable& c2)
    {
        return ComplexVariable(c1.a+c2.a,c1.b+c2.b,c1.c+c2.c);
    }
    ComplexVariable operator- (const ComplexVariable& c1, const ComplexVariable& c2)
    {
        return ComplexVariable(c1.a-c2.a,c1.b-c2.b,c1.c-c2.c);
    }
    ComplexVariable operator* (const ComplexVariable& c1, const ComplexVariable& c2)
    {
        if((c1.a!=complex(0.0,0.0)&&c2.a!=complex(0.0,0.0))||(c1.a!=complex(0.0,0.0)&&c2.b!=complex(0.0,0.0))||(c2.a!=complex(0.0,0.0)&&c1.b!=complex(0.0,0.0)))
        {
            throw runtime_error("power is bigger then 2");
        }
        return ComplexVariable(c1.a*c2.c+c2.a*c1.c+c1.b*c2.b,c1.b*c2.c+c2.b*c1.c,c1.c*c2.c);
    }   
    
    ComplexVariable operator^ (const ComplexVariable& C, const complex<double> n) 
    {
        if(n.imag() != 0) 
            throw runtime_error("Complex power is not valid");
        if(n.real() == 2) 
            return ComplexVariable(complex(1.0,0.0),complex(0.0,0.0),complex(0.0,0.0));
        if (n.real() == 1)
            return ComplexVariable(complex(0.0,0.0),complex(1.0,0.0),complex(0.0,0.0));
        if (n.real() == 0) 
            return ComplexVariable(complex(0.0,0.0),complex(0.0,0.0),complex(1.0,0.0));
        throw runtime_error("Complex power is not valid");
    }
    ComplexVariable operator/ (const ComplexVariable& c1, const ComplexVariable& c2) 
    {
        ComplexVariable x;
        if(!c2)
            throw runtime_error("divide by zero");
        if(!(c1-c2)) 
            return ComplexVariable(0,0,complex(1.0,0.0));
        if(!(c2-c1)) 
            return ComplexVariable(0,0,complex(1.0,0.0));
        if(!(c1-(x^2))&&!(c2-x))
            return ComplexVariable(0,complex(1.0,0.0),0);
        if(c2.a==complex(0.0,0.0)&&c2.b==complex(0.0,0.0))
            return c1/c2.c;
        throw runtime_error("error in divide");
    }
    ComplexVariable operator== (const ComplexVariable& c1, const ComplexVariable& c2) 
    {
        return c1-c2;
    }
    ComplexVariable operator== (const complex<double> n, const ComplexVariable& c) 
    {
        return n-c;
    }
    ComplexVariable operator== (const ComplexVariable& c, const complex<double> n)
    {
        return c-n;
    }
    bool ComplexVariable::operator! () const 
    {
        if ((this->a.real()==0.0)&&(this->a.imag()==0.0)&&(this->b.real()==0.0)&&(this->b.imag()==0.0)&&(this->c.real()== 0.0)&&(this->c.imag()==0.0))
            return true;
        return false;
    }

// /////////////////////////////////////////////////////////////

    double solve(const RealVariable& r)
    {
        double a=r.getA();
        double b=r.getB();
        double c=r.getC();
        if(!r)
            return 0;
        if(a==0&&b==0&&c!=0) 
            throw runtime_error("There is no solution");
        if(a==0&&b!=0) 
            return (c/-b);
        if(a==0&&b==0) 
            return c;
        double formula=b*b-4*a*c;
        if (formula>=0) 
            return ((-b + sqrt(formula))/(2*a));
        throw runtime_error("There is no real solution");
    }
    complex<double> solve(const ComplexVariable& c1)
    {
        complex<double> a=c1.getA();
        complex<double> b=c1.getB();
        complex<double> c=c1.getC();
        if(!c1)
            return 0;
        if(a==complex(0.0,0.0))
        { 
            if ((b==complex(0.0,0.0))&&(c!=complex(0.0,0.0)))
                throw runtime_error("There is no solution for this equation");
            else return -c/b;
        }
        complex<double> formula=b*b-complex(4.0,0.0)*a*c;
        return ((-b + sqrt(formula)) / (complex(2.0,0.0)*a));
    }
    double solve(const bool t)
    {
        if(!t)
            throw runtime_error("There is no solution for this equation");
        return 0;
    }
}