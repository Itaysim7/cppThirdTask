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
    
    
    
    ComplexVariable operator* (const ComplexVariable& x, const complex<double> n) {
        return ComplexVariable(x.a*n, x.b*n, x.c*n);

    }

    ComplexVariable operator* (const complex<double> n , const ComplexVariable& x) {
        return ComplexVariable(x.a*n, x.b*n, x.c*n);

    }

    ComplexVariable operator* (const ComplexVariable& left, const ComplexVariable& right) {
        if(left.a != complex(0.0,0.0) && right.a != complex(0.0,0.0)) throw invalid_argument("power is too big");
        if(left.a != complex(0.0,0.0) && right.b != complex(0.0,0.0)) throw invalid_argument("power is too big");
        if(left.b != complex(0.0,0.0) && right.a != complex(0.0,0.0)) throw invalid_argument("power is too big");

        return ComplexVariable(left.a * right.c + right.a * left.c + left.b * right.b ,
                                 left.b * right.c + right.b * left.c,
                                     left.c * right.c);

    }   

    ComplexVariable operator+ (const ComplexVariable& x, const complex<double> n) {
        return ComplexVariable(x.a, x.b , x.c+n);

    }
    
    ComplexVariable operator+ (const complex<double> n , const ComplexVariable& x) {
        return ComplexVariable(x.a, x.b , x.c+n);
    }

    ComplexVariable operator+ (const ComplexVariable& left, const ComplexVariable& right) {
        return ComplexVariable(left.a + right.a , left.b + right.b , left.c + right.c);
    
    }

    ComplexVariable operator- (const ComplexVariable& x, const complex<double> n) {
        return ComplexVariable(x.a, x.b , x.c-n);
        
    }
    
    ComplexVariable operator- (const complex<double> n , const ComplexVariable& x) {
        return ComplexVariable(-x.a, -x.b , -x.c+n);
        
    }

    ComplexVariable operator- (const ComplexVariable& left, const ComplexVariable& right) {
        return ComplexVariable(left.a - right.a , left.b - right.b , left.c - right.c);

    }

    ComplexVariable operator^ (const ComplexVariable& x, const complex<double> n) {
        if(n.imag() != 0) throw invalid_argument("complex power is not valid");

        if(n.real() == 2) {
            return ComplexVariable(complex(1.0,0.0),complex(0.0,0.0),complex(0.0,0.0));
        }
        if (n.real() == 1) {
            return ComplexVariable(complex(0.0,0.0),complex(1.0,0.0),complex(0.0,0.0));
        }
        if (n.real() == 0) {
            return ComplexVariable(complex(0.0,0.0),complex(0.0,0.0),complex(1.0,0.0));
        }
        throw invalid_argument("power is not vaild: ");

    }

    ComplexVariable operator/ (const ComplexVariable& x, const complex<double> n) {
        if(n == complex(0.0,0.0)) throw invalid_argument("Division by zero is not valid");

        return ComplexVariable(x.a/n, x.b/n, x.c/n);

    }
    
    ComplexVariable operator/ (const ComplexVariable& left, const ComplexVariable& right) {
        ComplexVariable x;

        if(!right) throw invalid_argument("Division by zero is not valid");

        if(!(left - right)) return ComplexVariable(0,0,complex(1.0,0.0));
        if(!(left - (x^2)) && !(right - x)) return ComplexVariable(0,complex(1.0,0.0),0);
        if(right.a == complex(0.0,0.0) && right.b == complex(0.0,0.0)) return left/right.c;

        throw invalid_argument("divison is out of our assignment's scope");
    }


    ComplexVariable operator== (const ComplexVariable& r1, const ComplexVariable& r2) {
        return r1-r2;

    }

    ComplexVariable operator== (const complex<double> n, const ComplexVariable& r) 
    {
        return r-n;
    }
    
    ComplexVariable operator== (const ComplexVariable& r, const complex<double> n)
    {
        return r-n;
    }
     bool ComplexVariable::operator! () const 
     {
        return ((this->a.real() == 0.0) && (this->a.imag() == 0.0) && 
                (this->b.real() == 0.0) && (this->b.imag() == 0.0) &&
                (this->c.real() == 0.0) && (this->c.imag() == 0.0));
    }

    string operator+ (const ComplexVariable& x, const string str) {
        return "("+ to_string(x.a.real()) + "+" + to_string(x.a.imag())  + "i)*x^2 +" + "("+ to_string(x.b.real()) + "+" + to_string(x.b.imag())  + "i)*x +("+ to_string(x.c.real()) + "+" + to_string(x.c.imag())  + "i)" + str;
    }
    
    string operator+ (const string str , const ComplexVariable& x) {
        return str + "("+ to_string(x.a.real()) + "+" + to_string(x.a.imag())  + "i)*x^2 +" + "("+ to_string(x.b.real()) + "+" + to_string(x.b.imag())  + "i)*x +("+ to_string(x.c.real()) + "+" + to_string(x.c.imag())  + "i)";

    }
// /////////////////////////////////////////////////////////////

    double solve(const RealVariable& r)
    {
        double a=r.getA();
        double b=r.getB();
        double c=r.getC();
        if(!r)
            return 0;
        if(a==0&&b!=0) 
            return (c/-b);
        if(a==0&&b==0) 
            return c;
        double formula=b*b-4*a*c;
        if (formula>=0) 
            return ((-b + sqrt(formula))/(2*a));
        throw runtime_error("There is no real solution");
    }

    complex<double> solve(const ComplexVariable& equation)
    {
        complex<double> a = equation.get_a(), b = equation.get_b(), c = equation.get_c();

        if(!equation) return 0;

        if(a == complex(0.0,0.0)) { // The equation is linear
            if ((b == complex(0.0,0.0)) && (c != complex(0.0,0.0))) throw runtime_error("There is no solution for this equation");
            else return -c/b;
        }
        
        complex<double> x1, discriminant;

        discriminant = b*b - complex(4.0,0.0)*a*c;
    
        x1 = (-b + sqrt(discriminant)) / (complex(2.0,0.0)*a);

        return x1;

    }

    double solve(const bool t)
    {
        if(t)
            return 0;
        throw runtime_error("There is no solution for this equation");
    }

  


}