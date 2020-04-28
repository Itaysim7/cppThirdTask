#include "doctest.h"
#include <iostream>
#include <complex>
#include "solver.hpp"

using namespace std; 
using solver::solve, solver::RealVariable, solver::ComplexVariable;
TEST_CASE("linear equations for RealVarivale")
{
    RealVariable x;
    CHECK(solve(0==0)==0);
    CHECK(solve(3==3)==0);
    CHECK(solve(5*x-8==12)==4);
    CHECK(solve(x*5-4==16)==4);
    CHECK(solve(2*x-4==10)==7);
    CHECK(solve(x*2-4==10)==7);
    CHECK(solve(x-4==-10)==-6);
    CHECK(solve(-x-4==-10)==6);
    CHECK(solve(2*x-4==11)==7.5);
    CHECK(solve(x*2-4==11)==7.5);
    CHECK(solve(-4+x*2==11)==7.5); 
    CHECK(solve(11==-4+x*2)==7.5); 
    CHECK(solve(-4-2+2+x*2==11)==7.5); 
    CHECK(solve(-4-2+(2+x)*2==11)==6.5); 
    CHECK(solve(-4-2+(2+x)*2==11+7)==10); 
    CHECK(solve(5*x-8==5*x-8)==0);
    CHECK(solve(2*x-4==2*5)==7);
    CHECK(solve(x+x-4==10)==7);
    CHECK(solve(3*x-x-4==10)==7);
    CHECK(solve(3*x+4*x-4==10)==2);
    CHECK(solve(3*x-4*x-4==10)==-14);    
    CHECK(solve(2*x+3*x-4*x-4==10)==14);
    CHECK(solve(x-4=8+2*x)==-4);
    CHECK(solve(2*x/2-4==10-x)==7);
    CHECK(solve(3*x-2*x==0) == 0);
    CHECK(solve(2.0*x==x*3.0) == 0);
    CHECK(solve((2+2*x)/2-5==10)==14);
    CHECK(solve(-2*x-4.0==10.0)==-7);
    CHECK(solve(5*x/2-4==10-x)==4);
    CHECK(solve(-x-4.0==10.0)==-14);
    CHECK(solve(x==7)==7);
    CHECK(solve(x+2==2)==0);
    CHECK(solve(x-4.0 == 10.0) == 14);
}
TEST_CASE("polinomial equations for RealVariable") 
{
    RealVariable x;
    CHECK(((solve((x^2)==16)==4)||(solve((x^2)==16)==-4)));
    CHECK(((solve((x^2)==25)==5)||(solve((x^2)==25)==-5)));
    CHECK(((solve((x^2)==36)==6)||(solve((x^2)==36)==-6)));
    CHECK((solve((x^2)==0)==0));
    CHECK(((solve((x^2)==32+4)==6)||(solve((x^2)==32+4)==-6)));
    CHECK(((solve((x^2)+2==18)==4)||(solve((x^2)+2==18)==-4)));
    CHECK(((solve((x^2)-4==32)==6)||(solve((x^2)-4==32)==-6)));
    CHECK(((solve((x^2)==9*4)==6)||(solve((x^2)==9*4)==-6)));
    CHECK(((solve(2*(x^2)==9*8)==6)||(solve(2*(x^2)==9*8)==-6)));
    CHECK(((solve(2*(x^2)==32)==4)||(solve(2*(x^2) == 32)==-4)));
    CHECK(((solve((x^2)*2==32)==4)||(solve((x^2)*2 == 32)==-4)));
    CHECK(((solve(2*(x^2)-4*x==6)==-1)||(solve(2*(x^2)-4*x==6)==3)));
    CHECK(((solve((x^2)+2*x==8)==-4)||(solve((x^2)+2*x==8)==2)));
    CHECK(((solve((x^2)+8*x==-12)==-2)||(solve((x^2)+8*x==-12)==-6)));
    CHECK(((solve((x^2)+2*x/2==16+x)==4)||(solve((x^2)+2*x/2==16+x)==-4)));
    CHECK(((solve(x*x==36)==6)||(solve(36)==-6)));
    CHECK(((solve(x*x==100)==10)||(solve(100)==-10)));
    CHECK(((solve(x*x==16)==4)||(solve(x*x==16)==-4)));
    CHECK(((solve(x*x==225)==15)||(solve(x*x==225)==-15)));
    CHECK(((solve(x*x+(x+x)/2==16+x)==4)||(solve(x*x+(x+x)/2==16+x)==-4)));
    CHECK(((solve((x+1)*(x+1)==16+x*2+1)==4)||(solve((x+1)*(x+1)==16+x*2+1)==-4)));
    CHECK(((solve((x^2)==16.0)==4)||(solve((x^2)==16.0)==-4)));
    CHECK(((solve((x^2)==64.0)==8)||(solve((x^2)==64.0)==-8)));
    CHECK(((solve((x^2)+0*x==16)==4)||(solve((x^2)+0*x==16)==-4)));
    CHECK(((solve(0*(x^2)+(x^2)==16)==4)||(solve(0*(x^2)+(x^2)==16)==-4)));
    CHECK(((solve((x^2)/x==16)==16)));
    CHECK(((solve((x^2)/x==32)==32)));
    CHECK(((solve((x^2)/(2*x)==32)==64)));
    CHECK(((solve((x^2)/(x*2)==32)==64)));
    CHECK(((solve((x^2)/(x*2)==50.0)==100)));
    CHECK(((solve((x^2)/x+(x^2)==16+x)==4)||(solve((x^2)/x+(x^2)==16+x)==-4)));
    CHECK(((solve((x^2)*2/2==16)==4)||(solve((x^2)*2/2==16)==-4)));
    CHECK(((solve(2*(x^2)/2==16)==4)||(solve(2*(x^2)/2==16)==-4)));
    CHECK(((solve(((x^2)+x)/x==33)==32)));
    CHECK(((solve(((x^2)+x)/x==65)==64)));
    CHECK(((solve(((x^2)+x)/x==101)==100)));
    CHECK(((solve(((2*(x^2))+x)/x==101)==50)));
    CHECK(((solve(((2*(x^2))+x)/x==67)==33)));
    CHECK(((solve((x^2)/(x^2)==1)==0)));
    CHECK(((solve(x/x==1)==0)));
    CHECK(((solve(2*(x^2)+x+x+1==17+(x^2)+2*x)==4)||(solve((x^2)+x+x+1==17+(x^2)+2*x)==-4)));
    CHECK(((solve((x^2)-16==0)==4)||(solve((x^2)-16==0)==-4)));
    CHECK(((solve(0==(x^2)-16)==4)||(solve(0==(x^2)-16)==-4)));
    CHECK(((solve((x^2)==16)==4)||(solve((x^2)==16)==-4)));
    CHECK(((solve((x^1)*x==16)==4)||(solve((x^1)*x==16)==-4)));
    CHECK(((solve((x^0)+(x^2)==16)==4)||(solve(((x^0)+(x^2)==16)==-4))));
    CHECK(((solve((x^0)*(x^2)==16)==4)||(solve(((x^0)*(x^2)==16)==-4))));
    CHECK(((solve((x^0)*(x^2)==36)==6)||(solve(((x^0)*(x^2)==36)==-6))));
    CHECK(((solve((x^0)*(x^2)==121)==11)||(solve(((x^0)*(x^2)==121)==-11))));
    CHECK(((solve((x^0)*(x^2)==400)==20)||(solve(((x^0)*(x^2)==400)==-20))));
    CHECK(((solve((x^1)*(x^1)==16)==4)||(solve((x^1)*(x^1)==16)==-4)));
    CHECK(((solve((x^1)*(x^1)==49)==7)||(solve((x^1)*(x^1)==49)==-7)));
    CHECK(((solve((x^1)*(x^1)==81)==9)||(solve((x^1)*(x^1)==81)==-9)));
    CHECK(((solve((x^2)-x==16-x)==4)||(solve((x^2)-x==16-x)==-4)));
    CHECK((solve((x^2)==0)==0));
    CHECK(solve((x^2)/x==1)==1);
}
TEST_CASE("checks throws for real numbers"){
    RealVariable x;
    CHECK_THROWS(solve(x/0 ==16));
    CHECK_THROWS(solve((x^5)==0));
    CHECK_THROWS(solve((x^2)/(x-x)==9));
    CHECK_THROWS(solve((x^2)*x==0));
    CHECK_THROWS(solve((x^1)/x==0));
    CHECK_THROWS(solve((x^2)*x==0));
    CHECK_THROWS(solve((x^-2)==0));
    CHECK_THROWS(solve((x^2)*(x^2) == 0));
    CHECK_THROWS(solve(2==3));
    CHECK_THROWS(solve(x/(x^2)==0));
    CHECK_THROWS(solve(x-x==2));
    CHECK_THROWS(solve(x/x==3));
    CHECK_THROWS(solve(x/(x-1)==1));
    CHECK_THROWS(solve(x-x==1));
    CHECK_THROWS(solve(2*x+2==5+2*x));
    CHECK_THROWS(solve(2*x/0==2));
    CHECK_THROWS(solve(2*x/(2-2)==2));
    CHECK_THROWS(solve((x+1)*(x^2)==2));
    CHECK_THROWS(solve(x/(x^2+2)==2));

}
