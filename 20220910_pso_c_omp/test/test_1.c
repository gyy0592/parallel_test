#include<stdio.h>
#include<math.h>
int main()
{
    double func1(double x0,double x1,double x2);
 
    // test 1
    double x0=-7.490342,  x1=-7.069157,  x2=1.610027; // fitFuncVal=156.267535
    double out = func1(x0,x1,x2);
    printf("test 1, res = %f\n",out-156.267535);

    // test 2
    x0=1.122411,  x1=0.344321,  x2=0.925163;  //fitFuncVal=7.061207
    out = func1(x0,x1,x2);
    printf("test 2, res = %f\n",out-7.061207);

    // test 3
    //tid=0, lpParticles=4
    x0=5.385253,  x1=3.834366,  x2=1.144854; // fitFuncVal=26.036910
    out = func1(x0,x1,x2);
    printf("test 3, res = %f\n",out-26.036910);

    // test 4
    //Best Location found: 
    x0 = -0.839485,  x1= 0.621972, x2 = 4.110079;   //Best Fitness Value: 6.514941
    out = func1(x0,x1,x2);
    printf("test 4, res = %f\n",out-6.514941);

    // test 5
    //tid=2, lpParticles=14
    x0=-1.771898,  x1=-5.485202,  x2=-6.515470;  //fitFuncVal=154.255845
    out = func1(x0,x1,x2);
    printf("test 5, res = %f\n",out-154.255845);

    // test 6
    x0=6.818244,   x1=-5.164693,  x2=-4.284120;      //fitFuncVal=138.243184
    out = func1(x0,x1,x2);
    printf("test 6, res = %f\n",out-138.243184);

    // test 7
    //tid=4, lpParticles=24
    x0=-3.713765,  x1=-2.325667,  x2=5.275360;  //fitFuncVal=46.108238
    out = func1(x0,x1,x2);
    printf("test 7, res = %f\n",out-46.108238);

    // test 8
    x0=5.385253,  x1=3.834366,  x2=1.144854;  //fitFuncVal=26.036910
    out = func1(x0,x1,x2);
    printf("test 8, res = %f\n",out-26.036910);

}

double func1(double x0,double x1,double x2)
{
    double out = (x0-1.0)*(x0-1.0) + (x1-2.0)*(x1-2.0) + (x2-3.0)*(x2-3.0);
    return out;
}
