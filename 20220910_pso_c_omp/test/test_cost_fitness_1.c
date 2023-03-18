#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <omp.h>
#include<time.h>

int main()
{
double sum,start1, end1,duration;	
clock_t start,end;
double x0=1.0,x1=2.0,x2=3.0;

// ----------------------- likelihood 1 -------------------------//
start = clock();
start1 = omp_get_wtime();
//
sum=0.0;
for(int i=0;i<1000;i++)
{
    sum = sum + pow(x0-1.0,2.0) + pow(x1-2.0,2.0) + pow(x2-3.0,2.0);  
}
//
end = clock();
end1 = omp_get_wtime();
duration = (double)(end-start)/CLOCKS_PER_SEC;
printf("likelihood 1, elpase time is %f sec by time clock.\n", duration);
printf("likelihood 1, elpase time is %f sec by omp_get_wtime.\n", end1-start1);
printf("\n");


// ----------------------- likelihood 2 -------------------------//
start = clock();
start1 = omp_get_wtime();
//
sum=0.0;
for(int i=0;i<10000;i++)
{
    sum = sum + pow(x0-1.0,2.0) + pow(x1-2.0,2.0) + pow(x2-3.0,2.0);  
}
//
end = clock();
end1 = omp_get_wtime();
duration = (double)(end-start)/CLOCKS_PER_SEC;
printf("likelihood 2, elpase time is %f sec by time clock.\n", duration);
printf("likelihood 2, elpase time is %f sec by omp_get_wtime.\n", end1-start1);
printf("\n");

// ----------------------- likelihood 3 -------------------------//
start = clock();
start1 = omp_get_wtime();
//
sum=0.0;
for(int i=0;i<100000;i++)
{
    sum = sum + pow(x0-1.0,2.0) + pow(x1-2.0,2.0) + pow(x2-3.0,2.0);  
}
//
end = clock();
end1 = omp_get_wtime();
duration = (double)(end-start)/CLOCKS_PER_SEC;
printf("likelihood 3, elpase time is %f sec by time clock.\n", duration);
printf("likelihood 3, elpase time is %f sec by omp_get_wtime.\n", end1-start1);

}

