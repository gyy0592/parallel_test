#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>
int main()
{

float start,end;

start = omp_get_wtime();
for(int i=0;i<10000000;i++)
{
    #pragma omp parallel
    {

    }
}
end = omp_get_wtime();
printf("cost 1 is %f sec by omp_get_wtime.\n", end-start);


start = omp_get_wtime();
for(int i=0;i<10000000;i++)
{
    #pragma omp parallel
    {
        #pragma omp parallel for
        for(int j=0;j<1;j++)
        {

        }
    } 
}
end = omp_get_wtime();
printf("cost 2 is %f sec by omp_get_wtime.\n", end-start);


start = omp_get_wtime();
for(int i=0;i<10000000;i++)
{
    #pragma omp parallel for
    for(int j=0;j<1;j++)
    {

    }
}
end = omp_get_wtime();
printf("cost 3 is %f sec by omp_get_wtime.\n", end-start);

}
