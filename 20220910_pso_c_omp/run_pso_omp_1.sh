# test pso_omp
rm a.out
rm *.o

export OMP_NUM_THREADS=5
gcc -Wall  -c ptapso_omp.c -fopenmp
gcc -Wall  -c ptapsotestfunc.c
gcc -Wall  -c test_ptapso_omp.c -fopenmp
gcc ptapso_omp.o  ptapsotestfunc.o test_ptapso_omp.o  -lgsl -lgslcblas -lm -fopenmp

./a.out
