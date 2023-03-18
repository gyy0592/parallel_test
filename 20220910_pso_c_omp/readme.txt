# introduction
This is PSO c code extracted from Mohanty&HUST_WangYan's Rappter code which search SMBHBs in PTA.
Here parallel the fitness function of 40 partilces using openmp. 

# code
In ptapso_omp.c, here is pso algorithm (for minimization).
In ptapsotstfunc.c, given the user specified fitness function.
test_ptapso_omp.c is the main file.

# run the code
run_pso_omp_1.sh
 
# The pso optimization output
The pso optimization process are strored in test_ptapso_Dump.txt,  
# (Coord),(velocity), (Pbest), (SnrPbest,SnrCurr,Inertia), (LocalBest),X,FitEvals
at each row for each particle.

# this fitness function include vectors.
# 2022/09

#performance test: my fitnesss function in ptapsotestfunc.c 
pso set up: 40 particels, 100 steps, 1 pso run
serial                       6.828308 sec
export OMP_NUM_THREADS=1     7.183281 sec
export OMP_NUM_THREADS=2     3.815021 sec
export OMP_NUM_THREADS=4     2.543299 sec
export OMP_NUM_THREADS=5     3.240243 sec
export OMP_NUM_THREADS=8     2.450021 sec
# the pso_omp works.  
# 2022/09











