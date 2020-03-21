#define _DEFAULT_SOURCE


#include<unistd.h>

typedef double v4df __attribute__ ((vector_size (32)));



#define G 6.67408E-11
#define Minute 60.0
#define Hour Minute*60.0
#define Day Hour*24

#define FL __FILE__ __LINE__

#define RING_PASS_TAG 0
#define EMIT_POINT_TAG 1


//#define DEV

//#define HGHEST_PROCESS_DOES_COMPUTE

//____________________________ Starting Debuggin macros ____________________________________

#define ForceToPrint 0
#define MAX_FORCE 1e100

//____________________________ Ending Debuggin macros ____________________________________



#define idprint(f_ , ...) printf( "[%d] " f_ , (id), ##__VA_ARGS__ );\



#define newlines( N )\
for(int sduihf9852738=0; sduihf9852738 < N; sduihf9852738++){\
    printf("\n");\
    fflush_stdout\
}\

#define ordp(f_,...)\
\
for(int uuuu =0; uuuu < p;uuuu++){\
    MPI_Barrier(MPI_COMM_WORLD);\
    if(id == uuuu){\
        idprint(f_, ##__VA_ARGS__);\
        fflush_stdout\
        usleep(200);\
    }\
}\



#define fflush_stdout fflush(stdout);


#define pv_format(in)\
    "[%d] x %s =%7.3e\n[%d] y %s =%7.3e\n[%d] z %s=%7.3e\n", id, #in, in[0],id, #in, in[1], id, #in, in[2]

#define pv(in)\
printf(pv_format( in ));\
fflush_stdout\

#define pv_special_format(low,high ,in)\
 "<%d><-><%d> x %s =%7.3e\n<%d><-><%d> y %s =%7.3e\n<%d><-><%d> z %s=%7.3e\n", low, high, #in, in[0],low, high, #in, in[1], low, high, #in, in[2]

#define pv_special(low, high, in)\
printf(pv_special_format(low, high, in ));\
fflush_stdout\

#define filterp(req_id, fmt, ...)\
if(id == req_id){\
    printf(fmt, ##__VA_ARGS__);\
    fflush_stdout\
}\


#define MYMPI_print(_____MPI_proc_id, f_,...)\
MPI_Barrier(MPI_COMM_WORLD);\
if(id == _____MPI_proc_id){\
    printf((f_),__VA_ARGS__);\
    fflush(stdout);\
}\
usleep(100);\
MPI_Barrier(MPI_COMM_WORLD);\


typedef struct RawSimulationData{
    int NumBodies;
    double* PosVelVec;
    double* MassVec;
}RawSimulationData;


RawSimulationData parseBodyKinematicsFile(char* );
v4df *run_simulation_parallel(double*, double*);
void emitPoints(unsigned int, v4df*, double*, int);

