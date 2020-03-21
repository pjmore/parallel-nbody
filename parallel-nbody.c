//Required
#define _DEFAULT_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#ifndef __clang__
#include <mpi.h>
#endif 
#include <stdarg.h>



#define ordp(f_,...)\
for(int uuuu =0; uuuu < p;uuuu++){\
    MPI_Barrier(MPI_COMM_WORLD);\
    if(id == uuuu){\
        printf((f_), __VA_ARGS__);\
        fflush(stdout);\
    }\
}\




#define MPI_flush_barrier_ord(file)\
for(int MPI_id_sldkfjlksd = 0; MPI_id_sldkfjlksd <p;MPI_id_sldkfjlksd++){\
    MPI_Barrier(MPI_COMM_WORLD);\
    if(id == MPI_id_sldkfjlksd){\
        fflush(file);\
    }\
    sleep(1);\
    MPI_Barrier(MPI_COMM_WORLD);\
}\





//Benchmarking 
#include <time.h>
#include <unistd.h>


//My defines
#include "parallel-nbody.h"

#define p1(f_, ...) if(id == 0){\
printf((f_), __VA_ARGS__ );\
fflush(stdout);\
}




int minInt(int,int);
double Sumv4df(v4df vec);
v4df Normalize(v4df vec);
v4df CalculateForce(v4df, v4df, double, double);
 





// ID of the process
int id;
// Number of processes running
int p;

// ID of the that will be sent data - mod(id - 1, p)
int send_to_id;
// ID of the process that this process will receive data from mod( id+1, p)
int recv_from_id;

// The maximum number of jobs that any one process has
int max_jobs;

// The minimum number of jobs that any one process has
// Will be at the least max_jobs-1
int min_jobs;


//The number of jobs that cannot be evenly distributed to processes
//These will be given to the lower rank processes until there are none left
int extra_jobs;


//The number of jobs that this process has
int num_jobs;

//Number of bodies that will be simulated
int n_bodies;

//The global start index that this process is responsible for.
int global_start_index;

//Time that the simulation will end
unsigned int num_steps;

//How many steps before a point will be emitted
unsigned int emit_steps;

// Time delta of each simulation step
double dt = -1;


//How many jobs each process has - initially used to scatter Position and velocity to processes. Later used to gather
// a single component of position at a time. Must be divided by 6 before this use transition occurs.
//MUST ONLY BE USED BY PID=0
int *jobAllocations;
// The displacement from 0 that eachs jobs section starts at - initially used to scatter Position and velocity to processes. Later used to gather
// a single component of position at a time. Must be divided by 6 before this use transition occurs.
//MUST ONLY BE USED BY PID=0
int *jobDispl;

//A temporary vector that stores a component of the position for all tasks to emit it from process 0
//MUST ONLY BE USED BY PID=0
double *EmitPointsBuf;

FILE *output;

// Must pass in <file name>, <dt>, <end time> <outputfil pathe>
int main(int argc, char *argv[]){
    if(argc < 5){
        fprintf(stderr,"Must pass in the arguments <data file path> <time delta> <number of steps> <number of steps before emitting point>         [output file path] after invoking program\n");
        fflush(stderr);
        exit(1);
    }


    char *filename = argv[1];
    dt = strtod(argv[2],NULL);
    num_steps = strtol(argv[3], NULL, 10);
    emit_steps = strtol(argv[4], NULL, 10);


    

    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //Calculate the sending and receive id. Add number of processes so it is guaranteed to be positive so
    // I can just use C's modulo operator


    //He says send to p-1
    //Reciev from p+1
    // Sending id. One less than current id
    send_to_id = (id -1 + p)%p;
    // Receiving id. one more than current
    recv_from_id = (id+1)%p;

    //test

    
     //send_to_id = (id +1)%p;
    // Receiving id. one more than current
    //recv_from_id = (id+p-1)%p;
    
    
    //Container for initial state of system
    RawSimulationData initial_conditions;


    
    //Holds all of the masses for the bodies
    double* MassBuffer;



    //File IO
    if (id == 0){
        if(argc>5){
            output = fopen(argv[5], "w");
        }else{
            output=stdout;
        }
        
        initial_conditions = parseBodyKinematicsFile(filename);

        EmitPointsBuf = (double*)malloc(sizeof(double)*initial_conditions.NumBodies);

        jobAllocations = (int*)malloc(sizeof(int)*p);
        jobDispl = (int*)malloc(sizeof(int)*p);

        int jobs_remaining = initial_conditions.NumBodies % p;
        int even_partition_size = initial_conditions.NumBodies/p;

        int curr_start_idx = 0;
        for(int i = 0; i < p; i++){
            //Assume even job size and set accordingly
            jobAllocations[i] = even_partition_size;
            //If there are more remaining jobs partition one to the current ith process
            if (jobs_remaining >0){
                jobAllocations[i]++;
                jobs_remaining--;
            }
            // Multiply by the number of doubles per body
            jobAllocations[i] *= 6;
            //Check for maximum buffer size
            if (jobAllocations[i]/6 > max_jobs){
                max_jobs = jobAllocations[i]/6;
            }
            //Set start indexes of body partition
            jobDispl[i] = curr_start_idx;
            // Increment the current start index
            curr_start_idx += jobAllocations[i];
        }
        MassBuffer = initial_conditions.MassVec;
        n_bodies = initial_conditions.NumBodies;
        printf("Simulating %d bodies\n", n_bodies);
    }



    //Broadcast the number of bodies to all processes running so that we can initialize the correct length for masses, as well as the correct length of positions.
    MPI_Bcast((void*)&n_bodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(id != 0){
        MassBuffer = (double*)malloc(sizeof(double)*n_bodies);
    }

    num_jobs = n_bodies/p;
    extra_jobs = n_bodies%p;

    // The global start index can be calculated by multiplying the current id time the number of jobs and then
    // adding the extra jobs that are governed that process
    global_start_index = id*num_jobs + minInt(id, extra_jobs);

    min_jobs = num_jobs;
    max_jobs = num_jobs;
    if (extra_jobs != 0){
        if (extra_jobs > id){
            num_jobs += 1;
        }
        max_jobs++;
    }


    // Allocating buffer to store positions and velociteis for the jobs given to process
    double *PosVelBuf = (double*)malloc(sizeof(double)*(6*max_jobs));
    //Broadcast the masses to all buffers
    MPI_Bcast((void*)MassBuffer, n_bodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Scatter the jobs accross all processes
    //Scatter evenly and if there are left over jobs add one to all of the lower jobs until there are no more
    MPI_Scatterv((void*)initial_conditions.PosVelVec, jobAllocations, jobDispl, MPI_DOUBLE, PosVelBuf, max_jobs*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //Some cleanup/prep before the simulation starts
    if(id ==0){
        free(initial_conditions.PosVelVec);
        for(int i = 0; i < p; i++){
            jobAllocations[i] = jobAllocations[i]/6;
            jobDispl[i] = jobDispl[i]/6;
        }
    }

    run_simulation_parallel(MassBuffer, PosVelBuf);
    if(id == 0){
        if(argc >=6){
            fclose(output);
        }
        free(jobAllocations);
        free(jobDispl);
    }
    free(MassBuffer);
    free(PosVelBuf);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();


    if(id != 0){
        return 0;
    }
    sleep(1);
    return 0;
}



 int minInt(int i,int j){
    if (i < j){
        return i;
    }
    return j;
}


 int localToGlobalStartIdx(int sending_id){
    return sending_id*min_jobs + minInt(sending_id, extra_jobs);
}







v4df* run_simulation_parallel(double *Masses, double *PosVelBuffer){
    int iteration;
    double t;
    double dt = 0.1;
    // Kinematics data about the particles owned by this process
    v4df *Positions = (v4df*)aligned_alloc(sizeof(v4df), sizeof(v4df)*num_jobs);
    v4df *Velocities = (v4df*)aligned_alloc(sizeof(v4df),sizeof(v4df)*num_jobs);
    v4df *Force = (v4df*)aligned_alloc(sizeof(v4df),sizeof(v4df)*num_jobs);

    

    for(int i=0; i<num_jobs; i++){
        int pvb_idx = i*6;
        Positions[i] = (v4df){PosVelBuffer[pvb_idx],PosVelBuffer[pvb_idx + 1], PosVelBuffer[pvb_idx+2],0.0};
        Velocities[i] = (v4df){PosVelBuffer[pvb_idx+3],PosVelBuffer[pvb_idx + 4], PosVelBuffer[pvb_idx+5],0.0};
    }
        


    //These are buffers for the particles owned by other processes
    v4df *tempPos = aligned_alloc(sizeof(v4df),sizeof(v4df)*max_jobs*2);
    v4df *tempForces = tempPos+max_jobs;


    v4df force;

    //Start index of the particles owned by other processes in the main large array
    int temp_glbl_idx;
    //Number of jobs recieved from other processes
    int num_recv_jobs;
    MPI_Status recv_status;
    filterp(0, "The maximum number of steps is %d\n", num_steps);
    unsigned int steps_at_next_emit = emit_steps;
    //emitPoints(0, Positions, Masses, num_jobs);
    for(unsigned int curr_step = 0; curr_step < num_steps; curr_step++){
        num_recv_jobs = num_jobs;
        //Zero force vector        
        memset(Force, 0, sizeof(v4df)*num_jobs);
        //zero the temporary force vector including the last entry in the position section of the buffer just in case
        memset(tempPos + (max_jobs-1), 0, sizeof(v4df)*(max_jobs+1));
        //Copy the new positions to the tempPos vector
        memcpy(tempPos, Positions, sizeof(v4df)*num_jobs);

        // The first send will mean that the working set is from the send_id process
        int sending_process_id = recv_from_id;
        //Perform p-1 ring passes calculating forces for other bodies not owned by process
        for(int ring_pass=0; ring_pass < p-1; ring_pass++){
            //Send num_recv_jobs*4*2 doubles to next process. Four because 4 doubles fit in the vector. 2 becuase that will send both the tempPositions and temp forces in a 
            //single send
            MPI_Sendrecv_replace((double*)tempPos, max_jobs*4*2, MPI_DOUBLE, send_to_id, RING_PASS_TAG,recv_from_id , RING_PASS_TAG, MPI_COMM_WORLD, NULL);
            num_recv_jobs = min_jobs + minInt(1, extra_jobs-sending_process_id);
            temp_glbl_idx = localToGlobalStartIdx(sending_process_id);
            // Only calculate if the sending process id is greater than the current one
            if(id < sending_process_id){
                for(int local_task_idx=0; local_task_idx < num_jobs; local_task_idx++){
                    for(int temp_task_idx = 0; temp_task_idx < num_recv_jobs; temp_task_idx++){
                        force = CalculateForce(Positions[local_task_idx], tempPos[temp_task_idx], Masses[(global_start_index+local_task_idx)], Masses[(temp_glbl_idx+temp_task_idx)]);
                        Force[local_task_idx] += force;
                        tempForces[temp_task_idx] -= force;
                    } 
                }
            }

            // If this is defined the highest process, which does pretty much nothing the whole time. Will compute the pair wise interactions
            // Within each tempPos vector and add/subtract them as neccesary to tempForces
            // My basic benchmarking found it was slower
            #ifdef HGHEST_PROCESS_DOES_COMPUTE
            else if(id == p-1 && p != 1){
                for(int l_idx=0; l_idx<num_recv_jobs; l_idx++){
                    for(int u_idx=l_idx+1; u_idx< num_recv_jobs; u_idx++){
                        force = CalculateForce(tempPos[l_idx], tempPos[u_idx], Masses[l_idx+temp_glbl_idx], Masses[u_idx+temp_glbl_idx]);
                        tempForces[u_idx] += force;
                        tempForces[l_idx] -= force;
                    }
                }
            }
            #endif

            // decrement in a mod manner 
            sending_process_id = (sending_process_id +1)%p;
        }
        //Final send receive
    
        MPI_Sendrecv_replace((void*)tempPos, max_jobs*4*2, MPI_DOUBLE, send_to_id , 0, recv_from_id, 0, MPI_COMM_WORLD, &recv_status);
        // pairwise interactions within a job
        #ifndef HGHEST_PROCESS_DOES_COMPUTE
        for(int l_idx = 0; l_idx < num_jobs -1;l_idx++){
            for(int u_idx=l_idx+1; u_idx< num_jobs; u_idx++){
                force = CalculateForce(Positions[l_idx], Positions[u_idx], Masses[global_start_index + l_idx], Masses[global_start_index + u_idx]);
                Force[l_idx]+=force;
                Force[u_idx]-=force;
            }
        }
        #endif
        
        //Adding the forces computed by other processes 
        //If the highset process doesn't do computations no other process with touch id=0 tempForce vector so it should be skipped
        #ifndef HGHEST_PROCESS_DOES_COMPUTE
        if(id != 0){
        #endif


        for(int i=0; i< num_jobs; i++){
                Force[i] += tempForces[i];    
        }
        

        #ifndef HGHEST_PROCESS_DOES_COMPUTE
        }
        #endif
        //Calulate the acceleration for all bodies owned by this process
        for(int local_idx=0; local_idx < num_jobs; local_idx++){
            v4df vPortion = dt*Velocities[local_idx];
            v4df AccPortion = dt*dt*Force[local_idx]/Masses[global_start_index+local_idx];
            Positions[local_idx] += vPortion + AccPortion;
            Velocities[local_idx] += AccPortion/dt;
        }
        
        filterp(0,"Finished step %d\n", curr_step);
        if(steps_at_next_emit == curr_step){
            steps_at_next_emit += emit_steps;
            emitPoints(curr_step, Positions, Masses, num_jobs);   
        }
    }
    //emitPoints(num_steps, Positions, Masses, num_jobs);
    
}



 double Sumv4df(v4df vec){
    return vec[0]+vec[1]+vec[2];
}

 v4df Normalize(v4df vec){
    return vec/hypot(vec[0], hypot(vec[1],vec[2]));
}



 v4df CalculateForce(v4df PosLoc, v4df PosTemp, double MassLoc, double MassTemp){
    v4df PartialGForce = {0,0,0,0};
    v4df Difference =  PosTemp - PosLoc;
    double Distance2 = Sumv4df(Difference*Difference);
    // Ignore gravity if bodies are sub 1000m to avoid signularites
    if(Distance2<1000000.0){
        return PartialGForce;
    }
    v4df Unit = Normalize(Difference);
    v4df Output = Unit*((MassLoc*MassTemp*G)/Distance2);
    return Output;
}




void emitPoints(unsigned int step, v4df* Positions, double *Masses, int NumBodies){
    static bool fileHeaderPrinted = false;
    static uint64_t pointEmitted = 0;
    
    double *DataBuf = calloc(max_jobs, sizeof(double));
    
    fflush(output);
    if(!fileHeaderPrinted){
            

        if(id == 0){
                
            fprintf(output,"Total particles %d\n", n_bodies);
                

        }
        fileHeaderPrinted = true;
            

    }
        

    if(id == 0){
        fprintf(output,"STEP=%lu\n", step);
    }
    char propertyString[5];
        

    for(int j = 0; j < 4; j++){
        switch(j){
            case 0:
                strcpy(propertyString,"X\0");
                break;
            case 1:
                strcpy(propertyString, "Y\0");
                break;
            case 2:
                strcpy(propertyString,"Z\0");
                break;
            case 3:
               strcpy(propertyString,"mass\0");
                break;
        }
        if(id == 0){
            fprintf(output,"%s%ld=[\n",propertyString, pointEmitted);
        }
        if(j != 3){
            for(int databuf_idx=0; databuf_idx < NumBodies; databuf_idx++){
                DataBuf[databuf_idx] = Positions[databuf_idx][j];
            }
            MPI_Gatherv(DataBuf, NumBodies, MPI_DOUBLE, EmitPointsBuf, jobAllocations, jobDispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(id == 0){
                for(int k = 0; k < n_bodies-1; k++){
                    fprintf(output,"%12.6e, ", EmitPointsBuf[k]);
                }
                fprintf(output,"%12.6e ];\n", EmitPointsBuf[n_bodies-1]);
            }
        }else{
            if(id == 0){
                for(int k = 0; k < n_bodies-1; k++){
                    fprintf(output,"%12.6e, ", Masses[k]);
                }
                fprintf(output,"%12.6e ];\n\n",Masses[n_bodies-1]);
                fflush(output);
            }
        }
    }
    pointEmitted++;
}




RawSimulationData parseBodyKinematicsFile(char* filepath){
    uint32_t num_bodies = 0;

    RawSimulationData returnData;
    returnData.NumBodies = -1;



    char str_buff[256];
    uint8_t buff_idx=0;
    FILE *fp_body_positions = fopen(filepath, "r");
    if(fp_body_positions == NULL){
            fprintf(stderr,"Could not open the body position file\n");
            return returnData;
    }

    char c;
    bool in_number = false;
    // Iterate over first line and extract number of bodies
    while(1){
        c = fgetc(fp_body_positions);
        if(isdigit(c)){
                in_number = true;
                str_buff[buff_idx] = c;
                buff_idx++;
        }else if(in_number){
                str_buff[buff_idx] = '\0';
                num_bodies = strtol(str_buff, NULL, 10);
                buff_idx = 0;
                break;
        }
    }
    //Move to nex line
    while(c != '\n'){
        c = fgetc(fp_body_positions);
    }


    //Array that is split across all processes
    // if vector was of length 3:
    //[Px, Py, Pz, Vx, Vy, Vz,  Px, Py, Pz, Vx, Vy, Vz, Px, Py, Pz, Vx, Vy, Vz]
    double* PosVelVec = (double*)malloc(sizeof(double)*num_bodies*6);
    memset((void*)PosVelVec,0.0 , sizeof(double)*num_bodies*6);
    double* Mass = (double*)malloc(sizeof(double)*num_bodies);
    memset((void*)Mass,0,sizeof(double)*num_bodies);



    // Since each value is layed out in the array Px,Py,Pz, Vx, Vy, Vz, Mass
    // This contains the offset of the current element within a body
    // e.g 1 corresponds to Py
    uint32_t offset = 0;

    //body_idx is the index of the body currently being parsed. 
    uint32_t body_idx = 0;

    char *endptr = str_buff;
    char *startptr = str_buff;


    in_number = false;
    while(fgets(&str_buff[buff_idx], 256 - buff_idx, fp_body_positions) != NULL){
        if(body_idx >= num_bodies){
            break;
        }
        int char_index = 0;
        for(char_index = 0; char_index < 256; char_index++){
            c = str_buff[char_index];
            if(c == '\0'){
                break;
            }
            if(isdigit(c) || c == '.' || c == 'E' ||c == '+' || c == '-'){
                    in_number = true;
            }else if(in_number){
                in_number = false;
                double val = strtod(startptr, &endptr);
                startptr = endptr;
                if(offset < 3){
                    PosVelVec[body_idx*6 + offset] = val;
                }else if (offset> 3){
                    PosVelVec[body_idx*6 + offset-1] = val;
                }else{
                    Mass[body_idx] = val;
                }
                offset++;
                offset = offset%7;
                if(offset == 0){
                    body_idx++;
                }
            }
        }
        if(endptr != &str_buff[255]){
            str_buff[0] = '\0';
            strcat(str_buff, endptr);
            buff_idx = strlen(str_buff);
            
        }else{
            buff_idx = 0;
        }
        endptr = str_buff;
        startptr = endptr;
    }

    //Parse the rest of the string buffer. This is the 
    PosVelVec[body_idx*6 + offset] = strtod(str_buff, NULL);

    //Set the return data values
    returnData.PosVelVec = PosVelVec;
    returnData.MassVec = Mass;
    returnData.NumBodies = num_bodies;
    return returnData;
}
