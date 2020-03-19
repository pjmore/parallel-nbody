#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <time.h>

#define G 6.67408E-11
#define Minute 60.0
#define Hour Minute*60.0
#define Day Hour*24

typedef double v4df __attribute__ ((vector_size (32)));

typedef struct prog_data{
    int NumBodies;
    double* Masses;
    v4df* Velocities;
    v4df* Positions;
}SimulationData;


typedef struct SimulationConfigurationParameters{
    double timestep;
    double emitPointsInterval;
    double endTime;
    int* IdxBodiesOfInterest;
    int NumberOfBodiesOfInterest;
    double emitPointsBodiesOfInterestInterval;

}SimConfig;

 


SimConfig defaultSimulationConfiguration(){
    SimConfig newConfig;
    newConfig.timestep = 0.1;
    newConfig.emitPointsInterval = 50000.0;
    newConfig.endTime = 231.0*Day;
    newConfig.NumberOfBodiesOfInterest = 0;
    newConfig.emitPointsBodiesOfInterestInterval = DBL_MAX;
    return newConfig;
}



 
SimulationData parseBodyKinematicsFile(char* filepath);


 
static inline v4df CalculatePositionDelta(int, double,int, v4df*, v4df*, double*);
static inline v4df CalculateAcceleration(int i, int numBodies, v4df *Positions, double *Masses);
static inline double Sumv4df(v4df vec);
static inline v4df Normalize(v4df vec);




void runSimulation(SimulationData *StartingParameters, SimConfig *ConfigurationParameters);
void emitBodiesOfInterestPoints(double time, v4df* Positions, int NumBodies, int* BodiesOfInterest, int NumBodiesOfInterest);



int main(int argc, char *arv[]){
    char temp_filename[50] = "./test.txt";
    v4df *Positions;
    v4df *Velocities;
    double *Masses;
    
    SimulationData startPos = parseBodyKinematicsFile(temp_filename);
    if (startPos.NumBodies < 1){
        fprintf(stderr, "There was an error parsing the kinematics file\n");
        return EXIT_FAILURE;
    }
    SimConfig SimulationParams;
    SimulationParams = defaultSimulationConfiguration();
    double timeSpentRunning = 0.0;
    clock_t start = clock();
    runSimulation(&startPos, &SimulationParams);
    clock_t end = clock();
    timeSpentRunning = (double)(end-start)/CLOCKS_PER_SEC;
    printf("It took %f seconds to run\n", timeSpentRunning);

}





static inline double Sumv4df(v4df vec){
    return vec[0]+vec[1]+vec[2];
}

static inline v4df Normalize(v4df vec){
    v4df sq = vec*vec;
    return vec/sqrt(Sumv4df(sq));
}

static inline v4df CalculateAcceleration(int i, int numBodies, v4df *Positions, double *Masses){
    v4df Acceleration = {0,0,0,0};
    for(int j = 0; j<numBodies; j++){
        //printf("Acceleration[%d] On iteration[%d]: x: %f y: %f z: %f \n",i,j, Acceleration[0], Acceleration[1], Acceleration[2]);
        if(j == i){
            continue;
        }
        v4df Difference =  Positions[j]-Positions[i];
        double Distance2 = Sumv4df(Difference*Difference);
        //printf("SquaredDistance[%d] On iteration[%d]: %f\n",i,j, Distance2);
        if (Distance2 <1000000.0){
            continue;
        }
        v4df AccDelta = Normalize(Difference)*((G*Masses[j])/Distance2);
        //printf("AccelerationDelta[%d] On iteration[%d]: x: %f y: %f z: %f \n",i,j, AccDelta[0], AccDelta[1], AccDelta[2]);
        Acceleration = Acceleration + AccDelta;
        //printf("Acceleration[%d] On iteration[%d]: x: %f y: %f z: %f \n",i,j, Acceleration[0], Acceleration[1], Acceleration[2]);
    }
    return Acceleration;
}

static inline v4df CalculatePartialGravity(int i, int j, v4df *Positions, double *Masses){
    v4df PartialGForce = {0};
    v4df Difference = Positions[j] - Positions[i];
    double Distance2 = Sumv4df(Difference*Difference);
    if(Distance2<1000000.0){
        return PartialGForce;
    }
    return (Difference/sqrt(Distance2))*(G/Distance2);
}


static inline v4df CalculatePositionDelta(int i, double timestep,int numBodies, v4df *Positions, v4df *Velocities, double *Masses){
    v4df Acc = CalculateAcceleration(i, numBodies, Positions, Masses);
    //printf("Acceleration[%d]: x: %f y: %f z: %f \n",i, Acc[0], Acc[1], Acc[2]);
    v4df V = Velocities[i];
    Velocities[i] += timestep * Acc;
    return  timestep*V + (timestep*timestep)*Acc;
}

static inline void UpdatePositionsByDualPositionDeltaCalculation(int i, double timestep, int NumBodies, v4df *Positions,v4df* Acceleration , v4df *Velocities, double *Masses){
    //static int UpdateNumber = 0;
    v4df PartialGForce = {0};
    double dt2 = timestep*timestep;
    if(i==0){
        for(int k = 0; k<NumBodies; k++){
            Acceleration[k] = (v4df){0};
        }
    }
    for(int j = i+1; j < NumBodies; j++){
        PartialGForce = CalculatePartialGravity(i,j, Positions, Masses);
        Acceleration[i] = Acceleration[i] + Masses[j]*PartialGForce;
        Acceleration[j] = Acceleration[j] - Masses[i]*PartialGForce;
    }
    Positions[i] += Velocities[i]*timestep + Acceleration[i]*dt2;
    Velocities[i] += Acceleration[i]*timestep;
}

void emitPoints(double time, v4df* Positions, int NumBodies){
    static bool fileHeaderPrinted = false;
    static uint64_t pointEmitted = 0;
    if(!fileHeaderPrinted){
        printf("Total particles %d\n", NumBodies);
        fileHeaderPrinted = true;
    }
    printf("TIME%ld=%f\n", pointEmitted, time);
    char propertyString[5];
    for(int j = 0; j < 3; j++){
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
        printf("%s%ld=[\n",propertyString, pointEmitted);
        for(int i = 0; i<NumBodies; i++){
            printf("%f", Positions[i][j]);
            if(i != NumBodies-1){
                printf(", ");
            }
        }
        printf(" ];\n");
    }
    pointEmitted++;
    printf("\n");
}


void runSimulation(SimulationData *StartingParameters, SimConfig *ConfigurationParameters){
    double current_timestamp = 0.0;
    double dt = ConfigurationParameters->timestep;
    double dtSq = dt*dt;

    double pointsInterval = ConfigurationParameters->emitPointsInterval;
    double pointsIntervalTracker = 0.0;

    double endTime = ConfigurationParameters->endTime;
    int NumBodies = StartingParameters->NumBodies;
    v4df *Positions = StartingParameters->Positions;
    v4df *Velocities = StartingParameters->Velocities;
    double *Masses = StartingParameters->Masses;
    //v4df *PositionDelta = aligned_alloc(32, NumBodies*sizeof(v4df));
    v4df *Acceleration = aligned_alloc(32, NumBodies*sizeof(v4df));
    emitPoints(current_timestamp, Positions, NumBodies);
    while(current_timestamp < endTime){
        /*
        for(int i = 0; i<NumBodies; i++){
            PositionDelta[i] = CalculatePositionDelta(i, dt, NumBodies, Positions, Velocities, Masses);
        }
        for(int i = 0; i<NumBodies; i++){
            Positions[i] = Positions[i] + PositionDelta[i];
        }
        */
        for(int i =0; i<NumBodies; i++){
            UpdatePositionsByDualPositionDeltaCalculation(i, dt, NumBodies, Positions, Acceleration , Velocities, Masses);
        }
        pointsIntervalTracker += dt;
        if(pointsInterval <= pointsIntervalTracker){
            pointsIntervalTracker = dt;
            emitPoints(current_timestamp, Positions, NumBodies);
        }
        current_timestamp+=dt;
    }
    emitPoints(current_timestamp, Positions, NumBodies);
}



SimulationData parseBodyKinematicsFile(char* filepath){
    int num_bodies = -1;
    SimulationData returnData;
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
    while(c != '\n'){
        c = fgetc(fp_body_positions);
    }
    //v4df* pArr = (v4df*)malloc(sizeof(v4df)*num_bodies);
    //v4df* vArr = (v4df*)malloc(sizeof(v4df)*num_bodies);
    //double* mArr = (double*)malloc(sizeof(double)*num_bodies);

    v4df* pArr = (v4df*)aligned_alloc(32,sizeof(v4df)*num_bodies);
    v4df* vArr = (v4df*)aligned_alloc(32,sizeof(v4df)*num_bodies);
    double* mArr = (double*)malloc(sizeof(double)*num_bodies);


    uint32_t offset = 0;
    uint32_t body_idx = 0;
    char *endptr = str_buff;
    char *startptr = str_buff;
    in_number = false;
    while(fgets(&str_buff[buff_idx], 256 - buff_idx, fp_body_positions) != NULL){
        if(body_idx >= num_bodies){
            printf("there aren't that many bodies");
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
                if (offset < 3){
                    pArr[body_idx][offset] = val;
                }else if(offset == 3){
                    mArr[body_idx] = val;
                }else{
                    vArr[body_idx][offset-4] = val;
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
    vArr[body_idx][2] = strtod(str_buff, NULL);
    returnData.Positions = pArr;
    returnData.Velocities = vArr;
    returnData.Masses = mArr;
    returnData.NumBodies = num_bodies;
    return returnData;
}





 
void emitBodiesOfInterestPoints(double time, v4df* Positions, int NumBodies, int* BodiesOfInterest, int NumBodiesOfInterest){
    if(NumBodiesOfInterest == 0){
        return;
    }
    printf("time - %f\n", time);
    int idx = 0;
    for(int i = 0; i < NumBodiesOfInterest; i++){
        idx = BodiesOfInterest[i];
        if(idx >= NumBodies){
            continue;
        }
        printf("%d %f %f %f\n", idx, Positions[idx][0], Positions[idx][1], Positions[idx][2]);
    }
}
