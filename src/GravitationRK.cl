#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point not supported by OpenCL implementation."
#endif

typedef struct {
    double2 pos;
    double2 vel;
    double mass;
    bool collision;
} p_state;

double2 accelByDistance(double2 distance, double mass) {
    double2 distSqr = distance * distance;
    double absDeltaRSqr = distSqr.x + distSqr.y;
    double massDivRSqr = mass * GRAVITAIONAL_CONSTANT / absDeltaRSqr;

    return normalize(distance) * massDivRSqr;
}

double2 RungeKutta4(double2 distance, double mass, double timestep) {
    double2 F1 = timestep * accelByDistance(distance, mass);
    double2 F2 = timestep * accelByDistance(distance + F1 * .5, mass);
    double2 F3 = timestep * accelByDistance(distance + F2 * .5, mass);
    double2 F4 = timestep * accelByDistance(distance + F3, mass);

    return (F1 + 2. * (F2 + F3) + F4) / 6.;
}

double2 RungeKutta5(double2 distance, double mass, double timestep) {
    double2 F1 = timestep * accelByDistance(distance, mass);
    double2 F2 = timestep * accelByDistance(distance + F1 * .25, mass);
    double2 F3 = timestep * accelByDistance(distance + F1 * .125 + F2 * .125, mass);
    double2 F4 = timestep * accelByDistance(distance + F2 * (-.5) + F3, mass);
    double2 F5 = timestep * accelByDistance(distance + F1 * (3. / 16.) + F4 * (9. / 16.), mass);
    double2 F6 = timestep * accelByDistance(distance + F1 * (-3. / 7.) + F2 * (2. / 7.) + F3 * (12. / 7.) + F4 * (-12. / 7.) + F5 * (8. / 7.), mass);

    return (7. * (F1 + F6) + 32. * (F3 + F5) + 12. * F4) / 90.;
}

__kernel void GravitationRK(__global p_state *statesIn, __global p_state *statesOut, const int roundOdd_c) {
    int index = get_global_id(0);
    __global p_state *statesPointers[2] = {(roundOdd_c == 1) ? statesOut : statesIn, (roundOdd_c == 1) ? statesIn : statesOut};

    if (statesPointers[0][index].collision) {
        // just skip the calculations if already collided, only copy to output for consistency
        statesPointers[1][index] = statesPointers[0][index];
        return;
    }

    int r;
    for(r = 0; r < SIMULATION_ROUNDS; r++) {
        double2 vel = statesPointers[0][index].vel;
        double2 pos = statesPointers[0][index].pos;

        double2 finalVel = vel;
        int i;
        for(i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++) {
            if (i == index) continue;

            double2 deltaR = statesPointers[0][i].pos - pos;
            double2 rk = RungeKutta5(deltaR, statesPointers[0][i].mass, SIMULATION_TIME_PER_STEP);

            finalVel += rk;
        }

        statesPointers[1][index].pos = pos + SIMULATION_TIME_PER_STEP * finalVel;
        statesPointers[1][index].vel = finalVel;
        statesPointers[1][index].mass = statesPointers[0][index].mass;

        // change pointers
        __global p_state *temp = statesPointers[0];
        statesPointers[0] = statesPointers[1];
        statesPointers[1] = temp;
        barrier(CLK_GLOBAL_MEM_FENCE); //  barrier to sync threads
    }
}

__kernel void InitialVelocity(__global p_state *statesIn, __global p_state *statesOut) {
    int index = get_global_id(0);
    
    // for consistency, first copy all properties so unchanged ones are still present in the output
    statesOut[index] = statesIn[index];
    
    double coordNorm = length(statesIn[index].pos);
    
    if (ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION == 2) { // RungeKutta5 (after all points are generated)
        double2 gravForce = {.0, .0};
        int i;
        for(i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++) {
            if (i == index) continue;
            
            double2 deltaR = statesIn[i].pos - statesIn[index].pos;
            
            double2 rk = RungeKutta5(deltaR, statesIn[i].mass, 1.0);
            
            gravForce += rk;
        }
        double gravNorm = length(gravForce);
        double absVel = sqrt(coordNorm * gravNorm);
        
        statesOut[index].vel.x = absVel * (-gravForce.y) / gravNorm;
        statesOut[index].vel.y = absVel * (+gravForce.x) / gravNorm;
    } else if(ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION == 1) { // circular orbit sqrt(GM*radius)*2
        double absVel = sqrt(GRAVITAIONAL_CONSTANT * (ENVIRONMENT_SPAWN_PARTICLES_TOTAL *
                                   ENVIRONMENT_SPAWN_PARTICLE_MASS) * coordNorm) * 2;
        
        statesOut[index].vel.x = absVel * (+statesIn[index].pos.y) / coordNorm;
        statesOut[index].vel.y = absVel * (-statesIn[index].pos.x) / coordNorm;
    } else { // default: plain angular function (alpha * radius)
        double absVel = ENVIRONMENT_SPAWN_START_ANGULAR_VELO * coordNorm;
        
        statesOut[index].vel.x = absVel * (+statesIn[index].pos.y) / coordNorm;
        statesOut[index].vel.y = absVel * (-statesIn[index].pos.x) / coordNorm;
    }
}
