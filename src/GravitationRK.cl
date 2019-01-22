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
} p_state;

__constant double GRAVITAIONAL_CONSTANT = 6.67408E-11; // [m^3 / (kg * s^2)]
__constant double TIME_STEP = 1.; // [s]

double2 accelByDistance(double2 distance, double mass) {
    double2 distSqr = distance * distance;
    double absDeltaRSqr = distSqr.x + distSqr.y;
    double massDivRSqr = mass * GRAVITAIONAL_CONSTANT / absDeltaRSqr;

    double2 result = distance * massDivRSqr;
    return result;
}

double2 RungeKutta4(double2 distance, double mass) {
    double2 F1 = accelByDistance(distance, mass);
    double2 F2 = accelByDistance(distance + F1 * .5, mass);
    double2 F3 = accelByDistance(distance + F2 * .5, mass);
    double2 F4 = accelByDistance(distance + F3, mass);

    double2 result = (F1 + 2. * (F2 + F3) + F4) / 6.;
    return result;
}

double2 RungeKutta5(double2 distance, double mass) {
    double2 one = {1., 1.};
    
    double2 F1 = accelByDistance(distance, mass);
    double2 F2 = accelByDistance(distance + F1 * .25, mass);
    double2 F3 = accelByDistance(distance + F1 * .125 + F2 * .125, mass);
    double2 F4 = accelByDistance(distance + F2 * (-.5) + F3, mass);
    double2 F5 = accelByDistance(distance + F1 * (3. / 16.) + F4 * (9. / 16.), mass);
    double2 F6 = accelByDistance(distance + F1 * (-3. / 7.) + F2 * (2. / 7.) + F3 * (12. / 7.) + F4 * (-12. / 7.) + F5 * (8. / 7.), mass);

    double2 result = (7. * (F1 + F6) + 32. * (F3 + F5) + 12. * F4) / 90.;
    return result;
}

__kernel void GravitationRK(__global p_state *statesIn, __global p_state *statesOut, const long size, const int rounds) {
    int index = get_global_id(0);
    __global p_state *statesPointers[2] = {statesIn, statesOut};

    int swapCount = (rounds == 0) ? 1 : 2;
    int roundsCnt = (rounds == 0) ? 1 : rounds;

    int r;
    for(r = 0; r < roundsCnt; r++) {
        int r_r;
        for(r_r = 0; r_r < swapCount; r_r++) {
            double2 vel = statesPointers[0][index].vel;
            double2 pos = statesPointers[0][index].pos;

            double2 finalVel = vel;
            int i;
            for(i = 0; i < size; i++) {
                if (i == index) continue;

                double2 deltaR = statesPointers[0][i].pos - pos;
                double2 rk = RungeKutta5(deltaR, statesIn[i].mass);

                finalVel += TIME_STEP * rk;
            }

            statesPointers[1][index].pos = pos + TIME_STEP * finalVel;
            statesPointers[1][index].vel = finalVel;
            statesPointers[1][index].mass = statesPointers[0][index].mass;
            
            barrier(CLK_GLOBAL_MEM_FENCE); //  barrier to sync threads
            // change pointers
            __global p_state *temp = statesPointers[0];
            statesPointers[0] = statesPointers[1];
            statesPointers[1] = temp;
            barrier(CLK_GLOBAL_MEM_FENCE); //  barrier to sync threads
        }
    }
}
