typedef struct {
    double2 pos;
    double2 vel;
    double mass;
} p_state;

__constant double GRAVITAIONAL_CONSTANT = 6.67408E-11; // [m^3 / (kg * s^2)]
__constant double TIME_STEP = 1.; // [s]

double2 accelByDistance(double2 distance, double mass) {
    double2 distSqr = distance;//mad(distance, distance, one);
    double absDeltaRSqr = distSqr.x * distSqr.x + distSqr.y * distSqr.y;
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

__kernel void GravitationRK(__global p_state *statesIn, __global p_state *statesOut, const int size) {
    int index = get_global_id(0);
    double2 vel = statesIn[index].vel;
    double2 pos = statesIn[index].pos;
    
    double2 finalVel = vel;
    int i;
    for(i = 0; i < size; i++) {
        if (i == index) continue;
        
        double2 deltaR = statesIn[i].pos - pos;
        double2 rk = RungeKutta5(deltaR, statesIn[i].mass);

        finalVel += TIME_STEP * rk;
    }

    statesOut[index].pos = pos + TIME_STEP * finalVel;
    statesOut[index].vel = finalVel;
    statesOut[index].mass = statesIn[index].mass;
}
