/*
 * Particle.cpp
 *
 *  Created on: Jan 12, 2019
 */

#include <algorithm>
#include <cstring>
#include <stdio.h>

#include "Particle.h"

using namespace std;

vector<Particle *> *Particle::particleList = new vector<Particle *>{};

Particle::Particle(double mass, double position[2], double velocity[2]) {
    this->mass = mass;
    memcpy(this->position, position, 2 * sizeof(double));
    memcpy(this->positionPre, position, 2 * sizeof(double));
    memcpy(this->velocity, velocity, 2 * sizeof(double));
}

Particle::~Particle() {}

vector<double> Particle::accelByDistance(vector<double> distance, double mass) {
    double absDeltaRSqr = distance[0] * distance[0] + distance[1] * distance[1];
    double massDivRSqr = mass * GRAVITAIONAL_CONSTANT / absDeltaRSqr;

    vector<double> result = {distance[0] * massDivRSqr,
                             distance[1] * massDivRSqr};
    return result;
}

vector<double> Particle::addVec(vector<double> v1, vector<double> v2) {
    vector<double> result = {v1[0] + v2[0], v1[1] + v2[1]};
    return result;
}

vector<double> Particle::addVec(vector<double> v1, vector<double> v2,
                                vector<double> v3) {
    vector<double> result = {v1[0] + v2[0] + v3[0], v1[1] + v2[1] + v3[1]};
    return result;
}

vector<double> Particle::multVec(vector<double> v1, double scalar) {
    vector<double> result = {v1[0] * scalar, v1[1] * scalar};
    return result;
}

vector<double> Particle::RungeKutta4(vector<double> distance, double mass) {
    vector<double> F1 = accelByDistance(distance, mass);
    vector<double> F2 =
        accelByDistance(addVec(distance, multVec(F1, .5)), mass);
    vector<double> F3 =
        accelByDistance(addVec(distance, multVec(F2, .5)), mass);
    vector<double> F4 = accelByDistance(addVec(distance, F3), mass);

    vector<double> result = {(F1[0] + 2. * (F2[0] + F3[0]) + F4[0]) / 6.,
                             (F1[1] + 2. * (F2[1] + F3[1]) + F4[1]) / 6.};
    return result;
}

vector<double> Particle::RungeKutta5(vector<double> distance, double mass) {
    vector<double> F1 = accelByDistance(distance, mass);
    vector<double> F2 =
        accelByDistance(addVec(distance, multVec(F1, .25)), mass);
    vector<double> F3 = accelByDistance(
        addVec(distance, multVec(F1, .125), multVec(F2, .125)), mass);
    vector<double> F4 =
        accelByDistance(addVec(distance, multVec(F2, -.5), F3), mass);
    vector<double> F5 = accelByDistance(
        addVec(distance, multVec(F1, 3. / 16.), multVec(F4, 9. / 16.)), mass);
    vector<double> F6 = accelByDistance(
        addVec(addVec(distance, multVec(F1, -3. / 7.), multVec(F2, 2. / 7.)),
               addVec(multVec(F3, 12. / 7.), multVec(F4, -12. / 7.),
                      multVec(F5, 8. / 7.))),
        mass);

    vector<double> result = {
        (7. * (F1[0] + F6[0]) + 32. * (F3[0] + F5[0]) + 12. * F4[0]) / 90.,
        (7. * (F1[1] + F6[1]) + 32. * (F3[1] + F5[1]) + 12. * F4[1]) / 90.};
    return result;
}

// after one second
void Particle::updateVelocity() {
    if (this->fixed)
        return;

    for (Particle *p : *Particle::particleList) {
        if (p == this)
            continue;

        vector<double> deltaR = {p->getPosition(0) - this->position[0],
                                 p->getPosition(1) - this->position[1]};
        double massP = p->getMass();

        printf("%f, %f\n", deltaR[0], deltaR[0]);

        vector<double> rk = RungeKutta5(deltaR, massP);

        this->velocity[0] += rk[0];
        this->velocity[1] += rk[1];
    }
}

// after one second
void Particle::updatePosition() {
    if (this->fixed)
        return;

    memcpy(this->positionPre, this->position, 2 * sizeof(double));

    this->position[0] += this->velocity[0];
    this->position[1] += this->velocity[1];
}

void Particle::setFixed(bool fixed) {
    this->fixed = fixed;
    if (fixed)
        memset(this->velocity, 0, 2 * sizeof(double));
}
