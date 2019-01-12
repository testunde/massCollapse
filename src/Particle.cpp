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

    return vector<double>(distance[0] * massDivRSqr, distance[1] * massDivRSqr);
}

// after one second
void Particle::updateVelocity() {
    if (this->fixed)
        return;

    for (Particle *p : *Particle::particleList) {
        if (p == this)
            continue;

        double deltaR[2] = {p->getPosition(0) - this->position[0],
                            p->getPosition(1) - this->position[1]};
        double absDeltaRSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1];
        double massDivRSqr =
            p->getMass() * GRAVITAIONAL_CONSTANT / absDeltaRSqr;

        this->velocity[0] += massDivRSqr * deltaR[0];
        this->velocity[1] += massDivRSqr * deltaR[1];
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
