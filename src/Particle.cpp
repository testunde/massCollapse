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
    memcpy(this->velocity, velocity, 2 * sizeof(double));
}

Particle::~Particle() {}

// after one second
void Particle::updateVelocity() {}

// after one second
void Particle::updatePosition() {}
