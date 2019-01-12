/*
 * Particle.h
 *
 *  Created on: Jan 12, 2019
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "math.h"
#include <vector>

#include "Global.h"

class Particle {
  private:
    double position[2] = {0., 0.}; // ([m]; [m])
    double velocity[2] = {0., 0.}; // ([m]; [m])
    double mass;                   // [kg]

    void updateVelocity();
    void updatePosition();

  public:
    static std::vector<Particle *> *particleList;

    Particle(double mass, double position[], double velocity[]);
    virtual ~Particle();

    double getPosition(int axis) const { return this->position[axis]; }
    double getVelocity(int axis) const { return this->velocity[axis]; }
    double getMass() const { return this->mass; }
};

#endif /* PARTICLE_H_ */
