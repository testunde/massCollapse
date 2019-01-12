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
    double position[2] = {0., 0.};    // ([m]; [m])
    double positionPre[2] = {0., 0.}; // ([m]; [m])
    double velocity[2] = {0., 0.};    // ([m]; [m])
    double mass = 0.;                 // [kg]
    bool fixed = false;

    static std::vector<double> accelByDistance(std::vector<double> distance,
                                               double mass);
    static std::vector<double> addVec(std::vector<double> v1,
                                      std::vector<double> v2);
    static std::vector<double> multVec(std::vector<double> v1, double scalar);

  public:
    static std::vector<Particle *> *particleList;

    Particle(double mass, double position[], double velocity[]);
    virtual ~Particle();

    void updateVelocity();
    void updatePosition();

    double getPosition(int axis) const { return this->position[axis]; }
    double getVelocity(int axis) const { return this->velocity[axis]; }
    double getMass() const { return this->mass; }

    void setFixed(bool fixed);
};

#endif /* PARTICLE_H_ */
