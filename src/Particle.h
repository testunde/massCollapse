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
    static std::vector<double> addVec(std::vector<double> v1,
                                      std::vector<double> v2,
                                      std::vector<double> v3);
    static std::vector<double> multVec(std::vector<double> v1, double scalar);
    static std::vector<double> RungeKutta1(std::vector<double> distance,
                                           double mass, double timestep);
    static std::vector<double> RungeKutta4(std::vector<double> distance,
                                           double mass, double timestep);
    static std::vector<double> RungeKutta5(std::vector<double> distance,
                                           double mass, double timestep);

  public:
    static std::vector<Particle *> *particleList;

    Particle(double mass, double position[], double velocity[]);
    virtual ~Particle();

    void updateVelocity();
    void updatePosition();
    void setVelocity(std::vector<double> newVel) {
        this->velocity[0] = newVel[0];
        this->velocity[1] = newVel[1];
    }

    double getPosition(int axis) const { return this->position[axis]; }
    double getVelocity(int axis) const { return this->velocity[axis]; }
    double getMass() const { return this->mass; }
    std::vector<double> const getcurrentGravForce();

    void setFixed(bool fixed);
};

#endif /* PARTICLE_H_ */
