/*
 * Particle.h
 *
 *  Created on: Jan 12, 2019
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "math.h"
#include <vector>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp> // for struct definition

#include "Global.h"

typedef struct {
    cl_double2 pos;
    cl_double2 vel;
    cl_double mass;
} p_state;

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
                                           double mass);
    static std::vector<double> RungeKutta4(std::vector<double> distance,
                                           double mass);
    static std::vector<double> RungeKutta5(std::vector<double> distance,
                                           double mass);

  public:
    static std::vector<Particle *> *particleList;

    Particle(double mass, double position[], double velocity[]);
    virtual ~Particle();

    void updateVelocity();
    void updatePosition();
    void setCLStruct(p_state *st);

    double getPosition(int axis) const { return this->position[axis]; }
    double getVelocity(int axis) const { return this->velocity[axis]; }
    double getMass() const { return this->mass; }
    p_state getCLStruct() const {
        p_state st;
        st.pos = {this->position[0], this->position[1]};
        st.vel = {this->velocity[0], this->velocity[1]};
        st.mass = this->mass;
        return st;
    };

    void setFixed(bool fixed);
};

#endif /* PARTICLE_H_ */
