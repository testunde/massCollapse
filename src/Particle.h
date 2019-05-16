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
    cl_double radius;
    cl_bool collision;
    cl_int cluster_nr;
} p_state;

class Particle {
  private:
    static constexpr double radiusPreFactor =
        3. / (ENVIRONMENT_SPAWN_PARTICLE_DENSITY * M_PI * 4.);
    double position[2] = {0., 0.};    // ([m]; [m])
    double positionPre[2] = {0., 0.}; // ([m]; [m])
    double velocity[2] = {0., 0.};    // ([m]; [m])
    double mass = 0.;                 // [kg]
    double radius = 0.;               // [m]
    bool fixed = false;
    bool collision = false;
    int cluster_nr = 0;

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

    double calcRadius();

  public:
    static std::vector<Particle *> *particleList;

    Particle(double mass, double position[], double velocity[], int cluster_nr);
    virtual ~Particle();

    void updateVelocity();
    void updatePosition();
    void setVelocity(std::vector<double> newVel) {
        this->velocity[0] = newVel[0];
        this->velocity[1] = newVel[1];
    }
    void setCLStruct(p_state *st);

    double getPosition(int axis) const { return this->position[axis]; }
    double getVelocity(int axis) const { return this->velocity[axis]; }
    double getMass() const { return this->mass; }
    double getRadius() const { return this->radius; }
    bool getCollission() const { return this->collision; }
    std::vector<double> const getcurrentGravForce();
    p_state getCLStruct() const {
        p_state st;
        st.pos = {this->position[0], this->position[1]};
        st.vel = {this->velocity[0], this->velocity[1]};
        st.mass = this->mass;
        st.radius = this->radius;
        st.collision = this->collision;
        st.cluster_nr = this->cluster_nr;
        return st;
    };

    void setFixed(bool fixed);
};

#endif /* PARTICLE_H_ */
