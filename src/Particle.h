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
    double radius;                 // [m]
    double coord[2] = {0., 0.};    // ([m]; [m])
    double coordPre[2] = {0., 0.}; // ([m]; [m])
    double mass;                   // [kg]
    double velocity;               // [m/s]
    Particle *mergedInto = nullptr;

    double updateMass();
    double updateVelocity();

  public:
    static int disposedDrops;
    static std::vector<Particle *> *dropList;

    static void sortListWidth();
    static void sortListHeight();
    static void sortListSize();

    static int remainingDrops();

    static Particle *bigger_h, *smaller_h; // head, tail
    static Particle *above_h, *below_h;    // head, tail
    static Particle *left_h, *right_h;     // head, tail

    Particle *bigger = nullptr, *smaller = nullptr; // head, tail
    Particle *above = nullptr, *below = nullptr;    // head, tail
    Particle *left = nullptr, *right = nullptr;     // head, tail

    void setAbove(Particle *drp);
    void setBelow(Particle *drp);

    Particle(double radius, double coord[]);
    virtual ~Particle();
    void deleteInstance();

    double getCoord(int axis) const { return this->coord[axis]; }
    double getCoordPre(int axis) const { return this->coordPre[axis]; }

    double getRadius() const { return this->radius; }

    double getMass() const { return this->mass; }
    double getVelocity() const { return this->velocity; }

    Particle *getFinalMergred();

    double growCondensation();
    void merge(std::vector<Particle *> *list);
    void fallBy(double way);
    void setMergedInto(Particle *drp);
};

#endif /* PARTICLE_H_ */
