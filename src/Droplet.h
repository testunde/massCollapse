/*
 * Droplet.h
 *
 *  Created on: Dec 16, 2018
 */

#ifndef DROPLET_H_
#define DROPLET_H_

#include "math.h"
#include <list>
#include <vector>

#include "Global.h"

class Droplet {
  private:
    double radius;                 // [m]
    double coord[2] = {0., 0.};    // ([m]; [m])
    double coordPre[2] = {0., 0.}; // ([m]; [m])
    double mass;                   // [kg]
    double velocity;               // [m/s]
    Droplet *mergedInto = nullptr;

    double updateMass();
    double updateVelocity();

  public:
    static int disposedDrops;
    static std::vector<Droplet *> *dropList;

    static void sortListWidth();
    static void sortListHeight();
    static void sortListSize();

    static int remainingDrops();

    static Droplet *bigger_h, *smaller_h; // head, tail
    static Droplet *above_h, *below_h;    // head, tail
    static Droplet *left_h, *right_h;     // head, tail

    Droplet *bigger = nullptr, *smaller = nullptr; // head, tail
    Droplet *above = nullptr, *below = nullptr;    // head, tail
    Droplet *left = nullptr, *right = nullptr;     // head, tail

    void setAbove(Droplet *drp);
    void setBelow(Droplet *drp);

    Droplet(double radius, double coord[]);
    virtual ~Droplet();
    void deleteInstance();

    double getCoord(int axis) const { return this->coord[axis]; }
    double getCoordPre(int axis) const { return this->coordPre[axis]; }

    double getRadius() const { return this->radius; }

    double getMass() const { return this->mass; }
    double getVelocity() const { return this->velocity; }

    Droplet *getFinalMergred();

    double growCondensation();
    void merge(std::list<Droplet *> *list);
    void fallBy(double way);
    void setMergedInto(Droplet *drp);
};

#endif /* DROPLET_H_ */
