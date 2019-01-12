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

int Particle::disposedDrops = 0;
vector<Particle *> *Particle::dropList = new vector<Particle *>{};

Particle *Particle::bigger_h = nullptr;  // head
Particle *Particle::smaller_h = nullptr; // tail
Particle *Particle::above_h = nullptr;   // head
Particle *Particle::below_h = nullptr;   // tail
Particle *Particle::left_h = nullptr;    // head
Particle *Particle::right_h = nullptr;   // tail

Particle::Particle(double radius, double coord[2]) {
    this->radius = radius;
    memcpy(this->coord, coord, 2 * sizeof(double));
    memcpy(this->coordPre, coord, 2 * sizeof(double));

    updateMass();
    updateVelocity();
}

Particle::~Particle() {}

// CALL IT ONCE AND LAST, IF CREATED WITH "new"!!!
void Particle::deleteInstance() {
    // update list heads
    if (bigger_h == this)
        bigger_h = this->smaller;
    if (smaller_h == this)
        smaller_h = this->bigger;

    if (above_h == this)
        above_h = this->below;
    if (below_h == this)
        below_h = this->above;

    if (left_h == this)
        left_h = this->right;
    if (right_h == this)
        right_h = this->left;

    // update list linkings
    if (bigger != nullptr)
        bigger->smaller = this->smaller;
    if (smaller != nullptr)
        smaller->bigger = this->bigger;

    if (above != nullptr)
        above->setBelow(this->below);
    if (below != nullptr)
        below->setAbove(this->above);

    if (left != nullptr)
        left->right = this->right;
    if (right != nullptr)
        right->left = this->left;

    dropList->erase(remove(dropList->begin(), dropList->end(), this));

    disposedDrops++;

    delete this;
}

void Particle::sortListWidth() {
    sort(dropList->begin(), dropList->end(), [](Particle *a, Particle *b) {
        return a->getCoord(0) < b->getCoord(0);
    });
    Particle *preDrop = nullptr;
    for (Particle *d : *dropList) {
        d->left = preDrop;
        if (preDrop != nullptr)
            preDrop->right = d;
        else
            Particle::left_h = d;
        preDrop = d;
    }
    if (preDrop != nullptr)
        preDrop->right = nullptr;
    Particle::right_h = preDrop;
}

void Particle::sortListHeight() {
    sort(dropList->begin(), dropList->end(), [](Particle *a, Particle *b) {
        return a->getCoord(1) < b->getCoord(1);
    });
    Particle *preDrop = nullptr;
    for (Particle *d : *dropList) {
        d->setAbove(preDrop);
        if (preDrop != nullptr)
            preDrop->setBelow(d);
        else
            Particle::above_h = d;
        preDrop = d;
    }
    if (preDrop != nullptr)
        preDrop->below = nullptr;
    Particle::below_h = preDrop;
}

void Particle::sortListSize() {
    sort(dropList->begin(), dropList->end(), [](Particle *a, Particle *b) {
        return a->getRadius() > b->getRadius();
    });
    Particle *preDrop = nullptr;
    for (Particle *d : *dropList) {
        d->bigger = preDrop;
        if (preDrop != nullptr)
            preDrop->smaller = d;
        else
            Particle::bigger_h = d;
        preDrop = d;
    }
    if (preDrop != nullptr)
        preDrop->smaller = nullptr;
    Particle::smaller_h = preDrop;
}

int Particle::remainingDrops() {
    return ENVIRONMENT_SPAWN_DROPS_TOTAL - disposedDrops;
}

double Particle::updateMass() {
    this->mass = DENSITY_WATER * (this->radius * this->radius * this->radius) *
                 (M_PI * 4. / 3.);
    return this->mass;
}

// after one second
double Particle::updateVelocity() {
    double tempD =
        2. * (this->radius > 5.E-3
                  ? 5.E-3
                  : this->radius); // no change if greater than 10 [mm] diameter

    // input: [mm], output: [cm/2]
    tempD *= 1.E+3; // [m] -> [mm]
    double tmpLog = log(tempD) - 2.4;
    double tempVel = exp(tmpLog * tmpLog * (-3. / 18.) + 6.9);
    this->velocity = tempVel * 1.E-2; // [cm/s] -> [m/s]

    return this->velocity;
}

// after one second
double Particle::growCondensation() {
    double seedRadius = max(this->radius, DROP_GROWTH_CONDENSATION_LIMIT_SIZE);

    double dR = (DROP_GROWTH_CONDENSATION_S - 1.) /
                (seedRadius * DROP_GROWTH_CONDENSATION_F);

    this->radius += dR;

    updateMass();
    updateVelocity();

    // No need to check for bigger neighbors if growing begins from the biggest
    // drop in the environment. Only check bigger neighbor(s) when drop radius
    // smaller than 0.31623 [µm] (~3E-7 [m]), because until then the grow rate
    // is bigger than the delta-radius would be.

    return dR;
}

void Particle::merge(vector<Particle *> *list) {
    if (list->empty())
        return;

    double totalMass = this->mass;

    double newCoord[] = {this->coord[0] * this->mass,
                         this->coord[1]}; // no height-change as approximation
    for (Particle *drp : *list) {
        totalMass += drp->getMass();
        newCoord[0] += drp->getCoord(0) * drp->getMass();

        drp->setMergedInto(this);
    }
    newCoord[0] /= totalMass;

    memcpy(this->coord, newCoord, 2 * sizeof(double));
    this->mass = totalMass;
    this->radius =
        pow(this->mass / (DENSITY_WATER * (M_PI * 4. / 3.)), 1. / 3.);
    double tempVel = this->velocity;
    updateVelocity();
    if (this->velocity < tempVel)
        printf("lower vel?!\n");

    for (Particle *drp : *list) {
        drp->deleteInstance(); // delete now for later faster&easier link
                               // updates
    }
}

void Particle::fallBy(double way) {
    memcpy(this->coordPre, this->coord, 2 * sizeof(double));
    this->coord[1] += way; // [m]
}

void Particle::setMergedInto(Particle *drp) { this->mergedInto = drp; }

Particle *Particle::getFinalMergred() {
    if (this->mergedInto != nullptr)
        return this->mergedInto->getFinalMergred();
    else
        return this;
}

void Particle::setAbove(Particle *drp) {
    if (this == drp)
        return;
    this->above = drp;
}

void Particle::setBelow(Particle *drp) {
    if (this == drp)
        return;
    this->below = drp;
}
