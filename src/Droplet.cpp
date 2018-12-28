/*
 * Droplet.cpp
 *
 *  Created on: Dec 16, 2018
 */

#include <cstring>
#include <stdio.h>

#include "Droplet.h"

using namespace std;

int Droplet::disposedDrops = 0;
list<Droplet*> * Droplet::dropList = new list<Droplet*>();

Droplet *Droplet::bigger_h = nullptr; // head
Droplet *Droplet::smaller_h = nullptr; // tail
Droplet *Droplet::above_h = nullptr; // head
Droplet *Droplet::below_h = nullptr; // tail
Droplet *Droplet::left_h = nullptr; // head
Droplet *Droplet::right_h = nullptr; // tail

Droplet::Droplet(double radius, double coord[2]) {
	this->radius = radius;
	memcpy(this->coord, coord, 2 * sizeof(double));
	memcpy(this->coordPre, coord, 2 * sizeof(double));

	updateMass();
	updateVelocity();
}

Droplet::~Droplet() {
}

// CALL IT ONCE AND LAST, IF CREATED WITH "new"!!!
void Droplet::deleteInstance() {
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

	dropList->remove(this);

	disposedDrops++;

	delete this;
}

void Droplet::sortListWidth() {
	dropList->sort([](Droplet *a, Droplet *b) {
		return a->getCoord(0) < b->getCoord(0);
	});
	Droplet *preDrop = nullptr;
	for (Droplet *d : *dropList) {
		d->left = preDrop;
		if (preDrop != nullptr)
			preDrop->right = d;
		else
			Droplet::left_h = d;
		preDrop = d;
	}
	if (preDrop != nullptr)
		preDrop->right = nullptr;
	Droplet::right_h = preDrop;
}

void Droplet::sortListHeight() {
	dropList->sort([](Droplet *a, Droplet *b) {
		return a->getCoord(1) < b->getCoord(1);
	});
	Droplet *preDrop = nullptr;
	for (Droplet *d : *dropList) {
		d->setAbove(preDrop);
		if (preDrop != nullptr)
			preDrop->setBelow(d);
		else
			Droplet::above_h = d;
		preDrop = d;
	}
	if (preDrop != nullptr)
		preDrop->below = nullptr;
	Droplet::below_h = preDrop;
}

void Droplet::sortListSize() {
	dropList->sort([](Droplet *a, Droplet *b) {
		return a->getRadius() > b->getRadius();
	});
	Droplet *preDrop = nullptr;
	for (Droplet *d : *dropList) {
		d->bigger = preDrop;
		if (preDrop != nullptr)
			preDrop->smaller = d;
		else
			Droplet::bigger_h = d;
		preDrop = d;
	}
	if (preDrop != nullptr)
		preDrop->smaller = nullptr;
	Droplet::smaller_h = preDrop;
}

int Droplet::remainingDrops() {
	return ENVIRONMENT_SPAWN_DROPS_TOTAL - disposedDrops;
}

double Droplet::updateMass() {
	this->mass = DENSITY_WATER * pow(this->radius, 3.) * (M_PI * 4. / 3.);
	return this->mass;
}

// after one second
double Droplet::updateVelocity() {
	double tempD = 2. * (this->radius > 5.E-3 ? 5.E-3 : this->radius); // no change if greater than 10 [mm] diameter

	// input: [mm], output: [cm/2]
	tempD *= 1.E+3; // [m] -> [mm]
	double tempVel = exp(pow(log(tempD) - 2.4, 2.) * (-3. / 18.) + 6.9);
	this->velocity = tempVel * 1.E-2; // [cm/s] -> [m/s]

	return this->velocity;
}

// after one second
double Droplet::growCondensation() {
	double seedRadius =
			(this->radius < DROP_GROWTH_CONDENSATION_LIMIT_SIZE) ? DROP_GROWTH_CONDENSATION_LIMIT_SIZE : this->radius;
	double dR = (DROP_GROWTH_CONDENSATION_S - 1.) / (seedRadius * DROP_GROWTH_CONDENSATION_F);

	this->radius += dR;

	updateMass();
	updateVelocity();

	// No need to check for bigger neighbors if growing begins from the biggest drop in the environment.
	// Only check bigger neighbor(s) when drop radius smaller than 0.31623 [µm] (~3E-7 [m]),
	// because until then the grow rate is bigger than the delta-radius would be.

	return dR;
}

void Droplet::merge(list<Droplet *> *list) {
	if (list->empty())
		return;

	double totalMass = this->mass;

	double newCoord[] = { this->coord[0] * this->mass, this->coord[1] }; // no height-change as approximation
	for (Droplet * drp : *list) {
		totalMass += drp->getMass();
		newCoord[0] += drp->getCoord(0) * drp->getMass();

		drp->setMergedInto(this);
	}
	newCoord[0] /= totalMass;

	memcpy(this->coord, newCoord, 2 * sizeof(double));
	this->mass = totalMass;
	this->radius = pow(this->mass / (DENSITY_WATER * (M_PI * 4. / 3.)), 1. / 3.);
	double tempVel = this->velocity;
	updateVelocity();
	if (this->velocity < tempVel)
		printf("lower vel?!\n");

	for (Droplet * drp : *list) {
		drp->deleteInstance(); // delete now for later faster&easier link updates
	}
}

void Droplet::fallBy(double way) {
	memcpy(this->coordPre, this->coord, 2 * sizeof(double));
	this->coord[1] += way; // [m]
}

void Droplet::setMergedInto(Droplet *drp) {
	this->mergedInto = drp;
}

double Droplet::getRadius() {
	return this->radius;
}

double Droplet::getCoord(int axis) {
	return this->coord[axis];
}

double Droplet::getCoordPre(int axis) {
	return this->coordPre[axis];
}

double Droplet::getMass() {
	return this->mass;
}

double Droplet::getVelocity() {
	return this->velocity;
}

Droplet *Droplet::getFinalMergred() {
	if (this->mergedInto != nullptr)
		return this->mergedInto->getFinalMergred();
	else
		return this;
}

void Droplet::setAbove(Droplet * drp) {
	if (this == drp)
		return;
	this->above = drp;
}

void Droplet::setBelow(Droplet * drp) {
	if (this == drp)
		return;
	this->below = drp;
}
