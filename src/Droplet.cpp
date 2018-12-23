/*
 * Droplet.cpp
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#include <cstring>

#include "Droplet.h"

using namespace std;

int Droplet::mergeCount = 0;

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
		above->below = this->below;
	if (below != nullptr)
		below->above = this->above;

	if (left != nullptr)
		left->right = this->right;
	if (right != nullptr)
		right->left = this->left;
}

// CALL IT ONCE AND LAST, IF CREATED WITH "new"!!!
void Droplet::recursiveDeleteWidthList() {
	// going towards tail using width-list links
	if (this->right != nullptr)
		this->right->recursiveDeleteWidthList();
	delete this;
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
	double tempVel = exp(pow(log(tempD) - 2.4, 2.) * (-3./18.) + 6.9);
	this->velocity = tempVel * 1.E-2; // [cm/s] -> [m/s]

	return this->velocity;
}

// after one second
double Droplet::growCondensation() {
	double seedRadius = (this->radius < DROP_GROWTH_CONDENSATION_LIMIT_SIZE) ? DROP_GROWTH_CONDENSATION_LIMIT_SIZE : this->radius;
	double dR = (DROP_GROWTH_CONDENSATION_S - 1.) / (seedRadius * DROP_GROWTH_CONDENSATION_F);

	this->radius += dR;

	updateMass();
	updateVelocity();

	// No need to check for bigger neighbors if growing begins from the biggest drop in the environment.
	// Only check bigger neighbor(s) when drop radius smaller than 0.31623 [µm] (~3E-7 [m]),
	// because until then the grow rate is bigger than the delta-radius would be.

	return dR;
}

void Droplet::merge(list<Droplet *> &list) {
	double totalMass = this->mass;

	double newCoord[] = {this->coord[0] * this->mass, this->coord[1]}; // no height-change as approximation
	for (Droplet * drp : list) {
		totalMass += drp->getMass();
		newCoord[0] += drp->getCoord(0) * drp->getMass();

		drp->setMergedInto(this);
		mergeCount++;
	}
	newCoord[0] /= totalMass;

	memcpy(this->coord, newCoord, 2 * sizeof(double));
	this->mass = totalMass;
	this->radius = this->mass / pow(DENSITY_WATER * (M_PI * 4. / 3.), 1./3.);
	updateVelocity();

	for (Droplet * drp : list) {
		delete drp; // for later faster&easier link updates
	}
	// TODO: update width list links
	// TODO: update size list links
}

void Droplet::fallBy(double way) {
	memcpy(this->coordPre, this->coord, 2 * sizeof(double));
	this->coord[1] += way; // [m]

	// Update height height-list links.
	// No need to check droplets above 'this' when beginning from bottom and fall-way is positive.

	Droplet *aboveDrop = this->above, *belowDrop = this->below;
	while (belowDrop != nullptr && belowDrop->getCoord(1) < this->coord[1]) {
		aboveDrop = belowDrop;
		belowDrop = belowDrop->below;
	}
	if (aboveDrop == this->above)
		return; // no change, even if 'this->above' was 'nullptr'

	// change pointers from new-neighbors to 'this'
	if (belowDrop != nullptr) {
		belowDrop->above = this;
	}
	aboveDrop->below = this;

	// change pointers from prior-neighbors which pointed to 'this'
	if (this->above != nullptr)
		this->above->below = this->below;
	if (this->below != nullptr)
		this->below->above = this->above;

	// change 'this' to new-neighbors
	this->above = aboveDrop;
	this->below = belowDrop;
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
