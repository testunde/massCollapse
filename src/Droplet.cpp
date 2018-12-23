/*
 * Droplet.cpp
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#include <cstring>

#include "Droplet.h"

int Droplet::mergeCount = 0;

Droplet::Droplet(double radius, double coord[2], std::list<Droplet*> *list) {
	this->radius = radius;
	memcpy(this->coord, coord, 2 * sizeof(double));
	memcpy(this->coordPre, coord, 2 * sizeof(double));

	this->parentList = list;

	updateMass();
	updateVelocity();
}

Droplet::~Droplet() {
	if (parentList != nullptr)
		parentList->remove(this);
}

double Droplet::updateMass() {
	this->mass = DENSITY_WATER * pow(this->radius, 3.) * (M_PI * 4. / 3.);
	return this->mass;
}

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
	double dR = (DROP_GROWTH_CONDENSATION_S - 1.) / (this->radius * DROP_GROWTH_CONDENSATION_F);

	this->radius += dR;

	updateMass();
	updateVelocity();

	return dR;
}

void Droplet::merge(Droplet *drp) {
	double totalMass = this->mass + drp->getMass();
	double newCoord[] = {(this->coord[0] * this->mass + drp->getCoord(0) * drp->getMass()) / totalMass,
			this->coord[1]};// no height-change as approximation //(this->coord[1] * this->mass + drp.getCoord(1) * drp.getMass()) / totalMass};

	memcpy(this->coord, newCoord, 2 * sizeof(double));
	this->mass = totalMass;
	this->radius = this->mass / pow(DENSITY_WATER * (M_PI * 4. / 3.), 1./3.);
	updateVelocity();

	drp->setMerged(this);
	mergeCount++;
}

void Droplet::fallBy(double way) {
	memcpy(this->coordPre, this->coord, 2 * sizeof(double));
	this->coord[1] += way; // [m]
}

void Droplet::setMerged(Droplet *drp) {
	this->mergedInto = drp;
}

void Droplet::setParentList(std::list<Droplet*> *list) {
	this->parentList = list;
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

bool Droplet::isInXRangeWith(Droplet* dropOther) {
	double thisBounds[2] = {coord[0] - radius, coord[0] + radius};
	double otherBounds[2] = {dropOther->getCoord(0) - dropOther->getRadius(), dropOther->getCoord(0) + dropOther->getRadius()};

	bool result = thisBounds[0] < otherBounds[0] && thisBounds[0] < otherBounds[1] && thisBounds[1] < otherBounds[0] && thisBounds[1] < otherBounds[1]
					&& thisBounds[0] > otherBounds[0] && thisBounds[0] > otherBounds[1] && thisBounds[1] > otherBounds[0] && thisBounds[1] > otherBounds[1];
	return result;
}
