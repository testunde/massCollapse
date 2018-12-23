/*
 * Droplet.h
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#ifndef DROPLET_H_
#define DROPLET_H_

#include "math.h"
#include <list>

#include "Global.h"

class Droplet {
private:
	double radius; // [m]
	double coord[2] = {0., 0.}; // ([m]; [m])
	double coordPre[2] = {0., 0.}; // ([m]; [m])
	double mass; // [kg]
	double velocity; // [m/s]
	Droplet *mergedInto = nullptr;
	std::list<Droplet*> *parentList = nullptr;

	double updateMass();
	double updateVelocity();
public:
	static int mergeCount;

	Droplet(double radius, double coord[], std::list<Droplet*> *list);
	virtual ~Droplet();

	double getRadius();
	double getCoord(int axis);
	double getCoordPre(int axis);
	double getMass();
	double getVelocity();
	Droplet *getFinalMergred();

	double growCondensation();
	void merge(Droplet *drp);
	void fallBy(double way);
	void setMerged(Droplet *drp);
	void setParentList(std::list<Droplet*> *list);

	bool isInXRangeWith(Droplet *dropOther);
};

#endif /* DROPLET_H_ */
