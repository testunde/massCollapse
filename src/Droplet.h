/*
 * Droplet.h
 *
 *  Created on: Dec 16, 2018
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

	double updateMass();
	double updateVelocity();
public:
	static int disposedDrops;
	static std::list<Droplet*> *dropList;

	static void sortListWidth();
	static void sortListHeight();
	static void sortListSize();

	static int remainingDrops();

	static Droplet *bigger_h, *smaller_h; // head, tail
	static Droplet *above_h, *below_h; // head, tail
	static Droplet *left_h, *right_h; // head, tail

	Droplet *bigger = nullptr, *smaller = nullptr; // head, tail
	Droplet *above = nullptr, *below = nullptr; // head, tail
	Droplet *left = nullptr, *right = nullptr; // head, tail

	void setAbove(Droplet * drp);
	void setBelow(Droplet * drp);

	Droplet(double radius, double coord[]);
	virtual ~Droplet();
	void deleteInstance();

	double getRadius();
	double getCoord(int axis);
	double getCoordPre(int axis);
	double getMass();
	double getVelocity();
	Droplet *getFinalMergred();

	double growCondensation();
	void merge(std::list<Droplet*> *list);
	void fallBy(double way);
	void setMergedInto(Droplet *drp);
};

#endif /* DROPLET_H_ */
