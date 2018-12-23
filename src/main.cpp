/*
 * main.cpp
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#include <stdio.h>
#include <opencv4/opencv2/opencv.hpp>
#include <list>
#include <stdlib.h> // rand
#include <time.h> // seed for rand
#include <random> // normal_distribution

#include "Global.h"
#include "Droplet.h"

//using namespace cv;
using namespace std;

list<Droplet*> ** * environment;

bool isDifferentEnvCellHeight(double height1, double height2) {
	return floor(height1 / ENVIRONMENT_CELL_SIZE) != floor(height2 / ENVIRONMENT_CELL_SIZE);
}

void update() {
	// condensation growth
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		for (int h = 0; h < ENVIRONMENT_CELLS_HEIGHT; h++) {
			for (Droplet* drop : *environment[w][h]) {
				drop->growCondensation();
			}
		}
	}

	// (plain) falling
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		for (int h = ENVIRONMENT_CELLS_HEIGHT - 2; h >= 0; h--) {
			list<Droplet*> toRemove;
			for (Droplet* drop : *environment[w][h]) {
				double preHeight = drop->getCoord(1);
				drop->fallBy(abs(drop->getVelocity())); // also sets pre-coord
				if (isDifferentEnvCellHeight(preHeight, drop->getCoord(1))) {
					toRemove.push_back(drop);
					int hCoord = (int) floor(drop->getCoord(1) / ENVIRONMENT_CELL_SIZE);
					hCoord = (hCoord < ENVIRONMENT_CELLS_HEIGHT) ? ((hCoord < 0) ? 0 : hCoord) : (ENVIRONMENT_CELLS_HEIGHT - 1);
					environment[w][hCoord]->push_back(drop);
					drop->setParentList(environment[w][hCoord]);
				}
			}
			environment[w][h]->remove_if([&toRemove](Droplet *d){return find(toRemove.begin(), toRemove.end(), d) != toRemove.end();});
		}
	}

	// TODO: merging
	// maintain "toDeleteDrops" list (to the end, just delete them will remove them from parent env lists (through deconstructor)
	list<Droplet*> toDelete;
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		for (int h = 0; h <= ENVIRONMENT_CELLS_HEIGHT - 2; h++) {
			toDelete.push_
		}
	}
	for (Droplet *d : toDelete)
		delete d;
}

void clearEnvironment() {
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		for (int h = 0; h < ENVIRONMENT_CELLS_HEIGHT; h++) {
			for (Droplet *drop : *environment[w][h])
				delete drop;
			delete environment[w][h];
		}
		delete[] environment[w];
	}
	delete[] environment;
}

double getRandom() {
	return ((double) random()) / ((double) RAND_MAX);
}

int main(int, char**) {
	srand(time(nullptr));
	normal_distribution<double> distribution(ENVIRONMENT_SPAWN_DROP_SIZE, ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 * .5);
	default_random_engine rnd_gen;

	// initialize environment
	printf("Init env...\n");
	environment = new list<Droplet*>**[ENVIRONMENT_CELLS_WIDTH];
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		environment[w] = new list<Droplet*>*[ENVIRONMENT_CELLS_HEIGHT];
		for (int h = 0; h < ENVIRONMENT_CELLS_HEIGHT; h++) {
			environment[w][h] = new list<Droplet*>();
		}
	}

	// generate droplets
	double tSmax = DBL_MIN, tSmin = DBL_MAX;
	printf("Generate droplets...\n");
	for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
		for (int h = 0; h < ENVIRONMENT_CELLS_SPAWN_THICKNESS; h++) {
			for (int c = 0; c < ENVIRONMENT_SPAWN_DROPS_PER_CELL; c++) {
				double tempCoord[] = {(((double) w) + getRandom()) * ENVIRONMENT_CELL_SIZE,
						(((double) h) + getRandom()) * ENVIRONMENT_CELL_SIZE};
				double tempSize = 0.;
				while (abs(tempSize - ENVIRONMENT_SPAWN_DROP_SIZE) > ENVIRONMENT_SPAWN_DROP_SIZE_STD_2) {
					tempSize = distribution(rnd_gen);
				}

				if (tempSize > tSmax)
					tSmax = tempSize;
				if (tempSize < tSmin)
					tSmin = tempSize;

				Droplet * tempDrop = new Droplet(tempSize / 2., tempCoord, environment[w][h]);
				environment[w][h]->push_back(tempDrop);
			}

		}
	}
	printf("min: %f | max: %f\n", tSmin * 1.E6, tSmax * 1.E6);

	// init openCV
	int heightDivider = 10;
	cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
	cv::Mat envVisu(ENVIRONMENT_CELLS_HEIGHT / heightDivider, ENVIRONMENT_CELLS_WIDTH, CV_8UC3);

	// simulation start
	Droplet *statDrop = environment[0][0]->front();
	double meanTotalMassPerCell = DENSITY_WATER * pow(ENVIRONMENT_SPAWN_DROP_SIZE * .5, 3.) * (M_PI * 4. / 3.) * ENVIRONMENT_SPAWN_DROPS_PER_CELL;
	printf("starting simulation...\n");
	for (int t = 0; t <= SIMULATION_TIME_MAX; t++) {
		update();

		if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
			printf("time: %d [s] | drop size: %f [mm] | velocity: %f [cm/s] | height: %f [m]\n", t, statDrop->getRadius() * 2. * 1000., statDrop->getVelocity() * 100., statDrop->getCoord(1));
			fflush(stdout);

			// visualization
			for (int w = 0; w < ENVIRONMENT_CELLS_WIDTH; w++) {
				for (int h = 0; h < ENVIRONMENT_CELLS_HEIGHT / heightDivider; h++) {
					int dropCount = 0;
					double totalMass = .0;
					for (int c = 0; c < heightDivider; c++) {
						dropCount += environment[w][h * heightDivider + c]->size();
						for (Droplet *d : *environment[w][h * heightDivider + c]) {
							totalMass += d->getMass();
						}
					}
					// TODO: insert color variations depending on avgRadius distribution
					int color = (int) (.5 * 255. * ((double) dropCount) / (double) (ENVIRONMENT_SPAWN_DROPS_PER_CELL * heightDivider));
					int colorMass = (int) (.5 * 255. * totalMass / (meanTotalMassPerCell * (double) heightDivider));
					int colorFinal = (colorMass > 255) ? 255 : colorMass;
					envVisu.at<cv::Vec3b>(h, w) = cv::Vec3b(colorFinal, colorFinal, colorFinal);
				}
			}
			cv::imshow("env_simu", envVisu);
			if (cv::waitKey(30) >= 0)
				break;
		}
	}

	// clear simulation environment
	clearEnvironment();
}
