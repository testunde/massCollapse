/*
 * main.cpp
 *
 *  Created on: Jan 12, 2019
 */

#include <float.h> // DLB_MIN + DBL_MAX
#include <random>  // normal_distribution
#include <stdio.h>
#include <stdlib.h> // rand
#include <time.h>   // seed for rand
#include <vector>

#include "Global.h"
#include "Particle.h"


using namespace std;

long currentMicroSec() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1E6L + t.tv_nsec / 1E3L;
}

long simplePowBase2(int e) {
    long result = 1L;
    for (int i = e; i > 0; i--)
        result *= 2;
    return result;
}

void update() {
    // velocity
    for (Particle *p : *Particle::particleList)
        p->updateVelocity();

    // position
    for (Particle *p : *Particle::particleList)
        p->updatePosition();
}

void clearEnvironment() {
    printf("Clearing environment...\n");
    fflush(stdout);
    for (Particle *d : *Particle::particleList) {
        delete d;
    }
    Particle::particleList->clear();
}

// [0.0, 1.0]
double getRandom() { return ((double)random()) / ((double)RAND_MAX); }

template <class T> vector<vector<T>> init_matrix(int a, int b, T v) {
    vector<vector<T>> matrix;
    for (int i = 0; i < a; i++) {
        matrix.emplace_back(b, v);
    }
    return matrix;
}

double sign(double x) { return (double)((x > 0.) - (x < 0.)); }

void getSpiralPoint(double result[2], int c, int max,
                    default_random_engine &rnd_gen) {
    normal_distribution<double> distribution(0., .5);
    int cMax = max / GALAXY_SPIRAL_NUM;
    int cCupped = c % cMax;

    double theta =
        ((double)c / (double)cMax) * (M_PI * 2. / (double)GALAXY_SPIRAL_NUM);
    double theta2 = ((double)(c + 1) / (double)cMax) *
                    (M_PI * 2. / (double)GALAXY_SPIRAL_NUM);

    double radius = 1. - cos(M_PI_2 * (double)cCupped / (double)cMax);
    double radius2 = 1. - cos(M_PI_2 * (double)(cCupped + 1) / (double)cMax);

    double dx = radius * cos(theta) - radius2 * cos(theta2);
    double dy = radius * sin(theta) - radius2 * sin(theta2);
    double dNorm = sqrt(dx * dx + dy * dy);
    //    double shift = distribution(rnd_gen) * SPIRAL_WIDTH / dNorm;
    double r = getRandom() - .5;
    double shift = .5 * sign(r / 2.) * (1. - pow(cos(r * M_PI), 0.5)) *
                   GALAXY_SPIRAL_WIDTH / dNorm;

    radius = sqrt(radius) * .5;
    result[0] = (radius * ENVIRONMENT_SPAWN_WIDTH * cos(theta)) + (+dy) * shift;
    result[1] =
        (radius * ENVIRONMENT_SPAWN_HEIGHT * sin(theta)) + (-dx) * shift;
}

void generateParticles(const int form,
                       normal_distribution<double> &distributionMass,
                       default_random_engine &rnd_gen) {
    double maxVel = 0.;
    double tMmax = DBL_MIN, tMmin = DBL_MAX;
    for (int c = 0; c < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; c++) {
        // initial position by distribution and form
        double tempCoord[2] = {0., 0.};
        double theta, radius;
        vector<double> tCoord;
        switch (form) {
        case 3: // spiral galaxy
            getSpiralPoint(tempCoord, c, ENVIRONMENT_SPAWN_PARTICLES_TOTAL,
                           rnd_gen);
            break;
        case 2: // concentrated ellipse
            theta = getRandom() * 2. * M_PI;
            radius = 1 - pow(cos(getRandom() * M_PI_2), 0.5);
            radius = sqrt(radius) * .5;
            tempCoord[0] = radius * ENVIRONMENT_SPAWN_WIDTH * cos(theta);
            tempCoord[1] = radius * ENVIRONMENT_SPAWN_HEIGHT * sin(theta);
            break;
        case 1: // ellipse
            theta = getRandom() * 2. * M_PI;
            radius = getRandom();
            radius = sqrt(radius) * .5;
            tempCoord[0] = radius * ENVIRONMENT_SPAWN_WIDTH * cos(theta);
            tempCoord[1] = radius * ENVIRONMENT_SPAWN_HEIGHT * sin(theta);
            break;
        case 0: // square
        default:
            tempCoord[0] = (getRandom() - .5) * ENVIRONMENT_SPAWN_WIDTH;
            tempCoord[1] = (getRandom() - .5) * ENVIRONMENT_SPAWN_HEIGHT;
            break;
        }
        double coordNorm =
            sqrt(tempCoord[0] * tempCoord[0] + tempCoord[1] * tempCoord[1]);

        // initial velocity by function
        double absVel = 0.;
        switch (ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION) {
        case 2: // RungeKutta5 (after all points are generated)
            break;
        case 1: // circular orbit sqrt(GM*radius)*2
            absVel = sqrt(GRAVITAIONAL_CONSTANT *
                          (ENVIRONMENT_SPAWN_PARTICLES_TOTAL *
                           ENVIRONMENT_SPAWN_PARTICLE_MASS) *
                          coordNorm) *
                     2.;
            break;
        case 0: // plain angular function (alpha * radius)
        default:
            absVel = ENVIRONMENT_SPAWN_START_ANGULAR_VELO * coordNorm;
            break;
        }
        if (absVel > maxVel)
            maxVel = absVel;

        double tempVel[] = {absVel * (+tempCoord[1]) / coordNorm,
                            absVel * (-tempCoord[0]) / coordNorm};
        double tempMass = 0.;
        while (abs(tempMass - ENVIRONMENT_SPAWN_PARTICLE_MASS) >
               ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2) {
            tempMass = distributionMass(rnd_gen);
        }

        if (tempMass > tMmax)
            tMmax = tempMass;
        if (tempMass < tMmin)
            tMmin = tempMass;

        Particle *tempParticle = new Particle(tempMass, tempCoord, tempVel);
        Particle::particleList->push_back(tempParticle);
    }
    printf("min mass: %f | max mass: %f [kg]\n", tMmin, tMmax);

    if (ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION == 2) { // RK5
        for (Particle *p : *Particle::particleList) {
                vector<double> gravForce = p->getcurrentGravForce();
                // assuming force points to center, since average position of
                // generated particles should be the center
                double fNorm = sqrt(gravForce[0] * gravForce[0] +
                                    gravForce[1] * gravForce[1]);
                double pNorm = sqrt(p->getPosition(0) * p->getPosition(0) +
                                    p->getPosition(1) * p->getPosition(1));

                double absVel = sqrt(
                    fNorm * pNorm); // by circular orbit eq. and grav. force eq.

                if (absVel > maxVel)
                    maxVel = absVel;

                vector<double> tempVel = {absVel * (-gravForce[1]) / fNorm,
                                          absVel * (+gravForce[0]) / fNorm};

                p->setVelocity(tempVel);
        }
    }
    printf("max velocity: %f [m/s]\n", maxVel);
}

int main(int, char **) {
    unsigned int rndSeed = time(nullptr);
    srand(rndSeed);
    normal_distribution<double> distributionMass(
        ENVIRONMENT_SPAWN_PARTICLE_MASS,
        ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2 * .5);
    default_random_engine rnd_gen;
    rnd_gen.seed(rndSeed);


    // initialize environment
    printf("Initialize environment + generate particles...\n");
    printf("Using seed for random generators: %u\n", rndSeed);

    // generate particles
    if (ENVIRONMENT_SPAWN_FORM >= 0 && ENVIRONMENT_SPAWN_FORM <= 3)
        generateParticles(ENVIRONMENT_SPAWN_FORM, distributionMass, rnd_gen);

    // init openCV
    int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
    int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;

    // simulation start
    double meanTotalMassPerPx = ENVIRONMENT_SPAWN_PARTICLE_MASS *
                                ENVIRONMENT_SPAWN_PARTICLES_TOTAL /
                                (visu_width * visu_height);
    long startTimeStamp = currentMicroSec();
    for (int t = 0; t <= SIMULATION_TIME_MAX / SIMULATION_TIME_PER_STEP;
         t++) { // [s]
        update();

        if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
            // count mass per pixel
            auto total_mass = init_matrix<double>(visu_width, visu_height, 0.);
            auto count_particle = init_matrix<int>(visu_width, visu_height, 0);

            for (Particle *p : *Particle::particleList) {
                int idx[2] = {
                    (int)((p->getPosition(0) + ENVIRONMENT_WIDTH * .5) *
                          VISU_WIDTH_PX_PER_METER / OPENCV_VIDEO_SCALE),
                    (int)((p->getPosition(1) + ENVIRONMENT_HEIGHT * .5) *
                          VISU_HEIGHT_PX_PER_METER / OPENCV_VIDEO_SCALE)};
                idx[0] = max(min(idx[0], visu_height - 1), 0);
                idx[1] = max(min(idx[1], visu_height - 1), 0);
                total_mass[idx[0]][idx[1]] += p->getMass();
                count_particle[idx[0]][idx[1]]++;
            }
            // visualization
            double avgColor = 0.;

            if (t > 0) {
                long currentTimeStamp = currentMicroSec();
                long eta_seconds =
                    (long)(((SIMULATION_TIME_MAX / SIMULATION_TIME_PER_STEP) -
                            (long)t) *
                           (currentTimeStamp - startTimeStamp) / (1E6L * (long)t));
                printf("time: %f [s] | avgColor: %f [0-256[ (ETA: %ldm %lds)\r",
                       ((float)t) * SIMULATION_TIME_PER_STEP, avgColor,
                       eta_seconds / 60, eta_seconds % 60);
                fflush(stdout);
            }
        }
    }
    long endTimeStamp = currentMicroSec();
    printf("\nTook %f seconds to run.\n", ((float) (endTimeStamp - startTimeStamp)) / 1.0E6f);

    // clear simulation environment
    clearEnvironment();
}
