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

#ifdef USE_OPENCV
#include <opencv4/opencv2/opencv.hpp>
#endif

using namespace std;

long currentMicroSec() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1E6L + t.tv_nsec / 1E3L;
}

void update() {}

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

int main(int, char **) {
    srand(time(nullptr));
    normal_distribution<double> distribution(
        ENVIRONMENT_SPAWN_PARTICLE_MASS,
        ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2 * .5);
    default_random_engine rnd_gen;

    // initialize environment
    printf("Initialize environment + generate particles...\n");

    // generate particles
    int dCount = 0;
    double tMmax = DBL_MIN, tMmin = DBL_MAX;
    for (int c = 0; c < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; c++) {
        double tempCoord[] = {getRandom() * ENVIRONMENT_WIDTH,
                              getRandom() * ENVIRONMENT_HEIGHT};
        double tempVel[] = {0., 0.};
        double tempMass = 0.;
        while (abs(tempMass - ENVIRONMENT_SPAWN_PARTICLE_MASS) >
               ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2) {
            tempMass = distribution(rnd_gen);
        }

        if (tempMass > tMmax)
            tMmax = tempMass;
        if (tempMass < tMmin)
            tMmin = tempMass;

        Particle *tempParticle = new Particle(tempMass, tempCoord, tempVel);
        Particle::particleList->push_back(tempParticle);

        dCount++;
        if (dCount % ((int)ENVIRONMENT_SPAWN_PARTICLES_TOTAL / 10) == 0) {
            printf("%.0f%%\r",
                   100. * dCount / ENVIRONMENT_SPAWN_PARTICLES_TOTAL);
            fflush(stdout);
        }
    }
    printf("min: %f | max: %f [kg]\n", tMmin, tMmax);

    // init openCV
    int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
    int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;
#ifdef USE_OPENCV
    cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
    cv::Mat envVisu(visu_height, visu_width, CV_8UC3);
    cv::VideoWriter video(
        OPENCV_VIDEO_FILENAME, cv::VideoWriter::fourcc(OPENCV_VIDEO_FORMAT),
        OPENCV_VIDEO_SECONDS_PER_SECOND, cv::Size(visu_width, visu_height));
    visu_width /= OPENCV_VIDEO_SCALE;
    visu_height /= OPENCV_VIDEO_SCALE;
#endif

    // simulation start
    double meanTotalMassPerPx = ENVIRONMENT_SPAWN_PARTICLE_MASS *
                                ENVIRONMENT_SPAWN_PARTICLES_TOTAL /
                                (visu_width * visu_height);
    long startTimeStamp = currentMicroSec();
    for (int t = 0; t <= SIMULATION_TIME_MAX; t++) { // [s]
        update();

        if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
            // count mass per pixel
            auto total_mass = init_matrix<double>(visu_width, visu_height, 0.);
            auto count_particle = init_matrix<int>(visu_width, visu_height, 0);

            int count_particle_max = 0;
            for (Particle *p : *Particle::particleList) {
                int idx[2] = {
                    (int)(p->getPosition(0) * VISU_WIDTH_PX_PER_METER /
                          OPENCV_VIDEO_SCALE),
                    (int)(p->getPosition(1) * VISU_HEIGHT_PX_PER_METER /
                          OPENCV_VIDEO_SCALE)};
                idx[1] = (idx[1] >= visu_height) ? (visu_height - 1) : idx[1];
                total_mass[idx[0]][idx[1]] += p->getMass();
                count_particle[idx[0]][idx[1]]++;
            }

            // visualization
            double avgColor = 0.;
#ifdef USE_OPENCV
            envVisu.setTo(cv::Scalar(0, 0, 0));
            for (int w = 0; w < visu_width; w++) {
                for (int h = 0; h < visu_height; h++) {
                    double totalMassRatio =
                        total_mass[w][h] / meanTotalMassPerPx;
                    // totalMassRatio =
                    //    (log10(totalMassRatio + 0.1) + 1) / (log10(1.1)
                    //    + 1);
                    int colorTotalMass = (int)(.68 * 255. * totalMassRatio);
                    colorTotalMass =
                        (colorTotalMass > 255) ? 255 : colorTotalMass;

                    int colorCountParticle =
                        (int)(255. * ((double)count_particle[w][h]) /
                              (double)count_particle_max);

                    cv::Vec3b finalColor(colorCountParticle, colorTotalMass, 0);
                    for (int ww = 0; ww < OPENCV_VIDEO_SCALE; ww++)
                        for (int hh = 0; hh < OPENCV_VIDEO_SCALE; hh++)
                            envVisu.at<cv::Vec3b>(h * OPENCV_VIDEO_SCALE + hh,
                                                  w * OPENCV_VIDEO_SCALE + ww) =
                                finalColor;
                    avgColor +=
                        (finalColor[0] + finalColor[1] + finalColor[2]) / 3.;
                }
            }
            avgColor /= visu_width * visu_height;

            video.write(envVisu);
            cv::imshow("env_simu", envVisu);
            if (cv::waitKey(30) >= 0)
                break;
#endif

            long currentTimeStamp = currentMicroSec();
            long eta_seconds =
                (long)((SIMULATION_TIME_MAX - (long)t) *
                       (currentTimeStamp - startTimeStamp) / (1E6L * (long)t));
            printf("time: %d [s] | avgColor: %f [0-256[ (ETA: %ldm %lds)\n", t,
                   avgColor, eta_seconds / 60, eta_seconds % 60);
            fflush(stdout);
        }
    }

    // clear simulation environment
#ifdef USE_OPENCV
    video.release();
#endif
    clearEnvironment();
}
