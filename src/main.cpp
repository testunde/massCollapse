/*
 * main.cpp
 *
 *  Created on: Jan 12, 2019
 */

#include <exception>
#include <float.h> // DLB_MIN + DBL_MAX
#include <fstream>
#include <iostream>
#include <random> // normal_distribution
#include <stdio.h>
#include <stdlib.h> // rand
#include <time.h>   // seed for rand
#include <vector>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include "Global.h"
#include "Particle.h"

#ifdef USE_OPENCV
#include <opencv4/opencv2/opencv.hpp>
#endif

#ifdef USE_OPENMP
#include <parallel/algorithm> // __gnu_parallel::for_each()
#endif

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

void update(cl::CommandQueue &queue, cl::Kernel &kernel,
            cl::Buffer *buffers[2]) {
    static p_state states[(long)ENVIRONMENT_SPAWN_PARTICLES_TOTAL];
    static const bool roundOdd = SIMULATION_ROUNDS % 2 == 1;
    static bool roundOdd_c = false;

    cl_int er = kernel.setArg(0, *buffers[0]);
    if (er != CL_SUCCESS)
        throw cl::Error(er, "arg0");
    er = kernel.setArg(1, *buffers[1]);
    if (er != CL_SUCCESS)
        throw cl::Error(er, "arg1");
    er = kernel.setArg(2, roundOdd_c ? 1 : 0);
    if (er != CL_SUCCESS)
        throw cl::Error(er, "arg2");
    queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                               cl::NDRange(ENVIRONMENT_SPAWN_PARTICLES_TOTAL),
                               cl::NullRange);
    queue.finish();

    queue.enqueueReadBuffer(roundOdd_c ? *buffers[0] : *buffers[1], CL_TRUE, 0,
                            ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state),
                            states);

    for (int i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++)
        Particle::particleList->at(i)->setCLStruct(&states[i]);

    if (roundOdd)
        roundOdd_c = !roundOdd_c;
}

void clearEnvironment() {
    printf("Clearing environment...\n");
    fflush(stdout);
    for (Particle *d : *Particle::particleList) {
        delete d;
    }
    // Particle::particleList->clear();
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

string loadKernelSource() {
    string filename(OPENCL_KERNEL_FILE);
    // look in current dir - if missing, look in parent dir - so exe can be in
    // Release subdir of source, for instance
    ifstream file(filename);
    if (!file.is_open()) {
        file.open("../" + filename);
        if (!file.is_open())
            throw runtime_error("File '" + filename + "' not found");
    }

    return string((istreambuf_iterator<char>(file)),
                  istreambuf_iterator<char>());
}

void generateParticles(const int form,
                       normal_distribution<double> &distributionMass,
                       default_random_engine &rnd_gen, cl::CommandQueue &queue,
                       cl::Program &program, cl::Buffer *buffers[2]) {
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

        cl::Kernel kernel = cl::Kernel(program, "InitialVelocity");

        static p_state states[(long)ENVIRONMENT_SPAWN_PARTICLES_TOTAL];

        for (int i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++)
            states[i] = (*Particle::particleList)[i]->getCLStruct();
        queue.enqueueWriteBuffer(
            *buffers[0], CL_TRUE, 0,
            ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state), states);

        cl_int er = kernel.setArg(0, *buffers[0]);
        if (er != CL_SUCCESS)
            throw cl::Error(er, "arg0");
        er = kernel.setArg(1, *buffers[1]);
        if (er != CL_SUCCESS)
            throw cl::Error(er, "arg1");

        queue.enqueueNDRangeKernel(
            kernel, cl::NullRange,
            cl::NDRange(ENVIRONMENT_SPAWN_PARTICLES_TOTAL), cl::NullRange);
        queue.finish();

        queue.enqueueReadBuffer(
            *buffers[1], CL_TRUE, 0,
            ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state), states);

        for (int i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++) {
            vector<double> tempVel = {states[i].vel.x, states[i].vel.y};

            double tempVelAbs =
                sqrt(tempVel[0] * tempVel[0] + tempVel[1] * tempVel[1]);
            if (tempVelAbs > maxVel)
                maxVel = tempVelAbs;

            Particle::particleList->at(i)->setVelocity(tempVel);
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

#ifdef USE_OPENMP
    printf("Using parallelisation by OpenMP.\n");
#endif

    // initialize environment
    printf("Initialize environment + generate particles...\n");
    printf("Using seed for random generators: %u\n", rndSeed);

    // init openCV
    int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
    int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;
#ifdef USE_OPENCV
    cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
    cv::Mat envVisu(visu_height, visu_width, CV_8UC3);
    cv::VideoWriter video(OPENCV_VIDEO_FILENAME,
                          cv::VideoWriter::fourcc(OPENCV_VIDEO_FORMAT),
                          OPENCV_VIDEO_SECONDS_PER_SECOND / SIMULATION_ROUNDS,
                          cv::Size(visu_width, visu_height));
    visu_width /= OPENCV_VIDEO_SCALE;
    visu_height /= OPENCV_VIDEO_SCALE;
#endif

    // init openCL
    cl::Program program;
    vector<cl::Device> devices;
    try {
        // get a platform and device
        vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platforms.size() == 0)
            throw runtime_error("OpenCL not available.");

        // create context and queue
        cl_context_properties cprops[3] = {
            CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[0](), 0};
        cl::Context context = cl::Context(CL_DEVICE_TYPE_GPU, cprops);

        devices = context.getInfo<CL_CONTEXT_DEVICES>();
        if (devices.size() == 0)
            throw runtime_error("GPU device not available.");

        cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);

        // compile source, get kernel entry point
        string source = loadKernelSource();
        cl::Program::Sources sources(1,
                                     make_pair(source.c_str(), source.size()));
        program = cl::Program(context, sources);
        char macros[512] = {'\0'};
        if (512 <
            sprintf(
                macros,
                "-DENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION=%d "
                "-DGRAVITAIONAL_CONSTANT=%0.16f "
                "-DENVIRONMENT_SPAWN_PARTICLES_TOTAL=%ld "
                "-DENVIRONMENT_SPAWN_PARTICLE_MASS=%lf "
                "-DENVIRONMENT_SPAWN_START_ANGULAR_VELO=%lf "
                "-DSIMULATION_TIME_PER_STEP=%lf -DSIMULATION_ROUNDS=%d "
                "-DCOLLISION_DISTANCE=%lf",
                ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION,
                GRAVITAIONAL_CONSTANT, (long)ENVIRONMENT_SPAWN_PARTICLES_TOTAL,
                ENVIRONMENT_SPAWN_PARTICLE_MASS,
                ENVIRONMENT_SPAWN_START_ANGULAR_VELO, SIMULATION_TIME_PER_STEP,
                SIMULATION_ROUNDS, COLLISION_DISTANCE))
            throw runtime_error(
                "Stack smash by passing macros to kernel! (buffer overflow)");
        program.build(devices, macros);
        cl::Kernel kernel = cl::Kernel(program, "GravitationRK");

        cl::Buffer bufIn(context, CL_MEM_READ_WRITE,
                         ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state));
        cl::Buffer bufOut(context, CL_MEM_READ_WRITE,
                          ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state));
        cl::Buffer *buffers[2] = {&bufIn, &bufOut};

        // generate particles
        generateParticles(ENVIRONMENT_SPAWN_FORM, distributionMass, rnd_gen,
                          queue, program, buffers);

        {
            p_state states[(long)ENVIRONMENT_SPAWN_PARTICLES_TOTAL];

            for (int i = 0; i < ENVIRONMENT_SPAWN_PARTICLES_TOTAL; i++)
                states[i] = (*Particle::particleList)[i]->getCLStruct();
            queue.enqueueWriteBuffer(
                bufIn, CL_TRUE, 0,
                ENVIRONMENT_SPAWN_PARTICLES_TOTAL * sizeof(p_state), states);
        }

        // simulation start
        double meanTotalMassPerPx = ENVIRONMENT_SPAWN_PARTICLE_MASS *
                                    ENVIRONMENT_SPAWN_PARTICLES_TOTAL /
                                    (visu_width * visu_height);
        long startTimeStamp = currentMicroSec();

        for (int t = 0; t <= SIMULATION_TIME_MAX / SIMULATION_TIME_PER_STEP;
             t += SIMULATION_ROUNDS) { // [s]
            update(queue, kernel, buffers);

            if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
                // count mass per pixel
                auto total_mass =
                    init_matrix<double>(visu_width, visu_height, 0.);
                auto count_particle =
                    init_matrix<int>(visu_width, visu_height, 0);
                auto count_particle_collided =
                    init_matrix<int>(visu_width, visu_height, 0);

#ifdef USE_OPENMP
                __gnu_parallel::for_each(
                    Particle::particleList->begin(),
                    Particle::particleList->end(), [&](Particle *p) {
#else
                for (Particle *p : *Particle::particleList) {
#endif
                        int idx[2] = {
                            (int)((p->getPosition(0) + ENVIRONMENT_WIDTH * .5) *
                                  VISU_WIDTH_PX_PER_METER / OPENCV_VIDEO_SCALE),
                            (int)((p->getPosition(1) +
                                   ENVIRONMENT_HEIGHT * .5) *
                                  VISU_HEIGHT_PX_PER_METER /
                                  OPENCV_VIDEO_SCALE)};
                        idx[0] = max(min(idx[0], visu_height - 1), 0);
                        idx[1] = max(min(idx[1], visu_height - 1), 0);

                        if (p->getCollission()) {
                            count_particle_collided[idx[0]][idx[1]]++;
                        } else {
                            total_mass[idx[0]][idx[1]] += p->getMass();
                            count_particle[idx[0]][idx[1]]++;
                        }
#ifdef USE_OPENMP
                    });
#else
                }
#endif
                // visualization
                double avgColor = 0.;
#ifdef USE_OPENCV
                envVisu.setTo(cv::Scalar(0, 0, 0));
                for (int w = 0; w < visu_width; w++) {
                    for (int h = 0; h < visu_height; h++) {
                        double totalMassRatio =
                            total_mass[w][h] / meanTotalMassPerPx;
                        // totalMassRatio =
                        //    (log10(totalMassRatio + 0.1) + 1) /
                        //    (log10(1.1)
                        //    + 1);
                        int colorTotalMass =
                            min(255, (int)(.68 * 255. * totalMassRatio));
                        int colorMass =
                            (int)(255. *
                                  (1. - ENVIRONMENT_SPAWN_PARTICLE_MASS /
                                            (total_mass[w][h] +
                                             ENVIRONMENT_SPAWN_PARTICLE_MASS)));

                        int cP = count_particle[w][h];
                        int cP_c = count_particle_collided[w][h];
                        int colorCountParticle =
                            (int)(255. * ((double)cP) /
                                  (double)Particle::particleList->size());

                        int colorCountParticle2 =
                            (cP > 0)
                                ? (int)(255. *
                                        (1. - 1. / (1. * (double)simplePowBase2(
                                                             cP))))
                                : 0;

                        cv::Vec3b finalColor(0, colorCountParticle2, colorMass);
                        if (finalColor[1] == 0 && finalColor[2] == 0) {
                            // get blue-ish color fold collision-positions
                            finalColor[0] =
                                (cP_c > 0)
                                    ? (int)(255. *
                                            (1. -
                                             1. / (1. * (double)simplePowBase2(
                                                            cP_c))))
                                    : 0;
                        }
                        for (int ww = 0; ww < OPENCV_VIDEO_SCALE; ww++)
                            for (int hh = 0; hh < OPENCV_VIDEO_SCALE; hh++)
                                envVisu.at<cv::Vec3b>(
                                    h * OPENCV_VIDEO_SCALE + hh,
                                    w * OPENCV_VIDEO_SCALE + ww) = finalColor;
                        avgColor +=
                            (finalColor[0] + finalColor[1] + finalColor[2]) /
                            3.;
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
                    (long)(((SIMULATION_TIME_MAX / SIMULATION_TIME_PER_STEP) -
                            (long)t) *
                           (currentTimeStamp - startTimeStamp) /
                           (1E6L * (long)t));
                printf("time: %f [s] | avgColor: %f [0-256[ (ETA: %ldm %lds)\n",
                       ((float)t) * SIMULATION_TIME_PER_STEP, avgColor,
                       eta_seconds / 60, eta_seconds % 60);
                fflush(stdout);
            }
        }
    } catch (cl::Error &err) {
        printf("Error: %s (%d)\n", err.what(), err.err());
        /* if it was a compilation error */
        if (err.err() == CL_BUILD_PROGRAM_FAILURE)
            cout << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
                 << endl;
    } catch (std::exception &err) {
        printf("Error: %s\n", err.what());
    }

    // clear simulation environment
#ifdef USE_OPENCV
    video.release();
#endif
    clearEnvironment();
}
