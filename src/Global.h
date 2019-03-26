/*
 * Global.h
 *
 *  Created on: Jan 12, 2019
 */

#define GRAVITAIONAL_CONSTANT 6.67408E-11 // [m^3 / (kg * s^2)]

#define ENVIRONMENT_WIDTH 8          // [m]
#define ENVIRONMENT_HEIGHT 8         // [m]
#define ENVIRONMENT_SPAWN_WIDTH 3.8  // [m]
#define ENVIRONMENT_SPAWN_HEIGHT 1.8 // [m]
#define ENVIRONMENT_SPAWN_PARTICLES_TOTAL 1500
#define ENVIRONMENT_SPAWN_PARTICLE_MASS 3.E2        // [kg]
#define ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2 1.8E1 // [kg] (2*std)
#define ENVIRONMENT_SPAWN_PARTICLE_DENSITY                                     \
    3.E8 // [kg/m^3] (Osmium: 2.26E4 kg/m^3)
#define ENVIRONMENT_SPAWN_START_ANGULAR_VELO                                   \
    2.E-3 // [1/s] (orthogonal towards the center)
#define ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION                     \
    2 // 0: plain angular function (alpha * radius), 1: circular orbit
      // sqrt(GM*radius)*2, 2: RungeKutta5 (-> by grav. force)
#define ENVIRONMENT_SPAWN_FORM                                                 \
    2 // 0: square, 1: ellipse, 2: concentrated ellipse, 3: spiral galaxy

#define GALAXY_SPIRAL_NUM 3
#define GALAXY_SPIRAL_WIDTH 0.9 // [m]

#define SIMULATION_TIME_MAX 3600       // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [simulation_time_per_step's]
#define SIMULATION_TIME_PER_STEP 0.05  // [s]
#define SIMULATION_ROUNDS                                                      \
    6 // X rounds = X time steps (per GPU bulk calculation)
#define SIMULATION_ETA_STEPS_ACCUMULATION 20

#define VISU_WIDTH_PX_PER_METER 100  // [px/m]
#define VISU_HEIGHT_PX_PER_METER 100 // [px/m]

#ifdef USE_OPENCV
#define OPENCV_VIDEO_FILENAME "video.avi"
#define OPENCV_VIDEO_FORMAT 'H', 'F', 'Y', 'U' // raw; as for .avi format
#define OPENCV_VIDEO_STEPS_PER_SECOND                                          \
    60                       // [s/s] simulation steps per video second (FPS)
#define OPENCV_VIDEO_SCALE 1 // only positive int!
#endif

#define OPENCL_KERNEL_FILE "GravitationRK.cl"
