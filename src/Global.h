/*
 * Global.h
 *
 *  Created on: Jan 12, 2019
 */

#define GRAVITAIONAL_CONSTANT 6.67408E-11 // [m^3 / (kg * s^2)]

#define ENVIRONMENT_WIDTH 1         // [m]
#define ENVIRONMENT_HEIGHT 1        // [m]
#define ENVIRONMENT_SPAWN_WIDTH .6  // [m]
#define ENVIRONMENT_SPAWN_HEIGHT .6 // [m]
#define ENVIRONMENT_SPAWN_PARTICLES_TOTAL .5E3
#define ENVIRONMENT_SPAWN_PARTICLE_MASS 5.E2        // [kg]
#define ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2 1.8E1 // [kg] (2*std)
#define ENVIRONMENT_SPAWN_START_ANGULAR_VELO                                   \
    -2.E-3 // [m/2] (orthogonal towards the center)
#define ENVIRONMENT_SPAWN_FUNCTIONAL_ANGULAR_VELO_FUNCTION                     \
    2 // 0: plain angular function (alpha * radius), 1: circular orbit
      // sqrt(GM*radius)*2, 2: RungeKutta5 (-> by grav. force)
#define ENVIRONMENT_SPAWN_FORM                                                 \
    2 // 0: square, 1: ellipse, 2: concentrated ellipse, 3: spiral galaxy

#define GALAXY_SPIRAL_NUM 5
#define GALAXY_SPIRAL_WIDTH 0.1 // [m]

#define SIMULATION_TIME_MAX 36       // [s]
#define SIMULATION_TIME_PER_STEP 1.0   // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [simulation_time_per_step's]

#define VISU_WIDTH_PX_PER_METER 800  // [px/m]
#define VISU_HEIGHT_PX_PER_METER 800 // [px/m]

#define OPENCV_VIDEO_SCALE 1 // don't change!
