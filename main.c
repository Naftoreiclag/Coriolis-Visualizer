/*

Copyright 2017 James Fong

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define VID_W (1280)
#define VID_H (720)

#define PATH_QUALITY 256
#define PATH_TIME 0.5
#define PATH_LEN 96
#define OCTOPUS_ARMS 8

#define FRAME_START 0
#define FRAME_COUNT 300

#define RENDER_THICK 50

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795028841971694
#endif

typedef struct Vec3 {
    double x, y, z;
} Vec3;

Vec3 vec_new(double x, double y, double z) {
    Vec3 vec = {x, y, z};
    return vec;
}

typedef struct Color {
    unsigned char r, g, b;
} Color;

typedef struct Path {
    Vec3 points[PATH_LEN];
} Path;

typedef struct Walker {
    Vec3 loc;
    Vec3 vel;
} Walker;

Path g_octopus[OCTOPUS_ARMS];
Color g_octo_color[OCTOPUS_ARMS];

Color g_frame[VID_W * VID_H];
Color* g_plate;
int g_plate_w;
int g_plate_h;
Color* g_starmap;
int g_starmap_w;
int g_starmap_h;

Vec3 g_earth_loc = {0, 0, 0};
double g_earth_radius = 6371.008; // Kilometers
double g_earth_ang_spd = 0.26251614; // Radians per hour
Vec3 g_cam_loc  = {0, 0, -3 * 6371.008};
Vec3 g_cam_right    = {1, 0, 0};
Vec3 g_cam_up       = {0, 1, 0};
Vec3 g_cam_forward  = {0, 0, 1};

double g_cam_focal_len = 2.0;

int clamp(int val, int min, int max) {
    if (val < min) return min;
    if (val >= max) return max - 1;
    return val;
}

double clampd(double val, double min, double max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

Color debug_to_color(Vec3 vec) {
    Color ret = {vec.x * 255, vec.y * 255, vec.z * 255};
    return ret;
}

Color debug_to_color2(Vec3 vec) {
    vec.x = (clampd(vec.x, -1, 1) + 1.0) / 2.0;
    vec.y = (clampd(vec.y, -1, 1) + 1.0) / 2.0;
    vec.z = (clampd(vec.z, -1, 1) + 1.0) / 2.0;
    Color ret = {vec.x * 255, vec.y * 255, vec.z * 255};
    return ret;
}

Color debug_to_color3(Vec3 vec) {
    vec.x = clampd(vec.x, 0, 1);
    vec.y = clampd(vec.y, 0, 1);
    vec.z = clampd(vec.z, 0, 1);
    Color ret = {vec.x * 255, vec.y * 255, vec.z * 255};
    return ret;
}

int load_plate() {
    int n, k, nn;
    unsigned char* data = stbi_load("plate.png", 
            &g_plate_w, &g_plate_h, &n, sizeof(Color));
    if (data) {
        printf("Loaded equirectangular map successfully.\n");
        g_plate = (Color*) data;
    } else {
        fprintf(stderr, "Fatal error loading Earth map image!\n");
        return 0;
    }
    data = stbi_load("starmap.png", 
            &g_starmap_w, &g_starmap_h, &n, sizeof(Color));
    if (data) {
        printf("Loaded equirectangular starmap successfully.\n");
        g_starmap = (Color*) data;
    } else {
        fprintf(stderr, "Fatal error loading Starmap image!\n");
        return 0;
    }
    data = stbi_load("set.png", 
            &nn, &k, &n, sizeof(Color));
    if (data) {
        printf("Loaded rainbow map successfully.\n");
        Color* rainbow = (Color*) data;
        for (int i = 0; i < OCTOPUS_ARMS; ++i) {
            g_octo_color[i] = rainbow[i];
        }
        stbi_image_free(data);
    } else {
        fprintf(stderr, "Fatal error loading rainbow map image!\n");
        return 0;
    }
    return 1;
    
}

Color get_equirec_color(
        double rad_nort, double rad_east, 
        Color* image, int img_w, int img_h) {
    int pix_x = (0.5 + (rad_east / (M_PI * 2))) * img_w;
    int pix_y = (0.5 - (rad_nort / M_PI)) * img_h;
    
    // Clamp to valid range
    pix_x = clamp(pix_x, 0, img_w);
    pix_y = clamp(pix_y, 0, img_h);
    return image[pix_x + pix_y * img_w];
}

Color get_starmap_color(double rad_nort, double rad_east) {
    return get_equirec_color(rad_nort, rad_east, 
            g_starmap, g_starmap_w, g_starmap_h);
}

Color get_plate_color(double rad_nort, double rad_east) {
    return get_equirec_color(rad_nort, rad_east, g_plate, g_plate_w, g_plate_h);
}

void set_frame_color(int x, int y, Color c) {
    if (x < 0 || x >= VID_W) {
        fprintf(stderr, "Tried to write out of bounds, x: %d", x);
        return;
    }
    if (y < 0 || y >= VID_H) {
        fprintf(stderr, "Tried to write out of bounds, x: %d", x);
        return;
    }
    g_frame[x + y * VID_W] = c;
}

void unload_plate() {
    stbi_image_free((unsigned char*) g_plate);
    stbi_image_free((unsigned char*) g_starmap);
}

int export_frame(const char* fname) {
    int success = stbi_write_jpg(fname, VID_W, VID_H, sizeof(Color), 
            (void*) g_frame, VID_W * sizeof(Color));
    if (success) {
        printf("Wrote frame to file %s\n", fname);
        return 1;
    } else {
        fprintf(stderr, "Error writing to file %s\n", fname);
        return 0;
    }
}

double vec_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Up and forward gives right
Vec3 vec_cross(Vec3 a, Vec3 b) {
    Vec3 ret = {a.y * b.z - a.z * b.y, 
            a.z * b.x - a.x * b.z, 
            a.x * b.y - a.y * b.x};
    return ret;
}

double vec_magsq(Vec3 vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

double vec_mag(Vec3 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

int vec_is_zero(Vec3 vec) {
    return vec.x == 0.0 && vec.y == 0.0 && vec.z == 0.0;
}

Vec3 vec_add(Vec3 a, Vec3 b) {
    Vec3 ret = {a.x + b.x, a.y + b.y, a.z + b.z};
    return ret;
}

Vec3 vec_sub(Vec3 a, Vec3 b) {
    Vec3 ret = {a.x - b.x, a.y - b.y, a.z - b.z};
    return ret;
}

Vec3 vec_mul(Vec3 vec, double amnt) {
    Vec3 ret = {vec.x * amnt, vec.y * amnt, vec.z * amnt};
    return ret;
}

Vec3 vec_norm(Vec3 vec) {
    double mag = vec_mag(vec);
    Vec3 ret = {vec.x / mag, vec.y / mag, vec.z / mag};
    return ret;
}

Vec3 vec_resize(Vec3 vec, double len) {
    return vec_mul(vec_norm(vec), len);
}

Vec3 vec_project(Vec3 vec, Vec3 norm_dir) {
    return vec_mul(norm_dir, vec_dot(vec, norm_dir));
}

Vec3 vec_reject(Vec3 vec, Vec3 norm_dir) {
    return vec_sub(vec, vec_project(vec, norm_dir));
}

Vec3 vec_angle_axis_rot(Vec3 axis, double angle, Vec3 vec) {
    // Rodrigues' rotation formula
    return vec_add(vec_add(
            vec_mul(vec, cos(angle)), 
            vec_mul(vec_cross(axis, vec), sin(angle))), 
            vec_mul(axis, vec_dot(axis, vec) * (1 - cos(angle))));
}

double sq(double x) {
    return x * x;
}

void orthogonalize(Vec3 dir, Vec3 skyward, 
        Vec3* forward, Vec3* right, Vec3* up) {
    *forward = vec_norm(dir);
    *right = vec_norm(vec_cross(skyward, *forward));
    *up = vec_norm(vec_cross(*forward, *right));
}

void cam_lookat(Vec3 lookat, Vec3 sidereal_up) {
    orthogonalize(vec_sub(lookat, g_cam_loc), sidereal_up, 
            &g_cam_forward, &g_cam_right, &g_cam_up);
}

#define NO_HIT -1
int sphere_intersect_ray(Vec3 O, Vec3 L, double R, Vec3 C, Vec3* interp) {
    if (vec_magsq(vec_sub(O, C)) < sq(R)) {
        *interp = O;
        return 0;
    }
    double radicand = 
            sq(vec_dot(L, vec_sub(O, C))) - vec_magsq(vec_sub(O, C)) + sq(R);
    if (radicand < 0) {
        return NO_HIT;
    }
    double dist = -vec_dot(L, vec_sub(O, C)) - sqrt(radicand);
    if (dist < 0) {
        return NO_HIT;
    }
    *interp = vec_add(O, vec_mul(L, dist));
    return dist;
}

int intersect_earth(Vec3 src, Vec3 norm_dir, Vec3* interp) {
    return sphere_intersect_ray(
            src, norm_dir, g_earth_radius, g_earth_loc, interp);
}

Vec3 lat_long_to_dir(double rad_nort, double rad_east, double mag) {
    Vec3 ret;
    ret.y = sin(rad_nort) * mag;
    ret.x = cos(rad_east) * cos(rad_nort) * mag;
    ret.z = sin(rad_east) * cos(rad_nort) * mag;
    
    return ret;
}

Color sky_color(Vec3 norm_dir) {
    Vec3 dir_collap = norm_dir;
    dir_collap.y = 0;
    dir_collap = vec_norm(dir_collap);
    
    double rad_nort = asin(norm_dir.y);
    double rad_east = atan2(dir_collap.z, dir_collap.x);
    
    return get_starmap_color(rad_nort, rad_east);
}

Color earth_color(Vec3 surf_loc) {
    Vec3 dir = vec_norm(vec_sub(surf_loc, g_earth_loc));
    Vec3 dir_collap = dir;
    dir_collap.y = 0;
    dir_collap = vec_norm(dir_collap);
    
    double rad_nort = asin(dir.y);
    double rad_east = atan2(dir_collap.z, dir_collap.x);
    
    return get_plate_color(rad_nort, rad_east);
}

Color color_add(Color ca, Color cb) {
    int r = clamp(ca.r + cb.r, 0, 256);
    int g = clamp(ca.g + cb.g, 0, 256);
    int b = clamp(ca.b + cb.b, 0, 256);
    Color ret = {r, g, b};
    return ret;
}

Color color_mul(Color c, double d) {
    Color ret = {c.r * d, c.g * d, c.b * d};
    return ret;
}

double my_fmod(double x, double y) {
    return fmod(fmod(x, y) + y, y);
}

Color g_atmo_color = {0x00, 0xCC, 0xFF};
Color g_pulse_color = {0xFF, 0xFF, 0xFF};
Color raytrace(Vec3 src, Vec3 norm_dir, int pixel_id, double anim_time) {
    Vec3 hit;
    double depth = DBL_MAX;
    Color color;
    Vec3 other_hit;
    double other_depth;
    
    other_depth = intersect_earth(src, norm_dir, &other_hit);
    if (other_depth != NO_HIT && other_depth < depth) {
        hit = other_hit;
        depth = other_depth;
        color = earth_color(hit);
        
        double glow_str = 
                1.0 - vec_dot(
                        vec_mul(norm_dir, -1), 
                        vec_norm(vec_sub(hit, g_earth_loc)));
        glow_str = clampd(glow_str * glow_str * glow_str, 0.0, 1.0) * 0.5;
        
        color = color_add(color, color_mul(g_atmo_color, glow_str));
    }
    
    double dsteps_per_day = 24.0 / PATH_TIME;
    int steps_per_day = dsteps_per_day;
    int pulse_shift = (int) (anim_time * dsteps_per_day);
    
    if (0 && sphere_intersect_ray(src, norm_dir, 
            g_earth_radius + (RENDER_THICK / 2), 
            g_earth_loc, &other_hit) != NO_HIT) {
        for (int arm_idx_raw = 0; arm_idx_raw < OCTOPUS_ARMS; ++arm_idx_raw) {
			
			// FOR INTRODUCTION
            /*
            if (anim_time < 1) {
				break;
			}
            */
			
            int arm_idx = arm_idx_raw;//(arm_idx_raw + pixel_id) % OCTOPUS_ARMS;
            int hit_arm = 0;
            double hit_step_time = 0.0;
            for (int step_idx = 0; step_idx < PATH_LEN; ++step_idx) {
                double step_time = (((step_idx - pulse_shift) % steps_per_day) 
                        + steps_per_day) % steps_per_day;
                step_time /= dsteps_per_day;
				
				// FOR INTRODUCTION
                /*
                int delta = (((step_idx - pulse_shift) % steps_per_day)
                        + steps_per_day) % steps_per_day;
				if (anim_time < 16) {
					if (step_idx - pulse_shift > -1 * steps_per_day || 
                            delta != steps_per_day - 2) {
						continue;
					}
				} else if (anim_time < 18) {
					if (step_idx - pulse_shift > -16 * steps_per_day && 
                            delta != steps_per_day - 2) {
						continue;
					}
				}
                */
				
                other_depth = sphere_intersect_ray(src, norm_dir, 
                        RENDER_THICK, 
                        g_octopus[arm_idx].points[step_idx], &other_hit);
                
                if (other_depth != NO_HIT && other_depth < depth) {
                    hit = other_hit;
                    depth = other_depth;
                    color = g_octo_color[arm_idx];
                    if (step_time > hit_step_time || !hit_arm) {
                        hit_step_time = step_time;
                    }
                    hit_arm = 1;
                }
            }
            if (hit_arm) {
                double intense = hit_step_time;
                intense *= intense;
                color = color_mul(color, 0.5 + (intense / 2.0));
            }
        }
    }
    
    if (depth == DBL_MAX) {
        color = sky_color(norm_dir);
    }
    
    return color;
}

int render_frame(double anim_time) {
    double asp_rat = ((double) VID_W) / ((double) VID_H);
    Vec3 canvas_disp = vec_mul(g_cam_forward, g_cam_focal_len);
    #pragma omp parallel for schedule(dynamic)
    for (int y = 0; y < VID_H; ++y) {
        double c_y = (((double) y) / ((double) VID_H)) * 2.0 - 1.0;
        for (int x = 0; x < VID_W; ++x) {
            double c_x = (((double) x) / ((double) VID_W)) * 2.0 - 1.0;
            c_x *= asp_rat;
            
            Vec3 ray_dir = vec_norm(
                    vec_add(vec_add(
                            vec_mul(g_cam_up, -c_y), 
                            vec_mul(g_cam_right, c_x)), 
                            canvas_disp));
            set_frame_color(x, y, raytrace(
                    g_cam_loc, ray_dir, y * VID_W + x, anim_time));
        }
    }
    return 1;
}

Walker apply_velocity_on_earth(Walker walker, double delta, double coriol) {
    if (vec_is_zero(walker.vel)) return walker;
    if (coriol != 0.0) {
        Vec3 effect = vec_cross(
                vec_new(0, g_earth_ang_spd * -2, 0), 
                vec_reject(walker.vel, vec_new(0, 1, 0)));
        walker.vel = vec_add(walker.vel, vec_mul(effect, delta * -coriol));
    }
    Vec3 surface_loc = vec_sub(walker.loc, g_earth_loc);
    Vec3 axis = vec_norm(vec_cross(surface_loc, walker.vel));
    double angle = (vec_mag(walker.vel) / g_earth_radius) * delta;
    walker.loc = vec_add(
            g_earth_loc, 
            vec_resize(vec_angle_axis_rot(axis, angle, surface_loc), 
            g_earth_radius));
    walker.vel = vec_reject(
            vec_angle_axis_rot(axis, angle, walker.vel), vec_norm(surface_loc));
    return walker;
}

void generate_octopus(Vec3 abs_origin, double speed, double coriol
        , Vec3 sidereal_up) {
    Vec3 rel_origin = 
            vec_mul(vec_norm(vec_sub(abs_origin, g_earth_loc)), g_earth_radius);
    
    Vec3 ortho_forward;
    Vec3 ortho_right;
    Vec3 ortho_up;
    orthogonalize(rel_origin, sidereal_up, 
            &ortho_forward, &ortho_right, &ortho_up);
    ortho_right = vec_mul(ortho_right, -1);
    
	#pragma omp parallel for
    for (int arm_idx = 0; arm_idx < OCTOPUS_ARMS; ++arm_idx) {
        double angle = 
                (((double) arm_idx) / ((double) OCTOPUS_ARMS)) * 2 * M_PI;
        Walker walker = {rel_origin, vec_add(
                vec_mul(ortho_right, cos(angle) * speed),
                vec_mul(ortho_up, sin(angle) * speed))};
        for (int step_idx = 0; step_idx < PATH_LEN; ++step_idx) {
            for (int rep = 0; rep < PATH_QUALITY; ++rep) {
                walker = apply_velocity_on_earth(
                        walker, 
                        PATH_TIME / ((double) PATH_QUALITY),
                        coriol);
            }
            g_octopus[arm_idx].points[step_idx] = walker.loc;
        }
    }
}

// Render a sequence of images 
int main(int argc, char* argv[]) {
    printf("Coriolis visualizer\n");
    #ifdef _OPENMP
        double timer_start, timer_end;
        printf("Using OpenMP!\n");
    #else
        clock_t timer_start, timer_end;
    #endif
    fflush(stdout);
    
    if (!load_plate()) return -1;
    fflush(stdout);
    
    // idle
    /*
    g_cam_loc = lat_long_to_dir(0, M_PI / -2, 18000);
    generate_octopus(g_cam_loc, 100, 0.5, vec_new(0, 1, 0));
    cam_lookat(g_earth_loc, vec_new(0, 1, 0));
    */

    double time_elapse;// 30 * 80 = 2400
    for (int frame_idx = FRAME_START; 
            frame_idx < FRAME_START + FRAME_COUNT; ++frame_idx) {
        #ifdef _OPENMP
            timer_start = omp_get_wtime();
        #else
            timer_start = clock();
        #endif
        double anim_time = ((double) frame_idx) / 30.0;
        // Earth arrival (300 frames)
        g_cam_loc = lat_long_to_dir(0, 
                M_PI / -2, 
                pow(2.71828, -anim_time * 1.5424948470) * 5000000 + 18000);
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        // Apollo blue marble recreation
        /*
        g_cam_loc = lat_long_to_dir(
                -0.4994, 
                0.6592, 
                g_earth_radius + 29000);
        g_cam_focal_len = 4.0;
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        */
		// Introduction
        /*
        g_cam_loc = lat_long_to_dir(0, M_PI / -2, 18000);
        generate_octopus(
                g_cam_loc, 
                100, 
                clampd((anim_time - 20.0) / 10, 0.0, 0.5), 
                vec_new(0, 1, 0));
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        */
        /* // Oscillate between latitudes
        g_cam_loc = lat_long_to_dir(
                sin(anim_time * (M_PI / 2) * 0.05) * 0.7, 
                M_PI / -2, 
                18000);
        generate_octopus(g_cam_loc, 100, 0.5, vec_new(0, 1, 0));
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        */
        /* // Fixed camera of above
        g_cam_loc = lat_long_to_dir((M_PI / 180.0) * 15.0, M_PI / -2, 18000);
        generate_octopus(lat_long_to_dir(
                sin(anim_time * (M_PI / 2) * 0.05) * 0.7, M_PI / -2, 10000), 
                100, 0.5, vec_new(0, 1, 0));
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        */
        /*// Oscillate between latitudes while also rotating around the globe
        g_cam_loc = lat_long_to_dir(sin(anim_time * (M_PI / 2) * 0.05) * 0.7, 
                anim_time * (M_PI / 2) * 0.2 + (M_PI / -2), 18000);
        generate_octopus(g_cam_loc, 100, 0.5, vec_new(0, 1, 0));
        cam_lookat(g_earth_loc, vec_new(0, 1, 0));
        */
        // Polar orbit around the globe
        /*
        g_cam_loc = lat_long_to_dir(
                (anim_time * (2 * M_PI)) / 80.0, M_PI / -2, 18000);
        generate_octopus(g_cam_loc, 100, 0.5, vec_new(1, 0, 0));
        cam_lookat(g_earth_loc, vec_new(1, 0, 0));
        */
        // Polar flyby
        /*
        g_cam_focal_len = 4.0;
        g_cam_loc = lat_long_to_dir(
                0, (M_PI / -2) + (M_PI / 16), g_earth_radius + 4000);
        Vec3 swirl_loc = lat_long_to_dir(
                (anim_time * (M_PI)) / 80.0 - M_PI / 2, 
                M_PI / -2, g_earth_radius);
        generate_octopus(swirl_loc, 100, 0.5, vec_new(1, 0, 0));
        cam_lookat(swirl_loc, vec_norm(g_cam_loc));
        */
        // Camera rolling transition
        /*
		g_cam_loc = lat_long_to_dir(0, M_PI / -2, 18000);
        Vec3 rot_up = vec_new(sin(
                anim_time * (M_PI / 2) * 0.2), 
                cos(anim_time * (M_PI / 2) * 0.2), 0);
        generate_octopus(g_cam_loc, 100, 0.5, rot_up);
        cam_lookat(g_earth_loc, rot_up);
		*/
        
        
        if (!render_frame(anim_time)) return -1;
        char fname_buff[16];
        sprintf(fname_buff, "frame/%06d.jpg", frame_idx);
        if (!export_frame(fname_buff)) return -1;
        #ifdef _OPENMP
            timer_end = omp_get_wtime();
            time_elapse = timer_end - timer_start;
        #else
            timer_end = clock();
            time_elapse = (double)(timer_end - timer_start) / CLOCKS_PER_SEC;
        #endif
        printf("Time taken: %f seconds\n", time_elapse);
        fflush(stdout);
    }
    
    unload_plate();
    return 0;
}
