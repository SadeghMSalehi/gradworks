//
//  NBody.h
//  imageParticles
//
//  Created by Joohwi Lee on 10/29/12.
//
//

#ifndef __imageParticles__NBody__
#define __imageParticles__NBody__

#include <iostream>
#include "armadillo"

class Body {
public:
    double x, y, z, filler, vx, vy, vz, mass;
    Body() {
        x = y = z = filler = vx = vy = vz = mass = 0;
    }
};

class NBodySystem {
private: std::vector<Body> bodies;

public:
    NBodySystem() {
    }

    void advance(double dt) {
        const unsigned N = (bodies.size() - 1) * bodies.size() / 2;
        struct R {
            double dx,dy,dz,unused;
        };
        static R r[1000];
        static __attribute__((aligned(16))) double mag[1000];

        for (unsigned i = 0, k = 0; i < bodies.size() - 1; ++i) {
            Body& iBody = bodies[i];
            for(unsigned j = i + 1; j < bodies.size(); ++j, ++k) {
                r[k].dx = iBody.x - bodies[j].x;
                r[k].dy = iBody.y - bodies[j].y;
                r[k].dz = iBody.z - bodies[j].z;
            }
        }

        for (unsigned i = 0; i < N; i += 2) {
            __m128d dx,dy,dz;
            dx = _mm_loadl_pd(dx, &r[i].dx);
            dy = _mm_loadl_pd(dy, &r[i].dy);
            dz = _mm_loadl_pd(dz, &r[i].dz);

            dx = _mm_loadh_pd(dx, &r[i+1].dx);
            dy = _mm_loadh_pd(dy, &r[i+1].dy);
            dz = _mm_loadh_pd(dz, &r[i+1].dz);

            __m128d dSquared = dx*dx + dy*dy + dz*dz;
            __m128d distance = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(dSquared)));
            for (unsigned j = 0; j < 2; ++j) {
                distance = distance * _mm_set1_pd(1.5) - ((_mm_set1_pd(0.5) * dSquared) * distance) * (distance * distance);
            }

            __m128d dmag = _mm_set1_pd(dt) / (dSquared) * distance;
            _mm_store_pd(&mag[i], dmag);
        }

        for (unsigned i = 0, k = 0; i < bodies.size() - 1; ++i) {
            Body& iBody = bodies[i];
            for (unsigned j = i + 1; j < bodies.size(); ++j,++k) {
                iBody.vx -= r[k].dx * bodies[j].mass * mag[k];
                iBody.vy -= r[k].dy * bodies[j].mass * mag[k];
                iBody.vz -= r[k].dz * bodies[j].mass * mag[k];

                bodies[j].vx += r[k].dx * iBody.mass * mag[k];
                bodies[j].vy += r[k].dy * iBody.mass * mag[k];
                bodies[j].vz += r[k].dz * iBody.mass * mag[k];
            }
        }

        for (unsigned i = 0; i < bodies.size(); ++i) {
            bodies[i].x += dt * bodies[i].vx;
            bodies[i].y += dt * bodies[i].vy;
            bodies[i].z += dt * bodies[i].vz;
        }
    }

public:
    double energy() {
        double dx, dy, dz, distance;
        double e = 0.0;

        for (unsigned i = 0; i < bodies.size(); ++i) {
            Body& iBody = bodies[i];

            // Kinetic Energy
            e += 0.5 * iBody.mass * (iBody.vx * iBody.vx + iBody.vy * iBody.vy + iBody.vz * iBody.vz);

            for (unsigned j = i + 1; j < bodies.size(); ++j) {
                Body& jBody = bodies[j];
                dx = iBody.x - jBody.x;
                dy = iBody.y - jBody.y;
                dz = iBody.z - jBody.z;

                // Potential Energy
                distance = sqrt(dx*dx + dy*dy + dz*dz);
                if (distance == 0) {
                    e -= (iBody.mass * jBody.mass) / distance;
                }
            }
        }
        return e;
    }
};


#endif /* defined(__imageParticles__NBody__) */
