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
public:
	std::vector<Body> bodies;
	double _x0, _y0, _x1, _y1;

public:
	NBodySystem() {
	}

	void setBounds(double x0, double y0, double x1, double y1) {
		_x0 = x0;
		_y0 = y0;
		_x1 = x1;
		_y1 = y1;
	}

	void addBody(Body& b) {
		bodies.push_back(b);
	}

	void advance(double dt) {
		const unsigned N = (bodies.size() - 1) * bodies.size() / 2;
		struct R {
			double dx, dy, dz, unused;
		};

		for (unsigned i = 0, k = 0; i < bodies.size() - 1; ++i) {
			Body& iBody = bodies[i];
			for (unsigned j = i + 1; j < bodies.size(); ++j, ++k) {
                struct R r;

				r.dx = iBody.x - bodies[j].x;
				r.dy = iBody.y - bodies[j].y;
				r.dz = iBody.z - bodies[j].z;

	            double d2 = r.dx * r.dx + r.dy * r.dy + r.dz * r.dz;
	            double mag = dt / (d2 * sqrt(d2));

	            iBody.vx -= r.dx * bodies[j].mass * mag;
				iBody.vy -= r.dy * bodies[j].mass * mag;
				iBody.vz -= r.dz * bodies[j].mass * mag;

				bodies[j].vx += r.dx * iBody.mass * mag;
				bodies[j].vy += r.dy * iBody.mass * mag;
				bodies[j].vz += r.dz * iBody.mass * mag;
			}
		}

		for (unsigned i = 0; i < bodies.size(); ++i) {
			double bx = bodies[i].x + dt * bodies[i].vx;
			double by = bodies[i].y + dt * bodies[i].vy;
			double bz = bodies[i].z + dt * bodies[i].vz;
			if (bx < _x0) {
				bodies[i].x = _x0;
			} else if (bx >= _x1) {
				bodies[i].x = _x1;
			} else {
				bodies[i].x = bx;
			}
			if (by < _y0) {
				bodies[i].y = _y0;
			} else if (by >= _y1) {
				bodies[i].y = _y1 - 1;
			} else {
				bodies[i].y = by;
			}
		}
	}

public:
	double energy() {
		double dx, dy, dz, distance;
		double e = 0.0;

		for (unsigned i = 0; i < bodies.size(); ++i) {
			Body& iBody = bodies[i];

			// Kinetic Energy
			e += 0.5 * iBody.mass
					* (iBody.vx * iBody.vx + iBody.vy * iBody.vy
							+ iBody.vz * iBody.vz);

			for (unsigned j = i + 1; j < bodies.size(); ++j) {
				Body& jBody = bodies[j];
				dx = iBody.x - jBody.x;
				dy = iBody.y - jBody.y;
				dz = iBody.z - jBody.z;

				// Potential Energy
				distance = sqrt(dx * dx + dy * dy + dz * dz);
				if (distance == 0) {
					e -= (iBody.mass * jBody.mass) / distance;
				}
			}
		}
		return e;
	}
};

#endif /* defined(__imageParticles__NBody__) */
