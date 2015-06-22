#include <iostream>
#include <cmath>
#include <cstring>
#include "Part.hpp"

using namespace std;

Part::Part(double m, double x, double y, double z) {
	setPosition(x, y, z);
	speed = -1;
	setMass(m);
	setPotentialEnergy(-1);
}

Part::Part(double m, double x, double y, double z, double vx, double vy,
		double vz, double fx, double fy, double fz) {
	setPosition(x, y, z);
	setVelocity(vx, vy, vz);
	setForce(fx, fy, fz);
	setMass(m);
	setPotentialEnergy(-1);
}

Part Part::operator =(const Part &p) {
	memcpy((void*) this, (const void*) &p, sizeof(Part));
	return *this;
}

Part::Part() {}

Part::~Part() {}

double* Part::distanceXYZ(Part p) const{
	double *d = new double[3];
	for (int i = 0; i < 3; i++)
		d[i] = pos[i] - p.getPosition()[i];
	return d;
}

void Part::verletPositionStep(double time) {
	double x, y, z;
	x = pos[0] + time * vel[0] + time * time * force[0] / (2.0 * mass);
	y = pos[1] + time * vel[1] + time * time * force[1] / (2.0 * mass);
	z = pos[2] + time * vel[2] + time * time * force[2] / (2.0 * mass);
	setPosition(x, y, z);
}

void Part::verletVelocityStep(double time) {
	double vx, vy, vz;
	vx = vel[0] + time * (force[0] + prevForce[0]) / (2.0 * mass);
	vy = vel[1] + time * (force[1] + prevForce[1]) / (2.0 * mass);
	vz = vel[2] + time * (force[2] + prevForce[2]) / (2.0 * mass);
	setVelocity(vx, vy, vz);
}


std::ostream & operator <<(std::ostream & os, const Part &p) {
	os << "Au " << p.mass << "\t" << p.pos[0] << " " << p.pos[1] << " " << p.pos[2]
			<< " \t" << p.vel[0] << " " << p.vel[1] << " " << p.vel[2]
			<< " \t" << p.force[0] << " " << p.force[1] << " " << p.force[2];
	return os;
}
