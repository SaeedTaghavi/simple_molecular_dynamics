#include <iostream>
#include <cmath>
#include "Cont.hpp"
#include "Part.hpp"
#include "rand.hpp"
#include "constants.hpp"
#include <omp.h>

using namespace std;

Cont::Cont(double l,int p) :
		n(p),side(l)
{
	part = new Part[n];
	double q = side / 2, v = MAX_V;
	for (int j = 0; j < n; j++) {
		double x = drand(-q, q);
		double y = drand(-q, q);
		double z = drand(-q, q);
		double vx = drand(-v, v);
		double vy = drand(-v, v);
		double vz = drand(-v, v);
		part[j] = Part(MASS, x, y, z, vx, vy, vz, 0.0, 0.0, 0.0);
	}
}


Cont::~Cont() {
	delete[] part;
}


/*** LENNARD-JONES
 * Recalculate forces for all particles in Cont.
 * Recalculate Potential Energy for particles
 * and Total for the container.
 */
void Cont::lennard_jones() {
	cout.precision(15);
	double *d_xyz, d;
	double total_energy = 0.0, particle_energy = 0.0;
	double fdist;
	double Fx,Fy,Fz;
	double q = side / 2;


	#pragma omp parallel for \
		private (Fx,Fy,Fz,particle_energy,d_xyz,d,fdist)\
		schedule(dynamic)\
		reduction(+:total_energy)
	for (int i = 0; i < n; i++) {
		Fx = 0.0;
		Fy = 0.0;
		Fz = 0.0;
		for (int j = 0; j < n; j++) {
			if (i != j) {
				d_xyz = part[i].distanceXYZ(part[j]);

				/*Periodic Boundary Condition*/
				if (d_xyz[0] > q) {
					d_xyz[0] = d_xyz[0] - side;
				} else if (d_xyz[0] <= -q) {
					d_xyz[0] = d_xyz[0] + side;
				}

				if (d_xyz[1] > q) {
					d_xyz[1] = d_xyz[1] - side;
				} else if (d_xyz[1] <= -q) {
					d_xyz[1] = d_xyz[1] + side;
				}

				if (d_xyz[2] > q) {
					d_xyz[2] = d_xyz[2] - side;
				} else if (d_xyz[2] <= -q) {
					d_xyz[2] = d_xyz[2] + side;
				}

				d =	d_xyz[0] * d_xyz[0] + d_xyz[1] * d_xyz[1]
								+ d_xyz[2] * d_xyz[2];

				double ssd6 = SIGMA6 / (d * d * d );
				double e_pp = 2.0 * EPSILON * (ssd6 * ssd6 - ssd6); //considerado ij y ji
				particle_energy += e_pp;
				
				d = sqrt(d);
				double f = 48.0 * EPSILON
				  * (1.0 / d * (ssd6 * ssd6 - 0.5 * ssd6)); //derivada de LJ normal

				fdist = f / d;
				Fx += fdist * d_xyz[0];
				Fy += fdist * d_xyz[1];
				Fz += fdist * d_xyz[2];

				delete[] d_xyz;

			}
		}
		part[i].setForce(Fx, Fy, Fz);
		part[i].setPotentialEnergy(particle_energy);
		total_energy += particle_energy;
		particle_energy=0.0;
	}

	potentialEnergy = total_energy;
}

double Cont::getKineticEnergy() const{
	double k=0.0;
	for(int i=0;i<n;i++)
		k+=part[i].getKineticEnergy();
	return k;
}


void Cont::nextPositionStep(double tau) {
	for (int i = 0; i < n; i++) {
		part[i].verletPositionStep(tau);
	}
}

void Cont::nextVelocityStep(double tau) {
	for (int i = 0; i < n; i++) {
		part[i].verletVelocityStep(tau);
	}
}

/***
 * Periodic Boundary Condition
 * Relocate particles out of bounds of
 * container basses on periodic bounds.
 */
void Cont::periodicBoundaryCondition() {
	double q = side / 2;
	for (int i = 0; i < n; i++) {
		Part &p = part[i];
		const double *pos = p.getPosition();
		double newX=pos[0], newY=pos[1], newZ=pos[2];
		if (pos[0] > q) {
			newX = pos[0] - side;
		} else if (pos[0] < -q) {
			newX = pos[0] + side;
		}

		if (pos[1] > q) {
			newY = pos[1] - side;
		} else if (pos[1] < -q) {
			newY = pos[1] + side;
		}

		if (pos[2] > q) {
			newZ = pos[2] - side;
		} else if (pos[2] < -q) {
			newZ = pos[2] + side;
		}

		p.setPosition(newX, newY, newZ);
	}
}

/*** CONDICION DE BORDE CERRADO
		 for(int i=0;i<p;i++){
		 int q=20;

		 if(a[i].getx()>q){
		 a[i].setvx(-a[i].getvx());
		 a[i].setx(2*q-a[i].getx());}
		 else if(a[i].getx()<-q){
		 a[i].setvx(-a[i].getvx());
		 a[i].setx(-2*q-a[i].getx());}


		 if(a[i].gety()>q){
		 a[i].setvy(-a[i].getvy());
		 a[i].sety(2*q-a[i].gety());}
		 else if(a[i].gety()<-q){
		 a[i].setvy(-a[i].getvy());
		 a[i].sety(-2*q-a[i].gety());}

		 if(a[i].getz()>q){
		 a[i].setvz(-a[i].getvz());
		 a[i].setz(2*q-a[i].getz());}
		 else if(a[i].getz()<-q){
		 a[i].setvz(-a[i].getvz());
		 a[i].setz(-2*q-a[i].getz());}

		 volume.setparti(a[i],i);
		 }
*/

Cont & Cont::operator =(const Cont & c) {
	n = c.getN();
	delete[] part;
	side = c.getSide();
	part = new Part[c.getN()];
	potentialEnergy = -1;
	for (int j = 0; j < n; j++) {
		part[j] = c.getParticle(j);
	}
	return *this;
}

std::ostream & operator <<(std::ostream & os, const Cont &c) {
	int n = c.getN();
	for (int j = 0; j < n; j++) {
		os << c.getParticle(j) << endl;
	}
	return os;
}
