#ifndef _Cont_
#define _Cont_
#include "Part.hpp"
#include <iostream>

/* Class Cont
 * A cubic container of particles with volume
 * equal to side^3.
 * Center in (0,0,0).
 * Periodic Boundary condition.
 * Based on Lennard-Jones energy and
 * Verlet time steps.
 */
class Cont {
private:
	int n;
	double side;
	Part * part; //particle array
	double potentialEnergy;

public:
	/***
	 * Constructors and destructor
	 */
	Cont(double side, int part_count);
	~Cont();

	/***
	 * Getters and Setters
	 */
	int getN() const;
	double getSide() const;
	Part& getParticle(int i) const;
	double getPotentialEnergy() const;
	double getKineticEnergy() const;

	void setParticle(Part p, int i) const;

	/***
	 * Other Methods
	 */
	void lennard_jones();
	void nextPositionStep(double tau);
	void nextVelocityStep(double tau);

	void periodicBoundaryCondition();

	Cont & operator =(const Cont &);
	friend std::ostream & operator <<(std::ostream &,const Cont &);

};



/***
 * Inline Methods
 */

inline int Cont::getN() const {
	return n;
}

inline double Cont::getSide() const{
	return side;
}

inline void Cont::setParticle(Part p, int i) const {
	part[i] = p;
}

inline Part& Cont::getParticle(int i) const {
	return part[i];
}

inline double Cont::getPotentialEnergy() const{
	return potentialEnergy;
}

#endif
