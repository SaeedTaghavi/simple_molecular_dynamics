#ifndef _Part_
#define _Part_

#include <iostream>

class Part {
private:
	double pos[3];
	double vel[3];
	double speed;
	double prevForce[3];//used for verlet step
	double force[3];
	double potEnergy;
	double mass;

public:
	/***
	 * Constructors and destructor
	 */
	Part();
	Part(double m, double x, double y, double z);
	Part(double m, double x, double y, double z, double vx, double vy,
			double vz, double fx, double fy, double fz);
	~Part();


	/***
	 * Getters and Setters
	 */
	double* getPosition() ;
	double* getVelocity() ;
	double* getForce() ;
	double getMass() const ;
	double getPotentialEnergy() const ;
	double getKineticEnergy() const ;
	double getSpeed() const ;

	void setPosition(double x, double y, double z);
	void setVelocity(double x, double y, double z);
	void setForce(double x, double y, double z);
	void setMass(double m);
	void setPotentialEnergy (double e);


	/***
	 * Others methods
	 */
	double* distanceXYZ(Part p) const;

	//Advance particle a time step bassed on Verlet Integration
	void verletPositionStep(double time);
	void verletVelocityStep(double time);

	Part operator =(const Part &);
	friend std::ostream & operator <<(std::ostream &, const Part &);
};





/***
 * Inline methods
 */


inline double* Part::getPosition() {
	return pos;
}

inline double* Part::getVelocity() {
	return vel;
}

inline double* Part::getForce() {
	return force;
}

inline double Part::getMass() const {
	return mass;
}

inline double Part::getPotentialEnergy() const {
	return potEnergy;
}

inline double Part::getKineticEnergy() const {
	return mass*speed*speed/2;
}

inline double Part::getSpeed() const {
	return speed;
}

inline void Part::setPosition(double x, double y, double z){
	pos[0]=x;pos[1]=y;pos[2]=z;
}

inline void Part::setVelocity(double x, double y, double z){
	vel[0]=x;vel[1]=y;vel[2]=z;
	speed = sqrt(x*x +y*y+z*z);
}

inline void Part::setForce(double x, double y, double z){
	prevForce[0]=force[0];prevForce[1]=force[1];prevForce[2]=force[2];
	force[0]=x;force[1]=y;force[2]=z;
}

inline void Part::setMass(double m){
	mass=m;
}

inline void Part::setPotentialEnergy (double e){
	potEnergy=e;
}


#endif

