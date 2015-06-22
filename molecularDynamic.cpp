#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <list>
#include "constants.hpp"
#include "Part.hpp"
#include "Cont.hpp"
#include "rand.hpp"

using namespace std;


list<Part>* read_files();
int read_stdin();
void volume_snapshot(const Cont &volume,float t_sim);
void draw_box();


/* FILES */
ofstream *f_etotal, *f_ekinetic, *f_epotential, *f_forces, *f_positions;
ifstream *init_state;

int main() {
	/***
	 * IOFILES INIT
	 */
	 f_positions = new ofstream("out/positions.xyz");
	 f_etotal = new ofstream("out/etotal.dat");
	 f_ekinetic = new ofstream("out/ekinetic.dat");
	 f_epotential = new ofstream("out/epotential.dat");
	 f_forces = new ofstream("out/forces.dat");

	 init_state = new ifstream("init_state.xyz");

	/***
	 * INIT CONTAINER
	 */
	list<Part> *particles = read_files();
	int extra = read_stdin();

	//box of side L with file* + extra particles
	Cont volume(L, extra + particles->size());

	//init file particles
	int k = 0;
	for (list<Part>::iterator it = particles->begin(); it != particles->end();
			k++, it++) {
		it->setVelocity(drand(-MAX_V, MAX_V), drand(-MAX_V, MAX_V),
				drand(-MAX_V, MAX_V));
		volume.setParticle(*it, k);
	}
	delete particles;

	/***
	 * INIT SIMULATION
	 */
	time_t t_exe = time(NULL);
	double t_sim = 0;

	//calculate initial forces and energy
	volume.lennard_jones();
	double pot = volume.getPotentialEnergy();
	double kin = volume.getKineticEnergy();

	//print initial state of system
	cout.precision(4);
	cout << "pot = " << pot << "    kin =" << kin << "    total = " << pot + kin
			<< endl;
	cout << volume << endl;

	/***
	 * INTEGRATION
	 */
	for (int frame = 0; frame < T_STOP; frame++) {
		//Move Particles
		volume.nextPositionStep(TAU);
		//recalculate forces and energy
		volume.lennard_jones();
		//Accelerate particles
		volume.nextVelocityStep(TAU);
		//Container Boundaries
		volume.periodicBoundaryCondition();

		/***
		 * SNAPSHOT OF SYSTEM
		 */
		if (frame % SNAPSHOT_STEP == 0)
			volume_snapshot(volume, t_sim);

		//Progress
		if (frame % ((int) T_STOP / 10) == 0)
			cout << "Progress: " << frame / ((int) T_STOP / 10) * 10 << "%"
					<< endl;

		t_sim += TAU;

	}
	cout << "Progress: 100%" << endl;

	/***
	 * Final Output
	 */
	time_t ttotal = time(NULL) - t_exe;
	tm tm2 = *gmtime(&ttotal);
	cout << "Tiempo transcurrido en la simulacion: "
			<< difftime(time(NULL), t_exe) << " segundos" << endl;
	cout << "Los datos se han guardado en el archivo positions.xyz " << endl;
	cout << "Los datos de energia cinetica vs tiempo estan en ekinetic.dat "
			<< endl;
	cout << "Los datos de energia potencial vs tiempo estan en epotential.dat "
			<< endl;
	cout << "Los datos de energia total vs tiempo estan en etotal.dat "
			<< endl;
	                
	cout << "Este programa calculo en " << tm2.tm_hour << ":" << tm2.tm_min
			<< ":" << tm2.tm_sec << endl;

	/***
	 * Closing files
	 */
	f_positions->close();
	f_etotal->close();
	f_epotential->close();
	f_ekinetic->close();
	f_forces->close();

	delete f_positions;
	delete f_etotal;
	delete f_epotential;
	delete f_ekinetic;
	delete f_forces;
	delete init_state;
	return 0;
}
/*** END OF MAIN ***/

list<Part>* read_files() {
	list<Part> *particles = new list<Part>();

	//read particles position from input files
	while (!init_state->eof()) {
		double x, y, z;
		*init_state >> x >> y >> z;
		particles->push_back(Part(MASS, x, y, z));
	}

	cout << "Read " << particles->size() << " elements from file "
			<< "entrada.xyz" << endl;

	init_state->close();
	return particles;
}

int read_stdin() {
	int extras;

	cout << "Extra particles in dynamic: ";
	cin >> extras;
	return extras;
}

void volume_snapshot(const Cont &volume,float t_sim) {
	double pot = volume.getPotentialEnergy();
	double kin = volume.getKineticEnergy();
	double forces = 0;
	//calculate force of gravity center
	for (int m = 0; m < volume.getN(); m++) {
		double* f = volume.getParticle(m).getForce();
		forces = f[0] + f[1] + f[2];
	}

	//print on different files
	*f_etotal << t_sim << " " << kin + pot << endl;
	*f_ekinetic << t_sim << " " << kin << endl;
	*f_epotential << t_sim << " " << pot << endl;
	*f_forces << t_sim << " " << forces << endl; //con esto veo la fuerza del centro de masas

	//print on *.xyz file format
	*f_positions << volume.getN() + 8 << endl;
	*f_positions << "E=" << pot + kin << " t=" << t_sim << endl;
	for (int i = 0; i < volume.getN(); i++) {
		double* pos = volume.getParticle(i).getPosition();
		*f_positions << "Au " << pos[0] << " " << pos[1] << " " << pos[2]
				<< endl;
	}
	draw_box();
}

void draw_box(){
	*f_positions<<"Au "<<-L/2<<" "<<-L/2<<" "<<-L/2<<endl;
	*f_positions<<"Au "<<L/2<<" "<<L/2<<" "<<L/2<<endl;
	*f_positions<<"Au "<<-L/2<<" "<<L/2<<" "<<-L/2<<endl;
	*f_positions<<"Au "<<L/2<<" "<<-L/2<<" "<<L/2<<endl;
	*f_positions<<"Au "<<-L/2<<" "<<-L/2<<" "<<L/2<<endl;
	*f_positions<<"Au "<<L/2<<" "<<L/2<<" "<<-L/2<<endl;
	*f_positions<<"Au "<<-L/2<<" "<<L/2<<" "<<L/2<<endl;
	*f_positions<<"Au "<<L/2<<" "<<-L/2<<" "<<-L/2<<endl;

}



