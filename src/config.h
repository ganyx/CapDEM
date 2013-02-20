
class Cconfig
{
	void set_random_grain(int,int,int);					/**< set particle randomly in a cell (flowing-grains as well as wall-grains). Take as parameters the number of flowing grains, the number of grains on the bot wall and on the top wall*/
   
	void predictor();			/**< Integrate the particle motion accroding to their acceleration.*/
	void renew_contact_list();	/**< Get the new contact list.*/ 
	void sum_force();			/**< Sum force and moment for each particles.*/ 
	void corrector();			/**< Get the new acceleration of particles, correct their velocities.*/ 
	// void solve_mobility();

	public: 

	void sum_heat();		/**< Get the sum of heat flux for each particles.*/
	void update_contact();	/**< General function to get the new contact list and measure the contact force/torque.*/
	void energy(); 

 	double Etot,Eela,Eforce,Emoment, Ekin,Erot,Etrans;
   
	double t,dt, dt_on_2,dt2_on_2;

	Cparameter parameter;		/**< Physical parameters.*/
	Ccell cell;					/**< Bondary condition. */
	std::vector <Ccontact>  C;		/**< List of contact.*/
	std::vector <Ccontact>  CThread[NTHREADS]; /**< Temperaty contact lists for threads for MPI, YG */
	std::vector <Cparticle> P;		/**< List of particle.*/
//	QList <Caggregate> G;		/**< List of aggregate, YG.*/
//	int Number_P_in_G;
	
	bool simule_thermal_expansion;	/**< if true, the thermal expansion of particle will simulate.*/
	bool simule_thermal_production;	/**< if true, the motion of particle will simulate.*/
	bool simule_thermal_conduction;	/**< if true, the heat transfer through contact will simulate.*/
  
	void Evale_conductivity_tensor();	/**< Measure the conductivity tensor*/
  
	void iterate(double);			/**< Make the configuration evolve over a time step dt. */
	void update_particle();


	void create_random();		/**< Create a new random config. */
	void set_wall_grain(int, int);/**< Set the grain of the walls as well as the planes if any. */
	void set_radius(std::vector <double> &radius, int &Nf);/**< Create a list of radius with either a uniform distribution of fractale distribution. */
	void fread(Cin_out); 		/**<Read from a file.*/
//	void fprint(Cin_out);		/**<Print in a file.*/
	void fprint(Cin_out);		/**<Print in a file.*/
	
//	void renew_aggregate();		/**< Update the aggregates */
//	bool merge_aggregate(Ccontact &);
//	bool breakage_aggregate(Caggregate &);
	
	double heat_in, heat_out;
	double PN, PS, PT, PR;
	
	std::vector <voronoicell_neighbor> voro;	
	bool Voronoi_Update;
	void liquid_transfer();
	void drag_force();
	double cap_pressure;
	double saturation;
	double water_content;
	
	double cap_pressure_mid;
	double saturation_mid;
	double water_content_mid;
	
	bool flag_wetting;
	double MAX_SCAN;
	double MIN_SCAN;
	double GRAVITY;
};






