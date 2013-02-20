/**
\class Ccell
\brief This class contains all the info about the simulation domain and its time evolution.
Periodic dupication are applied virtualy in dupicating the cell in all direction.
*/
class Ccell
{
	
	public:  
  
	string boundary;/**< Type of boundary condition, value can be PERIODIC_SHEAR, WALL_INCLINED, WALL_SHEAR*/
	Cvector L; 		/**< Size of the cell*/ 

	Cparticle *plan_bottom, *plan_top;

	double mass;
	Cvector rigid_velocity;  
	double normal_stress_ext,shear_stress_ext;
	double normal_stress_in, shear_stress_in;
	double shear_rate,dilat_rate;
	double cumul_strain;
	
	double dstress, dstrain_int;
	
	double Ashear,Ashear_p, Adilat,Vshear,Vdilat;
	double DeltaT,gradT;
	double Xshift; 	/**< Sifh of the bottom and upper cells when shearing*/
	double Yshift;
	
	// work input control variables
	bool shear_work_control;
	double shear_work_input;
	double shear_acc, shear_stress_rate;
	double shear_stress_in_old;
	
	// stick-slip control
	bool stick_slip;
	double slip_velocity;

	double gravity;				/**<Guess what...*/
	double slope;				/**<slope angle on the cell, in degree.*/ 

	bool normal_stress_control, shear_stress_control, gradT_control;

	Cmatrix stress;

	Cmatrix K; 		/**< Effective conductivity of the sample>*/
	Cmatrix H; 		/**< Effective convectivity of the sample>*/
	double energy_kinetic;
	double production;
	double coordination;
 	double solid_fraction;
	Ccell();			/** <Set every one to zero as default value*/ 
	
	void predictor(double,double);		/**< Move the cell.*/
	void corrector(double);  			/**< Get the velocity of the cell.*/
	
	void rescale(Cvector &X);		/**< Rescale a vector positon within the periodic domain*/
	void rescale(double &,double);	/**<Rescale a scalar positon within the periodic domain*/
	void PRINT(); 					/**< Print some cell properties on screen.*/
	friend ofstream &operator<<(ofstream &,Ccell);
	friend ifstream & operator>>(ifstream &,Ccell &);
};



