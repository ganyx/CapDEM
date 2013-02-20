/*! \class Cparticle 
*  \brief This class contains all the data relative to a particle, as weel as the predictor-corrector functions to integrate its motionover a time step.
*/


void Intergre_EULER2(Cvector&, Cvector&, Cvector&, double );

class Cnode;
class Cparticle
{
public:
	int id;
	Cvector X;             /**< Particle position. */ 
	Cvector V;             /**< Particle velocity.*/ 
	Cvector A;             /**< Particle acceleration. */ 
	Cvector Ome;           /**< Particle angular velocity.*/
	Cvector OmeDot;        /**< Particle angular acceleration. */ 
	Cvector Fsum;          /**< Sum of forces the particle is subjected to.*/
	Cvector Gsum;          /**< Sum of moment the particle is subjected to.*/

	Cvector Fext;          /**< External forces the particle is subjected to.*/
	Cvector Gext;          /**< External moment the particle is subjected to.*/

	double R;              /**< Particle radius. */
	double m;              /**< Particle mass. */
	double J;              /**< Particle moment of inertia.*/
	double RHO;				// for no expansion
	
	double E;				/**<Young modulus.*/
	double k;				/**<Thermal conductivity.*/
	double c;				/**<Heat capacity.*/
	double e;				/**<Thermal expantion coefficient.*/
	
	double Tm;
	double Lm;
	
	double T;				/**<Particle temperature.*/
	double Tdot;			/**<Particle temperature rate.*/
	double phi; 			/**<Particle total heat rate.*/
	double phi_ext;			/**<Particle external heat rate.*/
	double production;		/**<Total heat production by contact dissipation.*/
	double L;				/**< Latent heat */
//	bool phase;				/**< Phase identifier, 0-solid, 1-liquid */
//	bool interphase;			/**< 1-solid>>liquid or liquid>>solid, 0-no change */

	std::vector <Cparticle *> neighbour;	/**<List of particle within the verlet discance*/
	std::vector <Ccontact *> contact;
	std::vector <Cbox *> my_box;
	int AM_I_BOUNDARY;		/**<0 if it's a flowing particle, -1 or 1 for bottom and top boundary particle. -2 2 for bottom and top walls*/

	Cparticle(void);
	void PRINT();
	void set_me_in_main_cell(Ccell &);			/**<set a particle back in the cell if there is periodic boundary*/
	void set_in_box(Cmesh &mesh,Ccell &);
	void remove_from_box();
	void get_neighbour(Ccell &cell) ;
	
	
	void predictor(double dt, double dt2_on_2);	/**< Integrate the velocity and the acceleration to get the new postion. \warning Depend on wether there is inertia or not.*/
	void corrector(double dt_on_2,Ccell &); 	/**< Get the new velocity and acceleration \warning Only with inertia. */
	void expand_radius(double dt);						/**<Expand the radius of a particle according to its temperature and its coef. of thermal expansion. */
	friend ofstream &operator<<(ofstream &,Cparticle);		/**< Print the particle's data into a file.*/ 
	friend ifstream & operator>>(ifstream &,Cparticle &);	/**< Read the particle's data from a file.  */ 
	//bool operator == (Cparticle &);
	
//	bool member;			/**< If this particle belongs to an aggregate, YG.*/
	
	double RS, dRS, RS_old;
	
	//voronoi tesselation informations
	vector<double> v_face_areas;
	vector<int> v_neighbors;
	double voronoi_volume; 
	int num_neighbors;
	
	// water content
	double grain_volume;
	double water_volume;
	double water_volume_old;
	double void_volume;
	double saturation;
	double water_pressure;
	double positive_pressure;
	double N_water_bridge;
	double sum_vij;
	Cvector mass_transfer;
};





