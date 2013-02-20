//======================================
// Class contact data array and functions
// Functions details are in "contact.c"
//======================================

class Ccontact
{
	Cvector uDotT,dOmeN,dOmeT,nADot,OmeMean;
	
	public: 

	int A,B;
	Cparticle *pA;			/**< particle at the origin of the contact*/
	Cparticle *pB;			/**< particle at the end of the contact*/
	Cparameter *parameter;
	Ccell *cell;
	double dT;				/**<Diference of temperature between the two particles*/
	
	Cvector dX;				/**< Center to center vector*/
	double dx;				/**<centr to center distance*/
	Cvector dV;				/**< Center to center relative velocity*/
	
	Cvector nA;		/**< Unit vector normal to the surface.*/
	
	Cvector RA;		/**<Vector from center of A to its surface along nA.*/
	Cvector RB;		/**<Vector from center of B to its surface along nA.*/

	Cmatrix alpha;	/**<Normal projector.*/ 
	Cmatrix I_m_alpha;/**<Tangential projector.*/

	bool  TANGENTIAL_SLIDE, TWIST_SLIDE,ROLLING_SLIDE;	
	double phi;				/**<Heat rate through the contact.*/
	double production;		/**<Heat production of the contact.*/
	double production_slide,production_normal,  production_rolling, production_twist;
	double conductivity;
	
	Cvector F; 				/**< Contact force*/
	Cvector Fn, Ft,Fvis,Fela;
		
	Cvector G,Gn,Gt;
	double fn,ft,gn,gt;

	double deltaN; /** Normal deflection.*/
	double a; /**< Radius of the interacting region.*/
	double age;/**< Age of the contact.*/
	
	double E;
	double ct,cr;
	double mu;
	double Reff;
	double meff;
	
	
	void EVALE_Geo();								/**< Evaluate the geometry of contact (unit normal vector and projectors).*/
	void increment_force(double dt);
	bool rescale_slide(Cvector &, double &, double);/**< Limit the norm of a vector to a givien maximum, used to compute apply sliding criteria. It returns "true" if the vector's norm was to large */
	void relative_velocity( );
	void EVALE_heat_flow(); 						/** Get the total flux of heat rate through the contact (Watt)*/
	bool AM_I_CONTACTING();							/**<Return true if there is a contact, false if not*/
	void set_me_in_main_cell();
	void PRINT();
	Ccontact();
	Ccontact(Cparticle *,Cparticle *, Ccell *, Cparameter *);
	
	bool member;	/** Member contact for an aggregate, YG<.*/
//	double tm;		/** Melting time for calculating F0, counting tm>0 as bonded contact */
	double F0;
	double aB; 		// bonding radius;
	double deltaNB;	// bonding overlap;
	
	bool Flag_Boundary;	
	
	
	//void operator = (Ccontact ); /**< Replace the contact by another contact para. */
	
	friend ofstream &operator<<(ofstream &,Ccontact);
	friend ifstream & operator>>(ifstream &,Ccontact &);
	
	double voronoi_area;
	double water_volume;
	double water_volume_old;
	double water_vij;
	double dot_water_volume;
	double dwater_volume;
	double fcap; // the magnitude of capillary force
	double CONTACT_ANGLE;
};
 
