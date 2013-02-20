
class Cparameter
{
public:  

	double MODULE_N;				/**< Scalar elastic modulus associated to normal deflection.*/
	double friction_coefficient;
	double tang_constant; 			/**< Numerical constant for tangantial deflection.*/        
	double roll_constant;			/**< Numerical constant for rolling and twist.*/    

	double specific_heat;           /**< Specific Heat of the bulk material (J/kg/K).*/
	double bulk_conductivity;       /**< Bulk conductivity (W/m/K).*/
	double thermal_expansion;		/**< Thermal expansion coefficient (1/K).*/
	string GSD;						/**< Type of grains size distribution, can be FRACTALE or UNIFORM.*/
	double fractal_dim;				/**< Fractal dimension power of the grains size distribution.*/

	double t_inertia;				/**< Inertial time P \f$ \sqrt{ \frac{m}{Pd} }\f$. */
	double t_collision;				/**< Collision time E \f$ \sqrt{ \frac{m}{Ed} }\f$ */
	double t_thermal;				/**< Thermal time \f$ \frac{mc} {d k} \f$. */
	double t_shear;					/**< Shear time \f$ \frac{1}{\dot \gamma} \f$*/
	
	double RHOmean;					/**< Mean density.*/
	double Dmean;					/**< Unit particle diameter, defined as lenght unit.*/
	double Dmin;					/**< Smaller particle diameter*/
	double Dmax;					/**< Higher particle  diameter*/
	double Mmean;					/**< Unit particle mass, defined as mass unit*/
	double Mmin;					/**< Smaller particle mass*/
	double Mmax;					/**< Higher particle mass*/
	
	double total_mass;
	double average_temperature;
	
	double K;						/**< Dimentionless stiffness: \f$ \kappa = \frac{P}{E^n}\f$. */
	double J;						/**< Dimentionless thermal time: \f$ J = thermal_time/inertial time. \f$. */
	double I;						/**< Dimentionless inertial number: \f$ I_i= \dot \gamma \sqrt{\frac{m}{PR}}\f$. */

	double COMP_FRACTION;
	double INITIAL_SATURATION;
	double FIXED_SATURATION;
	double SURFACE_TENSION;
	double LIQUID_DIFFUSION;
	double CONTACT_ANGLE;
	
	void dimensionless_number(Ccell &, std::vector <Cparticle> &); 
	friend ofstream &operator<<(ofstream &,Cparameter); 
	friend ifstream & operator>>(ifstream &,Cparameter &);
};
