
/**This class provide statistical info about an array of data. */



class Cstat
{
	public:
	std::vector <double> data;
	
	double mean, std, min,max;
	
	double MIN();           /**< Minimum value of the array.*/
	double MAX();           /**< Maximum value of the array.*/
	double MEAN();          /**< Mean value of the array.*/
	double STD();   /**< Standart deviation of the array.*/
 	//void PRINT(); 
};


class Cdistribution: public Cstat
{
public:
	double step;
	std::vector <double> distrib;
	std::vector <double> distrib_cumul;
	std::vector <double> value;	

	void get_distribution(double );
};



