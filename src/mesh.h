#ifndef MESH_H
#define MESH_H
//

//

typedef std::vector <Cbox> oneDarray;
typedef std::vector <oneDarray> twoDarray;
typedef std::vector <twoDarray> threeDarray;

class Cbox 
{
	public:

	bool am_I_bottom;
	bool am_I_top;
	std::vector <Cparticle *> part; 
	std::vector <Cbox *> contact_box;
};

class Cmesh  
{
public:
	Cvector step;  							/**<Size of each box*/
//	threeDarray box; 						/**<3D array of box*/ 
	class Cbox box[100][100][30];
	int N[DIM];								/**<number of box in each direction, for convinience only*/
	Cmesh(Cvector, double, Ccell);
//	void get_neigbour(Ccell&);
};


#endif
