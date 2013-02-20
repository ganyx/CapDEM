//#include "mesh.h"
//
Cmesh::Cmesh(Cvector L, double step_target, Ccell cell) 
{

//	box.clear();	
	
	for(int k=0;k<DIM;k++)	
		{
//			N[k] = (int)(L.x[k]/step_target);//get the number of box 
			N[k] = (int) floor(L.x[k]/step_target);//get the number of box
			step.x[k] = L.x[k]/ ((double)N[k]);//get the actual step with this number of box
		}
	N[1]++;	//+1 means that one box is added on the positive side of the y direction for periodic boundary conditions
			//bottom paerticle will be duplicate up there only for periodic shear 
		
//	if(N[0]<3 || N[1]<4 || N[2]<3  ) //check if there is enough box
//		STOP("mesh.cpp","get_neigbour(QList <Cparticle> P, Ccell cell)"," There is not enough boxes in some direction (N<4)");	
		
//	twoDarray plan;
//	oneDarray line; 
//	Cbox box_single;	
	
//	for(int z=0;z<N[2];z++) line.push_back(box_single);	
//	for(int y=0;y<N[1];y++) plan.push_back(line);
//	for(int x=0;x<N[0];x++) box.push_back(plan);	//the order of these three lines ensures access to a box with x y z indexes by box[x][y][z] 	
	
	// SEG FAULT: destructor of box[][][]?
	for(int z=0;z<N[2];z++)	for(int y=0;y<N[1];y++)	for(int x=0;x<N[0];x++){
		box[x][y][z].part.clear();
		box[x][y][z].contact_box.clear();
	}
	
	for(int x=0;x<N[0];x++)for(int y=0;y<N[1];y++)for(int z=0;z<N[2];z++)//tell the ones close to the bottom and top
	{
		if(y==0) box[x][y][z].am_I_bottom=true;	else box[x][y][z].am_I_bottom=false;
		if(y==N[1]-2)box[x][y][z].am_I_top=true; else  box[x][y][z].am_I_top=false;
	}
		
	
	for(int x=0; x< N[0];x++) //get contacting box foreach box
		for(int y=0; y<N[1]-1;y++) 
			for(int z=0; z<N[2];z++) 
				for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++)for(int k=-1;k<=1;k++)
				{
					int a = x+i;
					int b = y+j;
					int c = z+k;
					
					if(a<0)			a =N[0]-1;
					else if(a >= N[0]) 	a =0;
					if(c<0)			c =N[2]-1;
					else if(c >= N[2])	c=0;
					
					if(b>=0) box[x][y][z].contact_box.push_back(& box[a][b][c]);	//if not, do nothing, because the bottom plan has been duplicated on the top
				}
		
		if(cell.boundary!="PERIODIC_SHEAR") return;
		for(int x=0; x< N[0];x++)	
			for(int z=0; z<N[2];z++) 
				for(int i=-1;i<=1;i++)	for(int k=-1;k<=1;k++)
				{
					int y = N[1]-1;
					int a = x+i;
					int b = y-1;
					int c = z+k;
					
					if(a<0)			a =N[0]-1;
					else if(a >= N[0]) 	a =0;
					if(c<0)			c =N[2]-1;
					else if(c >= N[2])	c=0;
			 		
			 		box[x][y][z].contact_box.push_back(& box[a][b][c]);
				}
}

//Cmesh::~Cmesh(){}
/*
void Cmesh::get_neigbour(Ccell &cell)
{

	for(int x=0; x< N[0];x++) for(int y=0; y< N[1];y++) for(int z=0; z<N[2];z++) //get neigbours from each pair of box
		for(int A = 0; A<box[x][y][z].part.size();A++)//foreach(Cparticle *pA,  box[x][y][z].part)					//for each particle of the box
			foreach(Cbox *boxB, box[x][y][z].contact_box)	//for each contacting box
				for(int B = 0; B<boxB->part.size();B++)	//for each particle of the contacting box
				{
				Cparticle *pA= box[x][y][z].part[A];
				Cparticle *pB= boxB->part[B];
				if(pA->AM_I_BOUNDARY==0 || pB->AM_I_BOUNDARY==0 || (pA->AM_I_BOUNDARY!=pB->AM_I_BOUNDARY)) //don't care of two stuck particle of the same plan
				if(boxB==&box[x][y][z]){  if(pA<pB) pA->neighbour.push_back(pB);}//Search within the same box: don't count twice the contact
				else {pA->neighbour.push_back(pB); }////Search in two different box
				}
	if(cell.boundary=="PERIODIC_SHEAR") return;
	
	//if there wall(s), get neigbourhoud 
	for(int x=0; x< N[0]-1;x++) for(int z=0; z<N[2]-1;z++) 
		{	
		//for(int A = 0; A<box[x][0][z].part.size();A++) cell.plan_bottom->neighbour.push_back(box[x][0][z].part[A]->neighbour.push_back());
		//for(int A = 0; A<box[x][1][z].part.size();A++) cell.plan_bottom->neighbour.push_back(box[x][1][z].part[A]->neighbour.push_back());
		foreach(Cparticle *pA,  box[x][0][z].part) pA->neighbour.push_back(cell.plan_bottom); 	//add plan
		foreach(Cparticle *pA,  box[x][1][z].part) pA->neighbour.push_back(cell.plan_bottom); 	//add plan
		}
		
	if(cell.boundary=="WALL_INCLINED") return;
	//then it is WALL_PERIODIC, with two walls
	for(int x=0; x< N[0]-1;x++) for(int z=0; z<N[2]-1;z++) 
		{
		for(int A = 0; A<box[x][N[1]-1][z].part.size();A++) box[x][N[1]-1][z].part[A]->neighbour.push_back(cell.plan_top);
		for(int A = 0; A<box[x][N[1]-2][z].part.size();A++) box[x][N[1]-2][z].part[A]->neighbour.push_back(cell.plan_top);
		}
		
}*/
