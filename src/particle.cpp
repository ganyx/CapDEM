 /* This file is part of the Soft_Dynamics library.
    The Soft_Dynamics library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.
 */


//======================================
// Class particle function details
// Class definition is in "particle.h"
//======================================

void  Cparticle::PRINT()
{
	cout<<endl;
	cout<<"Particle id\t "<<id<<endl;
	cout<<"\tWall status\t"<<AM_I_BOUNDARY<<endl;
	cout<<"\tPostion\t"; X.PRINT();
	cout<<"\tVelocity\t"; V.PRINT();
	cout<<"\tAcceleration\t"; A.PRINT();
	cout<<"\tAngular velocity\t"; Ome.PRINT();
	cout<<"\tAngular acceleration\t"; OmeDot.PRINT();
	
	cout<<"\tNumber of contacts\t"<<contact.size()<<endl;
	for(int ic=0;ic<contact.size();ic++)contact[ic]->PRINT();
	cout<<endl;
}





void Cparticle::expand_radius(double dt){R += e*R*Tdot*dt; RHO = m /(4.0/3.0 *PI *R*R*R);}

 
void Cparticle::predictor(double dt,double dt2_on_2)
{
//	Tm = 21;
//	Lm = 1.0*c;	// set to a small fraction of the total latent heat of the particle for a thin shell.
	double dl = 0.0;
	
	if(dt==0) { // recover previous value of RS from L for the first step
		RS = R*R*R - 3.0/4.0 *L /PI/RHO/Lm; 
		RS = pow(RS, 1.0/3.0);	
		}
	X += (V*dt)+ (A*dt2_on_2);
	V += (A*dt);
	Ome+= (OmeDot*dt);
	
	// Melting/cooling condition here: T<Tm; T=Tm, 0<L<Lm; T=Tm, L=Lm;	
	if(T<Tm) 
	{ 
		T += Tdot*dt;
		if(T>Tm) 
		{
			dl= (T-Tm)*c*m;
			L += dl; 
			T=Tm;
	}}	
	else
	{		 
		dl = Tdot*dt*c*m;
		if(L + dl<0.0) {T += (L+dl)/(c*m); L=0.0;}
		else L += dl;
	}
	
	if(L >= Lm*m) L=Lm*m;
	
	if(L>0.0) T = Tm + 0.1* L/(Lm*m); // increase melted temperature, max dT = 0.1
	
	RS_old = RS;
	RS = R;
	dRS = RS - RS_old;
	//update RS, dRS
	if(L>0.0)
	{
		RS = R*R*R - 3.0/4.0 *L /PI/RHO/Lm;
		RS = pow(RS, 1.0/3.0);
		if(RS<0.25*R && dRS<=0.0) {
			RS=0.25*R;dRS=0.0;}
		dRS = RS - RS_old;
	}
	
	// Cap mod
	RS=R;
	dRS=0.0;
}

void Cparticle::corrector(double dt_on_2,Ccell &cell)
{
	if(AM_I_BOUNDARY==-1 || AM_I_BOUNDARY==-2 )
		{
		Ome*=0;
		A.x[0]=	-cell.Ashear; A.x[1]= -cell.Adilat; 
		V.x[0]= -cell.Vshear; V.x[1]= -cell.Vdilat; 
		X.x[1]= -cell.L.x[1]/2;
		return;	
		}
	if(AM_I_BOUNDARY==+1 || AM_I_BOUNDARY==+2)
		{
		Ome*=0;	
		A.x[0]=	+cell.Ashear; A.x[1]= +cell.Adilat; 
		V.x[0]= +cell.Vshear; V.x[1]= +cell.Vdilat; 
		X.x[1]= cell.L.x[1]/2;
		return;
		}

//flowing particle
	Cvector Ap,OmeDotp;
	Ap=A;
	OmeDotp=OmeDot;

// Update total mass inside the cell, grain mass + water mass.
	double msum = m + water_volume*LIQUID_DENSITY;
	A = Fsum/msum;
	
	V += (A-Ap)*dt_on_2;
	OmeDot= Gsum/J;
	Ome += (OmeDot-OmeDotp)*dt_on_2;
	
	cell.rigid_velocity += V *m;
}

void Cparticle::set_me_in_main_cell(Ccell &cell)
{
	//periodic along x and z, true for avery kind of cell.boundary (WALL_INCLINED,WALL_SHEAR or PERIODIC_SHEAR
	
	cell.rescale(X.x[0],cell.L.x[0]);//check if out fron the right/left sides
	cell.rescale(X.x[2],cell.L.x[2]);//check if out fron the front/back sides	
		
	if(cell.boundary=="PERIODIC_SHEAR")	//specificity of PERIODIC_SHEAR
	{
		if (X.x[1]>=cell.L.x[1]/2.)		//out fron the top side
			{
				X.x[1]-=cell.L.x[1];
				X.x[0]-=cell.Xshift;	//shift particle along the right axis		
				cell.rescale(X.x[0],cell.L.x[0]);//check if out fron the right/left sides
				
				V.x[0]-=cell.Vshear;	//increment the velocity  
				V.x[1]-=cell.Vdilat;
				T-=cell.DeltaT;			//increment the temperature
			}	
		
		else if (X.x[1]<-cell.L.x[1]/2.)//out fron the bottom side
			{
				X.x[1]+=cell.L.x[1];	//out fron the top side
				X.x[0]+=cell.Xshift;	//shift particle along the right axis		
				cell.rescale(X.x[0],cell.L.x[0]);//check if out fron the right/left sides
			
				V.x[0]+=cell.Vshear;	//increment the velocity  
				V.x[1]+=cell.Vdilat;
				T+=cell.DeltaT;			//increment the temperature
			}	
	}
}

ofstream & operator<<(ofstream &file,Cparticle p)
{
	file<<p.X<<p.V<<p.Ome; //vector
	file<<p.R<<"\t"<<p.m<<"\t"<<p.water_pressure<<"\t"<<p.saturation<<"\t"<<p.water_volume;//scalar	
//	file<<"\t"<<p.v_face_areas[0]<<"\t"<<p.v_face_areas[1]<<"\t"<<p.v_face_areas[2]<<"\t"<<p.num_neighbors;
	file<<"\t"<<p.positive_pressure;
	file<<endl;	//new line
 	return file;
}

ifstream & operator>>(ifstream &file,Cparticle &p)
{
	file>>p.X>>p.V>>p.Ome; 	//vector
	file>>p.R>>p.m>>p.water_pressure>>p.saturation>>p.water_volume;	//scalar
 	return file;
}

Cparticle::Cparticle(void)
{
	AM_I_BOUNDARY=0;
	R=0;              
	m=0;             
	J=0;
	E=0;
	k=0;
	c=0;
	e=0;
	T=0;
	Tdot=0;			/**<Particle temperature rate.*/
	phi=0; 			/**<Particle total heat rate.*/
	phi_ext=0;			/**<Particle external heat rate.*/
	production=0;	
	L = 0;
	water_volume = 0;
	water_pressure = 0;
	RS=0;
}


void Cparticle::remove_from_box()
{
	for(int k=0;k<my_box.size();k++)
		for(int i=0;i<my_box[k]->part.size();i++) 
			if(my_box[k]->part[i]->id ==id ){
//				my_box[k]->part.removeAt(i);
				if(i < my_box[k]->part.size()-1) 
				{
				my_box[k]->part[i]=my_box[k]->part[my_box[k]->part.size()-1]; 
				}
				my_box[k]->part.pop_back();
				break;
					}
	my_box.clear();
}



void Cparticle::set_in_box(Cmesh &mesh, Ccell &cell) 
{
	my_box.clear();
	if(AM_I_BOUNDARY==+2 || AM_I_BOUNDARY==-2 ) return;//we don't set the planes in the mesh
	
	int x,y,z;
	
	Cvector Xr = X+(cell.L/2.);//rescale in positive value only
	
//	x = (int) (Xr.x[0]/mesh.step.x[0]);
//	y = (int) (Xr.x[1]/mesh.step.x[1]);
//	z = (int) (Xr.x[2]/mesh.step.x[2]);
	x = (int) floor(Xr.x[0]/mesh.step.x[0]);
	y = (int) floor(Xr.x[1]/mesh.step.x[1]);
	z = (int) floor(Xr.x[2]/mesh.step.x[2]);
		
	if(x<0 || x>=mesh.N[0] || y<0 || y>=mesh.N[1]-1|| z<0 || z>=mesh.N[2])
	{
		PRINT();
		STOP("particle.cpp","Cparticle::set_in_box(Cmesh &mesh, int ip, Ccell &cell)","Particle out of the box");
	}
	
	mesh.box[x][y][z].part.push_back(this);						//this is a pointer to "this particle"
	my_box.push_back(&mesh.box[x][y][z]);
	
	
	if(cell.boundary=="PERIODIC_SHEAR" && y==0)//Ouch, this one is tricky: duplicate the bottom layer on the top, with  shift
		{				
		int a,b,c;//new coordinate on the top layer
		double new_x =  (Xr.x[0]+cell.Xshift);
		if(new_x < 0) new_x += cell.L.x[0];
		if(new_x >= cell.L.x[0]) new_x -= cell.L.x[0];
//		a = (int) ( new_x/mesh.step.x[0] );//a is x + shifted
		a = (int) floor(new_x/mesh.step.x[0]);//a is x + shifted
		if(a<0)a+=mesh.N[0];
		else if(a>=mesh.N[0])a-=mesh.N[0];
		
		b = mesh.N[1]-1;//b is the top layer
		c=z;
				
		mesh.box[a][b][c].part.push_back(this);
		my_box.push_back(&mesh.box[a][b][c]);
		}
}


void Cparticle::get_neighbour(Ccell &cell) 
{
	neighbour.clear();
	for(int b =0;b< my_box.size();b++)//	foreach(Cbox *box,my_box)
	for(int cb = 0;cb< my_box[b]->contact_box.size();cb++ )	//for each contacting box		
		for(int B = 0; B< my_box[b]->contact_box[cb]->part.size();B++)	//for each particle of the contacting box
		{
		Cparticle *pB= my_box[b]->contact_box[cb]->part[B]; 
		if(AM_I_BOUNDARY==0 || pB->AM_I_BOUNDARY==0 || (AM_I_BOUNDARY!=pB->AM_I_BOUNDARY)) //don't care of two stuck particle of the same plan
		if(id>pB->id) neighbour.push_back(pB);					//Search within the same box: don't count twice the contact by accepting neigbours only with A.id>B.id
		}
		
	if(cell.boundary=="PERIODIC_SHEAR") return;
	
	if (my_box[0]->am_I_bottom) neighbour.push_back(cell.plan_bottom); 
	if(cell.boundary=="WALL_INCLINED") return;	//then it is WALL_PERIODIC, with two walls
	if (my_box[0]->am_I_top) neighbour.push_back(cell.plan_top); 
}





