 /* This file is part of the Soft_Dynamics library.
    The Soft_Dynamics library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.
 */

//=====================================================================
//============================== CELL PROPERTIES ======================
//=====================================================================

Ccell::Ccell()//initialisation
{

	mass=0.0;  
	normal_stress_ext=0.0;
	shear_stress_ext=0.0;
	normal_stress_in=0.0; 
	shear_stress_in=0.0;
	shear_rate=0.0;
	dilat_rate=0.0;
	cumul_strain=0.0;
	
	Ashear=0.0;
	Ashear_p=0.0;
	Adilat=0.0;
	Vshear=0.0;
	Vdilat=0.0;
	DeltaT=0.0;
	gradT=0.0;
	Xshift=0.0; 
	
	dstress=0.0; dstrain_int=0.0;	

	gravity=0.0;				
	slope=0.0;				

	normal_stress_control=false;
	shear_stress_control=false;
	gradT_control=false;

	energy_kinetic=0.0;
	production=0.0;
	coordination=0.0;
}


void Ccell::predictor(double dt,double dt2_on_2)
{
	if(boundary=="WALL_INCLINED")return;//nothing to do 
	
	//for shear only
	if(normal_stress_control)//if the normal stress is controled
		{
			// PR impl. of normal stress control
			Yshift = Vdilat*dt + Adilat*dt2_on_2;
			double max_shift = DILAT_LIMIT*L.x[1] * dt;
			if(Yshift > max_shift) Yshift = max_shift;
			L.x[1] += Yshift;
			Vdilat +=  Adilat*dt;
		if( Vdilat < -DILAT_LIMIT*L.x[1]) Vdilat= -DILAT_LIMIT*L.x[1]; //limit cell velocity
		else if(Vdilat >  DILAT_LIMIT*L.x[1]) Vdilat=  DILAT_LIMIT*L.x[1]; //limit cell velocity
			dilat_rate=Vdilat/L.x[1];
		}
	else //if the volume is controled
		{
		Vdilat=0;
		dilat_rate=0;
		}
		
	// block for STICK_SLIP control
	if(stick_slip){
//		double Host_Depth = 1e7/2.0;
//		double Host_Shear_Modulus = 1000/2.0; // G=E/2(1+\nv)
//		double MassFactor = 1.0*L.x[1]/100.0;
//		double MassFactor = 1.0/1000.0; // mass surface density = \sigma_n /g
		shear_stress_control = 1; // overwrite the flag for shear_stress control
//		shear_stress_ext += Host_Shear_Modulus *(slip_velocity*dt - Vshear*dt-Ashear*dt2_on_2)/Host_Depth;	// overwrite the value of external shear stress
		shear_stress_ext += WALL_KD *(slip_velocity*dt - Vshear*dt-Ashear*dt2_on_2);
		if(dt == 0) Ashear_p = Ashear;
//		shear_stress_ext += - MassFactor * (Ashear-Ashear_p);	// acceleration term
	}
		
	if(shear_stress_control)//if the shear stress is controled
	{
		Xshift += Vshear*dt+Ashear*dt2_on_2;
		Vshear +=  Ashear*dt;
		
		if( Vshear < -SHEAR_LIMIT*L.x[1]) Vshear = -SHEAR_LIMIT*L.x[1]; //limit cell velocity
		else if(Vshear >  SHEAR_LIMIT*L.x[1]) Vshear =  SHEAR_LIMIT*L.x[1]; //limit cell velocity
		shear_rate=Vshear/L.x[1];
	}
	//if(shear_work_control)
	//{
		//shear_stress_rate = (shear_stress_in - shear_stress_in_old)/dt;
		//shear_stress_in_old = shear_stress_in;
		//
		//double shear_stress = shear_stress_in;
		//if(fabs(shear_stress)<1.e-10) shear_stress = copysign(1.0,shear_stress_in);
		//shear_acc = - shear_stress_rate*shear_rate/shear_stress;
		//shear_rate = shear_work_input/shear_stress + shear_acc*dt;
		//
		//Vshear = shear_rate*L.x[1];
		//Xshift += Vshear*dt;
//	}
	if(!shear_stress_control && !shear_work_control) //if the shear rate is controled
	{
			Vshear=shear_rate*L.x[1];
			Xshift+=Vshear*dt;
	}
// ????? BUG
//	if(Xshift<-L.x[0]/2)Xshift+=L.x[0]/2;	//rescale in the box x size
//	else if(Xshift>=L.x[0]/2)Xshift-=L.x[0]/2;

// ??? BUG check contact.cpp, Xshift~[0,L), while here Xshift~[-L/2,L/2), change the impl. in contact.cpp
	if(Xshift<-L.x[0]/2.0) Xshift += L.x[0];		//rescale in the box x size
	else if(Xshift>=L.x[0]/2.0) Xshift -= L.x[0];

	cumul_strain+= shear_rate*dt;
}

void Ccell::corrector(double dt_on_2) 
{
	if(boundary=="WALL_INCLINED")return;	//nothing to do 
	
	double Adilat_p;
			
	normal_stress_in = stress.x[1][1]; 
	shear_stress_in = stress.x[0][1]; 
	mass= 1.0;
	//normal_stress_int shear_stress_int  
	
/*	// YG code start
	double Tint=1.0, Tder=1.e-5, Estar=1000;
	double dstress_p, dstrain;
	dstress_p = dstress;
	dstress = (-normal_stress_ext-normal_stress_in);
	// YG code end
*/
	
	if(normal_stress_control)//if the normal stress is controled
	{
		// PR impl. of normal stress control
		Adilat_p = Adilat;
		Adilat   = (-normal_stress_ext-normal_stress_in)/mass;
		Vdilat  += (Adilat-Adilat_p)*dt_on_2;
		if( Vdilat < -DILAT_LIMIT*L.x[1]) Vdilat= -DILAT_LIMIT*L.x[1]; //limit cell velocity
		else if(Vdilat >  DILAT_LIMIT*L.x[1]) Vdilat=  DILAT_LIMIT*L.x[1]; //limit cell velocity
		
		// YG PID controller simplified version
/*		dstrain = 1.0/Estar *dstress;				// Proportional part
		dstrain_int += dt*dstress/(Estar*Tint);		// Integral part
		dstrain += dstrain_int;
		dstrain += Tder*(dstress-dstress_p)/(Estar*dt);	// Derivative part
		Vdilat = dstrain/dt *L.x[1];
		if( Vdilat < -0.025*L.x[1]) Vdilat= -0.025*L.x[1]; //limit cell velocity
		else if(Vdilat >  0.025*L.x[1]) Vdilat=  0.025*L.x[1]; //limit cell velocity
		// YG PID controller simplified version, end
*/
	}
	else Adilat=0;


	if(shear_stress_control)//if the shear stress is controled
	{
		Ashear_p =  Ashear;
		double DAMP = WALL_DAMP *Vshear;;
		Ashear 	 =  (shear_stress_ext - shear_stress_in - DAMP)/mass;
		Vshear 	+= (Ashear-Ashear_p)*dt_on_2;
		//exp	
		if( Vshear < -SHEAR_LIMIT*L.x[1]) Vshear = -SHEAR_LIMIT*L.x[1]; //limit cell velocity
		else if(Vshear >  SHEAR_LIMIT*L.x[1]) Vshear =  SHEAR_LIMIT*L.x[1]; //limit cell velocity
	}
	else	Ashear=0;

// exp
	gradT_control=false;
	if(gradT_control) DeltaT = gradT*L.x[1];//if a difference of temperature is imposed
}


void Ccell::rescale(Cvector &X)
{
	rescale(X.x[0],L.x[0]);//check if out fron the right/left sides
	rescale(X.x[2],L.x[2]);//check if out fron the front/back sides
	
	if(boundary=="PERIODIC_SHEAR")	//specificity of PERIODIC_SHEAR
	{
		if (X.x[1]>=L.x[1]/2.)		//out fron the top side
			{
				X.x[1]-=L.x[1];
				X.x[0]-=Xshift;	//shift particle along the right axis		
				rescale(X.x[0],L.x[0]);//re-check if out fron the right/left sides
			}	
		
		else if (X.x[1]<-L.x[1]/2.)//out fron the bottom side
			{
				X.x[1]+=L.x[1];	//out fron the top side
				X.x[0]+=Xshift;	//shift particle along the right axis		
				rescale(X.x[0],L.x[0]);  //re-check if out fron the right/left sides
			}	
	}

}

void Ccell::rescale(double &x, double l)
{
	if (x>=l/2.)				x-=l;//out fron the positive side
	else if (x<-l/2.)			x+=l;//out from the negative side
}


void  Ccell::PRINT(){

cout<<"Cell properties"<<endl;
cout<<"\tSystem size:\t"; L.PRINT();

if(boundary=="WALL_INCLINED") 
{
	cout<<"\tInclinded plane: only one bottom wall, x and z direction are periodic"<<endl;
	cout<<"\tBottom normal and shear stress:\t"<<normal_stress_in<<"\t"<<shear_stress_in<<endl;
	return;
}

if(boundary=="PERIODIC_SHEAR") cout<<"\tPlane shear without wall: x,y and z directions are periodic"<<endl;	
if(boundary=="WALL_SHEAR") cout<<"\tShear with wall: bottom and top walls, x and z directions are periodic"<<endl;

if(normal_stress_control) cout<<"\tNormal stress is controlled:\t"<<normal_stress_ext<<endl;	
else	cout<<"\tNormal stress is not controlled: constant volume"<<endl;
if(shear_stress_control)cout<<"\tShear stress is controlled:\t"<<shear_stress_ext<<endl;	
else cout<<"\tShear rate is controlled:\t"<<shear_rate<<endl;	
cout<<"\tNormal/Shear stress inside:\t"<<normal_stress_in<<"/"<<shear_stress_in<<endl;

cout<<"\tCumulative shear strain:\t"<<cumul_strain<<endl;
cout<<"\tShift of the top/bottom cell:\t"<<Xshift<<endl;
cout<<"\tShear/Dilatation velocity:\t"<<Vshear<<"/"<<Vdilat<<endl;
}


ofstream & operator<<(ofstream &file,Ccell c)
{
	
	file<<c.L;
	file<<c.boundary<<"\t"<<c.mass<<"\t";
	file<<c.normal_stress_ext<<"\t"<<c.shear_stress_ext<<"\t"<<c.normal_stress_in<<"\t"<<c.shear_stress_in<<"\t";
	file<<c.shear_rate<<"\t"<<c.dilat_rate<<"\t"<<c.cumul_strain<<"\t";
	file<<c.Ashear<<"\t"<<c.Adilat<<"\t"<<c.Vshear<<"\t"<<c.Vdilat<<"\t";
	file<<c.normal_stress_control<<"\t"<<c.shear_stress_control<<"\t";
	file<<c.DeltaT<<"\t"<<c.gradT<<"\t"<<c.Xshift<<"\t"<<c.gradT_control<<"\t";
	file<<c.solid_fraction<<"\t";
	file<<endl;
	return file;
}

ifstream & operator>>(ifstream &file,Ccell &c)
{
	file>>c.L;
	file>>c.boundary>>c.mass;
	file>>c.normal_stress_ext>>c.shear_stress_ext>>c.normal_stress_in>>c.shear_stress_in;
	file>>c.shear_rate>>c.dilat_rate>>c.cumul_strain;
	file>>c.Ashear>>c.Adilat>>c.Vshear>>c.Vdilat;
	file>>c.normal_stress_control>>c.shear_stress_control;
	file>>c.DeltaT>>c.gradT>>c.Xshift>>c.gradT_control;
	file>>c.solid_fraction;
	return file;
}


