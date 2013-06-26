

void  Cconfig::update_particle()
{
	int nloop = 100; // convert percentage of compositions
	int nint = nloop*parameter.COMP_FRACTION;
	int i=0;
	
	for(int ip=0;ip<P.size();ip++)
	{
//		if(ip < COMP_FRACTION *P.size()){ // Biotite
		if(i<nint){ // Biotite
				P[ip].Tm = 21.0;
				P[ip].Lm = 5.0;
				P[ip].k = parameter.bulk_conductivity;
		}
		else { // Quartz, Plagioclase, K-feldstar
			P[ip].Tm = 25.0;
			P[ip].Lm = 3.0;
			P[ip].k = 5.0*parameter.bulk_conductivity;
		}
		
		P[ip].c = parameter.specific_heat;
		P[ip].E = parameter.MODULE_N;
//		P[ip].k = parameter.bulk_conductivity;
		P[ip].e = parameter.thermal_expansion;
		P[ip].J = 2./5.*P[ip].m* P[ip].R*P[ip].R; 	
		
		P[ip].Tdot = 0.0;
		
		i++;
		if(i>=nloop) i-=nloop;
	}
}
	
void Cconfig::create_random()
{
	int Npart_flow,Npart_wall_bottom,Npart_wall_top;
	
	cout<<endl<<"Attempt to create a random configuration of grains with no overlaping"<<endl<<endl;	
	//System size	
	get_secure("Enter the size of the system","SIZE_CELL",cell.L);
	//Kind of boundary	
	get_secure("Enter the kind of boundary you want along the y direction", "PERIODIC_SHEAR","WALL_INCLINED","WALL_SHEAR",cell.boundary);

	if(cell.boundary=="PERIODIC_SHEAR")
	{
	Npart_wall_bottom=0;	
	Npart_wall_top=0;	
	}

	if(cell.boundary=="WALL_INCLINED")
	{
		get_secure("Enter the number of grains of the bottom wall, including the plane (it counts as a grain)", "NPART_WALL_BOTTOM",Npart_wall_bottom);
		Npart_wall_top=0;	
		if(Npart_wall_bottom<1){serror="The number of grains on the bottom wall must be >=1";STOP("cdinit.cpp", "create_random()",serror);}
	}
	
	if(cell.boundary=="WALL_SHEAR")
	{
		get_secure("Enter the number of grains of the bottom wall (>=1), including the plan (it counts as one grain)", "NPART_WALL_BOTTOM",Npart_wall_bottom);
		get_secure("Enter the number of grains of the top wall  (>=1), including the plan (it counts as one grain)", "NPART_WALL_TOP",Npart_wall_top);	
		if(Npart_wall_bottom<1 || Npart_wall_top<1){serror="The number of grains on the wall must be >=1";STOP("cdinit.cpp", "create_random()",serror);}
	}
	

	
	get_secure("Enter the type of grain size distribution", "FRACTAL", "UNIFORM","VOLUME",parameter.GSD);
	get_secure("Enter the minimum diameter", "DMIN", parameter.Dmin);
	get_secure("Enter the minimum diameter", "DMAX", parameter.Dmax);
	if(parameter.GSD=="UNIFORM") get_secure("Enter the number of flowing grains (do not include the wall-grains)","NPART_FLOW",Npart_flow);
	if(parameter.GSD=="VOLUME") get_secure("Enter the number of flowing grains (do not include the wall-grains)","NPART_FLOW",Npart_flow);
	if(parameter.GSD=="FRACTAL")
		{
			 get_secure("Enter the fractal dimention of the grain size distribution","FRACTAL_DIM", parameter.fractal_dim);
			 get_secure("Enter the solid fraction wanted","SOLID_FRACTION", cell.solid_fraction);
		 }

	set_random_grain(Npart_flow, Npart_wall_bottom, Npart_wall_top); 


	for(int ip=0;ip<P.size();ip++)//example of how to initiate things, mass is to be initiated here, it won't be change after
	{
		double RHO;
		RHO = 6./PI; //mass 1 for diameter 1, at temperature 20 C (may change with thermal expansion)
		
		if(P[ip].AM_I_BOUNDARY==0)//flowing particle
		{
			P[ip].m = 4./3. * PI * pow(P[ip].R,3)*RHO;
			P[ip].T = 20;
		}
		if(P[ip].AM_I_BOUNDARY==-1 ||P[ip].AM_I_BOUNDARY==-2 )//bottom walls
		{
			P[ip].m = 0;  //don't care  the mass
			P[ip].T = 20;    //can be change as wish
		}
		
		if(P[ip].AM_I_BOUNDARY==1 ||P[ip].AM_I_BOUNDARY== 2 )//bottom walls
		{
			P[ip].m = 0;  //don't care the mass
			P[ip].T = 20;    //can be change as wish
		}
	}
	cout<<"Random configuration created: SUCCESS"<<endl<<endl;	
}

double UNIFORM_RANDOM_DOUBLE(double min, double max)
{
  double data;
  if(min==max) return min;
  if(min>max) {sprintf(error,"FCT UNIFORM_RANDOM_DOUBLE in random.cpp: min>max");ERROR(error);}
  data = min+(max-min)*(rand()/((double)RAND_MAX+1));
  return data;
}

    
void Cconfig::set_wall_grain(int Nb, int Nt)
{
	Cparticle part;
	
	if(Nb==0)return;
	part.R  =0; 
	part.AM_I_BOUNDARY =-2;
	part.id=P.size();
	P.push_back(part);
	cell.plan_bottom=&P[0];
	for(int p=1;p<Nb;p++)//bottom wall
		{
			part.R = UNIFORM_RANDOM_DOUBLE( parameter.Dmin/2., parameter.Dmax/2.); part.AM_I_BOUNDARY =-1;
			part.X.x[0]=UNIFORM_RANDOM_DOUBLE(-cell.L.x[0]/2.,cell.L.x[0]/2.);
	       	part.X.x[1]=-cell.L.x[1]/2.;
	       	if(!PSEUDO_2D) UNIFORM_RANDOM_DOUBLE(-cell.L.x[2]/2.,cell.L.x[2]/2.);// if we are in 2D, this component of position is zero  
        	part.id=P.size();
        	P.push_back(part);
        }
        
    if(Nt==0)return;  
	part.R  =0; 
	part.AM_I_BOUNDARY =2;
	part.id=P.size();
	P.push_back(part);
	cell.plan_top=&P[Nb];
	
	for(int p=1;p<Nt;p++)//top wall
		{
			part.R = UNIFORM_RANDOM_DOUBLE( parameter.Dmin/2., parameter.Dmax/2.); part.AM_I_BOUNDARY =1;
			part.X.x[0]=UNIFORM_RANDOM_DOUBLE(-cell.L.x[0]/2.,cell.L.x[0]/2.);
	       	part.X.x[1]= cell.L.x[1]/2.;
	       	if(!PSEUDO_2D) part.X.x[2]=UNIFORM_RANDOM_DOUBLE(-cell.L.x[2]/2.,cell.L.x[2]/2.);// if we are in 2D, this component of position is zero  
        	part.id=P.size();
        	P.push_back(part);
        }
	cout<<"Wall are done: SUCCESS"<<endl;
}
    
void  Cconfig::set_radius(std::vector <double> &radius, int &Nf)
{
	double R, Vtemp;
	double Vmin = pow(parameter.Dmin/2.0,-3.0);
	double Vmax = pow(parameter.Dmax/2.0,-3.0);

	if(parameter.GSD=="UNIFORM")
		for(int i=0;i<Nf;i++) 
		{
		R=UNIFORM_RANDOM_DOUBLE( parameter.Dmin/2., parameter.Dmax/2.);//radius  
		radius.push_back(R);
		}
	
	if(parameter.GSD=="VOLUME")
		for(int i=0;i<Nf;i++) 
		{
//		Vtemp = UNIFORM_RANDOM_DOUBLE(Vmin, Vmax);//volume
		Vtemp = UNIFORM_RANDOM_DOUBLE(Vmax,Vmin);//volume
		R = pow(Vtemp,-1.0/3.0); 
		radius.push_back(R);
		}
	
	double a = parameter.fractal_dim;//for readability only
	double M_on_rho =  cell.solid_fraction*cell.L.x[0]*cell.L.x[1]*cell.L.x[2];//for readability only
	double C;
	C= M_on_rho*(3.-a)/a*6./PI;
	C/= (pow(parameter.Dmax,3.-a)- pow(parameter.Dmin,3.-a));

	double  dstep = (parameter.Dmax-parameter.Dmin)/100;// 1./(C*a*pow(parameter.Dmax,-1.-a));
	if(parameter.GSD=="FRACTAL")
		for(double d =  parameter.Dmin;d< parameter.Dmax;d+=dstep)
		{
			int Npart =(int) (C*a*pow(d,-1.-a)*dstep);
			for(int p=0;p<Npart;p++) 
			{
				R=UNIFORM_RANDOM_DOUBLE(d-dstep/2.,d+dstep/2.)/2.;
				radius.push_back(R);		
			}
		}
		
//	qSort(radius.begin(), radius.end(), qGreater<double>());// sort the list from larger particle to smaller ones, so that it's would be easier to set them randomly (set bigger first)
	std::sort(radius.begin(), radius.end());
	std::reverse(radius.begin(), radius.end());
	Nf = radius.size();
//	cout<<radius[0]<<radius[Nf-1];
}

void Cconfig::set_random_grain(int Nf, int Nb, int Nt)
{

	srand( (unsigned int) time(NULL)+getpid());//Init of random

	cout<<endl<<"Start the random setting of grains"<<endl;
	
	set_wall_grain(Nb,Nt);//set the grain of the wall and the planes (if there is any wall)
	
	std::vector <double> radius;
	set_radius(radius,Nf);
	cout<<"Number of flowing grains\t"<<radius.size()<<endl;
	
	Cmesh mesh(cell.L, MESH_SIZE*parameter.Dmax, cell); 
	
	for(int p=0;p<radius.size();p++)//set flowing grains
	{
		bool position_occupied;
		
		Cparticle part;
		part.R=radius[p];
		
		part.id=P.size();
		part.AM_I_BOUNDARY=0;
		P.push_back(part);//save the new particle	

		do
			{
			position_occupied=false;	
			
			for(int i=0;i<DIM;i++) P[P.size()-1].X.x[i]=UNIFORM_RANDOM_DOUBLE(-cell.L.x[i]/2.,cell.L.x[i]/2.);	
			if(PSEUDO_2D) P[P.size()-1].X.x[2]=0.; // if we are in 2D, this component of position is zero
			P[P.size()-1].set_in_box(mesh,cell);
			P[P.size()-1].get_neighbour(cell);

			for(int in=0; in<P[P.size()-1].neighbour.size();in++) 
				{
				Ccontact cont(&P[P.size()-1], P[P.size()-1].neighbour[in], &cell, &parameter);		
				if(cont.AM_I_CONTACTING())
					{	
					P[P.size()-1].remove_from_box(); 
					position_occupied=true; 
					break;
					}
		}	
			}while(position_occupied);
		
		if(p%100==0) 
			cout<<"Number of particles set: "<<p<<endl;//write on screen every when 200 more grains found there place 
	}

double frac_sol = 0;
for(int p1 = 0;p1<P.size()-1;p1++) frac_sol+=4./3.*PI*pow(P[p1].R,3.);
frac_sol/=(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);
cout<<"Actual solid fraction:\t"<<frac_sol <<endl;

/*
	for(int p1 = 0;p1<P.size()-1;p1++)
		for(int p2 = p1+1;p2<P.size();p2++) 
		{
			Ccontact cont( &P[p1],&P[p2], &cell, &parameter);		
			if(cont.AM_I_CONTACTING()){cout<<"There is contacts"<<endl;cont.PRINT();}
			
		}
	*/
	cout<<"All the grains have found their position: SUCCESS"<<endl;
}


