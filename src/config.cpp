#include "cdinit.cpp"
bool CAP_FORCE=true;
		
void Cconfig::iterate(double time_step)
{ 	
	dt=time_step;
	dt_on_2 = dt/2.;
	dt2_on_2=dt*dt/2.;
	
	predictor();			//motion integration
	
	update_contact();		//find contact and get the force and torque
	sum_force(); 			//sum the force moment of each contact on particle   
	if(simule_thermal_conduction) sum_heat(); 					//Heat transfer 
	
	if(LIQUID_TRANSFER) liquid_transfer();
	
	// water input, controlling water volume
	if(LIQUID_TRANSFER){

	if(dt==0) flag_wetting = true;

	for(int ip=0; ip< P.size(); ip++) 
	{
		// initial
		if(dt==0) {P[ip].water_volume = parameter.INITIAL_SATURATION * P[ip].void_volume;
			P[ip].water_volume_old = P[ip].water_volume; }

		
		if(P[ip].void_volume <= 1e-3 * P[ip].grain_volume) 
			P[ip].void_volume = 1e-3 * P[ip].grain_volume; // avoid negative Void volume
			
		P[ip].saturation = P[ip].water_volume/P[ip].void_volume; // update the saturation of each cell
		if(P[ip].saturation > MAX_SATURATION) P[ip].saturation = MAX_SATURATION; //exp
	}}
	
	if(!LIQUID_TRANSFER){ // for pre-packing stage
	for(int ip=0; ip< P.size(); ip++) 
	{
		P[ip].saturation = 0.1;
		P[ip].water_volume = 0.1 * P[ip].void_volume; }}

	cell.rigid_velocity *= 0.0;
	corrector();			//acceleration of particles according to the sum of force/moment they experience
	cell.rigid_velocity /= parameter.total_mass;
// Velocity offset by rigid motion
//#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size(); ip++) 
		P[ip].V -= cell.rigid_velocity;

	if(dt==0) {
		for(int ip=0; ip< P.size(); ip++) P[ip].V *= 0.0;}
			
	for(int ip=0; ip< P.size(); ip++) {
		P[ip].V *= (1.0 - GLOBAL_DAMPING*dt); // Global damping
		P[ip].Ome *= (1.0 - GLOBAL_DAMPING*dt);
		}
		
	// calculation of global and local water pressure.
	for(int ip=0; ip< P.size(); ip++) {
		P[ip].water_pressure = 0.0;}
		
	for(int ic=0;ic<C.size();ic++) if(C[ic].fcap >0) {
		P[C[ic].A].water_pressure -= C[ic].fcap * C[ic].dx * P[C[ic].A].R /(P[C[ic].A].R+P[C[ic].B].R);
		P[C[ic].B].water_pressure -= C[ic].fcap * C[ic].dx * P[C[ic].B].R /(P[C[ic].A].R+P[C[ic].B].R);
	}
	for(int ip=0; ip< P.size(); ip++) {
// Modified effective stress term with the degree of saturation (micro-scale).
		if(P[ip].voronoi_volume>1.0e-10) P[ip].water_pressure /= (3.0 *P[ip].voronoi_volume);
		P[ip].positive_pressure = 0.0;
		if(P[ip].saturation > MAX_SATURATION_AIR && P[ip].saturation < 1.0) {
			P[ip].positive_pressure = 
				AIR_K *(P[ip].saturation - MAX_SATURATION_AIR)/(1.0 - P[ip].saturation); // positive pressure
			P[ip].water_pressure += P[ip].positive_pressure*1.0; // air compression the whole cell experiencing the pressure
			}
		if(P[ip].saturation>=1.0e-10) P[ip].water_pressure /= P[ip].saturation;
	}
		
	//Overall saturation and pressure
	saturation = 0.0;
	cap_pressure = 0.0;
	double void_volume=0.0;
	double water_volume=0.0;
	double total_volume_mid = 0.0;
	
	cap_pressure_mid = 0.0;
	double void_volume_mid=0.0;
	double water_volume_mid=0.0;
	
	for(int ip=0; ip< P.size(); ip++){
		void_volume += P[ip].void_volume;
		water_volume += P[ip].water_volume;
// Modified effective stress term with the degree of saturation (macro-scale).
//		if(P[ip].saturation <= 1.0) cap_pressure -= P[ip].water_pressure * P[ip].water_volume; // modified cap pressure
//			else cap_pressure -= P[ip].water_pressure * P[ip].void_volume;
		cap_pressure -= P[ip].water_pressure * P[ip].saturation *P[ip].voronoi_volume;
		if(P[ip].X.x[1] <= 5.0 && P[ip].X.x[1] >= -5.0){
			void_volume_mid += P[ip].void_volume;
			water_volume_mid += P[ip].water_volume;
			cap_pressure_mid -= P[ip].water_pressure * P[ip].saturation *P[ip].voronoi_volume;
			total_volume_mid += P[ip].voronoi_volume;
		}
	}
	saturation = water_volume / void_volume;
	if(saturation < 1.0e-10) saturation = 1.e-10;
	double total_volume = cell.L.x[0]* cell.L.x[1]* cell.L.x[2];
	cap_pressure /= saturation * total_volume;
	water_content = water_volume /total_volume;
	
	saturation_mid = water_volume_mid / void_volume_mid;
	if(saturation_mid<1.0e-10) saturation_mid = 1.e-10;
	if(total_volume_mid > 0.0){
		cap_pressure_mid /= saturation_mid * total_volume_mid;
		water_content_mid = water_volume_mid /total_volume_mid;
	}
	
	if(saturation >= MAX_SCAN && t>1.0 && flag_wetting)	flag_wetting=false;
	if(saturation <= MIN_SCAN && t>1.0 && !flag_wetting) flag_wetting=true;

}

void Cconfig::predictor()
{ 
	cell.predictor(dt,dt2_on_2);//move the cell

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int  ip=0; ip< P.size();ip++) //move the particles
	{ 	
		P[ip].predictor(dt,dt2_on_2);
		P[ip].set_me_in_main_cell(cell);//set the particle back in the cell if it has got out  
	}
	
	
	if(simule_thermal_expansion){
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI		
		for(int  ip=0; ip< P.size();ip++)  P[ip].expand_radius(dt);}
}    

void Cconfig::corrector()
{ 
cell.corrector(dt_on_2);

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
for(int ip=0; ip< P.size();ip++)	
		P[ip].corrector(dt_on_2,cell); 
	 
}

void  Cconfig::sum_heat()
{
	
	double sum_T=0;
	for(int ip=0; ip< P.size();ip++) sum_T += P[ip].T*P[ip].m;
	parameter.average_temperature = sum_T / parameter.total_mass;
	
	PN=0.0; PS=0.0; PT=0.0; PR=0.0;
	for(int ic=0;ic<C.size();ic++)
	{
		PN += C[ic].production_normal;
		PS += C[ic].production_slide;
		PT += C[ic].production_twist;
		PR += C[ic].production_rolling;
	}
	
	heat_in = cell.shear_stress_in * cell.shear_rate + cell.normal_stress_in * cell.dilat_rate;
	heat_in *= cell.L.x[0]* cell.L.x[1]* cell.L.x[2];
	heat_out = 0.0;
	
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++)
	{ 
		P[ip].phi_ext=0;
		P[ip].production=0;
		P[ip].phi =  P[ip].phi_ext; //initialisation for each particles; 0 by default
	}

	for(int ic=0;ic<C.size();ic++) //sum over each interaction
	{
		P[C[ic].A].phi+= C[ic].phi ;
		P[C[ic].B].phi-= C[ic].phi ;
		
		if(C[ic].Flag_Boundary){
			double dPhi = -2.*C[ic].conductivity*C[ic].a *(parameter.average_temperature-20.0)*fabs(C[ic].dX.x[1])/cell.L.x[1]; //exp
			P[C[ic].A].phi += dPhi;
			P[C[ic].B].phi += dPhi;
			heat_out += 2.0*dPhi;
		}

		if(simule_thermal_production)
		{
			P[C[ic].A].production+= C[ic].production/2.; 
			P[C[ic].B].production+= C[ic].production/2.; 		
		}
	}
	
#pragma omp parallel for num_threads(NTHREADS)
	for(int ip=0; ip< P.size();ip++) P[ip].Tdot =  (P[ip].phi +P[ip].production) /(P[ip].c*P[ip].m);
}


void Cconfig::sum_force()
{
 	Cmatrix stress;

	Cvector g;//vector gravity
	if(cell.boundary=="WALL_INCLINED"){
	g.x[0]=  cell.gravity*sin(PI/180 * cell.slope );
	g.x[1]= -cell.gravity*cos(PI/180 * cell.slope );
	g.x[2]=0;
	}//	else g is null
	
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++) 
	{
		P[ip].Fsum = (g*P[ip].m);
		P[ip].Gsum*=0;
	}
	
	if(LIQUID_TRANSFER) drag_force();
	
// YG, no MPI
	for(int ic=0;ic<C.size();ic++)
	{
		Cvector dx;
		P[C[ic].A].Fsum += C[ic].F;
		P[C[ic].B].Fsum -= C[ic].F;
		P[C[ic].A].Gsum +=  (C[ic].RA^C[ic].F) +C[ic].G;
		P[C[ic].B].Gsum -=  (C[ic].RB^C[ic].F) +C[ic].G;
		
// Additional force due to positive water pressure, CHECK !!
		Cvector FWATER_A = C[ic].nA * (P[C[ic].B].positive_pressure)*C[ic].voronoi_area;
		Cvector FWATER_B = C[ic].nA * (P[C[ic].A].positive_pressure)*C[ic].voronoi_area;
		C[ic].fwater = (P[C[ic].A].positive_pressure + P[C[ic].B].positive_pressure)*C[ic].voronoi_area/2.0;
		P[C[ic].A].Fsum += FWATER_A;
		P[C[ic].B].Fsum -= FWATER_B;
		
		dx = C[ic].dX;

		stress+=   (C[ic].F| dx);
	}
	
	double positive_pressure =0.0;
	for(int ip=0; ip<P.size();ip++){
// positive pressure, contribute to the full domain, using voronoi_volume !
		positive_pressure += P[ip].positive_pressure * P[ip].voronoi_volume;
//		 if(P[ip].saturation <= 1.0) positive_pressure += P[ip].positive_pressure * P[ip].water_volume;
//		 else positive_pressure += P[ip].positive_pressure * P[ip].void_volume;
	}
	
	for(int ii=0;ii<3;ii++)
		stress.x[ii][ii] -= positive_pressure;

	if(PSEUDO_2D)	stress/=(cell.L.x[0]*cell.L.x[1]);  
	else	stress/=(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);  
	
	//cell.stress = stress.symetric();  
	cell.stress =stress;
	
	
	if(cell.boundary!="WALL_INCLINED") return;
	cell.normal_stress_in=0;
	cell.shear_stress_in=0;
	for(int ic=0;ic<C.size();ic++)
		{
			if(P[C[ic].A].AM_I_BOUNDARY==-1||P[C[ic].A].AM_I_BOUNDARY==-2){ cell.normal_stress_in+=C[ic].F.x[1]; cell.shear_stress_in+=C[ic].F.x[0];}
			if(P[C[ic].B].AM_I_BOUNDARY==-1||P[C[ic].B].AM_I_BOUNDARY==-2){ cell.normal_stress_in-=C[ic].F.x[1]; cell.shear_stress_in-=C[ic].F.x[0];}
		}
	cell.normal_stress_in/=(cell.L.x[0]*	cell.L.x[2]);cell.shear_stress_in/=(cell.L.x[0]*	cell.L.x[2]);
	cell.normal_stress_in*=-1;//get it positive in compression
}

void Cconfig::update_contact()
{
	renew_contact_list();

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI 
	for(int ic=0;ic<C.size();ic++)
	{
		C[ic].EVALE_Geo();
		C[ic].relative_velocity();
		C[ic].increment_force(dt);
	}
	
	if(simule_thermal_conduction){
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI	
	for(int ic=0;ic<C.size();ic++)	C[ic].EVALE_heat_flow();
	}
}


void Cconfig::renew_contact_list()//Find contact via the mesh method
{

if(Voronoi_Update){

	voronoicell_neighbor vcell;
	// Voronoi tesselation
	double step_target = MESH_SIZE*parameter.Dmax;
	int nx = (int) floor(cell.L.x[0]/step_target);//get the number of box
	int ny = (int) floor(cell.L.x[1]/step_target);//get the number of box
	int nz = (int) floor(cell.L.x[2]/step_target);//get the number of box
	if(nx < 1) nx=1;
	if(ny < 1) ny=1;
	if(nz < 1) nz=1;

	container_poly  con(-cell.L.x[0]/2.0,cell.L.x[0]/2.0,
	-cell.L.x[1]/2.0,cell.L.x[1]/2.0,
	-cell.L.x[2]/2.0,cell.L.x[2]/2.0,
	nx,ny,nz,true,true,true,8);
	
	for( int ip=0;ip<P.size();ip++){
		con.put(ip,P[ip].X.x[0],P[ip].X.x[1],P[ip].X.x[2],P[ip].R);
		}

	for(int ijk=0; ijk < nx*ny*nz; ijk++){
		for(int q=0;q<con.co[ijk];q++){
			con.compute_cell(vcell,ijk,q);
			int ivp=con.id[ijk][q];
			vcell.neighbors(P[ivp].v_neighbors);
			vcell.face_areas(P[ivp].v_face_areas);
			P[ivp].num_neighbors = P[ivp].v_neighbors.size();
			P[ivp].voronoi_volume = vcell.volume();
		}}
	}
	
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI	
	for( int ip=0;ip<P.size();ip++){
		P[ip].grain_volume = 4.0/3.0 * PI * pow(P[ip].R,3.0);
		P[ip].void_volume = P[ip].voronoi_volume - P[ip].grain_volume;
		if(P[ip].void_volume <= 1e-3 * P[ip].grain_volume) 
			P[ip].void_volume = 1e-3 * P[ip].grain_volume; // avoid negative Void volume
		P[ip].saturation = P[ip].water_volume/P[ip].void_volume; // update the saturation of each cell
		if(P[ip].saturation > MAX_SATURATION) P[ip].saturation = MAX_SATURATION; //exp
		}

// YG, no MPI here
//	for( int ip=0;ip<P.size();ip++)	P[ip].set_in_box(mesh,cell); //set particles in box
		
//#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
//	for( int ip=0;ip<P.size();ip++)	P[ip].get_neighbour(cell);//guess what
		
// YG, Be careful for the following section		
	for(int it=0; it<NTHREADS; it++){ CThread[it].clear();} // done below
		
	int in, ip, c, tid; 
	bool exist_contact;
	Ccontact cont;

//#pragma omp parallel for private(in, c, exist_contact,cont, tid) schedule(dynamic) num_threads(NTHREADS)	// YG, MPI testing
	for(ip=0;ip<P.size();ip++)	// Problem with this loop!! Inverse the loop, different results!!
		{
		for(in=0; in<P[ip].num_neighbors;in++)
			{
			if(P[ip].id > P[ip].v_neighbors[in] && P[ip].v_neighbors[in] >=0)
				{
				exist_contact=false;
				for(c=0; c<P[ip].contact.size();c++)
					if(P[ip].contact[c]->pB == &P[P[ip].v_neighbors[in]] || P[ip].contact[c]->pA == &P[P[ip].v_neighbors[in]] )
						{	exist_contact= true; 
							
							P[ip].contact[c]->voronoi_area = P[ip].v_face_areas[in]; // update the voronoi_area, even for existing contacts.
							
							break; }//check if the contact exists
				if(!exist_contact)//if not, create a new contact
				{ 	
					tid = omp_get_thread_num();
					Ccontact cont(&P[ip], &P[P[ip].v_neighbors[in]], &cell, &parameter);	
					if( cont.AM_I_CONTACTING() ) //don't do it if particles are not contacting 	
					{
						cont.voronoi_area = P[ip].v_face_areas[in];	
						CThread[tid].push_back(cont); 
					}											
				}
			}}
		}
	for(int it=0; it<NTHREADS; it++){
		std::vector <Ccontact> Ctemp = C;
		C.insert(C.end(), CThread[it].begin(),CThread[it].end());
//		C+=CThread[it]; 
		CThread[it].clear();
		}
	

// YG, no MPI here	
	for(int ic=0; ic < C.size(); ic++) //delete contact within the list
		if (!C[ic].AM_I_CONTACTING()) 
			{ 
				if(ic< C.size()-1 ) 
				{
				C[ic]=C[C.size()-1]; 
				ic--; 
				}
			C.pop_back(); 
			}

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI   
	for(int ic=0; ic < C.size(); ic++)  C[ic].age+=dt;//get contact older
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI		
	for(int ip=0;ip<P.size();ip++) P[ip].contact.clear();//rebuilt the pointer list
// YG, no MPI here
    for(int ic=0; ic < C.size(); ic++) {
    	C[ic].pA->contact.push_back(&C[ic]);
		}
		
}




 void Cconfig::energy(void)
{ 
 // double Etot,Eela,Eforce,Emoment, Ekin,Erot,Etrans;
   Etot=0;
   Eela=0;Eforce=0;Emoment=0;
   Ekin=0;Erot=0; Etrans=0;
   
 /* foreach(Ccontact cont, C)
  {
//    Eforce += 2./5. *cont.fn*cont.deltaN + 0.5 * cont.ft*cont.ft/(cont.a*parameter.MODULE_T);
    
   cont.gt=cont.Gt.NORM(); 
   cont.gn=con.Gn.NORM();

  //  Emoment+= 0.5/(cont.a*cont.a*cont.a)*(cont.gn*cont.gn/parameter.MODULE_N +  cont.gt*cont.gt/parameter.MODULE_T);
  }
   Eela=Eforce+Emoment;
  
  foreach(Cparticle p, P)
  {
    Cvector V;
    V = p.V;
    V.x[0]-=cell.shear_rate*p.X.x[1];
    Etrans+=0.5*p.m*(V*V);
    Erot +=0.5* 2./5.*p.m* p.R*p.R*(p.Ome*p.Ome);
 
  }
  */
  Ekin=Erot+Etrans;
  Etot = Ekin+Eela;
}



void Cconfig::fread(Cin_out where_to_read)
{
	where_to_read.set_file_name();
	where_to_read.check_file();

    P.clear(); C.clear();

	ifstream file;

	file.open(where_to_read.save_file[0].c_str());  // read particles
	while(1<2)
	{
		Cparticle part;
		file>>part;
		if(!file.eof()) P.push_back(part);
		else break;
	}
	file.close();
	
	file.open(where_to_read.save_file[1].c_str());  // read contact
	while(1<2)
	{
		Ccontact cont;
		file>> cont;
		if(!file.eof()) C.push_back(cont);
		else break;
	}
	file.close();


	file.open(where_to_read.save_file[2].c_str());//save parameter
	file>>t;
	file>>parameter; 
	file.close();

	file.open(where_to_read.save_file[3].c_str()); //save cell
	file>>cell; 
	file.close();
	
	
	//bluild the pointers
	for(int ip=0;ip<P.size();ip++)//cell knows which particle is a plan
	{
		if(P[ip].AM_I_BOUNDARY==-2 ) cell.plan_bottom=&P[ip];
		if(P[ip].AM_I_BOUNDARY==2 )  cell.plan_top=&P[ip];
		P[ip].id=ip;
	}

	for(int ic=0;ic<C.size();ic++)
		{ 
			C[ic].pA = &P[C[ic].A]; 
			C[ic].pB = &P[C[ic].B];
			C[ic].cell = &cell;
			C[ic].parameter = &parameter;
			P[C[ic].A].contact.push_back(&C[ic]);
			if(INIT_BOND != 0) C[ic].deltaNB = INIT_BOND;
		} 

	update_particle();
	parameter.dimensionless_number(cell,P);	
	iterate(0.0);//use here to recover all the data that haven't been saved, but which derive from saved data

}

void Cconfig::fprint(Cin_out where_to_save)
{
	ofstream file;
	
	where_to_save.set_file_name();
	where_to_save.check_file();

	file.open(where_to_save.save_file[0].c_str());  // save particles
//	foreach(Cparticle part, P) file<<part;
	for(int iter=0; iter<P.size();iter++){
		Cparticle part = P[iter];
		file<<part;
	}
	file.close();

	file.open(where_to_save.save_file[1].c_str()); // save contact
//	foreach(Ccontact cont, C) file<<cont;
	for(int iter=0; iter<C.size();iter++){
		Ccontact cont = C[iter];
		file<<cont;
	}
	file.close();

	file.open(where_to_save.save_file[2].c_str()); //save parameter
	file<<t<<"\t";
	file<<parameter;
	file.close();
	
	file.open(where_to_save.save_file[3].c_str()); //save cell
	file<<cell; 
	file.close();

	//YG
//	file.open(where_to_save.save_file[4].c_str(), std::ios_base::app); 
//	file<<t<<"\t"<<cell.Xshift<<"\t"<<cell.L.x[1]<<"\t"<<cell.stress.x[1][1]<<"\t"<<cell.stress.x[1][0]<<"\t"
//	<<parameter.average_temperature<<"\t"<<heat_in<<"\t"<<heat_out<<"\t"<<PN<<"\t"<<PS<<"\t"<<PT<<"\t"<<PR<<endl; 
//	file.close();
	
//int who= RUSAGE_SELF;
//struct rusage usage;
//struct rusage *puse=&usage;
//getrusage(who,puse);
	
//	file.open(where_to_save.save_file[5].c_str(), std::ios_base::app); 
//	file<<t<<"\t"<<puse->ru_maxrss<<"\t"<<puse->ru_ixrss<<"\t"<<puse->ru_idrss<<"\t"<<puse->ru_isrss<<"\t"
//	<<puse->ru_minflt<<"\t"<<puse->ru_majflt<<endl;
}


void  Cconfig::Evale_conductivity_tensor() /**< Measure the conductivity tensor */
{
	Cmatrix K,H; //effective conductivity, convectivity

	K*=0.;
/*	foreach(Ccontact *cont, C)    K+= cont->alpha*(2.*cont->a*(cont->dX*cont->dX)) ;
	if(PSEUDO_2D)cell.K =K/(cell.L.x[0]*cell.L.x[1]);
	else cell.K = K/(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);
*/
	H*=0.;
//	foreach(Cparticle part, P) 
	for(int iter; iter<P.size(); iter++)
	{ 
		Cparticle part = P[iter];
		Cvector deltaV;
	deltaV=part.V;
	deltaV.x[0]-= cell.shear_rate*part.X.x[1];
	H+= ( deltaV|part.X*part.c*part.m); //*
	}

	if(PSEUDO_2D) cell.H =H/(cell.L.x[0]*cell.L.x[1]);
	else cell.H = H/(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);

	cell.production=0;
/*	foreach(Ccontact *cont, C) cell.production+=cont->production;
	if(PSEUDO_2D)  cell.production/=(cell.L.x[0]*cell.L.x[1]);
	else cell.production/=(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);
*/
}

void  Cconfig::liquid_transfer(){	
	double vwmin = MIN_SATURATION *1.0;
				
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++) 
	{
		P[ip].sum_vij = 0.0;
		P[ip].water_volume_old = P[ip].water_volume;
	}
	
	for(int ic=0;ic<C.size();ic++){
		if(C[ic].water_volume > 0.0) // if liquid bridge is existing
		{
			double dwater_volume = 0.0;				
			// define mass transfer of liquid phase, depdending on pressure gradient
			double pa_ave=0.0, pb_ave=0.0;
			bool flag_boundary = true;
			Cparticle *ptop, *pbottom;

//			double Gravity = 2.0e-3;
			if(fabs(C[ic].pA->X.x[1] - C[ic].pB->X.x[1])<10.0){
				flag_boundary = false;
				pa_ave = C[ic].pA->water_pressure + GRAVITY * C[ic].pA->X.x[1];
				pb_ave = C[ic].pB->water_pressure + GRAVITY * C[ic].pB->X.x[1];
			}
			else{ // no water tranport cross X.x[1] boundary
				flag_boundary = true;
				if(C[ic].pA->X.x[1] > 0) { ptop = C[ic].pA; pbottom = C[ic].pB;}
				if(C[ic].pA->X.x[1] < 0) { ptop = C[ic].pB; pbottom = C[ic].pA;}
			}
//			pa_ave = C[ic].pA->water_pressure;
//			pb_ave = C[ic].pB->water_pressure;

		
//	Implementation using water volume
//			double area = min(pow(C[ic].water_volume, 2.0/3.0),C[ic].voronoi_area);
			double area0 = C[ic].voronoi_area;
			double area = pow(C[ic].water_volume, 2.0/3.0);
			if(area > area0) area = area0;
			
			if(C[ic].dx > 1.e-10) C[ic].dot_water_volume = - area * parameter.LIQUID_DIFFUSION * 
				(pa_ave - pb_ave)/C[ic].dx;
			else C[ic].dot_water_volume = 0.0;
				
			dwater_volume = dt* C[ic].dot_water_volume;
			
			// Limit the transport increment!
			double max_dwater = MAX_DWATER *dt;
			if(dwater_volume > max_dwater) dwater_volume = max_dwater;
			if(dwater_volume < -1.0* max_dwater) dwater_volume = -1.0* max_dwater;
			
// The extreme of dry case
			if(dwater_volume > 0 && C[ic].pB->water_volume -vwmin < dwater_volume) 
				dwater_volume = C[ic].pB->water_volume -vwmin;
			if(dwater_volume < 0 && C[ic].pA->water_volume -vwmin < - dwater_volume) 
				dwater_volume = -(C[ic].pA->water_volume -vwmin);
				
// The extreme of wet case
			double Max_A_Water = MAX_SATURATION* C[ic].pA->void_volume;
			double Max_B_Water = MAX_SATURATION* C[ic].pB->void_volume;
			
			if(dwater_volume > 0 && C[ic].pA->water_volume >= Max_A_Water - dwater_volume) 
				dwater_volume = Max_A_Water - C[ic].pA->water_volume;
			if(dwater_volume < 0 && C[ic].pB->water_volume >= Max_B_Water + dwater_volume) 
				dwater_volume = -(Max_B_Water - C[ic].pB->water_volume);
				
			// add mass transfer vector for grains for dragging force (later)
			C[ic].pA->water_volume += dwater_volume;
			C[ic].pB->water_volume -= dwater_volume;
			
			C[ic].dwater_volume = dwater_volume;
			
			
		if(flag_boundary){
			double dwater_boundary = dt * area * parameter.LIQUID_DIFFUSION * INLET_OUTLET_PRESSURE;
//			double dwater_boundary = dt * 1.0* parameter.LIQUID_DIFFUSION * INLET_OUTLET_PRESSURE;

// Limit inlet / outlet flow rate
			double max_dwater = MAX_DWATER *dt;
			if(dwater_boundary > max_dwater) dwater_boundary = max_dwater;
			
			if(!flag_wetting){ // drainage from bottom
				pbottom->water_volume -= dwater_boundary;
				if(pbottom->water_volume < vwmin) 
					pbottom->water_volume = vwmin;
				}
			if(flag_wetting){ // wetting from top
				ptop->water_volume += dwater_boundary;
				if(ptop->water_volume > MAX_SATURATION * ptop->void_volume) 
					ptop->water_volume = MAX_SATURATION * ptop->void_volume;
				}
		}
			//if(flag_boundary && !flag_wetting){
				//pbottom->water_volume -= dwater_boundary;
				//if(pbottom->water_volume < vwmin) 
					//pbottom->water_volume = vwmin;
//				}
			
		}
		// water re-distribution over the contacts
		C[ic].water_vij = C[ic].pA->R * C[ic].pB->R *(C[ic].pA->R + C[ic].pB->R);
		C[ic].pA->sum_vij += C[ic].water_vij;
		C[ic].pB->sum_vij += C[ic].water_vij;
	}
/*	
	// Homogenous loading
	double dwater_volume;
	if(flag_wetting) dwater_volume = WETTING_RATE*dt;
	if(!flag_wetting) dwater_volume = -1.0 *WETTING_RATE*dt;
	
	for(int ip=0;ip<P.size();ip++){
		P[ip].water_volume += dwater_volume;
//		if(P[ip].water_volume > MAX_SATURATION * P[ip].void_volume) P[ip].water_volume = MAX_SATURATION * P[ip].void_volume;
		if(P[ip].water_volume < vwmin) P[ip].water_volume = vwmin;
	}
*/
	
	// water distribution inside contacting Voronoi cells
	for(int ic=0;ic<C.size();ic++){
		double dV_pA = C[ic].water_vij *C[ic].pA->water_volume/C[ic].pA->sum_vij;
		double dV_pB = C[ic].water_vij *C[ic].pB->water_volume/C[ic].pB->sum_vij;
		double dV = dV_pA + dV_pB;
		C[ic].water_volume_old = C[ic].water_volume;
		C[ic].water_volume = dV;
	}

}

void Cconfig::drag_force(){
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++) 
	{
		Cvector Fdrag;
		Fdrag *= 0.0;
		P[ip].Fsum += Fdrag;
	}
}




