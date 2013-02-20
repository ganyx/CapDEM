//
ofstream & operator<<(ofstream &file,Cslice s)
{
	file<<s.Y<<"\t";
	file<<s.frac_solid<<"\t";
	file<<s.V;
	file<<s.T<<"\t";
	file<<endl;	//new line
 	return file;
}

void Cslice::average(void)
{
	double volume_total_in_slice=0;
	for(int ip=0;ip<P.size();ip++)
		{
			volume_total_in_slice+=volume_in_slice[ip];
			T+= P[ip]->T* volume_in_slice[ip];
			V+= (P[ip]->V* volume_in_slice[ip]);
		}
	frac_solid = volume_total_in_slice/vol_slice;
	T/=volume_total_in_slice;
	V/=volume_total_in_slice;
}




Cprofile::Cprofile(double step_target,Cconfig config) 
{
	double Nslice;
	
	//initiate the slices
	Nslice = (int)(config.cell.L.x[1]/step_target);
	step = config.cell.L.x[1]/Nslice; //trick to get same size of slice, close to step target
	
	Cslice temp;
	if(PSEUDO_2D) temp.vol_slice=step*config.cell.L.x[0]; //in that case, we assume that the z size is 1
	else  temp.vol_slice=step*config.cell.L.x[0]*config.cell.L.x[2];
	for(int i=0;i<Nslice;i++) 
	{
		temp.Y= i*step + step/2 - config.cell.L.x[1]/2; 
		slice.push_back(temp); //bluild the list of slice
	}
	//end initiate the slices

	//dispach_particle_in_slice();
	for(int ip=0; ip <config.P.size();ip++ )
		if(config.P[ip].AM_I_BOUNDARY==0)//only flowing particles
		{
			double y = config.P[ip].X.x[1] + config.cell.L.x[1]/2.;
			double R = config.P[ip].R;
	
			int imin = (int)((y-R)/step);
			int imax = (int)((y+R)/step);
		
			double a = -R;
			double b = (imin+1)*step-y;
			for(int islice=imin;islice<=imax;islice++)
				{
					int islice_periodic=islice;
					if(islice<0)					islice_periodic = islice+slice.size();
					else if(islice>=slice.size()) 	islice_periodic = islice-slice.size();
					
					if(islice_periodic<0|| islice_periodic>=slice.size())STOP("profile.cpp", "profile()","Bad slice index");
					
					slice[islice_periodic].P.push_back(&config.P[ip]);//get the particle into the slice
				
					if(islice==imax) b=R;	
					
					double vol = fabs(PI*(R*R*(b-a) - (b*b*b - a*a*a)/3.));
					
					slice[islice_periodic].volume_in_slice.push_back(vol);//get its volume belonging to the slice
					
					a=b;b+=step;
				}
				
		}
//end dispach_particle_in_slice();

	for(int is=0;is<slice.size();is++) slice[is].average();
}
//
