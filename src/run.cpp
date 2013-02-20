


void Crun::init_evolve(void)
{
	//set where to read the initial configuration
	cout<<"Enter the path where the initial configuration files are (without final /) and the number of the file to be read"<<endl;
	cin>>where_read.path>>where_read.current_file;
	
	where_read.check_path('r');		//check wheter the path already exist. If no, stop.
	config.fread(where_read);		//read the config, would stop if the files don't exit
	cout<<"\tNumber of particle: "<< config.P.size()<<endl;
	cout<<"\tNumber of contact: "<< config.C.size()<<endl;
	cout<<"\tSystem syze: ";
	config.cell.L.PRINT();
	cout<<endl<<endl;

	//set where to save the files during th simulation
	cout<<"Enter the path where to write the saving files (without final /)  and the first number of the file to be written"<<endl;
	cin>>where_save.path>>where_save.current_file;
	
	where_save.check_path('w');		//check wheter the path already exist. If yes, stop.
	config.fprint(where_save);		// save the initial configuration, would stop if the files can't be opened
	
	//Cell boundary input  
	get_secure("Enter the type of boundary","WALL_INCLINED","PERIODIC_SHEAR","WALL_SHEAR",config.cell.boundary);
	
	//what to simulate
	get_secure("Do you want to simule the thermal conduction", "CONDUCTION",	config.simule_thermal_conduction);
	get_secure("Do you want to simule the thermal production", "PRODUCTION",	config.simule_thermal_production);
	get_secure("Do you want to simule the thermal expansion", "EXPANSION",		config.simule_thermal_expansion);
		
	
	string choice;
	if(config.cell.boundary=="WALL_INCLINED")
	{//input for inclined plane 
		get_secure("Enter the gravity","GRAVITY",config.cell.gravity);
		get_secure("Enter the slope angle", "SLOPE",config.cell.slope);
	
		config.cell.normal_stress_control=false;
		config.cell.normal_stress_ext=0;
		config.cell.normal_stress_control=false;
		config.cell.Vdilat=0; config.cell.Adilat=0;
	}
	
	else
	{//input for plane shear
		config.cell.gravity=0;
		config.cell.slope=0;
		config.cell.shear_stress_control=0;
		config.cell.shear_work_control=0;
		config.cell.stick_slip=0;
		
		//get the normal stress or volume control;
		get_secure("Do you want to control the normal stress or the volume", "NORMAL_STRESS","VOLUME", choice);
		if(choice=="NORMAL_STRESS"){config.cell.normal_stress_control=true; cin>>config.cell.normal_stress_ext; }
		if(choice=="VOLUME"){config.cell.normal_stress_control=false; config.cell.Vdilat=0; config.cell.Adilat=0;}
	
		//get the shear stress or shear rate control;
		get_secure("Do you want to control the shear stress or the shear rate",
		"SHEAR_STRESS","SHEAR_RATE",
//		"SHEAR_WORK",
		"STICK_SLIP",
		choice);
		if(choice=="SHEAR_STRESS"){config.cell.shear_stress_control=true; cin>>config.cell.shear_stress_ext;}
		if(choice=="SHEAR_RATE"){config.cell.shear_stress_control=false; cin>>config.cell.shear_rate;}
//		if(choice=="SHEAR_WORK"){config.cell.shear_work_control=true; cin>>config.cell.shear_work_input;}
		if(choice=="STICK_SLIP"){config.cell.stick_slip=true; cin>>config.cell.slip_velocity;
			config.cell.shear_stress_ext = config.cell.shear_stress_in;}
	
		config.cell.gradT_control=false;//apply a gradiant of temperature
	}//end input for plane shear
	//end of cell properties input
	  
	// Start parameter input
	get_secure("Enter the Young's Modulus for normal contact","MODULE_N", config.parameter.MODULE_N);
	get_secure("Enter the coefficient of friction","FRICTION", config.parameter.friction_coefficient);
	get_secure("Enter the constant for tangential contact (usually 1)","TANG_CONSTANT", config.parameter.tang_constant);
	get_secure("Enter the constant of rolling-twising resistance (usually 1)","ROLL_CONSTANT" ,config.parameter.roll_constant);
	
	get_secure("Enter the conductivity of bulk grains","CONDUCTIVITY",config.parameter.bulk_conductivity);
	get_secure("Enter the specific_heat of bulk grains","SPECIFIC_HEAT",config.parameter.specific_heat);
	get_secure("Enter the thermal expansion of bulk grains","THER_EXPANSION",config.parameter.thermal_expansion);
	
	get_secure("Enter the composite fraction", "COMP_FRACTION", config.parameter.COMP_FRACTION);
	get_secure("WETTING or DRYING, NO_LIQUID", "WETTING","DRYING","NO_LIQUID", choice);
	if(choice=="WETTING"){
		LIQUID_TRANSFER = true;
		config.parameter.INITIAL_SATURATION = 0.5;
		config.parameter.FIXED_SATURATION = 0.5;
		config.parameter.CONTACT_ANGLE = PI/2.0;
		}
	if(choice=="DRYING"){
		LIQUID_TRANSFER = true;
		config.parameter.INITIAL_SATURATION = 0.5;
		config.parameter.FIXED_SATURATION = 0.5;
		config.parameter.CONTACT_ANGLE = PI/4.0;
		}
	if(choice=="NO_LIQUID"){
		LIQUID_TRANSFER = false;
		config.parameter.INITIAL_SATURATION = 0.5;
		config.parameter.FIXED_SATURATION = 0.5;
		config.parameter.CONTACT_ANGLE = PI/4.0;
		}
	if(choice!="NO_LIQUID")
	{
	get_secure("Enter the GRAVITY", "GRAVITY", config.GRAVITY);
	get_secure("Enter the MAX_SCAN_SATURATION", "MAX_SCAN", config.MAX_SCAN);
	get_secure("Enter the MIN_SCAN_SATURATION", "MIN_SCAN", config.MIN_SCAN);
	}
	get_secure("Enter the SURFACE_TENSION", "SURFACE_TENSION", config.parameter.SURFACE_TENSION);
//	get_secure("Enter the LIQUID_DIFFUSION", "LIQUID_DIFFUSION", config.parameter.LIQUID_DIFFUSION);
	config.parameter.LIQUID_DIFFUSION = WATER_CONDUCTION_RATIO/config.parameter.SURFACE_TENSION;
	//End of parameter input
		
		
	//time_step
	config.parameter.dimensionless_number(config.cell,config.P); //typical time scale and  dimensionless numbers
	
	dt = config.parameter.t_collision/10.; //some fraction of collision time
	
	if(config.simule_thermal_conduction) if(config.parameter.t_thermal<config.parameter.t_collision)dt = config.parameter.t_thermal/8;//some fraction of thermal time if smaller
	cout<<"Time step: "<< dt<<endl;

	//Particule update
	get_secure("Enter the initial time of the simulation (usually 0)", "T_INIT",tstart);
	get_secure("Enter the final time of the simulation", "T_END",tend);
	
	if(tend<=tstart) sprintf(error," File run.cpp, function init_time(), simulation ends (tend = %.3f) before starting (tstart = %.3f)",tend,tstart);
	cout<<"Total number of iteration: "<< (int) ((tend-tstart)/dt) <<endl;
	
	get_secure("Enter the time to start saving","SAVE_BEGIN",save.next);
	get_secure("Enter the period between two saving","SAVE_PERIOD",save.period);
			
	screen.period=save.period; 
	screen.next=save.next;
 
	if(save.next>tend)
		sprintf(error," File run.cpp, function init_time(), simulation will not save configuration:\n tend = %.3f < save.next =%d",tend,save.next);
	
	config.Voronoi_Update = true;	
	config.update_particle();
	config.iterate(0.0);

	cout<<"Initialisation of the parameter:\t SUCCESS"<<endl;  
}






bool Cevent::should_do (double time)
{
	if(time>=next){next+= period; return true;}  
//	if(time>=next){next+= 1; return true;}       
  return false;
}

void Crun::evolve()
{
cout<<endl<<endl<<"\t\t==== Start the evolution over time ====="<<endl<<endl;
config.cell.shear_stress_in_old = config.cell.shear_stress_in;
int vflag=0;
config.Voronoi_Update = true;

string his_file = where_save.path+"/ahistory";
for(config.t=tstart;config.t<tend;config.t+=dt)   //start time loop
	{
		config.iterate(dt); // make a evolution of configuration over a time step dt, with or without refreshing the neighbours
		
		config.Voronoi_Update = false;
		if(vflag%50 == 0) config.Voronoi_Update = true;
		vflag++;
		
		if( save.should_do(config.t) ) //save if asked
		{ 
			config.fprint(where_save);
			where_save.current_file++;
		}

		if( config.t>=screen.next)     //print on screen if asked
		{ 
			screen.next += 1;
			int nc=0, np=0;
			for(int ip=0;ip<config.P.size();ip++) {
				if(config.P[ip].water_volume >0.0) np++;}
					
			double angle_max = 0.0;
			double angle_min = 10.0;
			double angle_ave = 0.0;
			for(int ic=0;ic<config.C.size();ic++) {
				if(config.C[ic].CONTACT_ANGLE > angle_max) angle_max = config.C[ic].CONTACT_ANGLE;
				if(config.C[ic].CONTACT_ANGLE < angle_min) angle_min = config.C[ic].CONTACT_ANGLE;
				angle_ave += config.C[ic].CONTACT_ANGLE;
				}
				angle_ave /= config.C.size();
			
			cout<<endl<<endl<<"Case: "<<where_save.path.c_str()<<endl; 
			cout<<"Progress: "<< config.t*100./(tend-tstart)<<" \% done"<<"\tTime step: "<<dt<<endl;
			cout<<"Number of grains: "<<config.P.size()<<"\tNumber of contact: "<<config.C.size()<<endl;
			cout<<"MAX CONTACT_ANGLE:"<<180.0/PI*angle_max<<"\tMIN CONTACT_ANGLE:"<<180.0/PI*angle_min
			<<"\tAVE CONTACT_ANGLE:"<<180.0/PI*angle_ave<<endl;
//			cout<<"Number of bonded contacts: "<<nc<<"\tNumber of melted particles: "<<np<<endl;
			cout<<"Number of wetted grains: "<<np<<"\t Water Potential/Saturation: "<<config.cap_pressure<<"\t"<<config.saturation<<endl;
//			cout<<"System temperature: "<< config.parameter.average_temperature<<endl;
//			cout<<"Heat in/out of the system:\t"<<config.heat_in<<"\t"<<config.heat_out<<endl;
//			cout<<"Heat generation modes (N/S/T/R):"<<"\t"<<config.PN<<"\t"<<config.PS<<"\t"<<config.PT<<"\t"<<config.PR<<endl;
			if(config.cell.shear_work_control)
				cout<<"Work control: "<<"\t"<<config.cell.shear_work_input
				<<"\t The system: "<<config.cell.shear_stress_in*config.cell.shear_rate<<endl;		
			config.cell.PRINT();
			
			
		// history file
		ofstream file;	
		file.open(his_file.c_str(), std::ios_base::app); 
		file<<config.t<<"\t"<<config.cell.Xshift<<"\t"<<config.cell.L.x[1]<<"\t"
		<<config.cell.stress.x[1][1]<<"\t"<<config.cell.stress.x[1][0]
		<<"\t"<<config.cap_pressure_mid*1.0/config.parameter.SURFACE_TENSION<<"\t"<<config.saturation<<"\t"
		<<config.cap_pressure*1.0/config.parameter.SURFACE_TENSION<<"\t"<<config.water_content<<endl; 
		file.close();
		}
	}

}



	








