

// place your code here
void Crun_post_process::init(void)
{
	cout<<"Enter the path where the configuration files are (without final /) and the number of the first and last (included) files to be read"<<endl;
	cin>>where_read.path>>first>>last;	

	

	path_to_write = where_read.path+"/post_process";
	if (!opendir(path_to_write.c_str() ))  //if the directory does not exist
		if(mkdir(path_to_write.c_str(),0755)==1) STOP("run.cpp","init(...)","Can not create the directory");//and if we can't create it: error

	cout<<"Start reading the files from "<<first <<" to "<< last<<endl;
		

	where_read.current_file=first;
	for(int iconf = 0; iconf<=last-first;iconf++)
	{
		config.fread(where_read);
		where_read.current_file++;
//		list_config.push_back(config);
	}
	cout<<"Start reading the files = SUCCESS"<<endl<<endl;
}



void Crun_post_process::post_process()
{

	for(int iconf=0;iconf<list_config.size();iconf++)
	{
		
		Cdistribution distrib_mass, distrib_diameter;
	
		//EXAMPLE TO ACCESS CONTACT AND PARTICLES
		//	foreach(Ccontact cont, list_config[iconf].C)
//		foreach(Cparticle part,  list_config[iconf].P)  
		for(int iter; iter<list_config[iconf].P.size(); iter++)
			{
				Cparticle part = list_config[iconf].P[iter];
				distrib_mass.data.push_back(part.m);
				distrib_diameter.data.push_back(2.*part.R);
			}
		distrib_mass.get_distribution(0.01); //get the distribution of the list data, with a step close to 0.01
		distrib_diameter.get_distribution(0.01);
			
		
		//do the profiles + example of write in files
		Cprofile profile(0.3, list_config[iconf]);	//the first parameter is the size of the slice of averaging
		
		
		
		
		
		string file_name;//example to write data for each saving files. Yes, it's ma pain in the ass
		stringstream i_string;  i_string<<int(iconf+first);
		ofstream file;
		/*
		file_name = path_to_write+"/profile_"+i_string.str()  ;
		file.open(file_name.c_str());
		for(int is=0;is<profile.slice.size();is++)//for each slice, we print something
				file<<profile.slice[is]; 	//to know what's plot, see the definition of the operator in profile.cpp
		file.close();
	
	
		file_name = path_to_write+"/distrib_mass_"+i_string.str()  ;
		file.open(file_name.c_str());
		for(int i=0;i<distrib_mass.distrib.size();i++)//for each slice, we print something
				file<< distrib_mass.value[i] <<"\t"<<distrib_mass.distrib[i]<<"\t"<<distrib_mass.distrib_cumul[i]<<endl; 	//to know what's plot, see the definition of the operator in profile.cpp
		file.close();
		*/
	/*	file_name = path_to_write+"/distrib_diameter_"+i_string.str()  ;
		file.open(file_name.c_str());
		for(int i=0;i<distrib_diameter.distrib.size();i++)//for each slice, we print something
				file<< distrib_diameter.value[i] <<"\t"<<distrib_diameter.distrib[i]<<"\t"<<distrib_diameter.distrib_cumul[i]<<endl; 	//to know what's plot, see the definition of the operator in profile.cpp
		file.close();
		*/
	}
	
	
}
