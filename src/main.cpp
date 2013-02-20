#include "soft_dynamics.h" 

//int main(int argc, char **argv)
int main(void)
{
omp_set_num_threads(NTHREADS); // MPI

     clock_t start, end;
     time_t  t0, t1;
     double cpu_time_used;
     double wall_time_used;

	
cout<<endl<<endl<<"Verion 2012 for Capillary Water Transport"<<endl<<endl;
cout<<"Type the action you want to perform:"<<endl;
cout<<"\tCREATE: creates a new random configuration"<<endl;
cout<<"\tEVOLVE: reads an existing configuration and makes it evolves"<<endl;
cout<<"\tPOST_PROCESS: reads an existing configuration and measures all sort of averages"<<endl;
//cout<<"\tVISUALISATION: reads an existing configuration and draws it on an opengl basis"<<endl;

string action;
cin>> action;

if(action!="CREATE" && action!="EVOLVE" && action!="POST_PROCESS")
	{serror="The action '"+action+"' does not exit. Choices are  'CREATE', 'EVOLVE' or 'POST_PROCESS'";
	STOP("main.cpp", "main()",serror);}
	
cout<<"The selected action is: "<<action<<endl<<endl;


if(action=="CREATE")//not parallel
{	
	Crun_create run;

	cout<<"Enter the path where to write the saving files (without final /) and the number of the file to be written (usually 1)"<<endl;
	cin>>run.where_save.path>>run.where_save.current_file;
	run.where_save.check_path('w');//check wheter the path already exist. If yes, stop.	
	
	run.config.create_random();	//create
	run.config.fprint(run.where_save);//save
}


if(action=="EVOLVE")//can be parallel
	{
	Crun Lrun[1];//24 is the max number of processor
	
	int Ntask;
	for(Ntask=0;Ntask<1;Ntask++)
		{	
		Lrun[Ntask].init_evolve();
	    string one_more;
	   	get_secure("Do you want to add one more task?","ONE_MORE_TASK","NO_MORE_TASK", one_more);
//		if(one_more=="ONE_MORE_TASK")cout<<"Ready to get the new task"<<endl;
	    if(one_more=="NO_MORE_TASK"){cout<<"No more task, ready to simulate"<<endl;Ntask++;break;}
		}
//	cout<<"Number of task to be run in parallel:\t"<<Ntask<<endl;	
//	for(int n=0;n<Ntask;n++)Lrun[n].start();// call the function evolve() for each run
//	for(int n=0;n<Ntask;n++)Lrun[n].wait();// wait for every one to finish
	start = clock(); t0 = time(NULL);
	Lrun[0].evolve();
	
	end = clock();t1 = time(NULL);
//	cpu_time_used = difftime(end, start);
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC / NTHREADS;
	wall_time_used = t1-t0;
	cout<<"Number of CPUs:\t"<<NTHREADS<<endl;
	cout<<"CPU time (s):\t"<<cpu_time_used<<endl;
	cout<<"Wall time (s):\t"<<wall_time_used<<endl;
			
//	cout<<"Number of task to be run in parallel:\t"<<Ntask<<endl;	
//	for(int n=0;n<Ntask;n++)Lrun[n].start();// call the function evolve() for each run
//	for(int n=0;n<Ntask;n++)Lrun[n].wait();// wait for every one to finish
	}

if(action=="POST_PROCESS")
{
	Crun_post_process Lrun[1];//24 is the max number of processor

	int Ntask;
	for(Ntask=0;Ntask<1;Ntask++)
		{	
		Lrun[Ntask].init();
	    string one_more;
	   	get_secure("Do you want to add one more task?","ONE_MORE_TASK","NO_MORE_TASK", one_more);
		if(one_more=="ONE_MORE_TASK")cout<<"Ready to get the new task"<<endl;
	    if(one_more=="NO_MORE_TASK"){cout<<"No more task, ready to simulate"<<endl;Ntask++;break;}
		}
	cout<<"Number of task to be run in parallel:\t"<<Ntask<<endl;	
//	for(int n=0;n<Ntask;n++)Lrun[n].start();// call the function evolve() for each run
//	for(int n=0;n<Ntask;n++)Lrun[n].wait();// wait for every one to finish
}

	return(1);
}







      
     
