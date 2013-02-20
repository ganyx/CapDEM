
void Cin_out::check_path(char c) 
{
   mode = c;
   if(mode!='w' && mode!='r')STOP("in_out.cpp","init(...)"," Variable'mode' must be 'w' or 'r'"); 
      
   if(c=='r') cout<< "Attempt to read in ";
   if(c=='w') cout<< "Attempt to write in ";
   cout<<"directory: "<<path<<endl;
   if(c=='r'&& (!opendir(path.c_str() )) ) STOP("in_out.cpp","init(...)","The directory does not exist");
   if(c=='w'&& (opendir(path.c_str() )) ) STOP("in_out.cpp","init(...)","The directory already exists");
   else if (c=='w' && mkdir(path.c_str(),0755)==1) STOP("in_out.cpp","init(...)","Can not create the directory");

   if(c=='r') cout<< "Path '"<<path<<"' found: success"<<endl;
   if(c=='w') cout<< "Path '"<<path<<"' created: success"<<endl;
}

void Cin_out::set_file_name(void)
{
  save_file.clear();//remove the old list
  stringstream i_string;  i_string<<current_file; //make a string from integer 
  save_file.push_back( path+"/pos_"+i_string.str()  );
  save_file.push_back( path+"/cont_"+i_string.str() );
  save_file.push_back( path+"/para_"+i_string.str() );
  save_file.push_back( path+"/cell_"+i_string.str() );
  save_file.push_back( path+"/ahistory" ); //YG
//  save_file.push_back( path+"/atime");
//  save_file.push_back( path+"/ascreen" ); //YG
//  save_file.push_back( path+"/agg_"+i_string.str() );

}

bool Cin_out::check_file()
{
//	foreach(string name, save_file)
	for(int iter=0;iter<save_file.size();iter++)
	{
		string name = save_file[iter];
		fstream file;
		file.open(name.c_str(),fstream::in | fstream::out | fstream::app);
		if(file.is_open())	file.close();
	
		else{
		serror="Can not open the file: '"+name+"'\n";
		STOP("inout.cpp", "Cin_out::check_file()",serror);
		return (false);
		} 
	}
return(true);	
}
void Cin_out::PRINT()
{
	for(int iter=0;iter<save_file.size();iter++){
		string s = save_file[iter];
		cout<<"\t"<<s<<endl;
	}
//	 foreach(string s,save_file) cout<<"\t"<<s<<endl;
}
 
 
 

void get_secure(string message, string key, double &value)
{
string choice;

cout<<message<<" ("<<key<<" value)";
cin>>choice>>value;
if(choice!=key){serror="The value'"+choice+"' is not expected. It should be '"+key+"'";
	STOP("inout.cpp", "get_secure(string key, double &value)",serror);}
cout<<"\t="<<value<<endl;
}

void get_secure(string message, string key, int &value)
{
string choice;

cout<<message<<" ("<<key<<" value)";
cin>>choice>>value;
if(choice!=key){serror="The value'"+choice+"' is not expected. It should be '"+key+"'";
	STOP("inout.cpp", "get_secure(string key, double &value)",serror);}
cout<<"\t="<<value<<endl;
}


void get_secure(string message, string key1, string key2, string &choice)
{
	cout<<message<<" ("<<key1<<" or "<<key2<<")";
	cin>>choice;

	if(choice!= key1 && choice!=key2){serror="The value'"+choice+"' does not match the expected choices: '"+key1+"' or'"+key2;
	STOP("inout.cpp", "Crun::get_secure(string, string, string, string &)",serror);}
	cout<<"\t="<<choice<<endl;
}

void get_secure(string message, string key1, string key2, string key3, string &choice)
{
	cout<<message<<" ("<<key1<<" or "<<key2<<" or "<<key3<<")";
	cin>>choice;

	if(choice!= key1 && choice!=key2  && choice!=key3){serror="The value '"+choice+"' does not match the expected choices: '"+key1+"' or '"+key2+"' or '"+key3+"'";
	STOP("inout.cpp", "get_secure(string message, string, string, string, string &)",serror);}
	cout<<"\t="<<choice<<endl;
}

void get_secure(string message, string key, bool &data)
{
	string choice,key_test;
	cout<<message<<" ("<<key<<" YES/NO" <<")";
	cin>>key_test;
	if(key_test!=key) {serror="The value '"+key_test+"' does not match the expected choice: '"+key+"'";
	STOP("inout.cpp", "get_secure(string, string, bool &)",serror);}
	
	cin>>choice;
	if(choice=="YES"){data=true;cout<<"\tYES"<<endl;}
	if(choice=="NO") {data=false;cout<<"\tNO"<<endl;}
	if(choice!="YES" && choice!="NO"){serror="The value '"+choice+"' does not match the expected choices: 'YES' or 'NO' ";
	STOP("inout.cpp", "get_secure(string, string, bool &)",serror);}
}

void get_secure(string message, string key, Cvector &V)
{
	string choice;
	cout<<message<<" ("<<key<<")";
	cin>>choice;	
	if(choice==key) V.READ();
	else{serror="The value '"+choice+"' does not match the expected choices: ('"+key+"')";
	STOP("inout.cpp", "get_secure(string, string, bool &)",serror);}
	cout<<"\t=";
	V.PRINT();	
	cout<<endl;
}
 

