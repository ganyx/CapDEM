
#ifndef IN_OUT_H
#define IN_OUT_H


#define AFF(x) {cout<<"\t\033[32m"<< #x<<"=\t"<<x <<"\n \033[00m"<<endl;}//print on screen the name of the data and its value

char error[300];
void ERROR(char *er);            //Print error text in red
void ERROR(char *er)
{  cout <<"\n\e[41;30mERROR\n" <<er<< " \033[0m\n"<<endl; exit(0);}

string serror;
void STOP(string file, string function, string error){
   cout<< endl << endl;	
   cout<< "===================================="<< endl;	
   cout<<"\tError in file " <<file <<", function "<< function <<endl<<"\t"<<error<<endl;
   cout<< "===================================="<< endl<< endl;	
	exit(0);
};/*<<Global function to display error.*/

//
class Cin_out  
{
public:
    std::vector <string> save_file;
	void check_path(char c);
	void set_file_name();
	bool check_file(); /**< Return true if all the file in the list save_file can be opened, false if not*/ 
	void PRINT(void);
	
	char mode; // mode: read or write
	int  current_file;// current file number where it is saved
	string path;// directory of the saving

};
#endif
