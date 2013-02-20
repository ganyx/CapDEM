 /* This file is part of the Soft_Dynamics library.
    The Soft_Dynamics library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.
 */





double Cstat::MIN()
{
  if(data.size()==0) return 0.;
  min =data[0];
//  foreach(double d, data) if(d<min)min=d;
	for(int iter=0; iter<data.size(); iter++){
		double d = data[iter];
		if(d<min)min=d;
	}
  return min;
}

double Cstat::MAX()
{
  if(data.size()==0) return 0.;
  max =data[0];
//  foreach(double d, data) if(d>max)max=d;
	for(int iter=0; iter<data.size(); iter++){
		double d = data[iter];
		if(d>max)max=d;
	}
  return max;	
}

double Cstat::MEAN()
{
  if(data.size()==0) return 0.;
  mean=0;
//  foreach(double d, data) mean += d;
	for(int iter=0; iter<data.size(); iter++){
		double d = data[iter];
		mean += d;
	}
  mean /= (double)data.size();
  return mean;	
}

double Cstat::STD()
{
  if(data.size()==0) return 0.;
  MEAN();
  std=0;
//  foreach(double d, data) std += (d-mean)*(d-mean);
	for(int iter=0; iter<data.size(); iter++){
		double d = data[iter];
		std += (d-mean)*(d-mean);
	}
  std = sqrt( std/ ((double)data.size()) );
  return std;	
}






void Cdistribution::get_distribution(double d)
{
	min=MIN();
	max=MAX();
	if(min==max) STOP("stat.cpp","Cdistribution::get_distribution(double d)", "All the data have the same value");
	int N;
	N =  (int) ((max-min)/d);
	step = (max-min)/N;
	N++;
		
	for(int i=0;i<N;i++)
	{	
		value.push_back(min + i*step + step/.2);
		distrib.push_back(0.);
		distrib_cumul.push_back(0.);	
	}

	std::vector <int> data_i;
//	foreach(double x, data)
	for(int iter=0; iter<data.size();iter++)
	{
		double x = data[iter];
		 data_i.push_back( (int) ((x-min)/step) );
	if(x<min|| x>max)cout<<x<<"\t"<<min<<"\t"<<max<<endl;	
	}
	for(int k=0;k<data.size();k++) 
		{
			if(data_i[k]<0 ||data_i[k]>=distrib.size() ){
				cout<<data_i[k]<<"\t"<<distrib.size()<<endl ;
				STOP("stat.cpp","Cdistribution::get_distribution(double d)", "Can not set the data in the distributuin: index out of range");
				}
			distrib[data_i[k]]++;
		}
	for(int k=0;k<distrib.size();k++) distrib[k]/=data.size();
	
	distrib_cumul[0]=distrib[0];
	for(int k=1;k<distrib.size();k++) distrib_cumul[k]=(distrib_cumul[k-1]+distrib[k]);
}






 






