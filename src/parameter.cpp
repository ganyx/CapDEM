

 /* This file is part of the Soft_Dynamics library.
    The Soft_Dynamics library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.
 */
 



void Cparameter::dimensionless_number(Ccell &cell,std::vector <Cparticle> &P)
{ 
	Cstat mass, diameter,RHO;
	//get the min,mean,max of grains mass and diameter
	
	for(int ip=0;ip<P.size();ip++)
	{
		if(P[ip].AM_I_BOUNDARY==0)
		{
		mass.data.push_back(P[ip].m);
		diameter.data.push_back(2.0*P[ip].R);
		P[ip].RHO = P[ip].m*4.0/(3.0*PI* pow(P[ip].R,3.0));
		RHO.data.push_back(P[ip].RHO);
		}
	}
	Dmin = diameter.MIN();
	Dmean = diameter.MEAN();
	Dmax = diameter.MAX();
	Mmin = mass.MIN();
	Mmean = mass.MEAN();
	Mmax = mass.MAX();
	total_mass = Mmean *P.size();
	RHOmean = RHO.MEAN();	
	
	t_inertia = sqrt(Mmin/(cell.normal_stress_in*Dmin));
	t_collision = sqrt(Mmin/(MODULE_N*Dmin));
	t_thermal = Mmin/Dmin * specific_heat/bulk_conductivity;
//	t_thermal = Mmin/Dmin * specific_heat*bulk_conductivity;
	t_shear = 1./cell.shear_rate;
	
	I = t_inertia / t_shear;
	J = t_thermal / t_inertia;
	K = cell.normal_stress_in / MODULE_N;
}  

ofstream &operator<<(ofstream &file,Cparameter p)
{
	file << p.MODULE_N<<"\t"<< p.friction_coefficient<<"\t"<<p.tang_constant<<"\t"<< p.roll_constant<<"\t";
	file << p.bulk_conductivity<<"\t"<< p.specific_heat<<"\t"<< p.thermal_expansion<<"\t";  
	file << p.Dmean<<"\t"<< p.Mmean<<"\t"<< p.RHOmean<<"\t"<<p.Dmin<<"\t"<<p.Dmax<<"\t";
	file << p.GSD  <<"\t"<< p.fractal_dim <<"\t";
	file << p.J <<"\t"<<p.K <<"\t"<<p.I<<"\t";
	file <<endl;
	return file;
}

ifstream & operator>>(ifstream &file,Cparameter &p)
{
	file >> p.MODULE_N >>p.friction_coefficient>> p.tang_constant>>p.roll_constant;	//c1-4
	file >> p.bulk_conductivity>> p.specific_heat>>p.thermal_expansion;  			//c 5-7
	file >> p.Dmean>>  p.Mmean>> p.RHOmean>>p.Dmin>>p.Dmax;							//c8 -12; 
	file >> p.GSD  >> p.fractal_dim;											//c13 -14; 
	file >> p.J >> p.K >>p.I; 													//c15-17
	return file;
}



