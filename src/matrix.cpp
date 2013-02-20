 /* This file is part of the Soft_Dynamics library.
    The Soft_Dynamics library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.
 */

//====================================================
//     class Cvector functions and operations
//     class definition is in "matrix.h"
//====================================================

Cvector::Cvector(){int i; for(i=0;i<DIM;i++)x[i]=0.;}

ofstream & operator<<(ofstream &file,Cvector v )
{
	for(int i=0;i<DIM;i++) file<<v.x[i]<<"\t"; return file;
}
ifstream & operator>>(ifstream &file,Cvector & v)
{
	for(int i=0;i<DIM;i++) file>>v.x[i];return file;
}
 
void Cvector::PRINT()
{  printf("%f\t%f\t%f\n",x[0],x[1],x[2]);}


void Cvector::READ()
{  cin>>x[0]>>x[1]>>x[2];}
	

/** Assert a vector value. If U and V are two vetors: \code V=U; \endcode each component of V is equal to that of U.*/
void Cvector::operator =(Cvector para)
{for(int i=0;i<DIM;i++) x[i]=para.x[i]; }

/** Sum of two vector. If U,V,W are three vectors: \code U=V+W; \endcode set in U the sum of V and W.   */

Cvector Cvector::operator + (Cvector para)
{Cvector res;  for(int i=0;i<DIM;i++)  res.x[i]=x[i]+para.x[i]; return res;}

/** Difference between two vectors. If U,V,W are three vectors: \code U=V-W; \endcode set in U the difference of V and W.   */
Cvector Cvector::operator - (Cvector para)
{Cvector res;  for(int i=0;i<DIM;i++)  res.x[i]=x[i]-para.x[i]; return res;}

/** Vector time a double. If U and V are two vectors and d is a double: \code  U = V*d;  \endcode for each component of U is equal to that of V times d; 
\attention U = d*V; is not valid. */
Cvector Cvector::operator *(double para )
{Cvector res;  for(int i=0;i<DIM;i++) res.x[i]=para*x[i];  return res;}
/** Vector divided by double. If U and V are two vectors and d is a double: \code U = V/d; \endcode for each component of U is equal to that of V divided by d; \attention U = V/0; is not safe */
Cvector Cvector::operator /(double para )
{Cvector res;  for( int i=0;i<DIM;i++) res.x[i]=x[i]/para;  return res;}

/** Scalar product. If U and V are two vectors, s a double:  \code s = U*V; \endcode set \f$ s = \sum_i U_iV_i \f$*/
double Cvector::operator *(Cvector para )
{double scalar=0.; for(int i=0;i<DIM;i++) scalar+=x[i]*para.x[i]; return scalar;}
/** Matrix time vector. If U and V are two vectors, and M is a matrix:  \code U = M*V; \endcode set \f$ U_i = \sum_j M_{ij} V_j \f$*/
Cvector Cvector::operator *(Cmatrix para)
{Cvector res; 
for(int i=0;i<DIM;i++) 
  for(int k=0;k<DIM;k++)res.x[i]+=x[k]*para.x[k][i]; return res;}
/** Cross product. If U, V and Z are three vectors: \code W = U^V; \endcode set in W the cross product of U and W.\attention The priority of this operation is not included: W = U + U^V will do W = (U+U)^V*/

Cvector Cvector::operator ^(Cvector para )
{Cvector res ; 
  res.x[0] = x[1]*para.x[2] - x[2]*para.x[1];
  res.x[1] = x[2]*para.x[0] - x[0]*para.x[2];
  res.x[2] = x[0]*para.x[1] - x[1]*para.x[0];
  return res;
}

/** Outer product. If U, V are two vectors and M a matrix: \code M = U|V; \endcode set \f$ M_{ij} = U_iV_j \f$*/
Cmatrix Cvector::operator | (Cvector para)
{ Cmatrix res; 
for(int i=0;i<DIM;i++) 
  for(int j=0;j<DIM;j++) res.x[i][j]=x[i]*para.x[j]; return res;}

 void Cvector::operator *=(double para) 
{  for(int i=0;i<DIM;i++) x[i]*=para;}
 void Cvector::operator /=(double para) 
{  for(int i=0;i<DIM;i++) x[i]/=para;}
 void Cvector::operator +=(Cvector para) 
{  for(int i=0;i<DIM;i++) x[i]+=para.x[i];}
void Cvector::operator -=(Cvector para) 
{  for(int i=0;i<DIM;i++) x[i]-=para.x[i];}



Cmatrix Cvector::Cross_Product_Matrix()
{
  Cmatrix res ; 
  res.x[0][0] =  0.;    res.x[0][1] = -x[2];  res.x[0][2] = x[1]; 
  res.x[1][0] =  x[2]; res.x[1][1] =  0.;     res.x[1][2] =-x[0];
  res.x[2][0] = -x[1]; res.x[2][1] =  x[0];  res.x[2][2] = 0.;
  return res;  
}

double Cvector::NORM( )
{double norm=0.; for(int i=0;i<DIM;i++) norm+=x[i]*x[i]; return sqrt(norm);}


Cvector Cvector::curl(Cvector U)
{
 Cvector res;
  
 res.x[0]= x[2]/U.x[1]- x[1]/U.x[2];
  res.x[1]=x[0]/U.x[2]- x[2]/U.x[0];
   res.x[2]=x[1]/U.x[0]- x[0]/U.x[1];
   return res;
}


Cmatrix::Cmatrix(){for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]=0.;}
Cmatrix::Cmatrix(double para){ for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) if(i==j)x[i][j]=para; else x[i][j]=0.;}

ofstream & operator << (ofstream &file, Cmatrix M)
{	for(int i=0;i<DIM;i++) 
	 for(int j=0;j<DIM;j++) 
	 	file << M.x[i][j]<<"\t"; 
	 	return file;
	 	}
	
	
ifstream & operator>>(ifstream &file, Cmatrix & M)
{	for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++)file>>M.x[i][j];return file; }
 
void Cmatrix::PRINT()
{
  printf("\n");
  for(int i=0;i<DIM;i++) {for(int j=0;j<DIM;j++) printf("%f\t",x[i][j]); printf("\n");}
  printf("\n");
}

void Cmatrix::READ(){
for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) scanf("%lf",&x[i][j]);}


Cmatrix Cmatrix::symetric()
{ Cmatrix res;  for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) res.x[i][j] = (x[i][j]+x[j][i])/2.;return res;}



void Cmatrix::add_to_trace(double para)
{ for(int i=0;i<DIM;i++) x[i][i]+=para;}

 double Cmatrix::DET()
{
 double det=0;
         
 det = x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]);
 det -= x[0][1]*(x[1][0]*x[2][2]-x[2][0]*x[1][2]);
 det += x[0][2]*(x[1][0]*x[2][1]-x[2][0]*x[1][1]);
	return det;
}

Cmatrix Cmatrix::INVERSE()
{
	Cmatrix res, B; double det;
	det = x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]);
 	det -= x[0][1]*(x[1][0]*x[2][2]-x[2][0]*x[1][2]);
 	det += x[0][2]*(x[1][0]*x[2][1]-x[2][0]*x[1][1]);
 	
 	if(det!=0){
 		for(int i=0;i<3;i++)
          for(int j=0;j<3;j++)    
               B.x[i][j]=x[j][i];
               
     res.x[0][0]=B.x[1][1]*B.x[2][2]-(B.x[2][1]*B.x[1][2]);
     res.x[0][1]=(-1)*(B.x[1][0]*B.x[2][2]-(B.x[2][0]*B.x[1][2]));
     res.x[0][2]=B.x[1][0]*B.x[2][1]-(B.x[2][0]*B.x[1][1]);
     res.x[1][0]=(-1)*(B.x[0][1]*B.x[2][2]-B.x[2][1]*B.x[0][2]);
     res.x[1][1]=B.x[0][0]*B.x[2][2]-B.x[2][0]*B.x[0][2];
     res.x[1][2]=(-1)*(B.x[0][0]*B.x[2][1]-B.x[2][0]*B.x[0][1]);
     res.x[2][0]=B.x[0][1]*B.x[1][2]-B.x[1][1]*B.x[0][2];
     res.x[2][1]=(-1)*(B.x[0][0]*B.x[1][2]-B.x[1][0]*B.x[0][2]);
     res.x[2][2]=B.x[0][0]*B.x[1][1]-B.x[1][0]*B.x[0][1];

 		res /= det;
		}
	  return res;
}
       
       
void Cmatrix::operator =(Cmatrix para )
{for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]=para.x[i][j];}

Cmatrix Cmatrix::operator +(Cmatrix para )
{ Cmatrix res;
  for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) 
      res.x[i][j]=x[i][j]+para.x[i][j];
  return res;
}
Cmatrix Cmatrix::operator -(Cmatrix para )
{Cmatrix res;
  for(int i=0;i<DIM;i++)for(int j=0;j<DIM;j++) res.x[i][j]=x[i][j]-para.x[i][j];
  return res;}

Cmatrix Cmatrix::operator *(double para )
{ Cmatrix res;
  for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) res.x[i][j]=x[i][j]*para;
  return res;
}
Cmatrix Cmatrix::operator /(double para )
{Cmatrix res;
  for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) res.x[i][j]=x[i][j]/para;
  return res;
}


Cmatrix Cmatrix::operator *(Cmatrix para )
{Cmatrix res;
  for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) for(int k=0;k<DIM;k++) res.x[i][j]+= x[i][k]* para.x[k][j];
  return res;
}

Cvector  Cmatrix::operator *(Cvector para )
{Cvector res;
  for(int i=0;i<DIM;i++) for(int k=0;k<DIM;k++) res.x[i]+=x[i][k]*para.x[k];
  return res;
}

 void Cmatrix::operator +=(Cmatrix para)
{for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]+=para.x[i][j];}
 
 void Cmatrix::operator -=(Cmatrix para)
{ for( int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]-=para.x[i][j];}


 void Cmatrix::operator +=(double para)
{for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]+=para;}
 
 void Cmatrix::operator -=(double para)
{for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]-=para;}


void Cmatrix::operator *=(double para)
{for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j]*=para;}
 
 void Cmatrix::operator /=(double para)
{for(int  i=0;i<DIM;i++) for(int  j=0;j<DIM;j++) x[i][j]/=para;}

