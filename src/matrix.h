class Cvector{
public:
  double x[DIM];
  
  Cvector();                     //intialisation to 0
  //double norm;                 //vector norm
 
  void operator =(Cvector);
  Cvector operator + (Cvector);  // add vector 
  Cvector operator - (Cvector);  // min vector  
  Cvector operator * (Cmatrix);  // transposed(vector) time matrix
  double  operator * (Cvector);  // scalar product
  Cvector operator ^ (Cvector);  // vectorial product 
  Cmatrix operator | (Cvector);  // outer product
  Cvector operator *(double);    // const time vector 
  Cvector operator /(double);    // const div vector
  void operator *=(double);   // time cte to vector
  void operator /=(double);   // div cte to vector
  void operator +=(Cvector);  // add to vector
  void operator -=(Cvector);  // min to vector
 
 // void init();                  //set 0 value
  Cvector curl(Cvector);
 
  friend ofstream &operator<<(ofstream &,Cvector );
  friend ifstream &operator>>(ifstream &,Cvector& );
  void PRINT();
  
  void READ();

  double NORM(void);
  Cmatrix Cross_Product_Matrix();// Build the equi matrix M such as x^u=M.u
};

class Cmatrix{
public:
  double x[DIM][DIM];
 
  Cmatrix();                     //intialisation to 0
 explicit  Cmatrix(double);  
  void operator =(Cmatrix );     //set matrix in matrix
  Cmatrix operator +(Cmatrix );  //add matrix
  Cmatrix operator -(Cmatrix );  //min matrix
  Cmatrix operator *(Cmatrix );  //matrix time matrix 
  Cvector operator *(Cvector );  //matrix time vector
  Cmatrix operator *(double );   //matrix time const
  Cmatrix operator /(double );   //matrix divided const
  void  operator +=(Cmatrix );   //add matrix to matrix 
  void  operator -=(Cmatrix );   //min matrix to matrix 
  void operator +=(double );  //add cte to matrix 
  void  operator -=(double ); //min cte to matrix 
  void  operator *=(double ); //time cte to matrix
  void  operator /=(double ); //div cte to matrix

  friend ofstream &operator<<(ofstream &,Cmatrix );
  friend ifstream &operator>>(ifstream &,Cmatrix & );

  Cmatrix symetric();
  void add_to_trace(double);// add to the diagonal terms
 // void init();                  //set 0 value
  
  double DET();
  Cmatrix INVERSE();
  
  void PRINT();
  void READ();

}MI(1.);








