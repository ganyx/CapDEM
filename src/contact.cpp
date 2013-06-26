

Ccontact::Ccontact(){}


void Ccontact::PRINT()
{	
	cout<<endl;
	cout<<"Contact between particle "<<A<<" and particle "<<B<<endl;
	cout<<"\tPosition of A\t"; pA->X.PRINT();
	cout<<"\tPosition of B\t"; pB->X.PRINT();
	cout<<"\tCenter to center vector (AB)\t"; dX.PRINT();
	cout<<"\tUnit normal vector (from A to B)\t"; nA.PRINT();
	cout<<"\tNormal penetration\t"<<deltaN<<endl;
	cout<<"\tElastic force (that A experiences)\t"; F.PRINT();
	cout<<"\tNormal elastic force\t"; Fn.PRINT();
	cout<<"\tTangential elastic force\t"; Ft.PRINT();
	cout<<"\tViscous force \t"; Fvis.PRINT();
	cout<<"\tEffective radius of particles\t"<<Reff<<endl;
	cout<<"\tEffective Young's Modulus of particles\t"<<E<<endl;
	cout<<"\tFriction coeeficient\t"<<mu<<endl;
	cout<<"\tRadius of contact area\t"<<a;
	cout<<endl;
}
Ccontact::Ccontact(Cparticle *partA, Cparticle *partB, Ccell *c, Cparameter *para)
{
	age=0;
	aB = 0.0; // initiation of bonding contact radius
	deltaNB = 0.0;
	fn=0;
	pA=partA;
	pB=partB;
	A = partA->id;
	B = partB->id;
	cell= c;
	parameter= para;
	Flag_Boundary=0;
	water_volume=0.0;
	dot_water_volume=0.0;
	CONTACT_ANGLE = CONTACT_ANGLE_MIN;
}

bool Ccontact::AM_I_CONTACTING()
{
	if(pA->AM_I_BOUNDARY==-1 && pB->AM_I_BOUNDARY==-1) return false;//no contact between suck particles of the same plane
	if(pA->AM_I_BOUNDARY==-1 && pB->AM_I_BOUNDARY==-2) return false;
	if(pA->AM_I_BOUNDARY==-2 && pB->AM_I_BOUNDARY==-1) return false;
	if(pA->AM_I_BOUNDARY== 1 && pB->AM_I_BOUNDARY== 1) return false;
	if(pA->AM_I_BOUNDARY== 1 && pB->AM_I_BOUNDARY== 2) return false;
	if(pA->AM_I_BOUNDARY== 2 && pB->AM_I_BOUNDARY== 1) return false;
	if(pA == pB) return false;
	set_me_in_main_cell(); //reset dX and dV according to the periodic boundary conditions	
	dx=dX.NORM();
	
	deltaN = dx-pA->R-pB->R;
	if(water_volume == 0.0) { // initial contact with sufficient water supply and smaller enough gap.
		if(pA->water_volume + pB->water_volume >= 0.1 
			&& deltaN <= 0.2*min(pA->R, pB->R) && LIQUID_TRANSFER) 
			return true;
	}
	if((deltaN>=0 && deltaNB <= 0 && !LIQUID_TRANSFER) // Solid contact without water
		|| (deltaN>=0 && deltaNB <= 0 && water_volume <= 1e-30) // non-existing capillary bridge
		|| (water_volume>1e-30 && deltaN >= min(pow(water_volume,0.3333),MAX_CAP_LENGTH) ) ) // existing capillary bridge, rupture
		return false; //there is no contact, and no pre-existing bond, and no capilary force (YG)
	else return true;
}

void Ccontact::EVALE_Geo()
{	
	nA= dX/dx;
	RA = nA*pA->R;    
	RB = nA*(-pB->R); 

	alpha = nA|nA; 
	I_m_alpha = MI-alpha; // MI is the unit diagonal matrix 9 flops 
}

void Ccontact::relative_velocity()// 41 fp
{
	Cvector d = dV - (RB^pB->Ome) + (RA^pA->Ome);//Relative velocity at the contact,  dV is the relative velocity of the grain centers, previously evaluated 
	Cvector dOme = pB->Ome - pA->Ome;
	Cvector	nADot= (I_m_alpha)*(dV/dx);

	uDotT = I_m_alpha*d;
	dOmeN = alpha*(dOme);
	dOmeT = I_m_alpha*(dOme);

	OmeMean = alpha*((pA->Ome+pB->Ome)/2); 
	OmeMean+= nA^nADot; //solid rotation of contact zone
}


void Ccontact::increment_force(double dt)
{
	double fnold,fntest;
	double aS = 0.0;
	double deltaNS, RSeff;
//	double dNB = 0.0;
	double h0 = 0.1; //exp
	double h = h0;
//	double eta = 1.0; // \eta for viscoucity
	
	Cvector Vn;
	Vn = alpha*dV;
	
	//intermediate name, just because shorter name are easier to pick up in the formula

	if(pA->AM_I_BOUNDARY==-2 || pA->AM_I_BOUNDARY==2) 
		{ Reff = pB->R;//plane
		  meff = pB->m;			// YG: effective mass for visco force	
		}
	else 
		if(pB->AM_I_BOUNDARY==-2 || pB->AM_I_BOUNDARY==2) 
			{ Reff=pA->R;//plane
			  meff = pA->m;		// YG: effective mass for visco force		
			}
		else {
			Reff = 2.0*pA->R*pB->R/(pA->R+pB->R);//two particles
			meff = 2.0*pA->m*pB->m/(pA->m+pB->m); // YG: effective mass for visco force
			}
	E = 2.0*( pA->E*pB->E)/(pA->E+pB->E);//geometric mean
			
	deltaNS = deltaN + (pA->R - pA->RS) + (pB->R - pB->RS); 
	RSeff = 2.0*pA->RS*pB->RS/(pA->RS+pB->RS); // effective solid contact radius
	
	// with solid contact, otherwise only molten contact	
	if(deltaNS < 0.0) {
		aS = sqrt(1.0*RSeff*fabs(deltaNS)); // solid contact radius
		h = h0;
	}
	else if(deltaNS >= 0.0) {
		h = deltaNS;
		deltaNS = 0.0;
		if(h < h0) h = h0;
	}
		
	// update deltaNB, aB	
	//if(pA->dRS != 0 || pB->dRS !=0)
	//{		
		//if(fabs(deltaNS) - pA->dRS - pB->dRS <0) // previous step with a solid gap
		//{
			//dNB = fabs(deltaNS) - (nA*Vn) * dt;
			//if(dNB<0) dNB = 0;
		//}
		//else dNB = pA->dRS + pA->dRS;	
		//deltaNB += dNB;
//	}
//	if(deltaNB < 0.0) 
		deltaNB = 0.0;
//	aB = sqrt(RSeff *deltaNB);
	aB = 0.0;
 	if(aS < aB) aS = aB;
 	
	mu= parameter->friction_coefficient;
	ct= parameter->tang_constant;
	cr= parameter->roll_constant;
	
	//exp
	double dwater = 0.0;
//	double KTHETA = 25.0*(CONTACT_ANGLE_MAX - CONTACT_ANGLE_MIN); //testing
	double KTHETA = 10.0/pow(Reff,3.0); //testing

	dwater = pA->water_volume + pB->water_volume - pA->water_volume_old - pB->water_volume_old;
//	dwater = 1.0 - water_volume_old/water_volume;
//	dwater = water_volume - water_volume_old;
//	if(dwater > 0.0)	 dwater = 1.0; 		// WETTING LIMIT, MAX CONTACT ANGLE
//	if(dwater < 0.0)	 dwater = -1.0; 	// DRAINAGE LIMIT, MIN CONTACT ANGLE
	
	CONTACT_ANGLE += dwater * KTHETA; 	
	if(CONTACT_ANGLE > CONTACT_ANGLE_MAX) CONTACT_ANGLE = CONTACT_ANGLE_MAX;
	if(CONTACT_ANGLE < CONTACT_ANGLE_MIN) CONTACT_ANGLE = CONTACT_ANGLE_MIN;
	
// For cases with almost full-saturated state: Contact angle is the maximum value.
//	if(pA->saturation >= MAX_SATURATION_AIR || pB->saturation >= MAX_SATURATION_AIR) 
//		CONTACT_ANGLE = CONTACT_ANGLE_MAX;
	
	if(dt == 0) CONTACT_ANGLE = CONTACT_ANGLE_MIN; //exp
// ????? BUG
//	a = sqrt(2*Reff*fabs(deltaN)); // PR's implementation
	a = sqrt(1.0*Reff*fabs(deltaN)); // exact contact radius
	

	//vector, forces and moments
/*	Fn = nA*(deltaN*E*a);							//normal force
	Ft += (uDotT*E* ct *a*dt)  + (OmeMean^Ft)*dt; 		//increment norm and rotation for tangential force
	Gn += (dOmeN*E*a*a*a*dt) + (OmeMean^Gn)*dt ; 		//rolling moment cr is the numerical constant, dimensionless, for the rolling and twist 
	Gt += (dOmeT*E*a*a*a*dt) + (OmeMean^Gt)*dt ;		//twist moment
*/
	fcap = 0.0;
	if(LIQUID_TRANSFER && water_volume >= 1e-10){
		double R2,D_bridge;
		if(pA->R > pB->R) R2 = pA->R;
		else R2 = pB->R;
		if(deltaNS > 0.0) D_bridge = deltaNS;
		else D_bridge = 0.0;
		
		double factor = PI*parameter->SURFACE_TENSION* sqrt(pA->R *pB->R);
		double V_R3 = water_volume / pow(R2,3.0);
		if(water_volume < 1e-6) V_R3 = 1e-6 / pow(R2,3.0);
		double Log_V_R3 = log(V_R3);
		double coefficient_a = -1.1*pow(V_R3,-0.53);
		double coefficient_b = (-0.148*Log_V_R3-0.96)*CONTACT_ANGLE*CONTACT_ANGLE
			-0.0082*Log_V_R3 + 0.48;
		double coefficient_c = 0.078 + 0.0018*Log_V_R3;

		fcap = factor *(coefficient_c+ exp(coefficient_a*D_bridge/R2 + coefficient_b));
	}

//	double FBond = deltaNB*E*aB;
//	FBond = 0.0;
// for solid contact, deltaNS < 0, deltaNB>=0!, Fn*nA negtive for compression.	
//	Fn = nA*(deltaNS*E*aS - deltaNB*E*aB);				//normal force, imp. before 15.11.2010
	fnold = fn; //renember the old value of normal force norm
	fn = deltaNS*E*aS;
	Fn = nA*(fn + fcap);						//normal force
	Ft += (uDotT*E* ct *aS*dt)  + (OmeMean^Ft)*dt; 		//increment norm and rotation for tangential force
	Gn += (dOmeN*E*aS*aS*aS*dt) + (OmeMean^Gn)*dt ; 	//rolling moment cr is the numerical constant, dimensionless, for the rolling and twist 
	Gt += (dOmeT*E*aS*aS*aS*dt) + (OmeMean^Gt)*dt ;		//twist moment
	
	
	//norm, forces and moments
//	fn = Fn.NORM();
	fn *= -1.0;
	ft = Ft.NORM();
//	fn = nA*Fn;

	gn = Gn.NORM();
	gt = Gt.NORM();

	if(fnold>fn) fntest = fn;//if unloading, (fnold>fn), take the smaller of these two for the sliding criteria 
	else fntest = fnold;	 //if loading, take the smaller of these two for the sliding criteria as well
	if(fntest <0.0) fntest=0.0; // capillary tensile, 0.0 for friction/twisting/rolling

	//if(deltaNB > 0.0) 
	//{ 
		//double aB_aS = aB/aS;
		//double aB_aS_cubic = pow(aB_aS,3.0);
		//double DForce = ft + 2.0*gt/aS * aB_aS_cubic + mu*( nA*Fn + 4.0*gn/aS *aB_aS_cubic);
		//double A_b = PI *aB*aB;
		//double BForce = BOND_STRENGTH * A_b;
		//double DDamage;
		//if(DForce >= BForce)
		//{	// Debonding by either normal stress or shear stress
			// Damage model!
			//DDamage = 10.0* DForce/BForce *dt;
			//if(DDamage >= 1.0) DDamage=1.0;
			//deltaNB *= 1.0 - DDamage;
		//}
		//if(deltaNB > 0.0){
			//cr *= 100.0; mu *= 100.0;
		//}
	//}
	//if(deltaNB <= 0.0 && (pA->RS < pA->R || pB->RS < pB->R)){ //exp, melting, but no bonding
		//cr *= 0.0; mu *= 0.0;
//	}
	
	TANGENTIAL_SLIDE 	= rescale_slide(Ft,ft, mu*fntest); 
	ROLLING_SLIDE 		= rescale_slide(Gn,gn, cr*aS*fntest);//cr is the empirical numerical constant, usually equal to 1
	TWIST_SLIDE 		= rescale_slide(Gt,gt, cr*mu*aS*fntest);

	//if(deltaNB >0.0){ //exp, bonding
	//TANGENTIAL_SLIDE 	= false;
	//ROLLING_SLIDE 		= false;
	//TWIST_SLIDE 		= false;
	//}

//	debonding due to sliding/rolling/twisting
//	if(deltaNB >0.0 && (TANGENTIAL_SLIDE || ROLLING_SLIDE || TWIST_SLIDE)) deltaNB = 0.0;

	//if(pA->RS < pA->R || pB->RS < pB->R){ //exp, melting
	//TANGENTIAL_SLIDE 	= false;
	//ROLLING_SLIDE 		= false;
	//TWIST_SLIDE 		= false;
	//}

	double xi;
	xi = 0.1;
	Fvis = Vn *xi;

	Fela=Fn+Ft;

	F = Fela+Fvis;	//sum of force
	G = Gn+Gt;		//sum of moments

	//mechanical dissipation
	if(deltaNB >0.0) production_normal = 0;
	else production_normal = fabs(Vn*Fvis); 
	if(TANGENTIAL_SLIDE) production_slide = fabs(uDotT*Ft);
	else production_slide=0;
	if(ROLLING_SLIDE)  production_rolling = fabs(dOmeT*Gt);
	else production_rolling =0;
	if(TWIST_SLIDE) production_twist = fabs(dOmeN*Gn);
	else production_twist =0;

	production= production_normal+production_slide +production_twist+production_rolling;
}

bool Ccontact::rescale_slide(Cvector &F, double &F_norm, double F_norm_max)
{
	if (F_norm< F_norm_max) return false;
	if(F_norm==0)return false;
	F*=  F_norm_max/F_norm;
	F_norm =  F_norm_max;
	return true;
}

void Ccontact::EVALE_heat_flow()
{
	conductivity = 2* pA->k*pB->k/( pA->k+pB->k); //geo average
	phi = 2.*conductivity*a*dT;
}	


void Ccontact::set_me_in_main_cell()
{
	Flag_Boundary = 0;
	dX =	pB->X - pA->X;  
	dV =	pB->V - pA->V;
	dT =	pB->T - pA->T;
	
	
	if( pA->AM_I_BOUNDARY == 2 || pA->AM_I_BOUNDARY == -2 || pB->AM_I_BOUNDARY == 2 || pB->AM_I_BOUNDARY == -2 ) { dX.x[0]=0; dX.x[2]=0;}
	
	//periodic along x and z, true for avery kind of cell.boundary (WALL_INCLINED,WALL_SHEAR or PERIODIC_SHEAR
	cell->rescale(dX.x[0],cell->L.x[0]);//check if out fron the right/left sides
	cell->rescale(dX.x[2],cell->L.x[2]);//check if out fron the front/back sides	
	
	if(cell->boundary=="PERIODIC_SHEAR")	//specificity of PERIODIC_SHEAR
	{
		if (dX.x[1]>=cell->L.x[1]/2.)		//out fron the top side
			{
				dX.x[1]-=cell->L.x[1];
				dX.x[0]-=cell->Xshift;	//shift particle along the right axis
				
// BUG !, after the following check, dX.x[0] might be larger than cell->L.x[0]/2., since Xshift~[-L/2,L/2)
// use function cell.rescale()	
//				if (dX.x[0]<-cell->L.x[0]/2.)	dX.x[0]+=cell->L.x[0];	//re-check out from the left side
				cell->rescale(dX.x[0], cell->L.x[0]);  //YG
				
				dV.x[0]-=cell->Vshear;	//increment the velocity  
				dV.x[1]-=cell->Vdilat;
				dT-=cell->DeltaT;			//increment the temperature
				Flag_Boundary = 1;
			}	

		else if (dX.x[1]<-cell->L.x[1]/2.)//out from the bottom side
			{
				dX.x[1]+=cell->L.x[1];	//out fron the top side
				dX.x[0]+=cell->Xshift;	//shift particle along the right axis
				
// BUG !, after the following check, dX.x[0] might be smaller than cell->L.x[0]/2., since Xshift~[-L/2,L/2)
// use function cell.rescale()			
//				if (dX.x[0]>=cell->L.x[0]/2.)	dX.x[0]-=cell->L.x[0];	//re-check if out from the right side
				cell->rescale(dX.x[0], cell->L.x[0]);  //YG
			
				dV.x[0]+=cell->Vshear;	//increment the velocity  
				dV.x[1]+=cell->Vdilat;
				dT+=cell->DeltaT;			//increment the temperature
				Flag_Boundary = 1;
			}	
	}
}



 ofstream & operator<<(ofstream &file,Ccontact c)
 {
	file<<c.A<<"\t"<<c.B<<"\t"; //c 1-2
	file<<c.Ft<<c.Gn<<c.Gt;
	file<<c.fn<<"\t";
	file<<c.age<<"\t";
	file<<c.fcap<<"\t";
	file<<c.fwater<<"\t";
	file<<c.dot_water_volume<<"\t";
	
// test output
//	file<<c.voronoi_area<<"\t";
//	file<<c.water_volume<<"\t";
	file<<endl;
	return file;
}
	
ifstream & operator>>(ifstream &file,Ccontact &c)
{
	file>>c.A;
	file>>c.B;
	file>>c.Ft>>c.Gn>>c.Gt;
	file>>c.fn;
	file>>c.age;
	file>>c.fcap;
	file>>c.dot_water_volume;
	return file;
 }



/*void Ccontact::operator = (Ccontact para)
{
	pA=para.pA;
	pB=para.pB;
	A=para.A;
	B=para.B;
	dx=para.dx;
	dX=para.dX;
	deltaN = para.deltaN;

	fn=para.fn;
	Ft=para.Ft;
	Gn=para.Gn;
	Gt=para.Gt;
	
	age = para.age;
}*/
/*
	dX = para.dX;
	dV = para.dV;
	RA = para.RA;
	RB = para.RB;
	nA = para.nA;
	a = para.a;
	deltaN=para.deltaN;

	dT=para.dT;
	production= para.production;
	production_normal=para.production_normal;
	production_slide=para.production_slide;
	production_rolling=para.production_rolling;
	production_twist=para.production_twist;

	fn=para.fn;
	ft=para.ft;
	Fn=para.Fn;
	Ft=para.Ft;
	Gn=para.Gn;
	Gt=para.Gt;
	G = para.G;
  
  Fvis=para.Fvis;
  F = para.F;*/
//}

