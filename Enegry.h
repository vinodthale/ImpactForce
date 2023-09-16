// Author: vinod thale  15 sep 2023 
double de, dee, se;
  double vd = 0, VD = 0;
	double ke1 = 0, ke2 = 0;
	double ke = 0;
	double area = 0;
	//double gpe = 0;
	double kn = (12./((1 - cfdbv.bubblediameter * cfdbv.bubblediameter * cfdbv.bubblediameter)*pi));
	double PreFactor = 2*pi; // 2*pi; for theta integration okay 
	foreach (reduction(+:ke1) reduction(+:ke2) reduction(+:ke) reduction(+:area) reduction(+:vd) reduction(+:VD))
		//kinetic energy of both gas and liquid
		ke1 +=  PreFactor*0.5*rho1*dv()*f[]*sq(u.x[]) + PreFactor*0.5*rho2*dv()*(1 - f[])*sq(u.x[]); // KE-X direction 
		ke2 +=  PreFactor*0.5*rho1*dv()*f[]*sq(u.y[]) + PreFactor*0.5*rho2*dv()*(1 - f[])*sq(u.y[]); // KE-Y direction 
		//the sum of kinetic energies of liquid and gas on vertical and horizontal components 
                 //kinetic energy contrast(only liquid)
		ke +=  PreFactor*0.5*rho1*(sq(u.x[]) + sq(u.y[]))*dv()*f[];
		if (f[] > 1e-6 && f[] < 1. - 1e-6)
		{
			coord p;
			coord n = mycs (point, f);
			double alpha = plane_alpha (f[], n);
			double s = plane_area_center (n, alpha, &p);
			//area = interface_area (f);
			area +=  PreFactor*s*dv()/Delta;
		}
		;
	    //calculate viscous dissipation of the liquid
		vd +=  PreFactor*f[]*mu1*dv()*((2.*(sq(u.x[1] - u.x[-1]) + sq(u.y[0, 1] - u.y[0, -1])) +
			sq(u.y[1] - u.y[-1] + u.x[0, 1] - u.x[0, -1])) / sq(2.*Delta) -
			(2./3.)* (u.y[0, 1] - u.y[0, -1] + u.x[1] - u.x[-1])/(2.*Delta) + 2*sq(u.y[]/y) - (2./3.)*u.y[]/y);
		//calculate viscous dissipation of the liquid and gas
		VD +=  PreFactor*f[]*mu1*dv()* ((2.*(sq(u.x[1] - u.x[-1]) + sq(u.y[0, 1] - u.y[0, -1])) +
			sq(u.y[1] - u.y[-1] + u.x[0, 1] - u.x[0, -1])) / sq(2.* Delta) -
			(2./3.)*(u.y[0, 1] - u.y[0, -1] + u.x[1] - u.x[-1])/(2.*Delta) + 2*sq(u.y[]/y) - (2./3.)*u.y[]/y) +
			 PreFactor*(1 - f[]) * mu2 * dv()*((2.*(sq(u.x[1] - u.x[-1]) + sq(u.y[0, 1] - u.y[0, -1])) +
				sq(u.y[1] - u.y[-1] + u.x[0, 1] - u.x[0, -1])) / sq(2. * Delta) -
				(2./3.)*(u.y[0, 1] - u.y[0, -1] + u.x[1] - u.x[-1])/(2.*Delta) + 2*sq(u.y[]/y) - (2./3.)*u.y[]/y);
	}
	//calculate total kinetic energy(both) 
	double KE1 = kn*ke1;
	double KE2 = kn*ke2;
	double KE = KE1+KE2;
	//calculate total kinetic energy(only liquid)	
	double KEE = kn*ke;
	//calcualte surface energy 
	se = area*cfdbv.Sigma;
	double SE = kn*se;
    //calculate dissipation energy(only liquid)
	de += vd*dt;
	double DE = kn*de;
	//calculate dissipation energy(both)
	dee += VD*dt;
	double DEE = kn*dee;
	double TENG = KE+SE+DEE ; //the total energy no gravity


