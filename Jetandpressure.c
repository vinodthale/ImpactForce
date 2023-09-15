
// Author: vinod thale  15 sep 2023 
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED            // Smear density and viscosity jumps
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "curvature.h"


#define VOFFOLDER "FVDb0.80V4.00"
scalar pressure[];
//1] qcc -O2 -Wall -D_FORTIFY_SOURCE=0 Jetandpressure.c -o Jetandpressure -lm
//2] ./Jetandpressure
double cfdbvbubblediameter = 0.00;
//struct CFDValues cfdbv;
//double interface_area_Vat(scalar c);
double PreFactor = 2*pi; // 2*pi; for theta integration okay
int main (int argc, char **argv)
{
  //numericalmainvalues(argv, argc, &cfdbv);
  run();
}

event init (t = 0)
{
    //scalar ppool[];
    //scalar viscstressdrop[];
    //scalar viscstresspool[];
    double PressureDropMaxima;
    //double PressurePoolMaxima;
    //double ViscStressPoolMaxima;
   // double ViscStressDropMaxima;
    double timebgn = 0.00;
    double timestp = 0.01;
    double timeend = 2.04;
    ;
    char namefile[500];
    char name[100];
    double timeload;
    FILE *ForceonLfet;
    FILE *SeVAt;
    FILE *Seegr;
    FILE *StresPre;
    char folder[500];
    strcpy(folder, "mkdir ");
    strcat(folder, VOFFOLDER);
    system(folder);
    for (timeload = timebgn; timeload <= timeend; timeload += timestp)
   {
    sprintf (namefile, "intermediate/snapshot-%5.4f", timeload);
    printf ("load the file %s!\r\n", namefile);
    restore (file = namefile);
    sprintf (name, "%s/VOFDb0.80V4.00-%.2f.gnu", VOFFOLDER, timeload);
    FILE *ip = fopen (name, "w");
    output_facets (f, ip);        // ################################################ tracer 
    fclose (ip);
    static int nfb = 0;
      sprintf (name, "%s/PmaxminDb0.80V4.00.txt", VOFFOLDER);
    if (!nfb)
      SeVAt = fopen(name, "w");
    else
    SeVAt = fopen(name, "a");
    //double XYD_Max[4];
    double XYD_Max[4];
    XYD_Max[0]= -1.0e20;
    XYD_Max[1]= 0.00;
    XYD_Max[2]= 0.00;
    XYD_Max[3]= 0.00;
    //foreach_boundary(left)
    foreach() 
    {
      if (XYD_Max[0] < pressure[])
      {
          XYD_Max[0] = pressure[];
          XYD_Max[1] = x;
          XYD_Max[2] = y;
          XYD_Max[3] = Delta;
      }               
    }
    fprintf (SeVAt, "%f %.10f %.10f %.10f %.10f\r\n", timeload, XYD_Max[0], XYD_Max[1], XYD_Max[2], XYD_Max[3]);
    fclose (SeVAt);
    nfb++;
   // Central difference derivates of both velocity components (u and v)
    // in both directions (x and y which correspond to r and z)
    static int nfe = 0;
      sprintf (name, "%s/StressDb0.80V4.00.txt", VOFFOLDER);
    if (!nfe)
      StresPre = fopen(name, "w");
    else
    StresPre = fopen(name, "a"); 
    foreach ()
    {                                                // ################################################ tracer 
    // Pressures in the droplet  find their maxima
    p[]=pressure[]*f[];                                                  // ################################################ tracer  
    // Viscous stresses in the droplet and pool
    //viscstressdrop[]=sqrt(f[]*(dudx[]*dudx[]+dudy[]*dudy[]+dvdx[]*dvdx[]+dvdy[]*dvdy[]))*mymu[];    // ################################################ tracer 
    }
      // Find the maxima of the pressure and viscous stress in the droplet and pool
   PressureDropMaxima=statsf(p).max;
   //ViscStressDropMaxima=statsf(viscstressdrop).max;
   fprintf (StresPre, "%f  %.10f \r\n", timeload, PressureDropMaxima);
   fclose (StresPre);
   nfe++;
   static int nfbc = 0;
      sprintf (name, "%s/SurafaceareapostDb0.80V4.00.txt", VOFFOLDER);
    if (!nfbc)
      Seegr = fopen(name, "w");
    else
    Seegr = fopen(name, "a"); 
    double kn = (12./((1 - cfdbvbubblediameter * cfdbvbubblediameter * cfdbvbubblediameter)*pi));
    double  se;
    double area = 0;
    foreach (reduction (+:area))
    {
      if (f[] > 1e-6 && f[] < 1. - 1e-6)
      {
        //coord p;
        //coord n = mycs (point, f);
        coord p, n = interface_normal (point, f);
        double alpha = plane_alpha (f[], n);
        double s = plane_area_center (n, alpha, &p);
        //area = interface_area (f);
        area +=  PreFactor * s * dv()/Delta;
        //len = line_length_center(n, alpha, &p);
        //se1 += 2.*pi*( y + p.y*Delta )*(len*Delta); // 2*pi*\int_l (r_c)dl
      }
    }
    //Sarea = interface_area_Vat(f);
    //calcualte surface energy 
    //se = Sarea *0.50; 
    //double SE = kn * se;
    //calcualte surface energy 
    se = area *0.50;
    double SE = kn * se;
    fprintf (Seegr, "%f  %.10f %.10f\r\n", timeload, area, SE);
    fclose (Seegr);
    nfbc++;
    //calculate the force on the substrate
    double pleft = 0.;
    double pForce  = 0.;
    static int nff = 0;
      sprintf (name, "%s/ForceDb0.80V4.00.txt", VOFFOLDER);
    if (!nff)
      ForceonLfet = fopen(name, "w");
    else
    ForceonLfet = fopen(name, "a");
    double pdatum = 0, wt = 0;
    foreach_boundary(top){
    pdatum += 2*pi*y*pressure[]*(Delta);
    wt += 2*pi*y*(Delta);
    } 
    if (wt >0){
    pdatum /= wt;
    }
    foreach_boundary(left)
    {
    pForce += 2*pi*y*(Delta)*(pressure[]-pdatum);
    pleft += pressure[];
    }
    boundary((scalar *){f, u.x, u.y, pressure});
    ;
    fprintf (ForceonLfet, "%f  %.10f %.10f\r\n", timeload, pForce, pleft);
    fclose (ForceonLfet);
    nff++;
 }
}

//This function returns the surface area of the interface as estimated using its VOF reconstruction.
//http://basilisk.dalembert.upmc.fr/src/fractions.h#interfacial-area
/*double interface_area_Vat (scalar c)
{
  double area = 0.;
  foreach (reduction(+:area))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) 
    {
      coord p;
      coord n = interface_normal(point, c);
      double alpha = plane_alpha (c[], n); 
      double s = plane_area_center (n, alpha, &p); 
      area += PreFactor*s*dv()/Delta;
    }
  return area;
}*/

event end(t = 0.0)
{
    printf("\r\n-------\r\nEND!\r\n");
}



