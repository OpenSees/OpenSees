// Date: 12/08/2004
//
// File name: SoilFootingSection2d.cpp
//
// Coded by: Sivapalan Gajan <sgajan@ucdavis.edu>
//
// Description: This file contains the members and methods for 
// SoilFootingSection2d class.


#include <math.h>
#include <stdlib.h>

#include <SoilFootingSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>

#include <classTags.h>

ID SoilFootingSection2d::code(3);



// default constructor

SoilFootingSection2d::SoilFootingSection2d(void)
   :SectionForceDeformation(0, SEC_TAG_SoilFooting2d),
    e(3), s(3),eCommit(3), sCommit(3), deModel(3), ks(3,3), ksE(3,3), ini_size(3)
{
   code(0) = SECTION_RESPONSE_P;   // vertical load
   code(1) = SECTION_RESPONSE_VY;  // shear
   code(2) = SECTION_RESPONSE_MZ;  // moment
}


// constructor

SoilFootingSection2d::SoilFootingSection2d
   (int tag, double fs, double vult, double l, double kv, double kh, double rv, double deltaL)
   :SectionForceDeformation(tag, SEC_TAG_SoilFooting2d),
   e(3), s(3), eCommit(3), sCommit(3), deModel(3), ks(3,3), ksE(3,3), ini_size(3), 
   FS(fs), Vult(vult), L(l), Kv(kv), Kh(kh), Rv(rv), dL(deltaL)
{
   if (FS <= 1.0)
   {
      opserr <<"SoilFootingSection:: Invalid input for FS\n"
             <<"FS should satisfy: FS > 1.0\n";
   }

   code(0) = SECTION_RESPONSE_P;   // axial load
   code(1) = SECTION_RESPONSE_VY;  // shear
   code(2) = SECTION_RESPONSE_MZ;  // moment

   V = Vult / FS;
   qult = Vult/L;
   ecc = 0.0;
   dTheta = 0.0;
   dThetaPrev = 0.0;
   dVtemp = 0.0;
   tolerance = 0.0001 * (1.0/FS);
   tolerance = 0.01;
   Kt = Kv*pow(L, 3.0)/12.0;

   incr = 0;
   noNodes = (int)(ceil(L/dL))+1;

   c1 = 0;
   c1T = 0;
   c2 = noNodes;
   c2T = noNodes;

   c1Commit = c1;
   c1TCommit = c1T;
   c2Commit = c2;
   c2TCommit = c2T;

   eccCommit = ecc;

   hPrev = -10.0;
   hCurr = -10.0;
   dVtemp = 0.0;
   Mmaxpast = 0.0;
   Mult = 0.0;
   Melastic = 0.0;


   initializeBoundingSurface();
   initializeInternalVariables();


   isOver = 0;
   isdV = 0;
   isElastic = 1;

}



// destructor


SoilFootingSection2d::~SoilFootingSection2d(void)
{
   if (foot != 0)
      for (int i = 0; i <= noNodes; i++)
         delete [] foot[i];
   if (soilMin != 0)
      for (int i = 0; i <= noNodes; i++)
         delete [] soilMin[i];
   if (soilMax != 0)
      for (int i = 0; i <= noNodes; i++)
         delete [] soilMax[i];
   if (pressure != 0)
      for (int i = 0; i <= noNodes; i++)
         delete [] pressure[i];
   if (pressMax != 0)
      for (int i = 0; i <= noNodes; i++)
         delete [] pressMax[i];


}    


// initializing the bounding surface

void
SoilFootingSection2d::initializeBoundingSurface (void)
{
   a = 0.32;
   b = 0.37;
   ccc = 0.25;
   d = 0.55;   
   eee = 0.8;
   f = 0.8;   

   Fv = V/Vult;  
   A = a * pow(Fv, ccc) * pow(1-Fv, d);
   B = b * pow(Fv, eee) * pow(1-Fv, f);
   
   beta = (A * h) / (pow(A*A*h*h + B*B*L*L, 0.5));
   
   if (beta < 0.0)
      beta *= (-1.0);
   
   alpha = Fv / (1 - beta*(1-Fv));
   pult = alpha;
   
   qult *= alpha;


}


// initializing internal variables

void
SoilFootingSection2d::initializeInternalVariables (void)
{


   foot = new double* [noNodes+1];
   soilMin = new double* [noNodes+1];
   soilMax = new double* [noNodes+1];
   pressure = new double* [noNodes+1];
   pressMax = new double* [noNodes+1];
      
   for (int j = 0; j <= noNodes; j++)
   {
      foot[j] = new double [ini_size];
      soilMin[j] = new double [ini_size];
      soilMax[j] = new double [ini_size];
      pressure[j] = new double [ini_size];
      pressMax[j] = new double [ini_size];
   }


   for (int i = 0; i <= noNodes; i++)
      for (int j = 0; j < ini_size; j++)
      {
         foot[i][j] = V/Kv;
         soilMin[i][j] = V/Kv;
         soilMax[i][j] = V/Kv;
         pressure[i][j] = 1.0/FS;
         pressMax[i][j] = 1.0/FS;
      }


   e.Zero();
   s.Zero();
   eCommit.Zero();
   sCommit.Zero();
   ks.Zero();
   ksE.Zero();

   ks(0,0) = Kv;
   ks(1,1) = Kh;
   ks(2,2) = Kt;
   ksE = ks;

   dTh = 0.0;
   dThP = 0.0;

   Mlimit = V*L/6.0;
   thetaPlus = Mlimit / (Kv * pow(L, 2.0) / 12.0);
   thetaMinus = -1.0 * Mlimit / (Kv * pow(L, 2.0) / 12.0);

   thetaRange = 2 * thetaPlus;



}


void
SoilFootingSection2d::tempFunction (void)
{

}



//  updating committed forces and displacements and internal variables

int
SoilFootingSection2d::commitState(void)
{


   incr++;  


//         cout<<"d/din = "<<Mmaxpast<<endl; 


         if (fabs(s(2)) > Mmaxpast)
            Mmaxpast = fabs(s(2)); 

         if (Mmaxpast > Melastic)
            isElastic = 0;


      thetaPlusPrev = thetaPlus;
      thetaMinusPrev = thetaMinus;


   double e_2 = e(2);


           if (e(2) > thetaPlus)
           {
              thetaPlus = e_2;
              thetaMinus = thetaPlus - thetaRange;
           }

           if (e(2) < thetaMinus)
           {
              thetaMinus = e_2;
              thetaPlus = thetaMinus + thetaRange;
           }


   HPrevCommit = sCommit(1);
   MPrevCommit = sCommit(2);


   eCommit = e;
   sCommit = s;
   ksE = ks;
   dThetaPrev = dTheta;

   c1Commit = c1;
   c1TCommit = c1T;
   c2Commit = c2;
   c2TCommit = c2T;

   eccCommit = ecc;

   hPrev = hCurr;

   for (int i = 0; i <= noNodes; i++)
      for (int j = 2; j > 0; j--)
      {     
         foot[i][j] = foot[i][j-1];
         soilMin[i][j] = soilMin[i][j-1];
         soilMax[i][j] = soilMax[i][j-1];
         pressure[i][j] = pressure[i][j-1];
         pressMax[i][j] = pressMax[i][j-1];
      }

   tolerance = 0.0000000000001 * (1/FS);


   isOver = 1;
   isdV = 0;
   return 0;
}


// going back to the last committed values

int
SoilFootingSection2d::revertToLastCommit(void)
{
  
   thetaPlus = thetaPlusPrev;
   thetaMinus = thetaMinusPrev;

   e = eCommit;
   s = sCommit;
   ks = ksE;
   dTheta = dThetaPrev;

   c1 = c1Commit;
   c1T = c1TCommit;
   c2 = c2Commit;
   c2T = c2TCommit;
   ecc = eccCommit;

   hCurr = hPrev;

   for (int i = 0; i <= noNodes; i++)
   {
      foot[i][1] = foot[i][2];
      soilMin[i][1] = soilMin[i][2];
      soilMax[i][1] = soilMax[i][2];
      pressure[i][1] = pressure[i][2];
      pressMax[i][1] = pressMax[i][2];
   }

   return 0;
}


// going back to the initial conditions 

int
SoilFootingSection2d::revertToStart(void)
{

   eCommit.Zero();
   sCommit.Zero();

   c1 = 0;
   c1T = 0;
   c2 = noNodes;
   c2T = noNodes;

   c1Commit = c1;
   c1TCommit = c1T;
   c2Commit = c2;
   c2TCommit = c2T;
   eccCommit = ecc;

   dTheta = 0.0;
   dThetaPrev = 0.0;

   for (int i = 0; i <= noNodes; i++)
      for (int j = 0; j < ini_size; j++)
      {
         foot[i][j] = V/Kv;
         soilMin[i][j] = V/Kv;
         soilMax[i][j] = V/Kv;
         pressure[i][j] = 1/FS;
         pressMax[i][j] = 1/FS;
      }


   return 0;
}




SectionForceDeformation* SoilFootingSection2d::getCopy()
{

   SoilFootingSection2d *theCopy =
   new SoilFootingSection2d (this->getTag(), FS, Vult, L, Kv, Kh, Rv, dL);

   return theCopy;
}   


// most important method in the class !!
// this function is called at the beginning of and at the middle of iterations

int
SoilFootingSection2d::setTrialSectionDeformation (const Vector &def)
{

   int temp = 0;
   Vector de(3), ds(3);
   double epsilon = pow(10.0, -20.0);



   e = def;
   de = e - eCommit;



   if (fabs(de(0)) < epsilon) de(0) = 0.0;
   if (fabs(de(1)) < epsilon) de(1) = 0.0;
   if (fabs(de(2)) < epsilon) de(2) = 0.0;



   deModel.Zero();

   dThP = dTh;
   dTh = de(2);

   if (de(0) == 0.0 && de(1) == 0.0 && de(2) == 0.0)
   {

   // give out the previous ks rather than ks_elastic

   }
   else
      temp = applyLoading(de);


   ds = ks * deModel;

   if (fabs(ds(0)) < epsilon) ds(0) = 0.0;
   if (fabs(ds(1)) < epsilon) ds(1) = 0.0;
   if (fabs(ds(2)) < epsilon) ds(2) = 0.0;


   s = sCommit + ds;

   return 0;
}




// applyLoading method

int
SoilFootingSection2d::applyLoading(Vector de)
{

   int c, nn = 0, switch1;
   double area1, area2, area1Prev = 0.0;
   double s_recover = Rv;
   double q_recover = 1.0 / L;
   double LcOverL;



   soilFree = 0.0;

   double ds, du, dTheta1, dss, theta;
   double dVt, dHt, dMt, dHt1;
   double Vinit, dMcal;
   double epsilon = pow(10.0, -20.0);
   char tempKey;
   double detKs;
   double expo = 0;
   double n_load, n_unload;
   double e_2;
 
   if (fabs(de(0)) < epsilon) de(0) = 0.0;
   if (fabs(de(1)) < epsilon) de(1) = 0.0;
   if (fabs(de(2)) < epsilon) de(2) = 0.0;


   du = de(1);
   dTheta = de(2);
   dTheta1 = de(2);
   theta = eCommit(2);
   dss = 0.0;


   c1 = c1Commit;
   c1T = c1TCommit;
   c2 = c2Commit;
   c2T = c2TCommit;
   ecc = eccCommit;

  

 
   double *footTemp = new double[noNodes+1];
   double *soilMinTemp = new double[noNodes+1];
   double *soilMaxTemp = new double[noNodes+1];
   double *pressureTemp = new double[noNodes+1];
   double *pressMaxTemp = new double[noNodes+1];

   double *ddH = new double[200];

   // for shear sliding model

   double deltaIn, delta;
   double a11, b11, c11, gr;
   double Fh1, Fh2, Fm1, Fm2;
   double ptFh1, ptFh2, ptFm1, ptFm2;
   double dist1, dist2, factor;
   double hNew, hNew1, uH, uM, dVt1;
   int ii;


   Vector tempLoad(3);

   // extracting the forces - 
   // should be close enough to the actual applied external loads

   tempLoad = ks * de;

   dVt = tempLoad(0);
   dHt = tempLoad(1);
   dMt = tempLoad(2);


   if (fabs(dVt) < epsilon) dVt = 0.0;
   if (fabs(dHt) < epsilon) dHt = 0.0;
   if (fabs(dMt) < epsilon) dMt = 0.0;


// set-up the state of switch depending on external loading

   switch1 = -1;

   if (dVt == 0.0)
      if (dHt == 0.0)
         switch1 = (dMt == 0) ? 0 : 1; 
      else
         switch1 = (dMt == 0) ? 2 : 3; 
   else
      if (dHt == 0.0)
         switch1 = (dMt == 0) ? 4 : 5; 
      else
         switch1 = (dMt == 0) ? 6 : 7; 



   if (isOver == 0)
      switch1 = 7;


   ds = 0.0;
   deModel.Zero();
   ks.Zero();


      for (int i = 0; i <= noNodes; i++)
      {
         footTemp[i] = foot[i][1];
         soilMinTemp[i] = soilMin[i][1];
         soilMaxTemp[i] = soilMax[i][1];
         pressureTemp[i] = pressure[i][1];
         pressMaxTemp[i] = pressMax[i][1];
      }

      for (int j = 0; j < 190; j++)
         ddH[j] = 0.0;


      Vinit = sCommit(0);




   if (switch1 == 0)
   {
      // ALL ZERO !

      ks(0,0) = epsilon;
      ks(1,1) = epsilon;
      ks(2,2) = epsilon;

      ks(0,0) = Kv;
      ks(1,1) = Kh;
      ks(2,2) = Kt;


      ks(0,0) = 1.0;
      ks(1,1) = 1.0;
      ks(2,2) = 1.0;


      deModel(0) = de(0);
      deModel(1) = de(1);
      deModel(2) = de(2);

      return 0;
   }

   else if (switch1 == 4)
   {
      // ONLY dVt
            
      deModel(0) = dVt / Kv;
      
      ks(0,0) = Kv;
      ks(1,1) = Kh;
      ks(2,2) = Kt;

      return 0;
   }

   else
   {
      // APPLY dVt in any case


   if (dVt == 0.0)
      if (dHt == 0.0)
         switch1 = (dMt == 0) ? 0 : 1; 
      else
         switch1 = (dMt == 0) ? 2 : 3; 

   if (isOver == 0)
      switch1 = 7;

   

         if (isOver != 0)
            isOver++;


     dVt = 0.0;

//     do {
     
      dVt1 = dVt;
      Vinit = sCommit(0)+dVt1;


     if (sCommit(1) != 0.0)
        hCurr = sCommit(2)/sCommit(1);
//     else if (dHt != 0.0)
//        hCurr = dMt / dHt;
     else
        hCurr = hPrev;




         if (hCurr < 0.01 && hCurr >= 0.0)
            hCurr = 0.01;
         if (hCurr > -0.01 && hCurr <= 0.0)
            hCurr = -0.01;



         FS = Vult/Vinit;
         Fv = Vinit/Vult;


     LcOverL = 1.0 / (FS * alpha);
     s_recover = Rv * (1.0 - LcOverL);




         A = a * pow(Fv, ccc) * pow(1-Fv, d);
         B = b * pow(Fv, eee) * pow(1-Fv, f);
         
         beta = (A * hCurr) / (pow(A*A*hCurr*hCurr + B*B*L*L, 0.5));
         
         if (beta < 0.0)
            beta *= (-1.0);
          
         // if no shear
         if (switch1 == 1 || switch1 == 5)
            beta = 1.0;  

         alpha = Fv / (1 - beta*(1-Fv));   
         pult = alpha;

         Mult = (Vinit*L/2.0) * beta * (1.0 - Fv);

         Melastic = Mult / (3.0 * (1.0 - Fv));


            // apply M-Theta model
        

         
            if (dTheta > 0.0)
            {
               c = c1;


               if (c > noNodes)
               c = noNodes;


               do // check for Vinit
               {
				   int i;
                  for (i = 0; i <= noNodes; i++)
                  {
                     ds = (i-c)*(L/noNodes)*tan(dTheta);

                     foot[i][0] = footTemp[i] + ds;

                     if (foot[i][0] >= soilMaxTemp[i])
                        soilMax[i][0] = foot[i][0];
                     else
                        soilMax[i][0] = soilMaxTemp[i];
            
                     if (foot[i][0] >= soilMinTemp[i])
                        soilMin[i][0] = foot[i][0];
                     else
                        soilMin[i][0] = soilMinTemp[i];
                  }
            
                  c1 = c2 = -1;

                  for (i = 0; i <= noNodes; i++)
                     if (foot[i][0] >= soilMax[i][0])
                     {
                        c1 = i;
                        break;
                     }
              
                  for (i = noNodes; i >= 0; i--)
                     if (foot[i][0] >= soilMax[i][0])  
                     {
                        c2 = i;
                        break;
                     }



                  if ((c1 == -1) && (c2 == -1))
                  {

                     if ((c >= 0) && (c <= noNodes))
                        c1 = c2 = c;
                     else
                        c1 = c2 = c1T;
                  }


                  for (i = 0; i <= c1; i++)
                  {
                     soilMin[i][0] = soilMax[i][0]
                                   - (soilMax[i][0] - soilFree)*s_recover;
                     if (foot[i][0] >= soilMin[i][0])
                        soilMin[i][0] = foot[i][0];
                  }
 
                  for (i = 0; i <= noNodes; i++) 
                     if (foot[i][0] >= soilMin[i][0])
                     {
                        c1T = i;
                        break;
                     }
              
                  for (i = noNodes; i >= 0; i--)
                  if (foot[i][0] >= soilMin[i][0])
                  {
                     c2T = i;
                     break;
                  }
                             
                  if (c1 == 0)
                     c1T = 0;
                  if (c2 == noNodes)
                     c2T = noNodes;
                 
                  if (c1T != 0) 
                     for (int i = 0; i <= c1T; i++)
                        pressure[i][0] = 0.0;

                  for (i = c1; i <= c2; i++)
                     pressure[i][0] = pressureTemp[i] +
                               (soilMax[i][0]-soilMin[i][1])*q_recover* Kv/ Vult;


                  expo = ((10.0 - 0.5)/(0.66*noNodes))*c1T + 0.5;

           n_load = 0.5;
           n_unload = 10.0;
                 
           e_2 = eCommit(2);
           n_load = 9.5 * pow((thetaPlus-e_2)/thetaRange, 1.0) + 0.5;
           n_unload = 9.5 * pow((e_2-thetaMinus)/thetaRange, 1.0) + 0.5;



                             
                  if (c1 != c1T)
                     for (int i = c1T; i <= c1; i++)
                        pressure[i][0] = pow((double)(i-c1T)/(c1-c1T), n_unload)
                                       * pressure[c1][0];  


                  if (c2 != c2T)
                     for (int i = c2; i <= c2T; i++)
                        pressure[i][0] = pow((double)(c2T-i)/(c2T-c2), n_load)
                                       * pressure[c2][0];
                 
                  for (i = c2T+1; i <= noNodes; i++)
                     pressure[i][0] = 0.0;

				  int a;

                  for (a = 0; a <= noNodes; a++)
                  {
                     if (pressure[a][0] > pult)
                        pressure[a][0] = pult;
                     if (pressure[a][0] < 0.0)
                        pressure[a][0] = 0.0;
                  }
                             
                  area1 = 0.0;
                  for (a = 0; a <= noNodes; a++)
                     area1 += pressure[a][0];
                  area1 *= (L/noNodes);
           
                  area2 = 0.0;
                  for (a = 0; a <= noNodes; a++)
                     area2 += (pressure[a][0] * (L/noNodes)
                           * (a - noNodes/2)*(L/noNodes));
  
                  for (a = 0; a <= noNodes; a++)
                     if (pressure[a][0] > pressMax[a][1])
                        pressMax[a][0] = pressure[a][0];
                     else
                        pressMax[a][0] = pressMax[a][1];
                
                  ecc = area2/area1;
                  area1 *= cos(theta);
                     
                  if ((area1/L) > (Vinit/Vult))
                     c += 1;
                  else
                     nn = -2; 

                    
                  if (c >= noNodes)
                     nn = -2;
                 if (c <= 0)
                     nn = -2;


                  nn++;
                  if ((area1/L-Vinit/Vult < 0.0) && (area1Prev/L-Vinit/Vult > 0.0))
                  {
                   tempKey = 'y'; 
                    if (tempKey == 'y' || tempKey == 'Y')
                    {

                          tolerance *= 2.0;
                    }
                    else
                    {
                    }  
                  }

                  area1Prev = area1;

               } while ((fabs(area1/L - Vinit/Vult) > tolerance) && (nn != -1));

            }

            else if (dTheta < 0.0)
            {  


               c = c2;


               if (c < 0)
                  c = 0;


               do // check for Vinit
               {  int i;            
                  for (i = 0; i <= noNodes; i++)
                  {
                     ds = (i-c)*(L/noNodes)*tan(dTheta);

        
                     foot[i][0] = footTemp[i] + ds;

                     if (foot[i][0] >= soilMaxTemp[i])
                        soilMax[i][0] = foot[i][0];
                     else
                        soilMax[i][0] = soilMaxTemp[i];
                 
                     if (foot[i][0] >= soilMinTemp[i])
                        soilMin[i][0] = foot[i][0];
                     else
                        soilMin[i][0] = soilMinTemp[i];
                  }
  
                  c1 = c2 = -1;
         
                  for (i = 0; i <= noNodes; i++)
                     if (foot[i][0] >= soilMax[i][0])
                     {  
                        c1 = i;
                        break;
                     }

                  for (i = noNodes; i >= 0; i--)   
                     if (foot[i][0] >= soilMax[i][0])
                     {   
                        c2 = i;
                        break;
                     }  
              
                  if ((c1 == -1) && (c2 == -1))
                  {

                     if ((c >= 0) && (c <= noNodes))
                        c1 = c2 = c;
                     else
                        c1 = c2 = c2T;
                  }  


                  for (i = c2; i <= noNodes; i++)
                  {
                     soilMin[i][0] = soilMax[i][0]
                                   - (soilMax[i][0] - soilFree)*s_recover;
                     if (foot[i][0] >= soilMin[i][0])
                        soilMin[i][0] = foot[i][0];
                  }

                  for (i = 0; i <= noNodes; i++)
                     if (foot[i][0] >= soilMin[i][0])
                     {
                        c1T = i;
                        break;
                     }
               
                  for (i = noNodes; i >= 0; i--)
                     if (foot[i][0] >= soilMin[i][0])
                     {
                        c2T = i;
                        break;
                     }
              
                  if (c1 <= 0)
                     c1T = 0;
                  if (c2 >= noNodes)
                     c2T = noNodes;

                 
                  if (c1T != 0)
                     for (int i = 0; i <= c1T; i++) 
                        pressure[i][0] = 0.0;

                  for (i = c1; i <= c2; i++)
                     pressure[i][0] = pressureTemp[i] +
                                    (soilMax[i][0]-soilMin[i][1])*q_recover* Kv/ Vult;


                  expo = ((10.0 - 0.5)/(0.66*noNodes))*(noNodes - c2T) + 0.5;


           n_load = 0.5;
           n_unload = 10.0;

           e_2 = eCommit(2);            
           n_unload = 9.5 * pow((thetaPlus-e_2)/thetaRange, 1.0) + 0.5;
           n_load = 9.5 * pow((e_2-thetaMinus)/thetaRange, 1.0) + 0.5;



                  if (c2 != c2T)
                     for (int i = c2; i <= c2T; i++)
                        pressure[i][0] = pow((double)(c2T-i)/(c2T-c2), n_unload)
                                         * pressure[c2][0];

                  if (c1 != c1T)
                     for (int i = c1T; i <= c1; i++)
                        pressure[i][0] = pow((double)(i-c1T)/(c1-c1T), n_load)
                                       * pressure[c1][0];  
                          
                  for (i = c2T+1; i <= noNodes; i++)
                     pressure[i][0] = 0.0;	
				   int a;
                  for (a = 0; a <= noNodes; a++)
                  {
                     if (pressure[a][0] > pult)
                        pressure[a][0] = pult;
                     if (pressure[a][0] < 0.0)
                        pressure[a][0] = 0.0;
                  }
                             
                  area1 = 0.0;
                  for (a = 0; a <= noNodes; a++)
                     area1 += pressure[a][0];
                  area1 *= (L/noNodes);
           
                  area2 = 0.0;
                  for (a = 0; a <= noNodes; a++)
                     area2 += (pressure[a][0] * (L/noNodes)
                           * (a - noNodes/2)*(L/noNodes));


                  for (a = 0; a <= noNodes; a++)
                     if (pressure[a][0] > pressMax[a][1])
                        pressMax[a][0] = pressure[a][0];
                     else      
                        pressMax[a][0] = pressMax[a][1];
                       
                  ecc = area2/area1;
                  area1 *= cos(theta);
        
                  if ((area1/L) > (Vinit/Vult))
                     c -= 1;
                  else
                     nn = -2;

               
                  if (c <= 0)
                     nn = -2; 
                  if (c >= noNodes)
                     nn = -2;




                  nn++;
                  if ((area1/L-Vinit/Vult < 0.0) && (area1Prev/L-Vinit/Vult > 0.0))
                  {
                    tempKey = 'y'; 

                     if (tempKey == 'y' || tempKey == 'Y')
                     {

                          tolerance *= 2.0;
                     }
                     else
                     {
                     }   
                  }

                  area1Prev = area1;


               } while ((fabs(area1/L - Vinit/Vult) > tolerance) && (nn != -1));

            }

            else
            { 
               dMcal = 0.0;
            }

            // Update Moment

            M = Vinit * ecc;
            dMcal = M - sCommit(2); 

            if (dTheta == 0.0) 
               dMcal = 0.0;


            dMt = dMcal;
            ds = foot[noNodes/2][0] - foot[noNodes/2][1];


/*
         
            if (((dTheta >= 0.0) && (eCommit(2)+dTheta <= thetaPlus)) ||
                ((dTheta < 0.0) && (eCommit(2)+dTheta >= thetaMinus)))
                   ds = 0.0;
                
*/


        if (isOver == 0)
           ds = 0.0;

        dVt = Kv * (de(0) - ds);
//        dVt = 1.0 * (de(0) - ds);


        if (sCommit(0)+dVt > Vult)
           dVt = 0.99*Vult - sCommit(0);
        if (sCommit(0)+dVt < 0.0)
           dVt = -1.0*sCommit(0) + 0.01*Vult;


//      } while ((fabs(dVt - dVt1) > 1.0) && (isOver != 0));



///////////////////////////////

      Fm2 = (sCommit(2)+dMt)/Vult/L;

      if (Fm2 > B)
      {
            dMt = 0.0;
//            exit(0);
      }


   for (int i = 0; i <= 180; i++)
   {        


      Fm2 = (sCommit(2)+dMt)/Vult/L;


      if (A*A*B*B - Fm2*Fm2*A*A > 0.0)
      {
         hNew = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult);
         dHt = -1.0 * hNew + (2.0*hNew/180.0)*i - sCommit(1);
      } 

      else
      {
/*
         if (isOver != 0)
         {
            hNew = sCommit(1);
            cout <<"TOP\n";
            cout <<"A, B, Fm2, dMt  = "<<A<<" "<<B<<" "<<Fm2<<" "<<dMt<<endl;
            cout <<"A*A*B*B - Fm2*Fm2*A*A = "<<A*A*B*B - Fm2*Fm2*A*A<<endl;
            exit(0);
         }
*/
       }
      
//      dHt = -1.0 * hNew + (2.0*hNew/180.0)*i - sCommit(1);

         Fh1 = sCommit(1)/Vult;
         Fm1 = sCommit(2)/Vult/L;
         Fh2 = (sCommit(1)+dHt)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;


         if (Fh1 != Fh2)
         {
            gr = (Fm2 - Fm1)/(Fh2 - Fh1);
    
            a11 = B*B + A*A*gr*gr;
            b11 = 2.0*gr*A*A * (Fm1 - gr*Fh1);
            c11 = A*A * (Fh1*Fh1*gr*gr + Fm1*Fm1 - 2.0*gr*Fh1*Fm1 - B*B);
      
            factor = b11*b11 - 4.0*a11*c11;
         
            if (factor < 0.0)
               factor = 0.0;


            if (factor >= 0.0)
            {
               ptFh1 = (-b11 + pow(factor, 0.5)) / 2.0 / a11;
               ptFh2 = (-b11 - pow(factor, 0.5)) / 2.0 / a11;
               ptFm1 = gr*(ptFh1 - Fh1) + Fm1;
               ptFm2 = gr*(ptFh2 - Fh1) + Fm1;
            }
            else
            {
          
            }
         }
         else
         {
      
            ptFh1 = fabs(Fh1);
            ptFh2 = -1.0 * fabs(Fh1);  
            if ((A*A*B*B-B*B*Fh1*Fh1)/(A*A) < 0.0)
               ptFm1 = ptFm2 = 0.0;
            else
            {  
               ptFm1 = pow((A*A*B*B-B*B*Fh1*Fh1)/(A*A), 0.5);
               ptFm2 = -1.0 * pow((A*A*B*B-B*B*Fh1*Fh1)/(A*A), 0.5);
            }
         }
   
     
         if (fabs(ptFh1) > 0.25 || fabs(ptFh2) > 0.25 ||
             fabs(ptFm1) > 0.15 || fabs(ptFm2) > 0.15 )
         {
         }


////////////////////////////////////////////////////////////////////////////
// One way to make sure that d <= din is to make Fh1 <= ptFh1 .... and so on
/////////////////////////////////////////////////////////////////////////////
      
         if (ptFh1 > ptFh2)
         {
            if (Fh1 > ptFh1)
               Fh1 = ptFh1;
            if (Fh2 > ptFh1)
               Fh2 = ptFh1;
            if (Fh1 < ptFh2)
               Fh1 = ptFh2;
            if (Fh2 < ptFh2)
               Fh2 = ptFh2;
         }
         else
         {  
            if (Fh1 < ptFh1)
               Fh1 = ptFh1;
            if (Fh2 < ptFh1)
               Fh2 = ptFh1;
            if (Fh1 > ptFh2)
               Fh1 = ptFh2;
            if (Fh2 > ptFh2)
               Fh2 = ptFh2;
         }

         if (ptFm1 > ptFm2)
         {  
            if (Fm1 > ptFm1)
               Fm1 = ptFm1;
            if (Fm2 > ptFm1)
               Fm2 = ptFm1; 
            if (Fm1 < ptFm2)
               Fm1 = ptFm2; 
            if (Fm2 < ptFm2)
               Fm2 = ptFm2;
         }
         else
         {
            if (Fm1 < ptFm1)
               Fm1 = ptFm1;
            if (Fm2 < ptFm1)
               Fm2 = ptFm1;
            if (Fm1 > ptFm2)
               Fm1 = ptFm2;
            if (Fm2 > ptFm2)
               Fm2 = ptFm2;
         }


         deltaIn = pow(pow(ptFm2-ptFm1, 2.0) + pow(ptFh2-ptFh1, 2.0), 0.5);

         deltaIn = fabs(2*A);
      
         dist1 = pow(pow(ptFm1-Fm1, 2.0) + pow(ptFh1-Fh1, 2.0), 0.5);
         dist2 = pow(pow(ptFm1-Fm2, 2.0) + pow(ptFh1-Fh2, 2.0), 0.5);
      
      
         if (dist2 < dist1)
            delta = dist2;
         else
            delta = pow(pow(ptFm2-Fm2, 2.0) + pow(ptFh2-Fh2, 2.0), 0.5);
      
         if (dist1 == dist2)
            delta = 0.0;
         
         if (deltaIn - delta <= 0.0)
         {
//            exit (0); 

            if (deltaIn != delta)
            {
            }

            deltaIn = delta + 0.0001;
         }

/*
         if (dist2 < dist1)
            hNew = ptFm1*L/ptFh1;
         else
            hNew = ptFm2*L/ptFh2;
*/


      if (delta/deltaIn > 0.05)
      {
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0))/((delta/(deltaIn-delta))/0.05);
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0)) * pow(1.0 - delta/deltaIn, 2.0);
//         uM = 0.0;
         uH = de(1) - uM;
         dHt1 = Kh * uH * pow(delta/(deltaIn), 1.0);           
         dHt1 = Kh * uH;           
      }
      else
      {
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0));
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0)) * pow(1.0 - delta/deltaIn, 2.0);
//         uM = 0.0;
         uH = de(1) - uM;
         dHt1 = Kh * uH * pow(delta/(deltaIn), 1.0);           
         dHt1 = Kh * uH;           
      }


/*
      if (delta/deltaIn > 0.02)
      {
         uH = dHt / (Kh * pow(delta/(deltaIn-delta), 1.0));
         ddH[i] = fabs(de(1) - uH);
      }
      else
      {
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0));
         ddH[i] = fabs(de(1) - uM);
      }

*/

      
     ddH[i] = fabs(dHt - dHt1);



   } // for loop for finding dHt ends here


   dHt1 = 100000000000.0;

   for (int j = 0; j <= 180; j++)
   {
      if (ddH[j] <= dHt1)
      {
         dHt1 = ddH[j];
         ii = j;         
      }
   }




//////////////////////////


      if (A*A*B*B - Fm2*Fm2*A*A > 0.0)
      {
         hNew = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult);
         dHt = -1.0 * hNew + (2.0*hNew/180.0)*ii - sCommit(1);
      }

      else
      {
/*
         if (isOver != 0)
         {
            hNew = sCommit(1);
            cout <<"BOT\n";
            exit(0);
         }
*/
      }

//         dHt = -1.0 * hNew + (2.0*hNew/180.0)*ii - sCommit(1);

     if (sCommit(1) != 0.0)
        hCurr = sCommit(2)/sCommit(1);
//     else if (dHt != 0.0)
//        hCurr = dMt / dHt;
     else
        hCurr = hPrev;


         
         if (hCurr < 0.01 && hCurr >= 0.0)
            hCurr = 0.01;
         if (hCurr > -0.01 && hCurr <= 0.0)
            hCurr = -0.01;



         hNew = 0.0;

         int gate = 1;
         dHt = dMt / ((hCurr+hPrev)/2.0);
         dHt1 = dMt / ((hCurr+hPrev)/2.0);


         do {

//            dHt = (dHt + dHt1)/2.0;
//            dHt = dHt/2.0;
//            dHt = dHt1;


        hNew = dHt - dHt1;

/*
        if (dHt1 > dHt)
           dHt += 10.0;
        else if (dHt1 < dHt)
           dHt -= 10.0;
*/

           dHt = (dHt + dHt1)/2.0;

         dHt = dMt / hCurr;


//         cout <<"dHt assumed = "<<dHt<<endl;


         Fh1 = sCommit(1)/Vult;
         Fm1 = sCommit(2)/Vult/L;
         Fh2 = (sCommit(1)+dHt)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;


         dHt1 = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5));
         if (Fh2 > dHt1)
            dHt = dHt1*Vult - sCommit(1);         
         if (Fh2 < -1.0*dHt1)
            dHt = -1.0*dHt1*Vult - sCommit(1);         




         if (Fh2*Fh2/A/A + Fm2*Fm2/B/B - 1.0 > 0.0)
         {

/*
//            if (dHt > 0.0)
           if (Fh2 > 0.0)
               dHt = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult)-sCommit(1);
//            if (dHt < 0.0)
           else if (Fh2 < 0.0)
               dHt = -1.0 * fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult)-sCommit(1);



           if (dMt > 0.0)
           dMt = pow(pow(A*B*hCurr, 2.0)/(A*A*hCurr*hCurr+B*B*L*L), 0.5)*Vult*L - sCommit(2);
           else if (dMt < 0.0)
           dMt = -1.0*pow(pow(A*B*hCurr, 2.0)/(A*A*hCurr*hCurr+B*B*L*L), 0.5)*Vult*L - sCommit(2);
           dHt = dMt/hCurr;


            cout<<"In top\n";
            cout <<"sCommit(1) and (2) = "<<sCommit(1)<<" "<<sCommit(2)<<endl;
            cout <<"dHt and dMt = "<<dHt<<" "<<dMt<<endl;
            cout <<"sCommit(0) and dVt = "<<sCommit(0)<<" "<<dVt<<endl;
            cout <<"Fh2 A Fm2 B = "<<Fh2<<" "<<A<<" "<<Fm2<<" "<<B<<endl;
            exit(0); 
*/

         }


         Fh2 = (sCommit(1)+dHt)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;


      
         if (Fh1 != Fh2)
         {
            gr = (Fm2 - Fm1)/(Fh2 - Fh1);
    
            a11 = B*B + A*A*gr*gr;
            b11 = 2.0*gr*A*A * (Fm1 - gr*Fh1);
            c11 = A*A * (Fh1*Fh1*gr*gr + Fm1*Fm1 - 2.0*gr*Fh1*Fm1 - B*B);
      
            factor = b11*b11 - 4.0*a11*c11;
         
            if (factor < 0.0)
               factor = 0.0;


            if (factor >= 0.0)
            {
               ptFh1 = (-b11 + pow(factor, 0.5)) / 2.0 / a11;
               ptFh2 = (-b11 - pow(factor, 0.5)) / 2.0 / a11;
               ptFm1 = gr*(ptFh1 - Fh1) + Fm1;
               ptFm2 = gr*(ptFh2 - Fh1) + Fm1;
            }
            else
            {
          
            }
         }
         else
         {
      
            ptFh1 = fabs(Fh1);
            ptFh2 = -1.0 * fabs(Fh1);  
            if ((A*A*B*B-B*B*Fh1*Fh1)/(A*A) < 0.0)
               ptFm1 = ptFm2 = 0.0;
            else
            {  
               ptFm1 = pow((A*A*B*B-B*B*Fh1*Fh1)/(A*A), 0.5);
               ptFm2 = -1.0 * pow((A*A*B*B-B*B*Fh1*Fh1)/(A*A), 0.5);
            }
         }
   
     
         if (fabs(ptFh1) > 0.25 || fabs(ptFh2) > 0.25 ||
             fabs(ptFm1) > 0.15 || fabs(ptFm2) > 0.15 )
         {
         }


////////////////////////////////////////////////////////////////////////////
// One way to make sure that d <= din is to make Fh1 <= ptFh1 .... and so on
/////////////////////////////////////////////////////////////////////////////
      
         if (ptFh1 > ptFh2)
         {
            if (Fh1 > ptFh1)
               Fh1 = ptFh1;
            if (Fh2 > ptFh1)
               Fh2 = ptFh1;
            if (Fh1 < ptFh2)
               Fh1 = ptFh2;
            if (Fh2 < ptFh2)
               Fh2 = ptFh2;
         }
         else
         {  
            if (Fh1 < ptFh1)
               Fh1 = ptFh1;
            if (Fh2 < ptFh1)
               Fh2 = ptFh1;
            if (Fh1 > ptFh2)
               Fh1 = ptFh2;
            if (Fh2 > ptFh2)
               Fh2 = ptFh2;
         }

         if (ptFm1 > ptFm2)
         {  
            if (Fm1 > ptFm1)
               Fm1 = ptFm1;
            if (Fm2 > ptFm1)
               Fm2 = ptFm1; 
            if (Fm1 < ptFm2)
               Fm1 = ptFm2; 
            if (Fm2 < ptFm2)
               Fm2 = ptFm2;
         }
         else
         {
            if (Fm1 < ptFm1)
               Fm1 = ptFm1;
            if (Fm2 < ptFm1)
               Fm2 = ptFm1;
            if (Fm1 > ptFm2)
               Fm1 = ptFm2;
            if (Fm2 > ptFm2)
               Fm2 = ptFm2;
         }



         deltaIn = pow(pow(ptFm2-ptFm1, 2.0) + pow(ptFh2-ptFh1, 2.0), 0.5);
         deltaIn = fabs(2.0*A);
//         deltaIn = 0.184;
      
         
         dist1 = pow(pow(ptFm1-Fm1, 2.0) + pow(ptFh1-Fh1, 2.0), 0.5);
         dist2 = pow(pow(ptFm1-Fm2, 2.0) + pow(ptFh1-Fh2, 2.0), 0.5);

      
         if (dist2 < dist1)
         {
//            delta = dist2;
//            delta = (dist1+dist2)/2.0;
            delta = dist1;

         delta = pow(pow(ptFh1-Fh1, 2.0), 0.5);


         }
         else
         {
//            delta = pow(pow(ptFm2-Fm2, 2.0) + pow(ptFh2-Fh2, 2.0), 0.5);
            dist1 = pow(pow(ptFm2-Fm1, 2.0) + pow(ptFh2-Fh1, 2.0), 0.5);
            dist2 = pow(pow(ptFm2-Fm2, 2.0) + pow(ptFh2-Fh2, 2.0), 0.5);
//            delta = (dist1+dist2)/2.0;
            delta = dist1;

            delta = pow(pow(ptFh2-Fh1, 2.0), 0.5);


         }     

//         if (dist1 == dist2)
//            delta = 0.0;
//            delta = deltaIn;
  


/*

         Mult = Vult*L*(Fv/2.0)*(hCurr/L)*(1.0 - Fv) / pow(pow(B/A, 2.0) + pow(hCurr/L, 2.0), 0.5);
         deltaIn = fabs(2.0 * Mult);

         if (dMt >= 0.0)
            delta = fabs(-1.0*Mult - (sCommit(2)+dMt));
         else
            delta = fabs(Mult - (sCommit(2)+dMt));

*/

       
         if (deltaIn - delta <= 0.0)
         {

//           exit(0);

            if (deltaIn != delta)
            {
            }

            deltaIn = delta + 0.0001;
         }

/*
       if (incr > 910)
       {
          cout <<"incr = "<<incr<<" "<<delta<<" "<<deltaIn<<endl;
    //      cout <<"delta = "<<delta<<endl;
      //    cout <<"deltaIn = "<<deltaIn<<endl;
       } 
*/


/*

      if (delta/deltaIn > 0.05)
      {
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0))/((delta/(deltaIn-delta))/0.05);
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0)) * pow(1.0 - delta/deltaIn, 10.0);
//         uM = 0.0;
         uH = de(1) - uM;
         dHt1 = Kh * uH;
//         cout <<delta/deltaIn<<endl;
//         cout <<"TOP\n";
      }

      else
      {
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0));
         uM = (2.0 * L * (dTheta) * (L/hCurr)*pow(B/A, 2.0)) * pow(1.0 - delta/deltaIn, 10.0);
//         uM = 0.0;
         uH = de(1) - uM;
         dHt1 = Kh * uH;
//         cout <<delta/deltaIn<<endl;
//         cout <<"BOT\n";
      }

*/      


//      if (delta/deltaIn < 0.1) 
//         delta = 0.1*deltaIn;

//      if (delta/deltaIn > 0.5) 
//         delta = deltaIn;


      Mmaxpast = delta/deltaIn;
      
      if (delta/deltaIn > 1.0)
         uM = 0.0;
      else
         uM = (L * (dTheta-dMt/Kt) * (L/hCurr)*pow(B/A, 2.0)) 
//                       * (pow(1.0 - delta/deltaIn, 4.0));
                       * (pow(1.0 - delta/deltaIn, 10.0));



//         uM = 0.0; 

      uH = de(1) - uM;
//      dHt1 = Kh * uH * pow((delta+0.0001)/(deltaIn - delta), 1.0);
//      dHt1 = Kh * uH * pow(0.0 + delta/deltaIn, 0.2);
//      dHt1 = Kh * uH * pow(epsilon + delta/deltaIn, 2.0);
//      dHt1 = Kh * uH * pow(delta/deltaIn, 2.0);
//      dHt1 = Kh * uH * pow(delta/(deltaIn-delta), 2.0);
      dHt1 = Kh * uH * (delta/deltaIn);
      dHt1 = Kh * uH;
      dHt1 = Kh* ((double)(c2-c1)/noNodes) * uH;


//      if (uH == 0.0) exit(0);
      

//      cout <<"dMt and dHt = "<<dMt<<" "<<dHt<<endl; 
//      cout <<"delta and deltaIn = "<<delta<<" "<<deltaIn<<endl; 
      

//       if (dHt1 == 0.0) 
//          dHt1 = de(1);  


         Fh2 = (sCommit(1)+dHt1)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;

//         if (Fh2*Fh2/A/A + Fm2*Fm2/B/B - 1.0 >= 0.0)
         if (Fh2*Fh2/A/A + Fm2*Fm2/B/B - 1.0 > 0.01)
         {
//            gate = 1;
/*
            cout <<"dHt1 = "<<dHt1<<endl;
            cout <<"dHt = "<<dHt<<endl;
            cout <<"dMt = "<<dMt<<endl;
            cout <<"Fh2 = "<<Fh2<<endl;
            cout <<"Fm2 = "<<Fm2<<endl;
            cout <<"A = "<<A<<endl;
            cout <<"B = "<<B<<endl;
            cout <<"Fh2*Fh2/A/A + Fm2*Fm2/B/B = "<<Fh2*Fh2/A/A + Fm2*Fm2/B/B<<endl;
*/

//          dHt1 = dMt / ((hCurr+hPrev)/2.0);
//          dHt1 = uH;

/*
            cout <<"delta = "<<delta<<endl;
            cout <<"deltaIn = "<<deltaIn<<endl;
            cout <<"dist1 = "<<dist1<<endl;
            cout <<"dist2 = "<<dist2<<endl;
            cout <<"dHt1 = "<<dHt1<<endl;
            cout <<"dHt = "<<dHt<<endl;
            cout <<"dMt = "<<dMt<<endl;
            cout <<"hCurr = "<<hCurr<<endl;
            cout <<"hPrev = "<<hPrev<<endl;
*/


/*
            if (delta/deltaIn < 0.05)
               dHt1 = uH;
            else
               dHt1 = dHt;
*/

/*
          if (-1.0 * pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5) < Fh2)
             dHt1 = -1.0 * pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult - sCommit(1);
          else if (1.0 * pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5) > Fh2)
             dHt1 = pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult - sCommit(1);
          else
             exit (0);
*/


//           dMt = 0.0;
//           dHt1 = 0.0;



//          gate = 0;
         }
         else
            gate = 101;


//            cout <<"gate = "<<gate<<endl;


/*
       if ((fabs(dHt1) < epsilon) || (fabs(dHt) < epsilon))
       {
          dHt1 = dMt / ((hCurr+hPrev)/2.0);
          gate = 0;
       } 
*/

/*
         cout <<"dHt = "<<dHt<<endl;
         cout <<"dHt1 = "<<dHt1<<endl;
         cout <<"delta = "<<delta<<endl;
*/


/*
        if (dHt1 > dHt)
           dHt += 1.0;
        else if (dHt1 < dHt)
           dHt -= 1.0;
*/


   
       if ((dHt - dHt1 > 0.0) && (hNew < 0.0))
            gate = 0;




/*
         Fh2 = (sCommit(1)+dHt1)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;

         if (Fh2*Fh2/A/A + Fm2*Fm2/B/B - 1.0 > 0.0)
         {
            dHt1 = 0.0;
         }

*/


/*
         Fm2 = (sCommit(2)+dMt)/Vult/L;
         dHt = pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult;

         if ((sCommit(1)+dHt1 > dHt) && (dHt1 > 0.0))
//            dHt1 = dHt - sCommit(1);
//              dHt1 = dMt/hCurr;
              dHt1 = 0.0;  
         if ((sCommit(1)+dHt1 < -1.0*dHt) && (dHt1 < 0.0))
//            dHt1 = -1.0*dHt - sCommit(1); 
//              dHt1 = dMt/hCurr;  
              dHt1 = 0.0;  
*/


//         if (fabs(dHt) > 100000) gate = 0;
       
     
//      } while ((fabs(dHt-dHt1) > 100.0) && (gate));
      } while (0);

//       if (dHt1 != 0.0)
          dHt = dHt1;  




         Fh2 = (sCommit(1)+dHt)/Vult;
         Fm2 = (sCommit(2)+dMt)/Vult/L;




         dHt1 = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5));
         if (Fh2 > dHt1)
            dHt = dHt1*Vult - sCommit(1);         
         if (Fh2 < -1.0*dHt1)
            dHt = -1.0*dHt1*Vult - sCommit(1);         


/*
         dHt1 = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5));
         if (Fh2 > dHt1)
            dHt = (Fv/2.0)*(1.0-Fv)/pow(B*B/A/A+hCurr*hCurr/L/L, 0.5)*Vult - sCommit(1);         
         if (Fh2 < -1.0*dHt1)
            dHt = -1.0*(Fv/2.0)*(1.0-Fv)/pow(B*B/A/A+hCurr*hCurr/L/L, 0.5)*Vult - sCommit(1);         
*/




         if (du == 0.0)
            dHt = 0.0;






         if (Fh2*Fh2/A/A + Fm2*Fm2/B/B - 1.0 > 0.0)
         {

/*

//            if (dHt > 0.0)
            if (Fh2 > 0.0)
               dHt = fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult)-sCommit(1);
//            if (dHt < 0.0)
             else if (Fh2 < 0.0)
               dHt = -1.0 * fabs(pow((A*A*B*B - Fm2*Fm2*A*A)/(B*B), 0.5)*Vult)-sCommit(1);




           if (dMt > 0.0)
           dMt = pow(pow(A*B*hCurr, 2.0)/(A*A*hCurr*hCurr+B*B*L*L), 0.5)*Vult*L - sCommit(2);
           else if (dMt < 0.0)
           dMt = -1.0*pow(pow(A*B*hCurr, 2.0)/(A*A*hCurr*hCurr+B*B*L*L), 0.5)*Vult*L - sCommit(2);
           dHt = dMt/hCurr;


            cout<<"In bottom\n";
            cout <<"sCommit(1) and (2) = "<<sCommit(1)<<" "<<sCommit(2)<<endl;
            cout <<"dHt and dMt = "<<dHt<<" "<<dMt<<endl;
            cout <<"sCommit(0) and dVt = "<<sCommit(0)<<" "<<dVt<<endl;
//            exit(0);

*/ 

         }





////////////////////////



/*
        if (isOver == 0)
           ds = 0.0;

        dVt = Kv * (de(0) - ds);


        if (sCommit(0)+dVt > Vult)
           dVt = 0.99*Vult - sCommit(0);
        if (sCommit(0)+dVt < 0.0)
           dVt = -1.0*sCommit(0) + 0.01*Vult;
*/


        if (isOver == 0)

               for (int i = 0; i <= noNodes; i++)
               {

                  foot[i][0] = de(0);
                  soilMin[i][0] = de(0);
                  soilMax[i][0] = de(0);
                  pressure[i][0] = dVt/Vult;
                  pressMax[i][0] = dVt/Vult;
               }

    }





// setup the stiffness matrix      

    deModel = de;
    ks.Zero();


    if (de(0)-ds != 0.0)
       ks(0,0) = dVt/(de(0)-ds);
    else
       ks(0,0) = Kv;
 
//    if (ks(0,0) == 0.0)
//       ks(0,0) = 1.0;

    if ((dTheta != 0.0) && (de(0)-ds != 0.0))
       ks(0,2) = -1.0 * (ds/dTheta) * (dVt/(de(0)-ds));
    else if (dTheta != 0.0)
       ks(0,2) = -1.0 * Kv * (de(0)/dTheta);
    else
       ks(0,2) = 0.0;


    if (isOver != 0)
    {
       if (uH != 0.0)
          ks(1,1) = dHt/uH;
       else
       {
          ks(1,1) = Kh;
       }
    }
    else
    {
       if (de(1) != 0.0)
          ks(1,1) = dHt/de(1);
       else
       {
          ks(1,1) = Kh;
       }
    }
    


    if ((dTheta != 0.0) && (uH != 0.0))
       ks(1,2) = -1.0 * (uM/uH) * (dHt/dTheta);
    else if (dTheta != 0.0)
       ks(1,2) = -1.0 * Kh * (de(1)/dTheta);
    else
       ks(1,2) = 0.0;



    if (de(2) != 0.0)
       ks(2,2) = dMt/de(2);
    else
       ks(2,2) = Kt;




 // make sure det|ks| is non-zero

   detKs = ks(0,0) * (ks(1,1)*ks(2,2) - ks(1,2)*ks(2,1))
         + ks(0,1) * (ks(2,0)*ks(1,2) - ks(1,0)*ks(2,2))
         + ks(0,2) * (ks(1,0)*ks(2,1) - ks(2,0)*ks(1,1));





   delete [] footTemp;
   delete [] soilMinTemp;
   delete [] soilMaxTemp;
   delete [] pressureTemp;
   delete [] pressMaxTemp;

   delete [] ddH;

   return 0;           
}


const Vector &
SoilFootingSection2d::getSectionDeformation (void)
{
  return e;
}

const Vector &
SoilFootingSection2d::getStressResultant (void)
{
  return s;
}

const Matrix &
SoilFootingSection2d::getSectionTangent(void)
{
  return ks;
}

const Matrix &
SoilFootingSection2d::getSectionFlexibility (void)
{
  static Matrix fs(3,3);
//  invert2by2Matrix(ks, fs);
  
  return fs;
}

const ID&
SoilFootingSection2d::getType ()
{
  return code;
}

int
SoilFootingSection2d::getOrder () const
{
  return 3;
}


const Matrix& 
SoilFootingSection2d::getInitialTangent()
{
  static Matrix IniTan(3,3);
                 
  return IniTan;
}

int
SoilFootingSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
SoilFootingSection2d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
SoilFootingSection2d::Print(OPS_Stream &s, int flag)
{
  s << "YieldSurfaceSection2d, tag: " << this->getTag() << endln;
  s << "\tSection Force:" << sCommit;
  s << "\tSection Defom:" << eCommit;
}

