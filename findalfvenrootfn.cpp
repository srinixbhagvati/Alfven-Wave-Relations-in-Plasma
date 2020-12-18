#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include<fstream>
using namespace std;

double findcubicroots(double kyL, double kzL, double beta, double mr, double z[3][2]){
	

	const double PI = 4.0 * atan( 1.0 );
	double a, b, c, d;
	double kL = sqrt(kyL*kyL + kzL*kzL);

	a = ((kzL*kzL)/(kL*kL))*(1 + (kyL*kyL)/(mr*beta));
	b = -(((kzL*kzL)/(kL*kL))*(1 + kyL*kyL) + (1 + beta + (kyL*kyL)/(mr*beta)));
	c = 1 + 2*beta + kyL*kyL;
	d = - beta;

   if (a==0) {
   	cout << "The equation can be reduced to an equation whose degree is less than three\n";
   	return 0.0;
   }
   
   if (d == 0){
   	cout << "The equation can be reduced to an equation whose degree is less than three\n";
   	return 0.0;
   }
	// Reduced equation: X^3 - 3pX - 2q = 0, where X = x-b/(3a)
   double p = ( b * b - 3.0 * a * c ) / ( 9.0 * a * a );
   double q = ( 9.0 * a * b * c - 27.0 * a * a * d - 2.0 * b * b * b ) / ( 54.0 * a * a * a );
   double offset = b / ( 3.0 * a );


   // Discriminant
   double discriminant = p * p * p - q * q;


   if ( discriminant > 0 )           // set X = 2 sqrt(p) cos(theta) and compare 4 cos^3(theta)-3 cos(theta) = cos(3 theta)
   {
      double theta = acos( q / ( p * sqrt( p ) ) );
      double r = 2.0 * sqrt( p );
      for ( int n = 0; n < 3; n++ ) { 
        //cout << r * cos( ( theta + 2.0 * n * PI ) / 3.0 ) - offset << '\n';
      	z[n][0] = r * cos( ( theta + 2.0 * n * PI ) / 3.0 ) - offset;
      	z[n][1] = 0.0;

      } 
      return 1; 
   }
   else 
   {
      double gamma1 = cbrt( q + sqrt( -discriminant ) );
      double gamma2 = cbrt( q - sqrt( -discriminant ) );

      z[0][0] = gamma1 + gamma2 - offset;
      z[0][1] = 0.0;

      double re = -0.5 * ( gamma1 + gamma2 ) - offset;
      double im = ( gamma1 - gamma2 ) * sqrt( 3.0 ) / 2.0;
      if ( discriminant == 0.0 )                // Equal roots (hmmm, floating point ...)
      {
     
         z[1][0] = re; z[1][1] = 0.0;

         z[2][0] = re; z[2][1] = 0.0;
      }
      else
      {
      
         z[1][0] = re; z[1][1] = im;

        
         z[2][0] = re; z[2][1] = -im;
      }
      return 1;
      
   }
}

double findalfvenroot(double kyL, double kzL, double beta, double mr, double sqroots[3][2]){
	double dum;
	dum = findcubicroots(kyL, kzL, beta, mr, sqroots);
   	if (dum != 0){

   		for(int n=0; n<3; n++){
   			if (sqroots[n][1] != 0){
   				cout << "Two roots of Omega squared are imaginary. This code cannot proceed further. Please choose appropriate input values\n";
	   			cout << "Roots for Omega Squared obtained are:\n";
	   			for(int n=0; n<3; n++){
	   				cout << sqroots[n][0] << " + " << sqroots[n][1] << " i\n";
  	 			}
   				return 0;
   			}
  		} 
   		double roots[3];
		

		for(int n=0; n<3; n++){
  			roots[n] = sqrt(sqroots[n][0]);		//ONLY POSITIVE SQUARE ROOTS CONSIDERED
          
  		}
   	
   		sort(roots, roots + 3);
  		//cout << "Positive roots for Omega are: \n";
  		//cout << roots[0] << ", " << roots[1] << ", " << roots[2] << "\n";
 		
  		return roots[1];
  	}
}

int main() {

   double kyL, kzL, beta, mr, dum, alfroot ;
   int len = 0, num;
   double sqroots[3][2];
   
   double min, max, kyL_arr[100], alfven_arr[100];
   cout << "Give minimum, maximum and the total number of values (less than 100) for x-axis kyL separated by spaces: \n";
   cin >> min >> max >> num;
  

   cout << "Give values for kzL, beta and mr separated by spaces: \n";
   cin >> kzL >> beta >> mr;

   if(num ==0){
      cout << "Invalid input data. Give at least one value for kyL.\n";
      return 0;
   }

   for (int i =0; kyL_arr[i-1] < max; i++){
      if(num == 1){
        
        cout << "Since only one value can be taken, we assume kyL = " << min << "\n";
        kyL_arr[i] = min;
      }

      if (num > 1){
        kyL_arr[i] = min + ((max - min)/(num-1))*i;
        if (kyL_arr[i] > max) break;
      }
      
      alfroot = findalfvenroot(kyL_arr[i], kzL, beta, mr, sqroots);
      alfven_arr[i] = alfroot;
      len++;
   }
   
   
   cout << "Array of Alfven roots is\n";
   for (int i = 0; i < len; i++){

    cout << "For kyL = " << kyL_arr[i] << " root is " << alfven_arr[i] << "\n";
   }
	 
   char yncsv;
   cout << "Do you want to export kyL and root values? [y/n]\n";
   //cin >> yncsv;
   yncsv = 'y';
   if (yncsv == 'y'){
    fstream fout; // opens an existing csv file or creates a new file. 
    fout.open("outputvals.csv", ios::out | ios::trunc); 
    fout.close();
    fout.open("outputvals.csv", ios::out | ios::app);
    fout << "For beta = " << beta << ", kzL  = " << kzL << " and mass ratio = " << mr << "\n";
    for (int i =0 ; i < len; i++){
      fout << kyL_arr[i] << ", " << alfven_arr[i] << "\n";
    }
    cout << "Done. CSV file exported as outputvals.csv\n";
   }


}

//ADDITIONAL NOTES: WHY ARE WE CONSIDERING ONLY POSITIVE SQUARE ROOTS OF OMEGA SQUARED? 