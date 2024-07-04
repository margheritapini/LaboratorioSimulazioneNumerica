/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.initialize_different_T();
  SYS.block_reset(0);
  
  int t_equilibration = 2000;
  
  double initial_T = 2.5;
  double final_T = 0.5;
  
  int N_simulations = 1;
  double T_increase = (initial_T-final_T)/N_simulations;
  
  for (int k=0; k<= N_simulations; k++) {
  
  	for (int l=0; l<t_equilibration; l++) SYS.step(); // Equilibration of the system

  	for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
   	 for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
   	   SYS.step();
   	   SYS.measure();
   	   if(j%10 == 0){
// 	       SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
   	     nconf++;
   	   }
   	 }
   	 SYS.averages(i+1);
   	 SYS.block_reset(i+1);
 	 	}
 	SYS.finalize();
   //SYS.initialize_properties();
 	SYS.change_temp(T_increase);
  }

  return 0;
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
