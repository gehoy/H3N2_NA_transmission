#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame hhEpidemic(DataFrame H,
                     double alpha,
                     double beta,
                     double delta,
                     double muchild,
                     double musHIUpper,
                     double musNAUpper,
                     double musStalkUpper,
                     double muiHIUpper,
                     double muiNAUpper,
                     double muiStalkUpper,
                     double dt) //used to be double followUp
  {
  

  // Simulate epidemic within a household

 

  // Data info
  int hhsize = H.nrows();

  IntegerVector age = H["age"];
  IntegerVector sex = H["sex"];
  NumericVector symptomOnset = H["symptomOnset"];
  NumericVector studyPeriod = H["studyPeriod"];
  IntegerVector flha_titer_q = H["flha_titer_q"];
  IntegerVector stalk_titer_q = H["stalk_titer_q"];
  IntegerVector na_titer_q = H["na_titer_q"];
  IntegerVector haihk14_titer_q = H["haihk14_titer_q"];
  int mainHHSize=5;

  // Initialize transmission vectors
  NumericVector infection_onset = rep(10000.0, hhsize); //infection time
  IntegerVector infectionStatus = H["infectionStatus"];
  NumericVector incub_period = rep(10000.0, hhsize);


  // Susceptible and index cases at the start of the epidemic
  IntegerVector infectors(0);
  IntegerVector sus(0);
  double lastDate = 0;
  double firstInfection = 10000.0;
  double startFollowUp = 10000.0;
  

  for (int ind=0; ind<hhsize; ind++) {
    

  	if (infectionStatus[ind] >= 1) { // Index cases

  	infectors.push_back(ind);

    	if ( infectionStatus[ind] == 1 ) { // Symptomatic cases 
  
        incub_period[ind] = incubPeriod_g();

    	} else { 
    
        incub_period[ind] = runif(1, 2.0, 7.0)[0];        

    	}
    	symptomOnset[ind]=symptomOnset[ind]+ runif(1)[0];
      infection_onset[ind] = symptomOnset[ind] - incub_period[ind];
      startFollowUp = min(startFollowUp, symptomOnset[ind]);
    	firstInfection = min(firstInfection, infection_onset[ind]);

    } else { // Susceptible household members
    	sus.push_back(ind);
    }
  }

  lastDate = studyPeriod[0]; //all household members have identical follow-up time

  // Simulate Epidemic within household
  for (int n = 0; n <= (lastDate-firstInfection)/dt; ++n) { 

   // Current time
    double t = firstInfection + dt * n;  
    

    
    // Draw new infections
    IntegerVector newInfected;
    NumericVector infected = runif(sus.size());
    // NumericVector severe_infected = runif(sus.size());
    
    // Loop on susceptibles 
    for (int s = 0; s < sus.size(); ++s) {
      
      // Individual
      int ind = sus[s];

      // Probability of getting infected
      int display=0;

      //START UPDATING HERE 11/3
      
      NumericVector FOIS = foi(
        t, 
        dt, 
        lastDate,
     
        incub_period[infectors],
        symptomOnset[infectors], 
        infection_onset[infectors],
                       
        age[infectors], 
        haihk14_titer_q[infectors],
        na_titer_q[infectors],
        stalk_titer_q[infectors], 
        infectionStatus[infectors], 
       
       
        alpha, 
        beta, 
        delta,
        muchild,
        musHIUpper,
        musNAUpper,
        musStalkUpper,
        muiHIUpper,
        muiNAUpper,
        muiStalkUpper,
        hhsize,
        mainHHSize,
        display
        );
      
      // Update foi according to age of contact 
      NumericVector fois = clone(FOIS);

      
      double totFoi=sum(fois);
      
      //Rcout << "hhsize " << hhsize << "\n";
      //Rcout << "FOIs  : " << totFoi << "\n";
      
      //age
      
      if (age[ind]==0) totFoi *= muchild;
      
      //sex
      
     // if (sex[ind]==1) totFoi *= musex;

      // Susceptiblity of HAI
      if (haihk14_titer_q[ind]==1 || haihk14_titer_q[ind]==2 || haihk14_titer_q[ind]==3) totFoi *= musHIUpper;
      
      // Susceptiblity of NA
      if (na_titer_q[ind]==1 || na_titer_q[ind]==2 || na_titer_q[ind]==3) totFoi *= musNAUpper;

   
      // Susceptibility of stalk
      
      if (stalk_titer_q[ind]==1 || stalk_titer_q[ind]==2 || stalk_titer_q[ind]==3) totFoi *= musStalkUpper;

      
      


      
      
      // Is ind infected at curr_time?
      double pInfection = 1-exp( -totFoi );

     // double pInfection_index = fois[1]/totFoi;

      if ( infected[s] < pInfection) {
        // Update ddi and dds upon infection
        infection_onset[ind] = t + dt*runif(1)[0]; //ddi infection_onset

        

        infectionStatus[ind] = 1;
        incub_period[ind] = incubPeriod_g();
        symptomOnset[ind] = infection_onset[ind] + incub_period[ind];
        
        
    //    if (severe_infected[s]<0.5){
      //    severity[ind]=1;
      //  }

        // Update new infected individuals
        newInfected.push_back(ind);
      }
    }
    
    // Update sus and infectors
    for (int index = 0; index < newInfected.size(); ++index) {
      
      infectors.push_back(newInfected[index]);

      for (int i = 0; i < sus.size(); ++i) {
        if (newInfected[index] == sus[i]) {
          sus.erase(i);
        }
      }
    }
  }
  
  // Output
  return DataFrame::create(
  	_("indid") = H["indid"],
  	_("hhid") = H["hhid"],
  	_("hhsize") = H["hhsize"],
  	_("symptomOnset") = symptomOnset,
  	_("infectionStatus") = infectionStatus,
  	_("studyPeriod") = rep(lastDate, hhsize),
  	_("age") = age,
  	_("sex") = sex,
  	_("flha_titer_q") = flha_titer_q,
  	_("haihk14_titer_q") = haihk14_titer_q,
  	_("na_titer_q") = na_titer_q,
  	_("stalk_titer_q") = stalk_titer_q
  	);
}




