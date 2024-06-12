#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double mIncub = 1; //log normal
double sdIncub = 1.2;
double maxPCRDetectability = 10.0;

double m_incub_g = 1; //gamma
double sd_incub_g = 1.2;
double shape_incub_g = pow(m_incub_g,2) / pow(sd_incub_g, 2);
double scale_incub_g = pow(sd_incub_g,2) / m_incub_g;

double mGamma = 4;
double vGamma = 2;
double shift = 25.6;
double shape_infectivity = pow(mGamma,2) / vGamma;
double scale_infectivity = vGamma / mGamma;

double shape_gt = 4.0;
double scale_gt = 5.0/4.0;

double mu_f = -4.0;
double sigma_f = 1.85; 
double alpha_f = 5.85;
double tau = exp(mIncub+0.5*pow(0.5, 2)); 

double cf = alpha_f / (sigma_f*(1- pow(1+exp((mu_f+tau)/sigma_f), -alpha_f)));
  

//////////////////////////////////////////////
// [[Rcpp::export]]
void test() {
	Rcout << "Mean incubation period: " << tau << endl;
	Rcout << "Cumulative tost: " << cf * sigma_f / alpha_f * pow(1+exp(-(40-mu_f)/sigma_f), -alpha_f) << endl;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double incubPeriod_g() {
	double d = rgamma(1, shape_incub_g, scale_incub_g)[0];
	//while(d<1) d = rgamma(1, shape_incub_g, scale_incub_g)[0];
	return d;
}



//////////////////////////////////////////////
// [[Rcpp::export]]
double cumuldensityPeriod_g(double x) {
  double d = R::pgamma(x, shape_infectivity, scale_infectivity,1,0);
  return d;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double incubPeriod() 
  {
  return rgamma(1, shape_incub_g, scale_incub_g)[0];;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector foi(
	double t, 
	double dt, 
	double lastDate,  


	NumericVector incub_periods,
	NumericVector symptomOnset, 
	NumericVector infection_onset,
	IntegerVector age, 
	IntegerVector haihk14_titer_q,
	IntegerVector na_titer_q,
	IntegerVector stalk_titer_q,
	IntegerVector infectionStatus, 
	
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

	double hhsize,
  double mainHHSize,
	int display
  ) {
  
  //Rcout  << alpha<< beta<< delta<<muCC<< muCM<<muCW<<hhsize<<mainHHSize<<"\n";
	// Force of infection from 
	// 	0: community
	// 	1 to #infected: infected household contacts
	NumericVector fois(symptomOnset.size() + 1);

	//Infection by the community
	fois[0] = alpha*dt;

  	// Infection by infected individuals within the same household
	for (int index = 0; index < symptomOnset.size(); ++index) {

		double k = 0.0;

		if ( infectionStatus[index] == 1 ) { // Symptomatic infector
			if ( t >= infection_onset[index] && infection_onset[index] < lastDate ) {
				k+= R::pgamma(t+dt-infection_onset[index], shape_infectivity, scale_infectivity,1,0) - R::pgamma(t-infection_onset[index], shape_infectivity, scale_infectivity,1,0);
			} 

		} 
		// Foi depends on household size - non parametric formulation (reference = father-mother)
		fois[index + 1] = beta * k / pow(hhsize / mainHHSize, delta);

		// Contact rate between infector - infectee pairs

		// Child-Mother
	//	if ( (age[index] == 0 && ageCurr == 1 && sexCurr ==1) || (age[index] == 1 && ageCurr == 0 && sex[index] ==1) ) fois[index + 1] *= muCM;

		// Child-Father
	//	if ( (age[index] == 0 && ageCurr == 1 && sexCurr ==2) || (age[index] == 1 && ageCurr == 0 && sex[index] ==2) ) fois[index + 1] *= muCW;
		
          // Child-child
	//	if ( (age[index] == 0 && ageCurr == 0)) fois[index + 1] *= muCC;
			
  //infectivity parameters?
  
			
			//HAIHK14 influence
			if (haihk14_titer_q[index]==1 || haihk14_titer_q[index]==2 || haihk14_titer_q[index]==3){
			  fois[index + 1] *= muiHIUpper;
			}
			
			//NA influence
			if (na_titer_q[index]==1 || na_titer_q[index]==2 || haihk14_titer_q[index]==3){
			  fois[index + 1] *= muiNAUpper;
			}
			
			//Stalk influence
			if (stalk_titer_q[index]==1 || stalk_titer_q[index]==2 || stalk_titer_q[index]==3){
			  fois[index + 1] *= muiStalkUpper;
			}
			
			//if (severity[index]==1){
			  //fois[index + 1] *= muIS;
			//}
			
			

	} 
	//Rcout << "FOIs  : " << fois << "\n";
  return fois;
}
