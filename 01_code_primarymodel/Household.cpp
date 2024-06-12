#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include "Household.h"
#include "utils.h"

using namespace std;

// Constant parameters
// double maxPCRDetectability = 15.0;
double mGamma = 4.0; //flu specific
double vGamma = 2.0; //flu specific
double shift = 25.6;
double shapeInf = pow(mGamma, 2) / vGamma;
double scaleInf = vGamma / mGamma;

double mIncub_g=1.0; //flu specific
double sdIncub_g=1.2; //flu specific
double shape_incub=pow(mIncub_g,2)/pow(sdIncub_g,2);
double scale_incub=pow(sdIncub_g,2)/mIncub_g;

double shape_gt = 4.0;
double scale_gt = 5.0/4.0;


//----------------------------------------------------------------------
// Infectivity profile
//----------------------------------------------------------------------
double infectivityProfile(
                          double origin,
                          double infectionDate,
                          int infectionStatus,
                          double tinf,
                          int studyPeriod,
                          int display
                          ) {

  // Density
  double out = 0.0;

  if (infectionStatus == 1) { // Symptomatic infector
    if (infectionDate < studyPeriod && infectionDate < tinf) {
      out = dgamma(tinf-infectionDate,shapeInf,scaleInf); 
    }

  } 
  return out;
}


double cumulativeInfectivity(
                             double origin,
                             double infectionDate,
                             int infectionStatus,
                             double tinf,
                             int studyPeriod,
                             int display
                             )
{

  // Cumulative infectivity profile
  double out = 0.0;

  // Cumulative density
  if (infectionStatus == 1) { // Symptomatic infector

    if (infectionDate < studyPeriod && infectionDate < tinf) {
      out = pgamma(tinf-infectionDate,shapeInf,scaleInf); 
    }

  }

  return out;
}

//----------------------------------------------------------------------
// Household class - Methods
//----------------------------------------------------------------------
Household::Household() : m_size(0), m_notInfected(0), m_startFollowUp(-1) {

  m_indid.resize(0);
  m_onsetTime.resize(0);
  m_augmentedOnsetTime.resize(0);
  m_infected.resize(0);
  m_confCase.resize(0);
  m_studyPeriod.resize(0);
  m_age.resize(0);
  m_sex.resize(0);
  m_infTime.resize(0);
  m_cumLambda.resize(0);
  m_instLambda.resize(0);
  m_flha_titer_q.resize(0);
  m_stalk_titer_q.resize(0);
  m_na_titer_q.resize(0);
  m_haihk14_titer_q.resize(0); 
  m_instcommunityinf.resize(0); 
	m_cumcommunityinf.resize(0); 
  m_incubPeriod.resize(0);
  

  m_contactPattern= "FALSE"; 
}



Household::Household(std::string contact_pattern,std::vector<double> instval,std::vector<double>cumval) {

  m_size=0;
  m_notInfected=0;
  m_startFollowUp=-1;

  m_indid.resize(0);
  m_onsetTime.resize(0);
  m_augmentedOnsetTime.resize(0);
  m_infected.resize(0);
  m_confCase.resize(0);
  m_studyPeriod.resize(0);
  m_age.resize(0);
  m_sex.resize(0);
  m_infTime.resize(0);
  m_cumLambda.resize(0);
  m_instLambda.resize(0); 
  m_flha_titer_q.resize(0);
  m_stalk_titer_q.resize(0);
  m_na_titer_q.resize(0);
  m_haihk14_titer_q.resize(0);
  m_incubPeriod.resize(0);
  m_instcommunityinf=instval;
  m_cumcommunityinf=cumval;
  m_contactPattern=contact_pattern; 
}




std::vector<int> Household::getSpInfected(std::vector<int> index) const{
    unsigned i, k;
    std::vector<int> output;
    output.resize(0);

    for(i=0; i<index.size();i++) {
        k = index[i];
        output.push_back(m_infected[k]);
    }

  return output;
}

size_t Household::nSymptomatic() const{
  size_t nSymptomatic = 0;
  for (size_t i=0; i<m_size; i++) {
    if (m_onsetTime[i] < 10000.0) nSymptomatic += 1;
  }

  return nSymptomatic;
}


std::string Household::getindid() const {
  std::string indid;
  for (size_t i=0; i<m_size;i++) {
    indid += std::to_string(m_indid[i]) + " "; 
  }

  return indid;
}


// Add individual to household
void Household::addIndividual(
  int indid, 
  double onsetTime, 
  int isCase, 
  int age, 
  int sex, 
  int startFollowUp,
  int studyPeriod, 
  int flha_titer_q,
  int stalk_titer_q,
  int na_titer_q,
  int haihk14_titer_q
  ) {

  m_size += 1; // Update size

  // Update vector attributes
  m_indid.push_back(indid);
  m_onsetTime.push_back(onsetTime);
  m_augmentedOnsetTime.push_back(onsetTime);
  m_infected.push_back(isCase);
  m_studyPeriod.push_back(studyPeriod);
  m_age.push_back(age);
  m_sex.push_back(sex);
  m_flha_titer_q.push_back(flha_titer_q);
  m_stalk_titer_q.push_back(stalk_titer_q);
  m_na_titer_q.push_back(na_titer_q);
  m_haihk14_titer_q.push_back(haihk14_titer_q);
  m_incubPeriod.push_back(0.0);
  
  m_startFollowUp = startFollowUp;

  // Add confirmed cases
    if (isCase > 0) { 
        m_confCase.push_back(m_size-1);
    } else {
        m_notInfected += 1;
    }

  // Initialize infection time
    m_infTime.push_back(10000.0);
}


// Change infection time or symptom onset
void Household::setOnsetTime(int index, double onsetTime) {
  m_augmentedOnsetTime[index] = onsetTime;
  m_infTime[index] = onsetTime-m_incubPeriod[index];
}

void Household::setInfTime(int index, double infTime) {
  m_infTime[index] = infTime;
  m_incubPeriod[index] = m_augmentedOnsetTime[index]-infTime;
  //m_incubPeriod;[index] = infTime; correcting for possibility of symptom onset being augmented before infection onset if using PCR dates
 // m_incubPeriod[index] = m_augmentedOnset_PCR[index-InfTime] correcting for possibility of symptom onset being augmented before infection onset
}


// Reinitialize object
void Household::newHousehold() {
  m_size = 0;
  m_notInfected = 0;  
  m_startFollowUp = -1;
  m_indid.clear();
  m_onsetTime.clear();
  m_augmentedOnsetTime.clear();
  m_infected.clear();
  m_confCase.clear();
  m_studyPeriod.clear();
  m_age.clear();
  m_sex.clear();
  m_infTime.clear();
  m_flha_titer_q.clear();
  m_stalk_titer_q.clear();
  m_na_titer_q.clear();
  m_haihk14_titer_q.clear();
  m_incubPeriod.clear();

 
}


// Initialize augmented data
void Household::initialOnsetTime(int index, std::mt19937_64& gen, int display) {

    if ( m_infected[index] > 0 ) {   // Infected household members 
      m_augmentedOnsetTime[index] = m_onsetTime[index] + runif(gen, 0.0, 1.0);
    } else {
      m_augmentedOnsetTime[index] = 10000.0;  // Unifected people
    }

    if (display) cout << "Observed onset time: " << m_onsetTime[index] << " - New onset time" << m_augmentedOnsetTime[index] << endl;
}




void Household::initialInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen, int display) {

    double incubPeriod(0.0);
    if ( m_infected[index] == 1 ) {   // Symptomatic case 
      m_incubPeriod[index] = rgamma(gen, shape_incub, scale_incub);
      m_infTime[index] = m_augmentedOnsetTime[index] - m_incubPeriod[index];

    } else {    // Covid-free household members or household members with unknown final outcome
      m_infTime[index] = 10000.0; //m_studyPeriod[index];
    }

    if (display) cout << "Inf time " << m_indid[index] << ": " << m_infTime[index] << endl;
}




// Compute the risk of infection
void Household::compute_infectivity_profiles(int display) {
  
  size_t infector, infectee;
  double t1;

  // Set m_cumLambda and m_instLambda size
  m_cumLambda.resize(m_size);
  m_instLambda.resize(m_size);
  for (infector=0; infector<m_size; infector++) {
    m_cumLambda[infector].resize(m_size);
    m_instLambda[infector].resize(m_size);
  }

  // Compute cumulative and instantaneous person-to-person transmission rates
  for (infector=0; infector<m_size; infector++) {
    for (infectee=0; infectee<m_size; infectee++) {

      if ( m_infected[infectee] == 0 ) {
        t1 = m_studyPeriod[infectee];
      }else{
        t1 = m_infTime[infectee];
      }

      if (display) cout << "Infector: " << infector << " " << m_infTime[infector]  << " " << m_augmentedOnsetTime[infector] << " Infectee: " << infectee << " " << t1 << " " << m_infected[infectee] << endl;

      // Cumulative transmission rate for household contacts
      if (m_infected[infector] > 0 && m_infTime[infector] < t1 && m_infected[infectee] >=0) {

        m_cumLambda[infector][infectee] += cumulativeInfectivity(
          m_augmentedOnsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector], 
          t1, 
          m_studyPeriod[infectee],
          display
        );

      }

      // Instantaneous transmission rate for secondary cases
      if (m_infected[infector] > 0 && m_infTime[infector] < m_infTime[infectee] && m_infected[infectee] > 0) {

        m_instLambda[infector][infectee] += infectivityProfile(
          m_augmentedOnsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector],
          m_infTime[infectee], 
          m_studyPeriod[infectee],
          display
        );
      }

    }
  }

   if (display) {
    // Cumulative transmission rate
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_cumLambda[i][j] << " "; 
      }
      cout << endl;
    }
    cout << endl;

    // Instantaneous transmission rate
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_instLambda[i][j] << " "; 
      }
      cout << endl;
    }
    cout << endl;

    // Infection time 
    for (size_t i=0; i<m_size; i++) {
      cout << m_infTime[i] << " ";
    }
    cout << endl;
  }

}

// Update m_cumLambda and m_instLambda matrices during data augmentation
void Household::update_infectivity_profiles(int ind, int display) {
  
  size_t infector, infectee;
  double t1;

  // Ind is the infector
  for (infectee=0; infectee<m_size; infectee++) {
    // Reinitialize
    m_cumLambda[ind][infectee] = 0;
    m_instLambda[ind][infectee] = 0;

    // Cumulative transmission rate for all household contacts
    if ( m_infected[infectee] == 0 ) {
      t1 = m_studyPeriod[infectee];
    }else{
      t1 = m_infTime[infectee];
    }

    if (m_infected[ind] > 0 && m_infTime[ind] < t1 && m_infected[infectee] >=0) {

      m_cumLambda[ind][infectee] = cumulativeInfectivity(
      m_augmentedOnsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind],
      t1, 
      m_studyPeriod[infectee]
      );
    }

    // Instantaneous transmission rate for secondary cases
    if (m_infected[ind] > 0 && m_infTime[ind] < m_infTime[infectee] && m_infected[infectee] > 0) {
      m_instLambda[ind][infectee] = infectivityProfile(
      m_augmentedOnsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind], 
      m_infTime[infectee], 
      m_studyPeriod[infectee]
      );
    }
  }


  // Ind is an infectee
  if (m_infected[ind] == 0) {
    t1 = m_studyPeriod[ind];
  } else {
    t1 = m_infTime[ind];    
  }

  for (infector=0; infector<m_size;infector++) {    
    // Reinitialize
    m_cumLambda[infector][ind] = 0;
    m_instLambda[infector][ind] = 0;

    // Cumulative transmission rate for all household contacts
    if (m_infected[infector] > 0 && m_infTime[infector] < t1 && m_infected[ind] >=0) {

      m_cumLambda[infector][ind] = cumulativeInfectivity(
      m_augmentedOnsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      t1, 
      m_studyPeriod[ind]
      );
    }

    // Instantaneous transmission rate for secondary cases
    if (m_infected[infector] > 0 && m_infTime[infector] < m_infTime[ind] && m_infected[ind] > 0) {
      m_instLambda[infector][ind] = infectivityProfile(
      m_augmentedOnsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      m_infTime[ind], 
      m_studyPeriod[ind]
      );
    }
  }

  if (display) {
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_cumLambda[i][j] << " "; 
      }
      cout << endl;
    }
    
    for (size_t i=0; i<m_size; i++) {
      cout << m_infTime[i] << " ";
    }
    cout << endl;
  }

}



// Propose new value for infection time and symptom onset time
double Household::newOnsetTime(int index, std::mt19937_64& gen) {

    double onsetTime(0.0);

    if ( m_infected[index] > 0 ) {           // Infected household member 
      onsetTime = m_onsetTime[index] + runif(gen, 0.0, 1.0);

    } else {    // Flu-free household members or household members with unknown final outcome
      onsetTime = 10000.0;
    }

    return onsetTime;
}



double Household::newInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen) {

    double infDate(0.0), incubPeriod(0.0);

    if ( m_infected[index] == 1 ) {           // Symptomatic case 

      incubPeriod = rgamma(gen, shape_incub, scale_incub);
      //incubPeriod = rlnorm(gen, mIncub, sdIncub);
      infDate += m_augmentedOnsetTime[index] - incubPeriod;

   
    } else {    // Covid-free household members or household members with unknown final outcome
      infDate = 10000.0;
    }

    return infDate;
}





double Household::log_dIncub(int index, double maxPCRDetectability, int display) {

    double out = 0.0;

    if ( m_infected[index] == 1 ) {       // Symptomatic case with known symptom onset 
      double incubPeriod = m_augmentedOnsetTime[index] - m_infTime[index];
      out += logdlnorm(incubPeriod, mIncub, sdIncub);

    } 

    if (display) cout << "log_dIncub: " << out << endl;

    return out;
}



double Household::log_dIncub_g(int index, double maxPCRDetectability, int display) {

    double out = 0.0;

    if ( m_infected[index] == 1 ) {       // Symptomatic case with known symptom onset 
      double incubPeriod = m_augmentedOnsetTime[index] - m_infTime[index];

      out += log(dgamma(incubPeriod, shape_incub, scale_incub));

    } 
    if (display) cout << "log_dIncub: " << out << endl;

    return out;
}







double Household::compute_log_lik(
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  double maxPCRDetectability,
  int display
  ) {

  // Identify index case
  size_t firstCase(0);
  auto it = std::min_element(m_infTime.begin(), m_infTime.end());
  double t0 = *it;
  for (size_t i =0; i < m_size; i++) {
    if ( m_infTime[i] == t0 ) firstCase = i;
  }

  // Log Likelihood
  double LL = 0.0;

  for (size_t i = 0; i < m_size; i++) { // All individuals

    if (display) {
      cout << i << " " << firstCase << " " << m_augmentedOnsetTime[i] << " " << m_infTime[i] << " " << m_infected[i] << endl; 
    }

    if ( i == firstCase ) {                                           		// Incubation period of 1st infected
      LL += log_dIncub_g(i, maxPCRDetectability, display); 

    } else if ( m_infected[i] > 0 && i != firstCase) {                 // Contribution of the secondary cases
      LL += log_S(i, t0, parameter, selectedParam, display);
    	LL += log_pInf(i, parameter, selectedParam, display);
    	LL += log_dIncub_g(i, maxPCRDetectability, display);

    } else {																// Non infected individuals tha contribute to the likelihood
      if (m_infected[i] ==0) LL += log_S(i, t0, parameter, selectedParam, display); 

    }
    if (display) cout << "end individual" << endl;

  }

  if (display) cout << "LL: " << LL << "\n\n";

  return LL;
}



// Probability of infection at t1
double Household::log_pInf(
  size_t curr, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  int display
  ){

    int tinf=(int) m_infTime[curr];

    //cout << m_instcommunityinf[tinf] << endl;
    // Instaneous risk of infection from community
    double alpha = parameter[0];

    // Instantaneous risk of infection from other household members
    double beta_i(0.0), beta(0.0);
    for (size_t ind=0; ind<m_size; ++ind) {
      if (curr != ind) {
        beta_i = m_instLambda[ind][curr];

        if (display) cout << "pInf: " << curr << " " << ind << " " << beta_i << endl;
        
    // Infectivity of children to do infectivity, use [ind] instead of [curr] and beta_i instead of beta
    //if (m_age[ind]==0 && selectedParam[3] ) beta_i *= parameter[3];

    //infectivity by NA titer
    if (m_na_titer_q[ind] == 1 || m_na_titer_q[ind] == 2 || m_na_titer_q[ind] == 3 && selectedParam[7] ) beta_i *= parameter[7];

    //infectivity by HI titer
    if (m_haihk14_titer_q[ind] == 1 || m_haihk14_titer_q[ind] == 2 || m_haihk14_titer_q[ind] == 3 && selectedParam[8] ) beta_i *= parameter[8];    

    //infectivity by Stalk titer
    if (m_stalk_titer_q[ind] == 1 || m_stalk_titer_q[ind] == 2 || m_stalk_titer_q[ind] == 3 && selectedParam[9] ) beta_i *= parameter[9];  


        beta += beta_i;
      }
    }


    // Instantaneous per capita transmission rate 
    beta *= parameter[1];

    // If transmission hazard depends on household size
    if (selectedParam[2] == 1) beta /= pow( m_size / 5.0, parameter[2] );

    beta=beta+alpha;

    // Suceptibility of children
    if (m_age[curr]==0 && selectedParam[3] ) beta *= parameter[3];




    // Susceptibility by stalk titer

    if (m_stalk_titer_q[curr] == 1 || m_stalk_titer_q[curr] == 2 || m_stalk_titer_q[curr] == 3 && selectedParam[4] ) beta *= parameter[4];  

    // Susceptibility by NA titer
    if (m_na_titer_q[curr] == 1 || m_na_titer_q[curr] == 2 || m_na_titer_q[curr] == 3 && selectedParam[5] ) beta *= parameter[5];

    // Susceptibility by HAI (HK14) titer

     if (m_haihk14_titer_q[curr] == 1 || m_haihk14_titer_q[curr] == 2 || m_haihk14_titer_q[curr] == 3 && selectedParam[6] ) beta *= parameter[6];

    if (display) cout << "log_pInf: " << log( (alpha+beta) ) << endl;

    return log( ( beta) );
}




// Survival from t0 to t1
double Household::log_S(
  size_t curr, 
  double t0, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  int display
  ) {

  // If the individual is not infected, it's infection
  double t1;
  if ( m_infected[curr] == 0 ) {
  	t1 = m_studyPeriod[curr];
  }else{
  	t1 = m_infTime[curr];
  }

  int tb = (int) t0;
  int te = (int) t1;
  
  // Cumulative risk of infection from community
  double alpha = parameter[0] * (t1 - t0);

  // Cumulative risk of infection from other household members
  double beta(0.0), beta_i(0.0);
	for (size_t ind=0; ind<m_size; ++ind) {
    if (curr != ind) {
      
      beta_i = m_cumLambda[ind][curr];

      if (display) cout << "pInf: " << curr << " " << ind << " " << beta_i << endl;

    
    // Infectivity of children to do infectivity, use [ind] instead of [curr] and beta_i instead of beta
    //if (m_age[ind]==0 && selectedParam[3] ) beta_i *= parameter[3];

    //infectivity by NA titer
    if (m_na_titer_q[ind] == 1 || m_na_titer_q[ind] == 2 || m_na_titer_q[ind] == 3 && selectedParam[7] ) beta_i *= parameter[7];

    //infectivity by HI titer
    if (m_haihk14_titer_q[ind] == 1 || m_haihk14_titer_q[ind] == 2 || m_haihk14_titer_q[ind] == 3 && selectedParam[8] ) beta_i *= parameter[8];    

    //infectivity by Stalk titer
    if (m_stalk_titer_q[ind] == 1 || m_stalk_titer_q[ind] == 2 || m_stalk_titer_q[ind] == 3 && selectedParam[9] ) beta_i *= parameter[9];  




      beta += beta_i;
    }
	}

  
    // Instantaneous per capita transmission rate 
    beta *= parameter[1];

    // If transmission hazard depends on household size
    if (selectedParam[2] == 1) beta /= pow( m_size / 5.0, parameter[2] );

    beta=beta+alpha;

    // Relative susceptibility of children
     if (m_age[curr]==0 && selectedParam[3] ) beta *= parameter[3];
 

    // Susceptibility by stalk titer

    if (m_stalk_titer_q[curr] == 1 || m_stalk_titer_q[curr] == 2 || m_stalk_titer_q[curr] == 3 && selectedParam[4] ) beta *= parameter[4];  

    // Susceptibility by NA titer
    if (m_na_titer_q[curr] == 1 || m_na_titer_q[curr] == 2 || m_na_titer_q[curr] == 3 && selectedParam[5] ) beta *= parameter[5];

    // Susceptibility by HAI (HK14) titer

     if (m_haihk14_titer_q[curr] == 1 || m_haihk14_titer_q[curr] == 2 || m_haihk14_titer_q[curr] == 3 && selectedParam[6] ) beta *= parameter[6];

  if (display) cout << beta << " " << alpha << endl;
  
  if (display) cout << "log_S: " << -(alpha+beta) << endl;

  return -( beta );
}

