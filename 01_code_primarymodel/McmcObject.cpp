#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "McmcObject.h"
#include "Household.h"
#include "utils.h"

using namespace std;

//----------------------------------------------------------------------
// MCMC class
//----------------------------------------------------------------------
// Constructor
McmcObject::McmcObject() {
    m_gen.seed(123);
    m_globalLogLik = 0.0;
    m_iterations = 0;
    m_iterTimeInfection = 0;
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
    m_rateForRandomWalk.resize(0);
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_accept.resize(0);
    m_parameter.resize(0);
    m_selectedParam.resize(0);
    m_hhLogLik.resize(0);
    std::vector<Household> m_data;
    m_data.resize(0);
    m_nHH = 0;
    m_lnormPrior.resize(0);
    m_sdLNormPrior.resize(0);
    m_maxPCRDetectability = 10.0;
    m_adaptive_steps = 0.0;
}

McmcObject::McmcObject(
                       unsigned seed,
                       int nIterations,
                       int nAdaptiveSteps,
                       std::vector<Household> data,
                       std::vector<double> parameter,
                       std::vector<int> selectedParameter,
                       std::vector<double> rateForRandomWalk,
                       std::vector<int> lnormPrior,
                       std::vector<double> sdLNormPrior,
                       int nIterTimeInfection,  
                       double maxPCRDetectability
                       ) {
    // Default values
    m_globalLogLik= 0.0;
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_accept.resize(parameter.size());
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;

    // Copy input vectors
    m_gen.seed(seed);
    m_hhLogLik.resize(data.size());
    m_iterations = nIterations;
    m_adaptive_steps = nAdaptiveSteps;
    m_iterTimeInfection = nIterTimeInfection;
    m_rateForRandomWalk = rateForRandomWalk;
    m_lnormPrior = lnormPrior;
    m_sdLNormPrior = sdLNormPrior;
    m_maxPCRDetectability = maxPCRDetectability;
    m_parameter = parameter;
    m_selectedParam = selectedParameter;
    m_data = data;
    m_nHH = data.size();
}

// Initialize parameter values
void McmcObject::initial_param_values() {
    
    for (size_t parID=0; parID < m_parameter.size(); parID++) {
        if ( m_selectedParam[parID] ) {
            
            switch (parID) {
            case 0: // alpha - Uniform(0,1)vaccinationStatus
                m_parameter[parID] = runif(m_gen); 
                //m_parameter[parID] = rexp(m_gen, 500); 
                break;

            case 1: // beta - Uniform(0,5)
                m_parameter[parID] = runif(m_gen, 0.0, 5.0); 
                break;

            case 2: // delta - Uniform(-3,3)
                m_parameter[parID] = runif(m_gen, -3.0, 3.0); 
                break;

            case 3:
                if ( m_lnormPrior[parID] ) { // child - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // child - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 4:
                if ( m_lnormPrior[parID] ) { // sFLHAQ2 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sFLHAQ2 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 5:
                if ( m_lnormPrior[parID] ) { // sFLHAQ3 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sFLHAQ3 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 6:
                if ( m_lnormPrior[parID] ) { // sFLHAQ4 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sFLHAQ4 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 7:
                if ( m_lnormPrior[parID] ) { // sStalkQ2 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sStalkQ2 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 8:
                if ( m_lnormPrior[parID] ) { // sStalkQ3 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sStalkQ3 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 9:
                if ( m_lnormPrior[parID] ) { // sStalkQ4 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sStalkQ4 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;
            
            case 10:
                if ( m_lnormPrior[parID] ) { // sNAQ2 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sNAQ2 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 11:
                if ( m_lnormPrior[parID] ) { // sNAQ3 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sNAQ3 - Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 12:
                if ( m_lnormPrior[parID] ) { // sNAQ4 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sNAQ4- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 13:
                if ( m_lnormPrior[parID] ) { // sHIQ2 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sHIQ2- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 14:
                if ( m_lnormPrior[parID] ) { // sHIQ3 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sHIQ3- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 15:
                if ( m_lnormPrior[parID] ) { // sHIQ4 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // sHIQ4- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 16:
                if ( m_lnormPrior[parID] ) { // iNAQ2 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // iNAQ3- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 17:
                if ( m_lnormPrior[parID] ) { // iNAQ3 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // iNAQ3- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 18:
                if ( m_lnormPrior[parID] ) { // iNAQ4 - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // iNAQ3- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 19:
                if ( m_lnormPrior[parID] ) { // iHIUpper - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // iHIUpper- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;

            case 20:
                if ( m_lnormPrior[parID] ) { // iStalkUpper - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormPrior[parID]);   
                } else {                // iStalkUpper- Uniform(0,5)
                    m_parameter[parID] = runif(m_gen, 0.0, 5.0);
                } 
                break;



            }
            
        }
    }
}

// Modification of attributes
void McmcObject::resetMoves() {
    m_numberOfMoveProposed.clear();
    m_numberOfMoveProposed.resize(m_parameter.size());
    m_numberOfMoveAccepted.clear();
    m_numberOfMoveAccepted.resize(m_parameter.size());
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
}

// Initialize infection time
void McmcObject::initialize_augmented_data() {
    
    int ind;
    size_t house;
    
    for (house=0; house<m_nHH; house++) {
        
        // Initialize infection time
        int nInd = m_data[house].getSize();
        
        for (ind=0; ind<nInd; ind++) {
            int display = 0;
            //if (house == 30) display=1;
            
            m_data[house].initialOnsetTime(ind, m_gen, display);
            m_data[house].initialInfTime(ind, m_maxPCRDetectability, m_gen, display);
            
        }

        // Compute person-to-person transmission rate within household
        int display = 0;
        //if (house == 2) display = 1;
        m_data[house].compute_infectivity_profiles(display);

    }
}

// Initialize infection time
void McmcObject::initial_log_lik() {
    
    size_t house;
    size_t nCases = 0;
    cout << m_parameter[1] << endl;
    for (house=0; house < m_nHH; house++) {
        int display = 0;
        //if (house == 0) display = 1;
        if (m_hhLogLik[house]< -10000000000) display = 1;
        m_hhLogLik[house] = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, display);
        //if (house == 0 || house  == 2) cout << "House " << house << ": " << m_hhLogLik[house] << endl;
        cout << "House " << house << ": " << m_hhLogLik[house] << endl;
        m_globalLogLik += m_hhLogLik[house];
        nCases += m_data[house].nSymptomatic();
    }

/*    cout << nCases << endl;
    cout << m_data.size() << endl;*/
}


// Update parameter value in the mcmc chain
void McmcObject::update_parameter(int parID, double step) {

    // Random walk
	double oldValue = m_parameter[parID];
	double newValue(0.0);
    double randomWalk(0.0);
    double logRatioProposal = 0.0;
    if (parID != 2) {
        randomWalk = m_rateForRandomWalk[parID] * rnorm(m_gen, 0.0, 1.0);
        newValue = oldValue * exp(randomWalk);
        logRatioProposal = randomWalk;
    } else {
        newValue = oldValue + rnorm(m_gen, 0.0, m_rateForRandomWalk[parID]);
    }

    // Log ratio of priors 
    double logRatioPrior = 0.0;
    switch (parID) {
    	case 0: // alpha - Uniform(0,1)
    		if (newValue > 1) logRatioPrior = log(0); 
            //logRatioPrior = logdexp(newValue, 500) - logdexp(oldValue, 500); 
    		break;

        case 1: // beta - Uniform(0,10)
        	if (newValue > 5) logRatioPrior = log(0); 
            //logRatioPrior = logdexp(newValue, 0.001) - logdexp(oldValue, 0.001); 
        	break;

        case 2: // delta - Uniform(-3,3)
            if (newValue > 3 || newValue < -3) logRatioPrior = log(0); 
            break;

        case 3:
        	if (m_lnormPrior[parID]) {  // rUC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rUC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;

        case 4:
        	if (m_lnormPrior[parID]) {  // rVC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rVC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;

        case 5:
        	if (m_lnormPrior[parID]) {  // r1VA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // r1VA - Uniform(0,5)
                if (newValue > 0.1 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 6:
        	if (m_lnormPrior[parID]) {  // r2VA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // r2VA - Uniform(0,5)
                if (newValue > 0.1 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 7:
        	if (m_lnormPrior[parID]) {  // rUA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rUA - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 8:
        	if (m_lnormPrior[parID]) {  // rBA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rBA - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 9:
        	if (m_lnormPrior[parID]) {  // rQ - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rQ - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;


    }


	// Compute new log likelihood
    m_parameter[parID] = newValue;
    std::vector<double> newLogLikHH(m_nHH, 0.0);
    double newLogLikGlobal(0.0);
    size_t house;
    for (house=0; house < m_nHH; house++) {
        int display = 0;
        //if (house == 0) display = 1;
        newLogLikHH[house] = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, display);        
        //cout << newValue << " - house " << house << ": " << newLogLikHH[house] << endl;
        //if (house == 0 || house  == 2) cout << "House " << house << ": " << newLogLikHH[house] << " - " << parID << ": " << newValue << endl;
        newLogLikGlobal += newLogLikHH[house];
    }

    //cout << m_parameter[parID] << endl;

	// Update log likelihood
	double Q = newLogLikGlobal - m_globalLogLik + logRatioProposal + logRatioPrior;

	if ( log(runif(m_gen)) < Q )
	{
		m_globalLogLik = newLogLikGlobal;
		m_hhLogLik = newLogLikHH;
        m_numberOfMoveAccepted[parID]++;
        m_accept[parID] = (1 + (step - 1) * m_accept[parID]) / step;
	}
	else {
		m_parameter[parID] = oldValue;
        m_accept[parID] = (0 + (step - 1) * m_accept[parID]) / step;
	}


    // Adaptative mcmc
    if (step <= m_adaptive_steps && m_adaptive_steps > 0) {
      double delta = 1.0 - step / m_adaptive_steps;
      double diff = m_accept[parID] - m_opt_acceptance;
      m_rateForRandomWalk[parID] *= 1.0 + delta * diff;
    }

	
    m_numberOfMoveProposed[parID]++;
}



// Update augmented infection times
void McmcObject::update_augmented_infection_times() {

    // Loop over houses
    size_t house, ind;
    for (house = 0; house<m_data.size(); ++house) {
        std::vector<int> infected = m_data[house].getInfectedIndex();

        // Loop over infected individuals
        for (ind = 0; ind < infected.size(); ind++) {

            int inf = infected[ind];
            int display = 0;

            // Independent sampler from the incubation period
            double oldValue = m_data[house].getSpInfTime(inf);
            double newValue = m_data[house].newInfTime(inf, m_maxPCRDetectability, m_gen);

            // Log ratio of the proposal (independent proposal) 
            double oldProposal = m_data[house].log_dIncub_g(inf, m_maxPCRDetectability, display);
            double newProposal = m_data[house].log_dIncub_g(inf, m_maxPCRDetectability, display);
            double logRatioProposal = oldProposal-newProposal;

            // Update infection time and person-to-person transmission rates
            m_data[house].setInfTime(inf, newValue);
            m_data[house].update_infectivity_profiles(inf, display);   
            
            // Log likelihood
/*            double differenceLogLik(0.0);
            double newLogLikOfHH(0.0);*/
            double newLogLikOfHH = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, display);
            double currentLogLik = m_hhLogLik[house];
            double differenceLogLik = newLogLikOfHH - currentLogLik;

            if( log(runif(m_gen)) < differenceLogLik+logRatioProposal ) {
                m_globalLogLik += differenceLogLik;
                m_hhLogLik[house] = newLogLikOfHH;
                m_acceptedMoveData++;

            } else {                
                m_data[house].setInfTime(inf, oldValue);
                m_data[house].update_infectivity_profiles(inf);
                
            }

            m_proposedMoveData++;
        }
    }
}





// Update augmented symptom onsets
void McmcObject::update_augmented_symptom_onset() {

    // Loop over houses
    size_t house, ind;
    for (house = 0; house<m_data.size(); ++house) {
        std::vector<int> infected = m_data[house].getInfectedIndex();

        // Loop over infected individuals
        for (ind = 0; ind < infected.size(); ind++) {

            int inf = infected[ind];
            int display = 0;

            // Independent sampler from uniform distribution [0,1]
            double oldValue = m_data[house].getSpOnsetTime(inf);
            double newValue = m_data[house].newOnsetTime(inf, m_gen);

            // Update infection time and person-to-person transmission rates
            m_data[house].setOnsetTime(inf, newValue);
            m_data[house].update_infectivity_profiles(inf, display);   
            
            // Log likelihood
            double newLogLikOfHH = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, display);
            double currentLogLik = m_hhLogLik[house];
            double differenceLogLik = newLogLikOfHH - currentLogLik;

            if( log(runif(m_gen)) < differenceLogLik) {
                m_globalLogLik += differenceLogLik;
                m_hhLogLik[house] = newLogLikOfHH;
                m_acceptedMoveData++;

            } else {                
                m_data[house].setOnsetTime(inf, oldValue);
                m_data[house].update_infectivity_profiles(inf);
                
            }

            m_proposedMoveData++;
        }
    }
}




// Return augmented data
std::string McmcObject::getAugmentedDataHH(int house) {

    std::string toprint;
    for (size_t i=0; i < m_data[house].getSize(); i++)
        toprint += std::to_string(m_data[house].getSpInfTime(i)) + " ";

    return toprint;
}

// Return ids of household members
std::string McmcObject::getIndId(int house){

    std::string toprint;
        for (size_t i=0; i < m_data[house].getSize(); i++)
        toprint += std::to_string(m_data[house].getSpId(i)) + " ";

    return toprint;
}
