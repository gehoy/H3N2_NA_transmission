#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "Household.h"
#include "McmcObject.h"

using namespace std;

//----------------------------------------------------------------------
// Load data
// datafile with columns : 1: indid (int), 2: currhhid (int), 3: hhsize (int), 4: onsetTime (double), 5: studyPeriod (double), 6: isCase (0 or 1), 7:  age (0 or 1), 8: sex (1 or 2), 9: flha_titer_q (int), 10: stalk_titer_q (int), 11: na_titer_q (int), 12: haihk14_titer_q (int)
//----------------------------------------------------------------------
std::vector<Household> buildData(std::string dataFile, std::string communityInfFile, std::string contact_pattern)
{

    std::vector<double> instinfcom;
    std::vector<double> cuminfcom;

    instinfcom.resize(0);
    cuminfcom.resize(0);

    std::ifstream infileinf(communityInfFile.c_str());
    if(infileinf)
    {
        cout << "=============INFILE===============" << endl;
        cout << "Read data file: " << dataFile << endl;

        int date(0);
        double instval = 0.0;
        double cumval = 0.0;

        // while(std::getline(infile, line))
        for ( std::string line; std::getline(infileinf, line); )
        {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> date >> instval >> cumval;

            //cout << currhhid << endl;

            instinfcom.push_back(instval);
            cuminfcom.push_back(cumval);
        }

    }

    // Variable with all the data
    std::vector<Household> output;
    output.resize(0);

    // Read file
    /*
    The data file should be a space-separated table with household ids
    sorted in ascending order and numbered from 0 to n-1 households
    */
    std::ifstream infile(dataFile.c_str());
    if(infile)
    {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << dataFile << endl;

        int indid(0), lasthhid(-1), currhhid(0), hhsize(0), isCase(0), startFollowUp(0), studyPeriod(0), age(0), sex(0),flha_titer_q(0),stalk_titer_q(0),na_titer_q(0),haihk14_titer_q(0),high_titer_sum(0);
        double onsetTime = 0.0;
        double ddi = 0.0;
        int numberOfCase(0), numberOfHousehold(0), numberOfSubject(0), numberOfDay(0);
        // datafile with columns : 1: indid (int), 2: currhhid (int), 3: hhsize (int), 4: onsetTime (double), 5: studyPeriod (double), 6: isCase (0 or 1), 7:  age (0 or 1), 8: sex (1 or 2), 9: flha_titer_q (int), 10: stalk_titer_q (int), 11: na_titer_q (int), 12: haihk14_titer_q (int)
        //int nSymptomatic(0);

        Household currHH(contact_pattern,instinfcom,cuminfcom);

        for ( std::string line; std::getline(infile, line); )
        {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> indid >> currhhid >> hhsize >> onsetTime >> studyPeriod >> isCase >> age >> sex >> flha_titer_q >> stalk_titer_q >> na_titer_q >> haihk14_titer_q >> high_titer_sum;

           // cout << recovered << endl;

            // Append previous household to the list of households
            if (currhhid != lasthhid)
                {
                    if (numberOfHousehold > 0) {
                        numberOfHousehold++;
                        output.push_back(currHH);
                        currHH.newHousehold();
                        lasthhid = currhhid;
                    } else {
                        numberOfHousehold++;
                        lasthhid = currhhid;
                    }
                }

            // Update household object and information parameters
            currHH.addIndividual(indid,
                                 onsetTime,
                                 isCase,
                                 age,
                                 sex, 
                                 startFollowUp,
                                 studyPeriod,
                                 flha_titer_q,
                                 stalk_titer_q,
                                 na_titer_q,
                                 haihk14_titer_q,
                                 high_titer_sum
                                 );
            numberOfSubject++;
            if (isCase != 0) numberOfCase++;
            numberOfDay = max(numberOfDay, studyPeriod);
        }

        // Add last household
        output.push_back(currHH);

        // General information
        cout << "Number of households: " << numberOfHousehold << endl;
        cout << "Number of individuals: " << numberOfSubject << endl;
        cout << "Number of cases: " << numberOfCase<< endl;
        cout << "Number of days: " << numberOfDay << "\n\n";

        return(output);

    }else{
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------
void runMCMC(McmcObject mcmc,
             std::string outputFile,
             std::string augmentedDataFile,
             int pas,
             /*int pasPrinting,*/
             std::vector<int> idOfSelectedParameter,
             std::vector<std::string> paramNames
)
{
	ofstream output(outputFile.c_str());
    //ofstream augmentedData(augmentedDataFile.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
    
	int nIterTimeInfection = mcmc.getNIterTimeInf();

	// Column names of the output file
	std::string colNames="iteration logLik ";
	for (size_t p=0; p < paramNames.size(); p++) colNames += paramNames[p] + " " + paramNames[p] + "_p " + paramNames[p] + "_a ";
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    // Initial state
    mcmc.initial_param_values();
	mcmc.initialize_augmented_data(); // Initialize infection time and symptom onset for all infected individuals
	mcmc.initial_log_lik(); 
	cout << "Initial log likelihood: " << mcmc.globalLogLik() << endl;

    output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
    for (size_t i = 0; i < mcmc.nParameters(); i++)
        output << mcmc.parameter(i) << " 0 0 ";
    output << "0 0" << endl;

    //for (size_t i = 0; i < mcmc.getNumbHH(); i++)
    //    augmentedData << mcmc.getIndId(i);
    //augmentedData << endl;

    //for (size_t i = 0; i < mcmc.getNumbHH(); i++)
    //    augmentedData << mcmc.getAugmentedDataHH(i);
    //augmentedData << endl;


    // Chain
	for (iteration = 0; iteration < numberOfIteration; iteration++)
	{
        cout << "iteration " << iteration << endl;
	    mcmc.resetMoves();

		for (iter = 0; iter < pas; iter++)
		{
			for (size_t selectedParameter = 0; selectedParameter < idOfSelectedParameter.size(); selectedParameter++)
			{
				parID = idOfSelectedParameter[selectedParameter];
				mcmc.update_parameter(parID, iteration*pas+iter+1);
			}
            //cout << "numberOfIteration " << iter << endl;

			// Augmented data
			for (int i=0; i < nIterTimeInfection; i++) {
                mcmc.update_augmented_infection_times(); // Update loglik at the household level
                mcmc.update_augmented_symptom_onset();
			}
            //cout << iteration * pas + iter << " " << mcmc.globalLogLik() << endl;
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood
        //cout << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood

		for (size_t i = 0; i < mcmc.nParameters(); i++)
            output << mcmc.parameter(i) << " " << mcmc.proposedMove(i) << " " << mcmc.acceptedMove(i) << " ";

        output << mcmc.proposedMoveData() << " " << mcmc.acceptedMoveData() << endl;


        // write the output file of augmented data
        //for (size_t i = 0; i <mcmc.getNumbHH(); i++)
        //    augmentedData << mcmc.getAugmentedDataHH(i);
        //augmentedData << endl;
	}

	output.close();
    //augmentedData.close();
}



//----------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------
int main(int argc, char **argv)
{
    
    // Arguments passed to main
    std::string inference_contact = argv[1];  // If we estimate contact pattern parameters
    int numberOfIteration = std::stoi(argv[2]);  // Number of iterations in the MCMC
    std::string name_output = argv[3];  // Name of the output file
    std::string name_data = argv[4];  // Name of the date file
    std::string path_data = argv[5];  // Path data file
    std::string path_to_save = argv[6];  // Path otput file

   double maxPCRDetectability = 10.0;     
    std::string chainID = "2";

    //==========Model parameters==========
    // Initial values will need to update
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "child", "sHigh_titer_one", "sHigh_titer_two", "sHigh_titer_three","iHigh_titer_one","iHigh_titer_two","iHigh_titer_three"};
    //                                            0        1      2        3           4         5          6          7          8          9     
    int numberOfParameters = parameterNames.size();
    //Use a lognorm prior or not for parameters (0 or 1) (useless for alpha, beta, delta)
    std::vector<int> lnormPrior(numberOfParameters);
    lnormPrior = {1,1,1,1,1,1,1,1,1,1};

    //Standard deviation of the lognormal prior
    std::vector<double> sdLNormPrior(numberOfParameters);
    sdLNormPrior = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    
    std::vector<double> parameter(numberOfParameters);

    // Parameters to infer according to model
    std::vector<int> selectedParameter(numberOfParameters, 0);
    
    //alpha
    parameter[0] = 0.001;
    selectedParameter[0] = 1;
    
    //beta
    parameter[1] = 0.6;
    selectedParameter[1] = 1;
    
    //delta 
    parameter[2] = 1.0;
    selectedParameter[2] = 1;

    //child
    parameter[3] = 1.0;
    selectedParameter[3] = 1;
    
    //sHigh_titer_one
    parameter[4] = 1.0;
    selectedParameter[4] = 1;

    //sHigh_titer_two
    parameter[5] = 1.0;
    selectedParameter[5] = 1;

    //SHigh_titer_three
    parameter[6] = 1.0;
    selectedParameter[6] = 1;

    //iHigh_titer_one
    parameter[7] = 1.0;
    selectedParameter[7] = 1;

    //iHigh_titer_two
    parameter[8] = 1.0;
    selectedParameter[8] = 1;

    //iHigh_titer_three
    parameter[9] = 1.0;
    selectedParameter[9] = 1;


	int parameterNumber;
	std::vector<int> idOfSelectedParameter(0);
    for(parameterNumber=0; parameterNumber<numberOfParameters; parameterNumber++)
    {
        if(selectedParameter[parameterNumber]==1) idOfSelectedParameter.push_back(parameterNumber);
    }
    
    //Print id of selected parameters
    cout << "ID of selected parameters: "; 
    for (auto i = idOfSelectedParameter.begin(); i != idOfSelectedParameter.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;


    //==========MCMC parameters==========
/*    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();*/
    unsigned seed;
    if ( chainID == "1" ) seed=20201208;
    if ( chainID == "2" ) seed=20210415;
    if ( chainID == "3" ) seed=20210601;
    cout << "Seed: " << seed << endl;

    int pas = 20; 
    //int nAdaptiveSteps = pas * numberOfIteration / 4;   // Doesn't work much...
    int nAdaptiveSteps = 0;
    int numberOfIterationTimeInfection = 1;

    // Variance of random walk
    std::vector<double> rateForRandomWalk(numberOfParameters); //higher the random walk, bigger variance of proposed moves, lower chance of accepting any given move
    rateForRandomWalk = { 1.5,    0.3,     0.3,   0.4,  0.75,  0.75,  0.75, 0.75, 0.75, 0.75};
    //                  "alpha", "beta", "delta", "child", "sFLHAQ2", "sFLHAQ3","sFLHA14","rFC"
    

    //==========Output files==========
    //Paths
    std::string condition;
    std::string dataFile, outputFile, augmentedDataFile,communityfile;


    
    //Data file
    // File structure :         0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: studyPeriod, 6: adult
    dataFile=path_data+name_data+".txt";         //Input file 
    outputFile=path_to_save + "mcmc_" + name_output + ".txt"; //Output file
    augmentedDataFile=path_to_save + name_output + ".txt"; //Output file
    communityfile=path_data+"community_infection_israel.txt";// communityfile

    
    cout << "Input file: " << dataFile << endl;
    cout << "Output file: " << outputFile << "\n\n";

    //==========Build data==========
    // Load data
    std::vector<Household> hhData = buildData(dataFile, communityfile ,inference_contact);

    // Initialize MCMC object
    McmcObject mcmc(
        seed, 
        numberOfIteration, 
        nAdaptiveSteps,
        hhData, 
        parameter, 
        selectedParameter, 
        rateForRandomWalk,
        lnormPrior,
        sdLNormPrior,
        numberOfIterationTimeInfection, 
        maxPCRDetectability          
        );

    //==========MCMC==========
    runMCMC(
        mcmc,
        outputFile,
        augmentedDataFile,
        pas,
        idOfSelectedParameter,
        parameterNames
        );

    return 0;
}
