#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H
#include <string>
#include <vector>
#include <random>


// Constant parameters
// extern double maxPCRDetectability;
extern double mGamma;
extern double vGamma;
extern double shift;
extern double shapeInf;
extern double scaleInf;


//----------------------------------------------------------------------
// Infectivity profile
//----------------------------------------------------------------------
double infectivityProfile(
	double origin,
	double infectionDate,
	int infectionStatus,
	double tinf,
	int studyPeriod,
	int display = 0
	);
double cumulativeInfectivity(
	double origin,
	double infectionDate,
	int infectionStatus,
	double tinf,
	int studyPeriod,
	int display = 0
	);


//----------------------------------------------------------------------
// Household class
//----------------------------------------------------------------------
class Household
{
public:
	Household();
	Household(std::string contact_pattern,std::vector<double> instval,std::vector<double>cumval);
	~Household() {};
	// int getStudyPeriod() const { return m_studyPeriod; };
	
	size_t getSize() const { return m_size; };
	size_t nInfected() const { return m_confCase.size(); };
	size_t nSymptomatic() const;
	std::string getindid() const;
	int getSpInfected(int index) const { return m_infected[index] ; };
	std::vector<int> getSpInfected(std::vector<int> index) const;
	std::vector<int> getAllInfected() const { return m_infected; };
	std::vector<int> getInfectedIndex() const { return m_confCase; };
	double getSpInfTime(int index) const { return m_infTime[index]; };
	double getSpOnsetTime(int index) const { return m_augmentedOnsetTime[index]; };
	std::vector<double> getAllInfTime() const { return m_infTime; };
	int getSpId(int index) const { return m_indid[index]; };

	void addIndividual(int indid, double onsetTime, int isCase, int age,int sex, int startFollowUp, int studyPeriod, int flha_titer_q, int stalk_titer_q, int na_titer_q, int haihk14_titer_q, int high_titer_sum);
	void setInfTime(int index, double infTime);
	void setOnsetTime(int index, double onsetTime);
	void newHousehold();
	void initialInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen, int display=0);
	void initialOnsetTime(int index, std::mt19937_64& gen, int display=0);
	void compute_infectivity_profiles(int display = 0);
	void update_infectivity_profiles(int ind, int display = 0);

	double newInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen);
	double newOnsetTime(int index, std::mt19937_64& gen);
	double log_dIncub(int index, double maxPCRDetectability, int display = 0);
	double log_dIncub_g(int index, double maxPCRDetectability, int display = 0);

	double compute_log_lik(std::vector<double> parameters, std::vector<int> selectedParam, double maxPCRDetectability, int display = 0);
	double log_pInf(size_t curr, std::vector<double> parameter, std::vector<int> selectedParam, int display = 0);
	double log_S(size_t curr, double t0, std::vector<double> parameter, std::vector<int> selectedParam, int display = 0);

private:
	size_t m_size;
	int m_notInfected;
	int m_startFollowUp;
	std::string m_contactPattern;
	std::vector<int> m_indid;
	std::vector<int> m_sex;
	std::vector<int> m_flha_titer_q;
	std::vector<int> m_stalk_titer_q;
	std::vector<int> m_na_titer_q;
	std::vector<int> m_haihk14_titer_q;
	std::vector<int> m_high_titer_sum;
	std::vector<double> m_onsetTime;
	std::vector<double> m_augmentedOnsetTime;
	std::vector<double> m_instcommunityinf;
	std::vector<double> m_cumcommunityinf;
	std::vector<double> m_incubPeriod;
	std::vector<int> m_infected;
	std::vector<int> m_confCase;
	std::vector<int> m_studyPeriod;
	std::vector<int> m_age;
	std::vector<double> m_infTime;
	std::vector<std::vector<double>> m_cumLambda;
	std::vector<std::vector<double>> m_instLambda;
};

#endif