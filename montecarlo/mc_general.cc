#include "quintcosmos.h"
#include "global.h"
#include <list>
#include <fstream>
#include <string> 
#include "cleanvector.h"
#include "safevector.h"
#include "lensing.h"
#include "controlpanel.h"
#include "cl.h"
#include "cmbcalc.h"
#include "spline.h"
#include "model.h"
#include "perturbation.h"
#include "analyzethis.h"
#include <unistd.h>
#ifndef NOWMAP7
#include "WMAP_7yr_options.h"
#endif

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "multigaussian.h"
#include "rollingaverage.h"
#include "sdsslrg.h"
#include "cbi2.h"

const string outputname = "/output/";

string base = ControlPanel::cmbeasyDir(outputname);
string errorname = base+"errorlog.txt";

const int  MAGIC_CHECK=12345678;

/*!
  Example Monte-Carlo driver. Most things are controlled via the variables that follow
  right below.

  If you would like to run your own models, go ahead and change the parts in slave
  that set up the cosmology (and evaluate the models). 
  
  Most propably, you won't have to change the master() routine

*/

/* 
   Experiments to use 
 */

const bool Astier = false;
const bool Riess06 = false;

const bool SDSS = true;
const bool TwodF = false;
const bool SDSS_BAO = false;  // Baryon acoustic peak 
const bool SDSSLRG = false; // Luminous reg galaxies

const bool LYA_MCDONALD=false;  // Lyman alpha code of Pat McDonald

#ifndef NOWMAP7
const bool WMAP7 = false; // WMAP 7-year yes or no ?
#else 
const bool WMAP7 = false;  // don't touch this
#endif
const bool BOOMERANG03 = false; // Boomerang from 2003 flight
const bool VSA = false; 
const bool ACBAR = false;  
const bool CBI=false;


/*! The number of chains; 1-9 chains are supported, if you want more,
   you'll have to slightly modify the code. No big problem, though. */
const unsigned int CHAINS=2;

/*! The number of parameters in our model */
unsigned int PARAMETERS=6;

/*! Other information (derived parameters such as sigma 8 etc) can be stored from this
position on.  MULTIPURPOSE_REQUIRED is the number of variables that you can
store in the chain on top of the parameters and likelihoods. 
*/
const int MULTIPURPOSE_REQUIRED = 1; 

/*! The position in the Taskarray where the likelihoods will be written
   LOGLIKEPOS is automatically set by the program
*/
int LOGLIKEPOS=0;

/*! The size of the array, i.e. all possible values (parameters, likelihoods, derived parameters) that we wish 
to write to a file. Automatically set by the program   */
unsigned int TASKARRAY_SIZE=0;  
const unsigned int MAX_TASKARRAY_SIZE = 100; // this is the maximum task array size
int MULTIPURPOSEPOS=0; // will be set automatically 

/*
  Needed to keep track of   parameters  needed for experiments
*/
int AstierParameterPos=0;
int WMAP7ParameterPos=0;
int SDSSLRGParameterPos=0;

/*! minimum number of points for R-statistic calculation*/
const  unsigned int RMIN_POINTS = 30;  // minimum number of points for R-statistic calculation


// FLAGS & VARIABLES FOR ADAPTIVE STEPSIZE ONLY
/*!Flags for the adaptive stepsize algorithm
 with or without adaptive stepsize? */
bool ADAPTIVE=true;

/*! If you want to stop adaptive stepsize after convergence has been reached;
freeze-in will be performed once Gelman-Rubin statistics indicate convergence
sets ADAPTIVE to false once converged.*/
bool FREEZE_IN=true;


/*! only if FREEZE_IN=true: If R-statistic for all parameters less than RBREAK and
 number of points more than MIN_SIZE_FOR_FREEZE_IN then freeze-in. 
 If you want to be really conservative, you can set RBREAK=1.1 */
const double RBREAK=1.2;
const unsigned int MIN_SIZE_FOR_FREEZE_IN=500;


/*!For the adaptive stepsize algoritm, variable step length; these numbers are
rather arbitrary, but seem to work well
Please note, that after FREEZE_IN, the covariance matrix without the 
step factor will be used. Hence, the average stay at a certain parameter point
is approx 3 */
const double INCREASE_STEP=1.15;
const double DECREASE_STEP=0.9;
const int  HIGH_STEP_BOUND=5;
const int  LOW_STEP_BOUND=3;

/*! For the covariance updater; recompute covariance matrix after UPDATE_TIME steps */
const int UPDATE_TIME=5;

/*! Start using the estimated covariance matrix for steps after BEGINCOVUPDATE points have
been computed in a chain */
const unsigned int BEGINCOVUPDATE=100;

/*! For the adaptive stepsize: Do not compute covariance Matrix with more than MAXPREVIOUSPOINTS*/
const unsigned int MAXPREVIOUSPOINTS=5000;


/*! A task is just a double array.  The benefit of having
  this structure anyhow is that we can very easily store tasks
  in lists and all that.
  Each double variable corresponds to either a parameter,
  a likelihood, some multipurpose stuff or just no information
  at all. 
  I find it convenient to leave [0..7] for parameters
  [8..14] is for likelihoods and after that come multipurpose 
  stuff. If you only need 4 parameters, there will be garbage
  in [5..7], but as you don't have to look at it, who cares, if 
  you don't ?
*/
struct Task {  
  double  f[MAX_TASKARRAY_SIZE];
  int Multiplicity; //!< the weigh of this point in the chain, i.e. how long does it rest at the point
  /*!
    The number of times a different model has been simulated until the step has been
    taken. This is Multiplicity minus the number of times the step proposal crossed the boundaries
    of parameter Space. We take this number as a meassure to increase or decrease step sizes
  */
  int ReallyInvestigated;
  Task() : Multiplicity(1), ReallyInvestigated(1)  {};
};


const char inttochar[10]={'0','1','2','3','4','5','6','7','8','9'};

/*! This function implements the transition Kernel. In our case, this is a Gaussian distribution with 
  variance sigma, which we call the "stepsize".
  \param old the parameters at present
  \param newtask a pointer with the new task. 
  \param gauss the multigaussian 
 */
bool throwDice(const Task& old, Task *newtask, MultiGaussian *gauss) {
  for (unsigned int l=PARAMETERS; l < TASKARRAY_SIZE; l++) newtask->f[l]=0; 
  bool success = gauss->throwDice(old.f);
  for (unsigned int l=0; l < PARAMETERS; l++) newtask->f[l]=gauss->getRandomValue(l);
  return success;
}

void setVariables() {	
    if (Astier) {
      AstierParameterPos = PARAMETERS;
      PARAMETERS += 2; // two more parameters for Astier
    }
    if (SDSSLRG) {
      SDSSLRGParameterPos = PARAMETERS;
      PARAMETERS +=1; // bias
    }

    LOGLIKEPOS = PARAMETERS + 1;  // leave at least one blank collumn
    // at least at position 15. That means that most of the time, the likelihoods will be at exactly the same
    // position in the file. That's more convenient
    LOGLIKEPOS = max(LOGLIKEPOS, 15);
    MULTIPURPOSEPOS = LOGLIKEPOS + 22; // leave one blank collumn between
    TASKARRAY_SIZE = MULTIPURPOSEPOS + MULTIPURPOSE_REQUIRED + 1;

    if (TASKARRAY_SIZE > MAX_TASKARRAY_SIZE) throw Bad_Error("TASKARRAY_SIZE exceeds limit");
}


// Some constants needed for multiprocessor application
const int GIMMETASK=1; //!< MPI flag indicating slave()'s wish to receive a task from master()
const int TAKERESULT=2;  //!< MPI flag indicating slave()'s request to master() to take result
const int TAKETASK=3;  //!< MPI flag indicating master()'s request to slave to take task.

/*! The master is responsible for sending tasks to the slaves and collecting the information from the slaves
after the calculation. In other words, this is the coordinator. The slaves don't know about each others 
existence. The master also monitors the output. There is only one master, but an arbitrary number of slaves.
 */
void master(bool Restart) {     
  try {
    // we also use a cmbcalc briefly to ask initjl to read in bessel functions
    // if there is no jlgen.dat in the resource directory, this will automatically generate
    // bessel functions up to l = 5000
    CmbCalc *cmbcalc = new CmbCalc();
    string cmbeasydir = ControlPanel::cmbeasyDir();
    //    filename = cmbeasydir + "/resources/jlgen.dat";  // bessel function filename
    cmbcalc->initjl(cmbeasydir + "/resources/jlgen.dat",2000);  // will trigger the generation if not already there
    delete cmbcalc; 

    vector< list<Task> > chain(CHAINS);
    unsigned int chainSize[CHAINS];   // chain size
    unsigned int chainTotalPerformed[CHAINS];   // total number of performed models
    //    int count[CHAINS];
    unsigned int checkSize=0;
  
    // for dynamical stepsize
    int steps_since_update[CHAINS];
    double EntireFactor[CHAINS];
    vector< RollingAverage* >  roll(CHAINS);   // rolling average of EntireFactor for after burn-in
    vector< vector<double> > sigma(CHAINS);
    vector<vector<double> > covMatrix[CHAINS];
    vector<Task*> next(CHAINS);

    double AverageMultiplicity[CHAINS]; // just for our information
    //initialize to 0 
    for (unsigned int i=0; i < CHAINS; i++) AverageMultiplicity[i]=0;

    MultiGaussian *mGauss[CHAINS]; // multivariate gaussian random number generator

    
    cout << "PARAMETERS: " << PARAMETERS << endl;
    cout << "LOGLIKEPOS: " << LOGLIKEPOS << endl;
    cout << "MULTIPURPOSEPOS: " << MULTIPURPOSEPOS << endl;

    vector<double> lowbound(PARAMETERS);
    vector<double> highbound(PARAMETERS);
  
    //set the (flat)  priors, i.e. parameter boundaries
    lowbound[0]=0.05 ; highbound[0]=0.24;  // omega_mh2;
    lowbound[1]=0.016 ; highbound[1]=0.03 ; // omega_bh2
    lowbound[2]=0.50; highbound[2]=0.85 ; // h
    lowbound[3]=0.0 ; highbound[3]=0.3; // optdlss
    lowbound[4]=0.8 ; highbound[4]=1.4;  // n
    lowbound[5]=2.5; highbound[5] = 3.2; // ln (10^10 A_s) - 2\tau 
    // add your parameters here
    

    if (Astier) {
      lowbound[AstierParameterPos] = 1; highbound[AstierParameterPos] = 2.2; // alpha Astier
      lowbound[AstierParameterPos+1]=1.1; highbound[AstierParameterPos+1]=2.4;  // beta  Astier
    }
    if (SDSSLRG) {
      lowbound[SDSSLRGParameterPos]=0.2;  highbound[SDSSLRGParameterPos]=5;  // bias 
    }
    
    
    
    //  Just keep on, if you have more parameters
   

    // this sets the starting points and stepsizes for each chain for the first step
    for (unsigned int i = 0; i < CHAINS; i++) {
      roll[i] = new RollingAverage(500);
      Task t;
      for (unsigned int k=0; k < TASKARRAY_SIZE; k++) t.f[k]=0; // for safety reasons, initialize (otherwise garbage in first line)
      mGauss[i]=new MultiGaussian(PARAMETERS);
      mGauss[i]->setBounds(lowbound,highbound);
    
      //Generate Starting points (inside the prior)
      for(unsigned int j=0; j<PARAMETERS; j++) { 
	t.f[j] = lowbound[j] + Miscmath::posRnd(highbound[j]-lowbound[j]);
      }

      // the sigmas; these need to be chosen with care, which is why we use dynamic stepsize
      // the sigmas only matter in the beginning (or if you de-select ADAPTIVE) 
      sigma[i].resize(TASKARRAY_SIZE);
      sigma[i][0] = 0.01; // omega_mh2
      sigma[i][1] = 0.001; // omega_bh2
      sigma[i][2] = 0.032; //h
      sigma[i][3] = 0.03; // optdlss
      sigma[i][4] = 0.02; // n
      sigma[i][5] = 0.03; //  ln (10^10 A_s) - 2\tau  
      // add your sigmas here

      if (Astier) {
	sigma[i][AstierParameterPos] = 0.15; // alpha
	sigma[i][AstierParameterPos+1] = 0.15; // beta
      }
      if (SDSSLRG) {
	sigma[i][SDSSLRGParameterPos]=0.4;  // bias 
      }
       

      // .. more sigmas, if you have more parameters
    
      covMatrix[i].resize(PARAMETERS);
      for(unsigned int k=0; k< PARAMETERS;k++) covMatrix[i][k].resize(PARAMETERS); // create covariance matrix
      for(unsigned int k=0; k< PARAMETERS; k++){                                 //fill with sigmas
	for(unsigned int j=0; j< PARAMETERS; j++){
	  if (j==k) covMatrix[i][j][j]=sigma[i][j]*sigma[i][j];
	  else covMatrix[i][j][k]=0.0;
	}
      }

      mGauss[i]->generateEigenvectors(covMatrix[i],1.0);
      mGauss[i]->printExtInfo(cout);


      steps_since_update[i]=0; 
      EntireFactor[i]  = 1.0;
    
      cout << "i: " << i << endl;
   
    
      next[i]=0;    // next model to be computed

      // we set the loglike for this starting point to be really small, this will wash
      // away in burn in
      t.f[LOGLIKEPOS] = -1e100;
      chain[i] = list<Task>();   // a list of tasks for every chain
      if (! Restart) chain[i].push_back(t); // this is our very first task
      chainSize[i] = 1; // list.size() is an order (N) operation
      chainTotalPerformed[i]=0;
    }  
  
    // The chain output will be written to these files
    ofstream data[CHAINS];
    ofstream head[CHAINS]; // for information: this is the model currently at the top of the chain
    ofstream investigated[CHAINS]; //all models, whether taken or not taken and their chi^2's
    //
    // Change the string "base" below, if you like to get your results into a different 
    // (existing) directory
    //
    string montecarlo_name="montecarlo_chain";
    string head_name= "head";
    string investigated_name="investigated";
    string covfile_name="covfile";
    ios::openmode Mode = ios::out;
    // If we have restarted, we need to read in and set some stuff

    if (Restart) { 
      Mode = ios::app;  // append mode when re-opening the files
      for (unsigned int k =0; k < CHAINS; k++) {	  
	ifstream readcov;
	readcov.open((base+covfile_name + inttochar[k+1]+ ".dat").c_str());
	Task t;
	for (unsigned int j=0; j < TASKARRAY_SIZE; j++) {
	  t.f[j] =  AnalyzeThis::read<double>(readcov);
	}
	if (readcov.fail()) throw Bad_Error("Covariance file corrupt or non-existent. Cannot restart.");
	t.Multiplicity =  AnalyzeThis::read<int>(readcov);
	t.ReallyInvestigated =  AnalyzeThis::read<int>(readcov);
	cout << "read in: " << t.Multiplicity << "  :: " << t.ReallyInvestigated << endl;
	chain[k].push_back(t);
	for(unsigned int j=0; j < PARAMETERS; j++){
	  for (unsigned int m=0; m < PARAMETERS; m++){
	    covMatrix[k][j][m] = AnalyzeThis::read<double>(readcov);
	  }
	}
	int check = AnalyzeThis::read<int>(readcov);
	if (check != MAGIC_CHECK) throw Bad_Error("MAGIC CHECK FAILED");
	double OptimalFactor = 2.4 / sqrt((double)PARAMETERS);
	mGauss[k]->generateEigenvectors(covMatrix[k],OptimalFactor*OptimalFactor);
	EntireFactor[k] = OptimalFactor;
	mGauss[k]->lock();
	ADAPTIVE=false;
	FREEZE_IN=false;
      }
    }
    // End of reading in to restart`
    for (unsigned int k=1; k<= CHAINS; k++) {
      data[k-1].open((base+montecarlo_name + inttochar[k]+ ".dat").c_str(),Mode); // full name is montecarlo_chain1.dat
      head[k-1].open((base+head_name + inttochar[k]+ ".txt").c_str(),Mode);
      investigated[k-1].open((base+investigated_name + inttochar[k]+ ".dat").c_str(),Mode);
    }
    for(unsigned int m=0; m< CHAINS ; m++) { 
      data[m].setf(ios::scientific);
      head[m].setf(ios::scientific);
      investigated[m].setf(ios::scientific);
    }
    
    // a summary of information will be written to this file 
    ofstream progress((base+"progress.txt").c_str(),Mode);
    if (Restart) progress << "\n\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: RESTARTED FROM HERE :::::::::::::::::::::::::::::::::::::\n\n";
    // Gelman-Rubin statistical output into this file
    ofstream gelmanRubin((base+"gelmanRubin.txt").c_str(),Mode);
    // covariance output for this file

    ofstream covarianceMatrices((base+"covMatrix.txt").c_str(),Mode);
    covarianceMatrices.setf(ios::scientific);
   
    // END OF INITIALIZATION

    bool work = true;  //as long as its value is true, the chains will run.  
    MPI_Status  status;
    double result[MAX_TASKARRAY_SIZE]; // should be larger than TASKARRAY_SIZE

    /*
    **
    ** Now that we have initialized everything, we go into an 
    ** eternal loop: We wait for a message from the slaves (either
    ** GIMMETASK or TAKERESULT
    ** and send out work and receive the   results.
    ** The master will also calculate the Gelman-Rubin statistic,
    ** which takes up most of it's size
    */
    do {    
      // these are messages from the slaves to the master
      MPI_Recv(&result, TASKARRAY_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  
      // we have a message, check what it sais; 
      if (status.MPI_TAG == GIMMETASK) {
	// the slave wants something to do, so give him something to do
	int i = status.MPI_SOURCE - 1; // which slave sent this?
	if (next[i]) delete next[i];  
	// generate a new set of parameters from the transition kernel and send to slave
	next[i] =new Task();
	for (;;) {  // forever (if not break)
	  if ( throwDice(chain[i].back(),next[i],mGauss[i]) ) break;  // if succesful step, quit
	  // so we have crossed the boundary. Hence we have to increase the weight of current point
	  chain[i].back().Multiplicity++;
	}
	// just for information: this file is always one ahead 
	for(unsigned int k=0; k< TASKARRAY_SIZE;k++) head[i] << next[i]->f[k] << "  ";
	head[i] << endl; head[i].flush();
	MPI_Send(&next[i]->f, TASKARRAY_SIZE, MPI_DOUBLE,status.MPI_SOURCE, TAKETASK, MPI_COMM_WORLD);
      }
    
      if (status.MPI_TAG == TAKERESULT) {
	// the slave sent us a result
	int i = status.MPI_SOURCE - 1; // which slave sent this?
	// cout << "TAKING RESULT FROM SOURCE: " << i << endl;
	// for (unsigned int k = 0;k < TASKARRAY_SIZE; k++) cout << "TOOK BACK " << k << "    :::: " << result[k] << endl;
	//	count[i]++;
      
	// First, keep track of the model, no matter whether we take or do not take the step
	// store this in a file called "investigatedxyz.dat"

	for(unsigned int k=0; k< TASKARRAY_SIZE;k++) investigated[i] << result[k] << "  ";
	investigated[i] << " 1" << endl; investigated[i].flush();

	chainTotalPerformed[i]++; 
	// we have the result, now we have to decide whether or not to take the step
	bool take = false;
	double loglike1 = chain[i].back().f[LOGLIKEPOS]; // get the Likelihood of the previous parameter set
	if (result[LOGLIKEPOS] <= loglike1) {    
	  // calculate likeli2 / likeli1. We subtract loglike2 first, doesn't change
	  // result, but helps prevent 0/0
	  double ratio = 1.0/exp( loglike1 - result[LOGLIKEPOS]);
	  double x = Miscmath::posRnd(1.0);
	  if (x < ratio) take=true;
	} else take = true; // always take if the likelihood is larger than the previous
      
	if( ADAPTIVE ) roll[i]->push(EntireFactor[i]);  // no matter if take or not, note EntireFactor for average

	// we know now whether or not to take the step. 
	if (take) {
	  // first, write the last one in the chain to data file
	  Task &last = chain[i].back();
	  for(unsigned int k=0; k< TASKARRAY_SIZE;k++) data[i] << last.f[k] << "  ";
	  data[i] << last.Multiplicity << endl; // Multiplicity, how often is this point in the chain?
	  data[i].flush();
	  
	  
	 
	  // ADAPTIVE: increase stepsize, if we take steps  too often
	  if( ADAPTIVE ){ 
	    if (chain[i].back().ReallyInvestigated < LOW_STEP_BOUND) {
	      if (EntireFactor[i] < 10) { 
		EntireFactor[i] *=  INCREASE_STEP;
		mGauss[i]->generateEigenvectors(covMatrix[i],EntireFactor[i]*EntireFactor[i]);
	      }
	    }
	  }
	  // end of ADAPTIVE increase
	  Task t;
	  for (unsigned int k=0; k < TASKARRAY_SIZE; k++) t.f[k] = result[k];
	    chain[i].push_back(t);
	    chainSize[i]++;		
	    steps_since_update[i]++;  // a taken step contributes towards steps_since_update counting
	}  else {
	  chain[i].back().Multiplicity++;  // enhance multiplicty of the old one again
	  chain[i].back().ReallyInvestigated++; // we did a simulation and we didn't take the step
	}
	
	/*
	**
	** Adaptivly setting stepsize according to covariance matrix
	** Here, we estimate the covariance
	**
	*/
	if(ADAPTIVE && steps_since_update[i] >= UPDATE_TIME && chainSize[i] >=  BEGINCOVUPDATE) {
	  double average[PARAMETERS];	
	  double cov[PARAMETERS][PARAMETERS];	  
	  int TotalSize =0;

	  // Here, we determine the number of distinct points in the chain
	  // that will be taken into account to estimate the covariance matrix
	  unsigned int covSize;  	 
	  if (chainSize[i] < 2*BEGINCOVUPDATE) { // early on
	    covSize = chainSize[i]/2;
	  } else {
	    covSize = chainSize[i] - BEGINCOVUPDATE;  // all but the first few begincovupdate points
	  }
	  covSize = min(covSize,MAXPREVIOUSPOINTS);  // at most MAXPREVIOUSPOINTS

  
	  for(unsigned int k=0; k < PARAMETERS; k++){ // initialize
	    average[k]=0;
	    for(unsigned int j=0; j<PARAMETERS; j++) cov[k][j]=0.0;
	  }
  
	  vector<list<Task>::iterator > chainIter(CHAINS); 
	  chainIter[i] = chain[i].begin();
	  // ignore the first part of the chain, i.e. the ones to be excluded 
	  for (unsigned int j=0; j < chainSize[i]-covSize; j++) chainIter[i]++;
  
	  for( unsigned int n=0; n < covSize; n++) {
	    for (unsigned int k=0; k< PARAMETERS; k++) average[k]+= chainIter[i]->f[k] *chainIter[i]->Multiplicity;
	    TotalSize += chainIter[i]->Multiplicity; // total weight of all points
	    chainIter[i]++; //next point
	  }
  
	  for(unsigned int k=0; k < PARAMETERS; k++){ //normalize
	    average[k]=average[k]/(TotalSize);
	  }
  
	  // reset iterator
	  chainIter[i]=chain[i].begin();
	  for (unsigned int j=0; j < chainSize[i]-covSize; j++) chainIter[i]++;
  
	  for (unsigned int n = 0; n  < covSize ; n++) { // run over all points
	    for (unsigned int k =0; k < PARAMETERS; k++) {
	      for(unsigned int j=0; j<PARAMETERS;j++){
		cov[k][j] += (chainIter[i]->f[k] - average[k])*
		  (chainIter[i]->f[j] - average[j]) * chainIter[i]->Multiplicity;
	      } 
	    }
	    chainIter[i]++;
	  }
  
	  for(unsigned int k=0; k < PARAMETERS; k++){ //normalize
	    for(unsigned int j=0; j < PARAMETERS; j++) {
	      cov[k][j]/=(TotalSize-1.0);
	    }
	  }
  
  
	  for(unsigned int k=0; k  < PARAMETERS; k++) {
	    for (unsigned int j=0; j< PARAMETERS; j++){
	      covMatrix[i][k][j]=cov[k][j];
	    }
	  }
	  mGauss[i]->generateEigenvectors(covMatrix[i],EntireFactor[i]*EntireFactor[i]);
	  mGauss[i]->printExtInfo(progress);
  
  
	  // output information, in text file as well as binary format for re-starting
	  covarianceMatrices << "Chain: " << i+1 << " Step: " << chainSize[i] << " Points used: " 
			     << covSize << endl;
	  ofstream covfile;
	  covfile.open((base+covfile_name + inttochar[i+1]+ ".dat").c_str());
	  for (unsigned int j=0; j < TASKARRAY_SIZE; j++) {
	    AnalyzeThis::write<double>(covfile,chain[i].back().f[j]);
	  }
	  AnalyzeThis::write<int>(covfile,chain[i].back().Multiplicity);
	  AnalyzeThis::write<int>(covfile,chain[i].back().ReallyInvestigated);
	  
	  for(unsigned int j=0; j < PARAMETERS; j++){
	    for (unsigned int k=0; k < PARAMETERS; k++){
	      AnalyzeThis::write<double>(covfile,cov[j][k]);
	      covarianceMatrices << cov[j][k] << "  ";
	    }
	    covarianceMatrices << endl;
	  }
	  AnalyzeThis::write<int>(covfile,MAGIC_CHECK);
	  covfile.close();
	  covarianceMatrices << "Entire factor: " << EntireFactor[i] << endl;
	  covarianceMatrices << "***********************************************************" << endl;
	  
	  steps_since_update[i]=0;
	}
	// END OF ADAPTIVE  COVARIANCE
      
	/*
	** In addition and to help the adaptive covariance,
	** we use an adaptive step size multiplicator, called EntireFactor
	** Whenever we take a step, we consider increasing the step size,
	** as frequent approval corresponds to small stepsizes. 
	** This is performed several lines above when storing the new step
	** in the chain. Here, we consider the opposite. If we hang around
	** for some time at the same spot, we decrease the step size..
	*/

	// ADAPTIVE: decrease if we didn't take for a long time
	if (ADAPTIVE) {
	  if (chain[i].back().ReallyInvestigated > HIGH_STEP_BOUND) {
	    if (EntireFactor[i] > 1e-1) { 
	      EntireFactor[i] *=  DECREASE_STEP;
	      mGauss[i]->generateEigenvectors(covMatrix[i],EntireFactor[i]*EntireFactor[i]);
	    }
	  }
	} // END ADATPIVE decrease

	// now, we re-evaluate 
	// first, let us find out how long the shortest chain is 
	unsigned int min_size=1000*1000;
	for (unsigned int i = 0; i < CHAINS; i++) {
	  if ( chainSize[i] < min_size) min_size = chainSize[i];	
	}
      
	// output of progress information
	progress << "performed/Chainsize: " ;
	for (unsigned int i=0; i<CHAINS; i++){   progress << chainTotalPerformed[i] << " / " << chainSize[i] << "  ";  }
	progress << "min_size: " << min_size   << "  Multiplicity [Really]: ";
	for (unsigned int i=0; i< CHAINS; i++){  progress  << chain[i].back().Multiplicity << " [" << chain[i].back().ReallyInvestigated << "]  ";}
	progress << endl;
      
      
	/*
	** Start with the Gelman and Rubin(1992) statistic;
	** compare variances within the chain with the variances between the chains
	** R[k] should be smaller than 1.1 for convergence
	** If you would like to run a synthetic distribution for checks, this statistics
	** will slow you down. You can speed things up if you replace the if
	** statement below by something like:
	**        if ( (min_size) > RMIN_POINTS && Miscmath::posRnd(1.0) > 0.95) { 
	** hence, only about every 20 times, the statistics is re-calculated.
	*/
	if (min_size > RMIN_POINTS) { 
	  unsigned int break_at = (min_size)/2; 
	  double N = min_size - break_at;
	  vector<list<Task>::iterator > iter(CHAINS); 
	  for (unsigned int i = 0; i < CHAINS; i++) {
	    iter[i] = chain[i].begin();
	    // ignore the first half of the minimal chain size
	    for (unsigned int k =1; k < break_at; k++,iter[i]++);
	  }
	
	  // now let us calculate the mean of each parameter
	  double y[CHAINS][PARAMETERS];
	  double dist_y[PARAMETERS]; // distribution mean
	  double B[PARAMETERS]; // variance between chains
	  double W[PARAMETERS]; // variance within chains
	  double R[PARAMETERS]; // monitoring parameter
	  for (unsigned int k =  0; k < PARAMETERS; k++) { // initialize
	    dist_y[k]= 0; 
	    B[k] = 0;
	    W[k]=0;
	    for (unsigned int i=0; i < CHAINS; i++) y[i][k]=0;
	  }
	  double M = CHAINS;
	  for (unsigned int i = 0; i < CHAINS; i++) {  // for all chains
	    int TotalMultiplicity = 0; 
	    for (unsigned int n = break_at; n  < min_size; n++) {  // traverse through the chains
	      for (unsigned int k = 0; k < PARAMETERS; k++) { // for all parameters
		// add up all parameters	      
		y[i][k] += iter[i]->f[k] *iter[i]->Multiplicity; 
	      }
	      TotalMultiplicity += iter[i]->Multiplicity;
	      iter[i]++;
	    }
	    //progress << "new eval: N: " << N <<  "  check_count: " << check_count << endl;
	    N = (double) TotalMultiplicity;
	    for (unsigned int k =0; k < PARAMETERS; k++){
	      y[i][k] /= (double) N;
	      dist_y[k] += y[i][k] / M;  // four chains hence mean is sum / CHAINS
	    }
	  }	
	
	  // variance between chains
	  for (unsigned int k=0; k < PARAMETERS; k++) {
	    for (unsigned int i = 0; i < CHAINS; i++) {
	      double t1 = y[i][k] - dist_y[k];
	      B[k] +=t1*t1 / (M-1.0);  // as defined in (22) (why 1/(M-1) and not M ???) 
	      // variance within a chain
	    }
	  }
	  // once again, go back to break_at in the chains
	  for (unsigned int i = 0; i < CHAINS; i++) {
	    iter[i] = chain[i].begin();
	    for (unsigned int k =1; k < break_at; k++, iter[i]++);
	  }
	
	  // get W, the variance within chains
	  // run over all points
	  for (unsigned int i = 0; i < CHAINS; i++) {  // run over all chains
	    for (unsigned int n =  break_at; n  < min_size; n++) { // run over all points
	      for (unsigned int k =0; k < PARAMETERS; k++) {
		double t2 = iter[i]->f[k] - y[i][k];
		W[k] += t2*t2*iter[i]->Multiplicity;   
	      }
	      iter[i]++; //next point
	    }
	  }
	  // normalize W and get R
	  for (unsigned int k =0; k < PARAMETERS; k++) {
	    W[k] /= (M*(N-1.0));
	    R[k] = (N-1.0)/N * W[k]  + B[k]*(1.0 + 1.0/N);
	    R[k] /= W[k];
	  }

	  // if we want FREEZE_IN, test if convergence has been reached:
	  if(FREEZE_IN && min_size > MIN_SIZE_FOR_FREEZE_IN ){
	    bool ConvergenceReached=true;
	    for(unsigned int m=0; m<PARAMETERS; m++) {
	      if(R[m] > RBREAK) ConvergenceReached=false;
	    }
	    if(ConvergenceReached) {
	      ADAPTIVE=false; //stop with adaptive stepsize
	      FREEZE_IN=false; //never check again (otherwise no real freeze-in)
	      ofstream final("freezeInEigenvectors.txt");
	      for (unsigned int j=0; j < CHAINS; j++) {
		Task keep = chain[j].back();
		// the next 2 lines write the flag to each montecarlo file 
		// indicating that all models before this have to be discarded
		// the flag is a zero weight 
		
	 
		for(unsigned int k=0; k < TASKARRAY_SIZE; k++) data[j] << keep.f[k] << " ";
		data[j] << 0 << endl;  // zero Multiplicity signals freeze-in
		
		// Now, we clear all chains, i.e. we set discard everything...
		chain[j].clear();
		// and we push back the last one...
		chain[j].push_back(keep);
		chainSize[j]=1;
		// USE Dunkley et. al. result for optimal sigma_t: 2.4 / sqrt(Parameters)
		double OptimalFactor = 2.4 / sqrt((double)PARAMETERS);
		OptimalFactor = roll[j]->average();
		EntireFactor[j]=OptimalFactor;   
		mGauss[j]->generateEigenvectors(covMatrix[j],OptimalFactor*OptimalFactor);
		mGauss[j]->lock();
		final << "::: CHAIN[" << j << "] Eingenvector and value info: " << endl;
		mGauss[j]->printInfo(final);
		final << endl << endl;		
	      }
	      covarianceMatrices << "Stopped adaptive stepsize at " 
				 << min_size << " for all chains (FREEZE_IN=true)" << endl; 
	    }
	  }

	  progress << "Statistics: " << endl;
	  for (unsigned int k = 0; k < PARAMETERS; k++) {
	    progress << k << "  dist_y: " << dist_y[k] << "  B: " << B[k] << "  W: " << W[k] << "    R[k]:  " << R[k] << endl; 
	  }
	  if (FREEZE_IN) progress <<  "Multiplicity (and EntireFactor): "; else 
	    progress << "Multiplicity (and frozen EntireFactor): ";
	  for (unsigned int i = 0; i < CHAINS; i++) progress  << chain[i].back().Multiplicity << " ("<< EntireFactor[i] << ")  ";
	  progress << endl;
	
	  if((min_size)!=checkSize){
	    gelmanRubin << break_at << "  " ;
	    for (unsigned int n=0 ; n < PARAMETERS; n++){
	      gelmanRubin << R[n] << "  ";
	    }
	    gelmanRubin << (min_size) << endl;
	    checkSize=min_size;
	  }
	} // end if min_size > old_size
	progress.flush();
	gelmanRubin.flush();
      } // if TAKE_RESULT
    }  while (work); // end of DO
  } catch (Bad_Error x) {
    cout << "\n\n******** BAD ERROR OCCURED IN MASTER ****\n\n";
    cout << x.s << "\n\n";
    cout << "******************************************\n";
  }
}   


/*! The slaves are responsible for computing  models and its associated likelihoods. They recive the
parameters from the Master, report back to him the results of the computation and request new jobs.  */
void slave(int rank,bool Restart) {
  try {
  bool virgin = true;
  
  CmbCalc cmbcalc;


  // set your gauge here
  //  cmbcalc.setGauge(Gauge::dmdSynchronous);
  cmbcalc.setGauge(Gauge::speedyInvariant);
  AnalyzeThis ai;

  int irt;
  bool accept=true;
  MPI_Status  status;


  string filename;  //! Used frequently for storing the names of various files
  string scalarFileName, tensorFileName, lensedFileName;
  string filejlens; // the filename of the lensing bessels
   
  map<int,string> transferFile;
  double task[MAX_TASKARRAY_SIZE];
  double cpy[MAX_TASKARRAY_SIZE];
  
  /*start an infinite loop;
   basically, a lot of this is just setting the parameters we have received from the master
   and running cmbeasy, obtaining the likelihoods and derived parameters (e.g. sigma8 ) 
   and sending the result to the master
  */ 
  int cnt = 0;
  for (;;) {
    try {
      Cosmos cosmos;    // use QuintCosmos for Quintessence
      ControlPanel control;
      accept=true;
   
      MPI_Send(&task,TASKARRAY_SIZE,MPI_DOUBLE,0,GIMMETASK,MPI_COMM_WORLD); // request a task   
      MPI_Recv(&task, TASKARRAY_SIZE,MPI_DOUBLE,0,TAKETASK,MPI_COMM_WORLD,&status); //master's answer
     
      control.cmb = WMAP7 || BOOMERANG03 || VSA || CBI || ACBAR;   // want cl - spectrum  ?
      control.power_cdm = SDSS || TwodF || SDSSLRG || control.cmb ; // want cdm power-spectrum ? We compute it if cmb required to get sigma_8 in multipurpose

      // always set hubble first, cause setting omega_b h^2 etc will depend on h()
      cosmos.seth(task[2]);   
      // fixed parameters (can be converted to variable parameters of course)
      cosmos.setT_cmb(2.725); //Temperature of the CMB
      cosmos.setY_he(0.24);  // Primordial He-Abundance
      cosmos.setNuR(3.04); // Number of (massless) relativistic neutrios
      cosmos.setNuNR(0);  // Number of non-relativistic neutrinos
      cosmos.setOmega_nuNR(0.00); // Contribution of non-relativistic neutrinos
      /* **********************************************************************************
      ** if you would like to run quintessence instead of a cosmological constant,
      ** replace Comsos by QuintCosmos, uncomment the next two lines
      ** and make sure to modify and uncomment the two lines below
      ** that set the quintessence parameters and call tuneQuintessence()
      ************************************************************************************ */
      // cosmos.setQuintessence(Quintessence::corasaniti);   // the model you like
      // cosmos.setOmega_quintessence(1 - task[0]/cosmos.h2()); // for quintessence 

      cosmos.setOmega_vacuum(1 - task[0]/cosmos.h2()); // for cosmological constant    
      cosmos.setOmegaH2_cdm(task[0] - task[1]);   
      cosmos.setOmegaH2_b(task[1]);
      cosmos.setOptDistanceLss(task[3]); 
     
     
      if (control.power_cdm) { 
	control.highPrecisionTransfer = true;  // if at all cdm, high precision ?
	control.transferMaxK=5*cosmos.h();  // maximal k for cdm  
	control.transferPerLog=5;  // k-values per log interval
	control.transferZ[0] = 0;   // cdm-power output redshift
	transferFile[0] = "trans.dat"; // cdm-power output filename
      }  
      if (control.cmb) {	
	control.scalar = true;  // scalar pert. 
	control.tensor = false; // tensor per.
	// some lines for the spectal indices
	if (control.scalar)  cosmos.setInitialPower(task[4]) ;	

	if (control.scalar && control.tensor) { 
	  //  set tensor spectral index to scalar index - 1
	  for (uint in = 0; in < cosmos.InitialPower.size(); ++in) {
	    cosmos.InitialTensorPower[in] = cosmos.InitialPower[in] - 1.0; 
	  }

	}
	if (control.tensor && (!control.scalar) ) {
	  int n =1;
	  cosmos.InitialTensorPower.clear();
	  for (int in = 0; in < n; ++in) {
	    cosmos.InitialTensorPower[in]=0;
	    //cosmos.TensorRatio[in] = 1.0;
	  }	 	
	}	

	if (control.scalar && control.power_cdm ) {
	  // Enter (0) unlensed Cls only
	  // Enter (1) lensed Cls, linear evolution
	  // Enter (2) lensed Cls, non-linear evolution
	  int lensflag = 0;
	  control.setLensing(lensflag);
	}
	
	if (control.isLensing()   && (!control.power_cdm)) {
	  cout << "You did not request the transfer function" << endl;
	  cout << "calculation needed to do the lensing" << endl;
	  cout << "you will have to start again" << endl;
	  throw Bad_Error("no transfer functions requested, needed for lensing");
	}	  
      
	if (control.scalar) {
	  // Enter output filename for SCALAR cl
	  if (control.isLensing()) {
	    // If lensing was requested this will
	    // be the unlensed power spectrum
	  }
	  scalarFileName = "plot.dat";  // unlensed output file name
 
	  if (control.isLensing()) {
	    // Enter output filename for LENSED SCALAR cl
	    lensedFileName = "lensed.dat"; // output file name for lensed spectra
	    filejlens = "jlens.dat";  // file with bessel functions for lensing
	  }
	}
	if (control.tensor) {
	  // Enter output filename for TENSOR cl
	  tensorFileName="plottensor.dat";
	}	
      } else {
	control.scalar=true; control.tensor=false;
	cosmos.InitialPower.clear();
	control.setLensing(0);
      }
      if  (control.scalar)
	control.setInitialConditions(ControlPanel::adiabatic);  // initial conditions
    
      if (fabs(cosmos.omega_k()) <= .001) { // only flat universes
	if (control.cmb) {
	  string cmbeasydir = ControlPanel::cmbeasyDir();
	  filename = cmbeasydir + "/resources/jlgen.dat";  // bessel function filename
	  if (virgin) cmbcalc.initjl(filename,2000);  // max nr of l's
	  virgin = false;
	}
      } else   throw Bad_Error("cmbeasy only supports flat universes");
   
      
      int n = cosmos.InitialPower.size();
      CL cl;
      cl.resize(n);
      cosmos.reset(); //reset is always good, epsecially, if you run a loop of models
    
      // For quintessence. You might also have to call tuneQuintessnce() depending on the model
      // cosmos.setQParameters(task[6],task[7],task[8],pow(10.0,task[9]));
      // control.setPhantomCrossing(true);  // is the model potentially crossing the w=-1 line ?

      string gaga("42"); // don't ask, just wonder :-) hint: has no purpose whatsoever
      try {
	if (control.cmb || control.power_cdm) cmbcalc.cmbflat(&cosmos,gaga, control, cl);
	else cosmos.history();
      } catch (Bad_Error x) {
	cout << "leandriver bad error:" << endl;
	cout << x.s << endl;
	ofstream ErrorFile;
	ErrorFile.open(errorname.c_str(),ios::app);
	ErrorFile<< "\n\n******** BAD ERROR OCCURED IN SLAVE Nr: " << rank << "  *****\n\n";
	ErrorFile<< x.s << "\n\n";
	ErrorFile << "Parameters have been: " << endl;
	for (unsigned int i =0; i < TASKARRAY_SIZE; i++) {
	  ErrorFile << "task["<<i<<"] = " << task[i] << endl;
	}
	ErrorFile << "******************************************\n";
	accept = false;
	cout << "Slave: " << rank << " error " << endl;
      }	  
      if (accept) {
	if (control.cmb) {
	  cl.ts[0]->setChildrensN(); 	
	  ai.scaleCls(cl,0,pow(2.725e6,2));	// analyzethis needs all  Cl's in units of muK^2	
	  vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
	  ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
	  // For our single spectral index, we choose A_s
	  A_s[0] =  exp(task[5] + 2*task[3])*1e-10;
	  // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
	  ai.applyInflationaryTensorRatio(cosmos,A_s,A_t); 
	  // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
	  // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
	  // automatically]
	  ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);

	  if (control.isLensing()) {		
	    Lensing lens(cosmos,control,cl,cmbcalc);
	    cl = *lens.lensedCls(); // get the lensed cl's
	    cl.ts[0]->arm(Spline::all); // arm the lensed ones
	  }
	}

	/* *************************************************************************************
	** Having computed the model, we can compare it with the data
	** using the rich number of functions available in the analyzeit-class
	************************************************************************************* */
#ifndef NOWMAP7
	if (WMAP7) {
	  AnalyzeThis::WMAP7Likelihood like = ai.WMAP7computeLikelihood(cl);       
	  /* *************************************************************************************
	  ** Here, we note the first  of several likelihoods. Please note, that
	  ** we start at task[LOGLIKEPOS+1], because the total likelihood of all 
	  ** experiments will be in LOGLIKEPOS
	  ************************************************************************************* */
	  task[LOGLIKEPOS+1] = -0.5*like.total;
	  task[LOGLIKEPOS+2] = -0.5*like.TTlike;
	  task[LOGLIKEPOS+3] = -0.5*like.TTlowllike;
	  task[LOGLIKEPOS+4] = -0.5*like.TTlowldet;
	  task[LOGLIKEPOS+5] = -0.5*like.Beamlike;
	  task[LOGLIKEPOS+6] = -0.5*like.TElike;
	  task[LOGLIKEPOS+7] = -0.5*like.TEdet;
	  task[LOGLIKEPOS+8] = -0.5*like.Lowllike;
	  task[LOGLIKEPOS+9] = -0.5*like.Lowldet;
	}
#endif
	
	if (control.cmb) cl.ts[0]->arm(Spline::all);      // arm all output splines
	/* *************************************************************************************
	** begin CBI and ACBAR, VSA,BOOMERANG03. First, we deselect several bands. Particularily
	** The very high l bands. Hence, it sufficies to calculate up to l=2000,
	** as the window functions would give almost zero weight to these high l's
	************************************************************************************* */
	//  CBI and ACBAR select; include only l < 1500 
	for (int i = 1; i < 5; i ++) ai.acbarBandDeselect[i] = true;
	for (int i = 12; i  < 15; i++)  ai.acbarBandDeselect[i] = true;

	for (int i = 1; i  < 11; i++)  ai.VSABandDeselect[i] = true;

	if (CBI) {
		CBI2 cbi2(cl);
		for (int i = 0; i < 5; i++) cbi2.bandDeselect[i] = true;
		for (int i = 11; i <= 15; i++) cbi2.bandDeselect[i] = true;
		task[LOGLIKEPOS+10] = -0.5*cbi2.chi2WithCalibration();
	}

	if (VSA)  task[LOGLIKEPOS+11] = -0.5*ai.VSAChi2WithCalibration(*cl.ts[0]);
	if (ACBAR) task[LOGLIKEPOS+12] = -0.5*ai.ACBARChi2WithCalibration(*cl.ts[0]);

	//B-Mode deselect:
	ai.BOOMERANGDeselect(false, false,false,true);
	
	//for combining with WMAP consider only l_eff > 924
	if (WMAP7) {
	  for ( int i = 1; i <= 17; ++i ) ai.BOOMERANGDeselect( i );
	}
	if (BOOMERANG03) {
	  ai.BOOMERANGInit( cl );
	  task[LOGLIKEPOS+13]  =  -0.5*ai.BOOMERANGChi2WithCalibration();  
	}

	
	/* *************************************************************************************
	** begin 2dF stuff.
	** So if you don't need the power_cdm in a real run, switch it off above !
	************************************************************************************* */
	if (TwodF) {
	  // obtain a powerspectrum at the redshift of the 2df-Survey (z=0.17)
	  Spline *pwr = cosmos.createPower(0,"hallo",cosmos.power_cdm(),0,cosmos.z2tau(0.17));
	  double bestbias = ai.TwoDF_bestBias(pwr);
	  *pwr *= bestbias * bestbias;
	  pwr->arm();
	  double loglike2df = -0.5*ai.TwoDF_convolutePowerSpectrum(pwr);
	  delete pwr; 
	  task[LOGLIKEPOS+14] = loglike2df;
	}
	
	
	/* *************************************************************************************
	** begin SDSS stuff
	** So if you don't need the power_cdm in a real run, switch it off above !
	************************************************************************************* */
	if (SDSS) {
	  // obtain a powerspectrum at the redshift of the SDSS survey (z=0.1)
	  Spline *sdsspwr =  cosmos.createPower(0,"SDSS",cosmos.power_cdm(),0,cosmos.z2tau(0.10));
	  sdsspwr->arm();
	  double sdssbias = ai.SDSS_bestBias(sdsspwr);
	  task[LOGLIKEPOS+15] =  -0.5*ai.SDSS_chiSquared(sdsspwr,sdssbias);
	  delete sdsspwr;
	}
		

	/* *************************************************************************************
	** For the Sne Ia, we first need a Spline containing the luminosity distance.
	** AnalyzThis will then calculate the chi^2 for us.
	** Please note that you have the choice between three different SNe Ia
	** compilations: Tonry et. al. , Knop et al. and Riess et. al etc.
	************************************************************************************* */
	if (Riess06) {
	  Spline lum(100, "Model::lum");
	  for (double z = 0; z <= 2; z+= 0.05) lum.set(z, cosmos.luminosityDistance(z)); 
	  lum.arm();
	  task[LOGLIKEPOS+16] = -0.5 *ai.Sn1aRiess06(lum);
	}

	if (Astier) { 
	  double astier = ai.SNIaAstier05(cosmos,task[AstierParameterPos],task[AstierParameterPos+1]);
	  task[LOGLIKEPOS+17] =  -0.5 *astier;
	}
      
	/*! ***********************************************************************************
	** SDSS Baryon coustic osciallation likelihood
	**  
	** *********************************************************************************** */
	if (SDSS_BAO) {
	  double ns = task[4], ode = 0;
	  // in case you have dark energy uncomment next line 
	  // ode = cosmos.omebar();
	  task[LOGLIKEPOS+18] = -0.5* ai.SDSS_BAP_A_chiSquared(cosmos, task[4],ode);
	}

	/* *************************************************************************************
	** SDSS LRG
	************************************************************************************* */
	if (SDSSLRG) {
	  Spline* cdm=cosmos.createPower(0, "cdm",  cosmos.power_cdm(), 0, cosmos.z2tau(0));       
	  cdm->arm();
	  SdssLrg s(cosmos,ControlPanel::cmbeasyDir("/resources/sdss_lrg/") );
	  double bias = task[SDSSLRGParameterPos];
	  task[LOGLIKEPOS+19] = -0.5*s.chi2(*cdm, bias);
	  delete cdm;
	}

	/* *************************************************************************************
	** Lyman Alpha Pat McDonald
	************************************************************************************* */
	if (LYA_MCDONALD) {
	  double chi2 = ai.lymanAlphaPatMcDonaldChi2(cosmos);	  
	  task[LOGLIKEPOS+20] = -0.5*chi2;
	}


	/* *****************************************************************************
	** Now, we add up all loglike's we would like to include.
	** master() will only compare the value of task[LOGLIKEPOS]
	** for determining the likelihood. 
	** Hence, you may compute as many likelihoods as you like, and
	** store them in task[]. For as long as you don't add a log likelihood
	** to task[LOGLIKEPOS], this will not influence the run of the MCMC
	******************************************************************************* */
	task[LOGLIKEPOS] = 0;
	if (WMAP7) task[LOGLIKEPOS] += task[LOGLIKEPOS+1]; //  wmap 7year
	if (CBI)  task[LOGLIKEPOS] += task[LOGLIKEPOS+10]; 
	if (VSA)  task[LOGLIKEPOS] += task[LOGLIKEPOS+11];  
	if (ACBAR)  task[LOGLIKEPOS] += task[LOGLIKEPOS+12]; 
	if (BOOMERANG03)  task[LOGLIKEPOS] += task[LOGLIKEPOS+13]; 

	if (TwodF)  task[LOGLIKEPOS] += task[LOGLIKEPOS+14 ]; 
	if (SDSS)  task[LOGLIKEPOS] += task[LOGLIKEPOS+15 ];  	
	
	if (Riess06)  task[LOGLIKEPOS] += task[LOGLIKEPOS+16];  
	if (Astier)  task[LOGLIKEPOS] += task[LOGLIKEPOS+17];  
	if (SDSS_BAO)  task[LOGLIKEPOS] += task[LOGLIKEPOS+18];  
	if (SDSSLRG)  task[LOGLIKEPOS] += task[LOGLIKEPOS+19];  
	if (LYA_MCDONALD)  task[LOGLIKEPOS] += task[LOGLIKEPOS+20];  

	/* **********************************************************************
	** As an example, we keep track of sigma8. Hence, you can
	** monitor and plot the expectation value of sigma_8
	** You can add as many multipurpose columns as you like.
	************************************************************************ */

	task[MULTIPURPOSEPOS] =  cosmos.sigma8[0];   // add more stuff if you need it


	/* ***********************************************************************
	** send the whole task - line back to master()+ last minute sanity check 
	************************************************************************ */
	for (unsigned int k = 0;k < TASKARRAY_SIZE; k++)  { 
	  cpy[k] = task[k];
	  if (isnan(task[k])) throw Bad_Error("isnan element");
	}
	MPI_Send(&cpy, TASKARRAY_SIZE, MPI_DOUBLE, 0, TAKERESULT,MPI_COMM_WORLD);  
      } 
    } catch (Bad_Error x) {
      cout << "\n\n******** BAD ERROR OCCURED IN SLAVE *****\n\n";
      cout << x.s << "\n\n";
      cout << "******************************************\n";
      ofstream ErrorFile;
      ErrorFile.open(errorname.c_str(),ios::app);
      ErrorFile<< "\n\n******** BAD ERROR OCCURED IN SLAVE Nr: " << rank << "  *****\n\n";
      ErrorFile<< x.s << "\n\n";
      ErrorFile << "Parameters have been: " << endl;
      for (unsigned int i =0; i < TASKARRAY_SIZE; i++) {
	ErrorFile << "task["<<i<<"] = " << task[i] << endl;
      }
      ErrorFile << "******************************************\n";
    } catch (SafeVectorOutOfRange) {
      cout << "\n\n*** SafeVector out of range *** \n\n";
    }
  }
  } catch (Bad_Error x)  {
    cout << "\n\n******** BAD ERROR LATE OCCURED IN SLAVE *****\n\n";
    cout << x.s << "\n\n";
    cout << "******************************************\n";
  }
}



int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  cout << myrank << endl;
  bool Restart = false;
  if (argc == 2) {
    string arg(argv[1]);
    if (arg == "-restart") Restart = true; 
    else {
      cerr << "only argument accepted is -restart";
      exit(0);
    }
  }
  setVariables();
  if (myrank) slave(myrank,Restart); else master(Restart);
  MPI_Finalize();
  return 0;
}
