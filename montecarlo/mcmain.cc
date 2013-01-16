// vim:ts=2:sw=2:et

#include "mcerror.h"
#include "mcrunner.h"
#include "nestedsampler.h"
#include "mclikelihoodcalculator.h"
#include "mcsettings.h"
#include "global.h" // for Bad_Error

#include "mpi.h"

#include <iostream>
#include <string>
#include <sstream>

struct McMPIError
{
  McMPIError(int c, std::string t): errorCode(c), errorText(t) {}

  std::string errorText;
  int         errorCode;
};

void mpiErrorhandler(MPI::Comm& comm, int* errorcode, ...)
{
  char* errCStr = new char[MPI_MAX_ERROR_STRING];
  int errorCStringLen;
  MPI::Get_error_string(*errorcode, errCStr, errorCStringLen);
  std::string errStr(errCStr);
  delete[] errCStr;
  throw McMPIError(*errorcode, errStr);
}

void runMonteCarlo(int argc, char* argv[], bool restart)
{
  int rank = -1;
  try {
    MPI::Init(argc, argv);
    // we don't use c++ exceptions by default because many MPI installations
    // apparently do not have exception support compiled in, even
    // though they have the C++ bindings
    // MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
    MPI::COMM_WORLD.Set_errhandler(MPI::Comm::Create_errhandler(mpiErrorhandler));

    McSettings *cfg = McSettings::self();
    cfg->setNumberOfChains(MPI::COMM_WORLD.Get_size()-1);

    rank = MPI::COMM_WORLD.Get_rank();

    McRunner* runner;
    McRunner::globalInit();

    if (rank==0) {
      runner = &cfg->master();
    } else {
      runner = &cfg->slave();
    }

    struct Cleanup {
      ~Cleanup() { delete delRunner; }
      McRunner* delRunner;
    } finally;
    finally.delRunner=runner;

    runner->setRank(rank);
    runner->start(restart);

  } catch (MPI::Exception e) {
    std::cerr << "From rank: " << rank << ": MPI error: "
      << e.Get_error_code() << " - "
      << e.Get_error_string() << std::endl;
    MPI::COMM_WORLD.Abort(-1);
  } catch (McMPIError e) {
    std::cerr << "From rank: " << rank << ": McMPI error: "
      << e.errorCode << " - "
      << e.errorText << std::endl;
    MPI::COMM_WORLD.Abort(-1);
  } catch (McError e) {
    std::cerr << "MonteCarlo error from rank: "
      << rank << ":  error: "
      << e.errorText << std::endl;
    MPI::COMM_WORLD.Abort(-1);
  } catch (Bad_Error e) {
    std::cerr << "uncaught Bad_Error error from rank: "
              << rank << ":  error: "
              << e.s << std::endl;
    MPI::COMM_WORLD.Abort(-1);
  } catch (std::exception e) {
    std::cerr << "an error occurred in rank " << rank << ":\n"
              << e.what() << "\n exiting..." << endl;
  }

  MPI::Finalize();
}

void produceNestedSamplerPosterior()
{
  try {
    McSettings *cfg = McSettings::self();
    cfg->setNumberOfChains(0);

    McRunner::globalInit();
    McRunner* runner = &cfg->master();
    NestedSampler* nestedSampler = dynamic_cast<NestedSampler*>(runner);
    if (!nestedSampler) {
      cout << "'-posterior' only makes sense when using nested sampling.\n";
      return;
    }

    struct Cleanup {
      ~Cleanup() { delete delRunner; }
      McRunner* delRunner;
    } finally;
    finally.delRunner=runner;

    nestedSampler->produceCurrentPosterior();

  } catch (Bad_Error e) {
    std::cerr << "caught Bad_Error error while checking a single model: "
              << ":  error: " << e.s << std::endl;
  } catch (McError e) {
    std::cerr << "MonteCarlo error: " << e.errorText << std::endl;
  }
}

void runSingleModel(const std::string& chainFile, const std::string& lineNoStr)
{
  try {
    McRunner::globalInit();
    McModel& model = McSettings::self()->model();
    McLikelihoodCalculator calcLike;

    McTaskInfo params;
    if (lineNoStr=="binary") { // reading from binary file, contains only one model
      ifstream in(chainFile.c_str());
      if (!in) {
        std::cerr << "could not open " << chainFile << endl;
        return;
      }
      params = McTaskInfo::fromBinaryStream(in);
    } else {
      std::stringstream conv(lineNoStr);
      unsigned int line; conv >> line;
      params = McTaskInfo::fromTextChainFile(chainFile, line);
    }

    std::cout << prettyPrint << "computing model for the following parameters: "
              << endl  << params << endl;;
    model.setParameters(params);
    model.compute();
    model.outputDebugInfo();
    calcLike.computeLogLike(model, params);
    std::cout << prettyPrint << "result is: " << endl << params << endl;;

  } catch (Bad_Error e) {
    std::cerr << "caught Bad_Error error while checking a single model: "
              << ":  error: " << e.s << std::endl;
  } catch (McError e) {
    std::cerr << "MonteCarlo error: " << e.errorText << std::endl;
  }
}

void usage(const std::string& progName)
{
  std::cerr << "unrecognized arguments. Usage is one of:\n";
  std::cerr << progName << "                                     #starts a new monte carlo run\n";
  std::cerr << progName << " -restart                            #resumes a previously halted monte carlo run\n";
  std::cerr << progName << " -posterior                          #only for nested sampling: manually produce current posterior file\n";
  std::cerr << progName << " -singlemodel chainfile lineNo'  #computes only the model in textfile 'chainfile' at line 'lineNo',";
  std::cerr <<" then exits\n";
  std::cerr << progName << " -lastmodel  filename                #computes the model in the binary file 'filename',";
  std::cerr << " which contains the last model of a chain" << std::endl;
}

int main(int argc, char* argv[])
{
  using namespace std;

  bool restart = false;
  vector<string> args(argv, argv+argc);

  if (argc==1) {
    runMonteCarlo(argc, argv, restart);
  } else if (argc==2 && args[1]=="-restart") {
    restart = true;
    runMonteCarlo(argc, argv, restart);
  } else if (argc==2 && args[1]=="-posterior") {
    produceNestedSamplerPosterior();
  } else if (argc==3 && args[1]=="-lastmodel") {
    runSingleModel(args[2], "binary");
  } else if (argc==4 && args[1]=="-singlemodel") {
    runSingleModel(args[2], args[3]);
  } else {
    usage(args[0]);
  }

  return 0;
}
