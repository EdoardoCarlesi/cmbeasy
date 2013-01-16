// vim:ts=2:sw=2:et
#include "mctaskinfo.h"

#include "mcerror.h"

#include "analyzethis.h"

#include <limits>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <iterator>
#include <fstream>

#include <iostream> //debug only


using namespace std;

vector<string>                McTaskInfo::mNames;
vector<McTaskInfo::EntryType> McTaskInfo::mTypes;
vector<string>                McTaskInfo::mParameterNames;
map<string, McTaskInfo::ParameterInfo> McTaskInfo::mParameterInfos;
map<string, McTaskInfo::ParameterInfo> McTaskInfo::mFixedParameterInfos;
vector<string>                McTaskInfo::mLogLikeNames;
vector<string>                McTaskInfo::mExtraInfoNames;

static const double not_a_number = numeric_limits<double>::quiet_NaN();

static int taskInfoFmtIndex()
{
  static int idx = ios::xalloc();
  return idx;
}

typedef bitset<sizeof(unsigned long int)> FormatFlags;
static const int prettyPrintIdx = 1;
static const int convertNanToZeroIdx = 2;
static const int debugOnIdx = 3;

static void setFlag(ostream& os, int idx, bool b)
{
  FormatFlags flags(os.iword(taskInfoFmtIndex()));
  flags.set(idx, b);
  os.iword(taskInfoFmtIndex()) = flags.to_ulong();
}

static bool formatFlag(ostream& os, int idx)
{
  FormatFlags flags(os.iword(taskInfoFmtIndex()));
  return flags.test(idx);
}

ostream& prettyPrint(ostream& os)
{
  setFlag(os, prettyPrintIdx, true);
  return os;
}

ostream& noPrettyPrint(ostream& os)
{
  setFlag(os, prettyPrintIdx, false);
  return os;
}

ostream& debugOn(ostream& os)
{
  setFlag(os, debugOnIdx, true);
  return os;
}

ostream& debugOff(ostream& os)
{
  setFlag(os, debugOnIdx, false);
  return os;
}

ostream& convertNanToZero(ostream& os)
{
  setFlag(os, convertNanToZeroIdx, true);
  return os;
}

ostream& keepNan(ostream& os)
{
  setFlag(os, convertNanToZeroIdx, false);
  return os;
}

ostream& prettyTaskInfoOutput(ostream& os, const McTaskInfo& t)
{
  struct PrettyPrint {
    PrettyPrint(const vector<string>& v, ostream& os, const McTaskInfo& t) {
      bool convertNans = formatFlag(os, convertNanToZeroIdx);
      vector<string>::const_iterator it, end;
      end = v.end();
      for (it = v.begin(); it != end; ++it) {
        os << "(" << *it << ": ";
        double val = t(*it);
        double val2 = t.value(*it);
        if (convertNans && isnan(val)) {
          os << "0";
        } else {
          os << val;
        }
        os << ")  ";
      }
    }
  };

  PrettyPrint(McTaskInfo::mParameterNames, os, t);
  PrettyPrint(McTaskInfo::mLogLikeNames, os, t);
  PrettyPrint(McTaskInfo::mExtraInfoNames, os, t);

  os << "(Multiplicity: " << t.Multiplicity << ")";

  return os;
}

ostream& debugTaskInfoOutput(ostream& os, const McTaskInfo& t)
{
  os << "\nTaskInfo debug output: ";
  os << t.mEntries.size() << " entries (" << t.mParameterNames.size() << " parameters, ";
  os << t.mLogLikeNames.size() << " loglikes, " << t.mExtraInfoNames.size() << " extra infos)\n";
  os << "parameters are: ";
  copy(t.mParameterNames.begin(), t.mParameterNames.end(), ostream_iterator<string>(os, " "));
  os << "\nloglikes are: ";
  copy(t.mLogLikeNames.begin(), t.mLogLikeNames.end(), ostream_iterator<string>(os, " "));
  os << "\nextra infos are: ";
  copy(t.mExtraInfoNames.begin(), t.mExtraInfoNames.end(), ostream_iterator<string>(os, " "));
  vector<double> pv = t.parameterVector();
  os << "\nparameterVector() returns: ";
  copy(pv.begin(), pv.end(), ostream_iterator<double>(os, " "));
  os << "\nentries as an array of doubles: ";
  double *arr = t.array();
  copy(arr, arr+t.entryCount(), ostream_iterator<double>(os, " "));
  os << "(Multiplicity: " << t.Multiplicity << ")\n";

  return os;
}

ostream& operator<<(ostream& os, const McTaskInfo& t)
{
  if (formatFlag(os, debugOnIdx)) {
    return debugTaskInfoOutput(os, t);
  }

  if (formatFlag(os, prettyPrintIdx)) {
    return prettyTaskInfoOutput(os, t);
  }

  // bare output, no fancy decoration
  struct Print {
    Print(const vector<string>& v, ostream& os, const McTaskInfo& t) {
      bool convertNans = formatFlag(os, convertNanToZeroIdx);
      vector<string>::const_iterator it, end;
      end = v.end();
      for (it = v.begin(); it != end; ++it) {
        double val = t(*it);
        if (convertNans && isnan(val)) {
          os << "0";
        } else {
          os << t(*it);
        }
        os << "    ";
      }
    }
  };

  Print(McTaskInfo::mParameterNames, os, t);
  Print(McTaskInfo::mLogLikeNames, os, t);
  Print(McTaskInfo::mExtraInfoNames, os, t);

  os << t.Multiplicity;

  return os;
}

istream& operator>>(istream& is, McTaskInfo& t)
{
  struct Read {
    Read(const vector<string>& v, istream& is, McTaskInfo& t) {
      vector<string>::const_iterator it, end = v.end();
      for (it = v.begin(); it != end; ++it) {
          if (!is) {
            //throw McError("Couldn't read " + *it);
            return;
          }
          double val;
          is >> val;
          if (is) {
            t(*it) = val;
          }
      }
    }
  };

  Read(McTaskInfo::mParameterNames, is, t);
  Read(McTaskInfo::mLogLikeNames, is, t);
  Read(McTaskInfo::mExtraInfoNames, is, t);

  if (is)
    is >> t.Multiplicity;

  return is;
}

McTaskInfo::McTaskInfo()
              : mArray(0)
{
  vector<string>::iterator it, end = mNames.end();
  for (it = mNames.begin(); it != end; ++it) {
    mEntries.insert(make_pair(*it, not_a_number));
  }

  ParameterInfoMap::const_iterator it2, end2 = mFixedParameterInfos.end();
  for (it2 = mFixedParameterInfos.begin(); it2 != end2; ++it2) {
      mEntries[it2->first] = it2->second.fixedValue;
  }

  Multiplicity = ReallyInvestigated = 1;
  skipFurtherLikelihoodComputations = false;
}

McTaskInfo::~McTaskInfo()
{
  delete[] mArray;
  mArray = 0;
}

McTaskInfo McTaskInfo::fromTextChainFile(const std::string& chainFileName, unsigned int lineNo)
{
    ifstream paramFile(chainFileName.c_str());
    if (!paramFile) {
      throw McError("could not open "+chainFileName);
    }
    unsigned int skippedLines = 0;
    while (++skippedLines < lineNo) {
      paramFile.ignore(65000, '\n');
      if (paramFile.eof() || !paramFile) {
        std::stringstream conv; conv << skippedLines;
        throw McError(" error while skipping line " + conv.str() + " - line number too big?");
      }
    }

    McTaskInfo newInfo;
    paramFile >> newInfo;
    return newInfo;
}

McTaskInfo McTaskInfo::fromBinaryStream(istream& is)
{
  McTaskInfo newInfo;
  for (unsigned int j=0; j<entryCount(); ++j) {
    const double val = AnalyzeThis::read<double>(is);
    newInfo.setEntryValue(j, val);
  }

  newInfo.Multiplicity = AnalyzeThis::read<int>(is);
  newInfo.ReallyInvestigated = AnalyzeThis::read<int>(is);

  return newInfo;
}

McTaskInfo McTaskInfo::fromArray(double* array, int size)
{
  if (size != entryCount()) {
    stringstream errInfo;
    errInfo << "Expected data of size " << entryCount()
            << "but got an array of size " << size << ".";
    throw McError("McTaskInfo::fromArray() - Size mismatch: "
                    + errInfo.str());
  }

  McTaskInfo newInfo;
  for (int i = 0; i < size; ++i) {
    newInfo.setEntryValue(i, array[i]);
  }
  return newInfo;
}

vector<double> McTaskInfo::parameterVector() const
{
  vector<double> params(paramCount());
  for (int i = 0; i < paramCount(); ++i) {
    params[i]=(*mEntries.find(mParameterNames[i])).second;
  }
  return params;
}


double* McTaskInfo::array() const
{
  if (mArray)
    delete[] mArray;
  double* array = new double[entryCount()];
  for (int i=0; i < entryCount(); ++i) {
    array[i] = indexToValue(i);
  }
  return array;
}

double McTaskInfo::indexToValue(const int i) const
{
  pair<string, double> entry;
  entry = *(mEntries.find(mNames[i]));
  return entry.second;
}


void McTaskInfo::addEntry(const EntryType t, const string& name)
{
  if (find(mNames.begin(), mNames.end(), name) != mNames.end()
      || name == "TotalLogLike") {
    throw McError("McTaskInfo::addEntry() - entry with name " + name + "already exists.");
  }

  // insert "TotalLogLike" when inserting the first loglike entry
  if (t==LogLike && mLogLikeNames.size()==0) {
    // not LogLike but ExtraInfo, since we sum all entries with
    // type LogLike to get the total likelihood
    mExtraInfoNames.push_back("TotalLogLike");
    mTypes.push_back(ExtraInfo);
    mNames.push_back("TotalLogLike");
  }

  switch (t) {
    case Parameter: mParameterNames.push_back(name);
                    break;
    case LogLike:   mLogLikeNames.push_back(name);
                    break;
    case ExtraInfo: mExtraInfoNames.push_back(name);
                    break;
    default:        throw McError("McTaskInfo::addEntry() - adding entry with unknown type");
  }

  mNames.push_back(name);
  mTypes.push_back(t);
}


void McTaskInfo::addMcParameter(const string& name, const double lBound,
                                const double uBound, const double sigma)
{
  addEntry(Parameter, name);
  mParameterInfos.insert(make_pair<string, ParameterInfo>(
                         name, ParameterInfo(lBound, uBound, sigma)));
}

void McTaskInfo::fixMcParameter(const std::string& name, const double value)
{
    ParameterInfoMap::iterator it = mParameterInfos.find(name);
  if (it==mParameterInfos.end()) {
    throw McError("McTaskInfo::parameterInfo() - parameter '"
                   + name + "' not found. Forgot to call McTaskInfo::addMcParameter()?");
  }
  mFixedParameterInfos.insert(*it);
  ParameterInfo info = it->second;
  mParameterInfos.erase(it);
  ParameterInfo& param = mFixedParameterInfos[name];
  param.isFixed = true;
  param.fixedValue = value;
  mParameterNames.erase(std::find(mParameterNames.begin(),
                                  mParameterNames.end(), name));

  unsigned int idx = -1;
  for (unsigned int i=0; i < entryCount(); ++i) {
      if (mNames[i]==name) {
          idx = i;
          break;
      }
  }
  mTypes[idx] = ExtraInfo;
  mExtraInfoNames.push_back(name);
}

McTaskInfo::ParameterInfo& McTaskInfo::parameterInfo(const std::string& name)
{
  ParameterInfoMap::iterator it;
  it = mParameterInfos.find(name);
  if (it==mParameterInfos.end()) {
    throw McError("McTaskInfo::parameterInfo() - parameter '"
                   + name + "' not found. Forgot to call McTaskInfo::addMcParameter()?");
  }
  return it->second;
}

void McTaskInfo::setEntryValue(const string& name, const double x)
{
  mEntries[name] = x;
}

void McTaskInfo::setEntryValue(unsigned int index, const double x)
{
  mEntries[mNames[index]] = x;
}

void McTaskInfo::setParameterValues(const vector<double>& values)
{
  for (int i=0; i < values.size(); ++i) {
    mEntries[mParameterNames[i]] = values[i];
  }
}

void McTaskInfo::setParameterValue(unsigned int index, const double x)
{
  mEntries[mParameterNames[index]] = x;
}

double& McTaskInfo::value(const string& name)
{
  map<string, double>::const_iterator it;
  it = mEntries.find(name);
  if (it==mEntries.end()) {
    throw McError("McTaskInfo::value() - no parameter named '"
                   + name
                   + "'; forgot to call McTaskInfo::addParameter()?");
  }
  return mEntries[name];
}

double McTaskInfo::value(const string& name) const
{
  map<string, double>::const_iterator it;
  it = mEntries.find(name);
  if (it==mEntries.end()) {
    throw McError("McTaskInfo::value() - no parameter named '"
                   + name
                   + "'; forgot to call McTaskInfo::addParameter()?");
  }
  return it->second;
}

bool McTaskInfo::hasEntry(const string& name)
{
  vector<string>::const_iterator it;
  it = find(mNames.begin(), mNames.end(), name);
  if (it==mNames.end()) {
    return false;
  }
  return true;
}

vector<double> McTaskInfo::lowerParameterBounds()
{
  if (mParameterInfos.size() <= 0)
    throw McError("McTaskInfo::lowerParameterBounds() - no parameters available.");
  vector<double> lBounds;
  for (int i = 0; i < paramCount(); ++i) {
    lBounds.push_back(mParameterInfos[mParameterNames[i]].lowerBound);
  }
  return lBounds;
}

vector<double> McTaskInfo::upperParameterBounds()
{
  if (mParameterInfos.size() <= 0)
    throw McError("McTaskInfo::upperParameterBounds() - no parameters available.");
  vector<double> uBounds;
  for (int i = 0; i < paramCount(); ++i) {
    uBounds.push_back(mParameterInfos[mParameterNames[i]].upperBound);
  }
  return uBounds;
}

vector<double> McTaskInfo::initialParameterSigmas()
{
  if (mParameterInfos.size() <= 0)
    throw McError("McTaskInfo::initialParameterSigmas() - no parameters available.");
  vector<double> sigmas;
  for (int i = 0; i < paramCount(); ++i) {
    sigmas.push_back(mParameterInfos[mParameterNames[i]].initialSigma);
  }
  return sigmas;
}

double McTaskInfo::totalLogLike() const
{
  vector<string>::const_iterator it, end;
  end = mLogLikeNames.end();
  double total=0;
  for (it = mLogLikeNames.begin(); it != end; ++it) {
    double current = value(*it);
    if (isnan(current)) {
      total = not_a_number;
      break;
    }
    total += current;
  }

  const_cast<McTaskInfo*>(this)->setEntryValue("TotalLogLike", total);
  return total;
}

void McTaskInfo::setErrorLogLike()
{
  vector<string>::const_iterator it, end;
  end = mLogLikeNames.end();
  for (it = mLogLikeNames.begin(); it != end; ++it) {
    setEntryValue(*it, not_a_number);
  }
  setEntryValue("TotalLogLike", not_a_number);
}

void McTaskInfo::writeParameterNamesFile()
{
  const string fileName = "automatic.parameterNames";
  ofstream nameFile(fileName.c_str());
  unsigned int colNr = 0;

  vector<string> names = orderedEntryNames();
  vector<string>::const_iterator it, end = names.end();
  for (it = names.begin(); it != end; ++it) {
    nameFile << colNr++ << "\t" << *it << "\n";
  }
}

