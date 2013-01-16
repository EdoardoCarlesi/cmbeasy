#ifndef MCTASKINFO_H
#define MCTASKINFO_H
// vim:ts=2:sw=2:et

#include <string>
#include <vector>
#include <map>
#include <ostream>
#include <istream>
#include <limits>

class McTaskInfo
{
  public:
    enum EntryType { Parameter, LogLike, ExtraInfo };
    struct ParameterInfo {
      ParameterInfo() : lowerBound(-1e100), upperBound(-1e100),
                        initialSigma(-1e100),
                        isFixed(false),
                        fixedValue(std::numeric_limits<double>::quiet_NaN()),
                        isFast(false)
      {}

      ParameterInfo(double lBound, double uBound, double sigma)
                      : lowerBound(lBound), upperBound(uBound),
                        initialSigma(sigma),
                        isFixed(false),
                        fixedValue(std::numeric_limits<double>::quiet_NaN()),
                        isFast(false)
      {}

      double lowerBound, upperBound;
      double initialSigma;
      bool   isFixed;
      double fixedValue;
      bool   isFast;
    };
    typedef  std::map<std::string, ParameterInfo> ParameterInfoMap;

    McTaskInfo();
    ~McTaskInfo();

    static McTaskInfo fromTextChainFile(const std::string& chainFileName, unsigned int lineNo);
    static McTaskInfo fromBinaryStream(std::istream& is);

    static McTaskInfo fromArray(double* array, int size);
    double* array() const;

    //double operator() (const unsigned int i) const;
    double& operator() (const std::string& name) { return value(name); }
    double  operator() (const std::string& name) const { return value(name); }
    std::vector<double> parameterVector() const;

    void setEntryValue(const std::string& name, const double x);
    void setEntryValue(const unsigned int index, const double x);

    void   setParameterValues(const std::vector<double>& values);
    void   setParameterValue(unsigned int index, const double x);
    double& value(const std::string& name);
    double value(const std::string& name) const;

    double totalLogLike() const;
    void   setErrorLogLike();

    static void addEntry(const EntryType t, const std::string& name);
    static void addMcParameter(const std::string& name, const double lBound,
                               const double uBound, const double sigma);
    static void fixMcParameter(const std::string& name, const double value);
    static bool hasEntry(const std::string& name);
    static unsigned int entryCount() { return mNames.size(); }
    static unsigned int paramCount() { return mParameterNames.size(); }
    static unsigned int logLikeCount() { return mLogLikeNames.size(); }
    static unsigned int extraInfoCount() { return mExtraInfoNames.size(); }
    static std::string  parameterName(int i)  { return mParameterNames[i]; }
    static void         writeParameterNamesFile();

    static std::vector<std::string> orderedEntryNames()  {
              std::vector<std::string> v(mParameterNames.begin(), mParameterNames.end());
              v.insert(v.end(), mLogLikeNames.begin(), mLogLikeNames.end());
              v.insert(v.end(), mExtraInfoNames.begin(), mExtraInfoNames.end());
              return v;
           };

    static std::vector<double> lowerParameterBounds();
    static std::vector<double> upperParameterBounds();
    static std::vector<double> initialParameterSigmas();
    static ParameterInfoMap&   parameterInfos() { return mParameterInfos; }
    static ParameterInfo&      parameterInfo(const std::string& name);


    unsigned int Multiplicity;
    unsigned int ReallyInvestigated;

    bool skipFurtherLikelihoodComputations;

    friend std::ostream& operator<<(std::ostream& os, const McTaskInfo& t);
    friend std::istream& operator>>(std::istream& is, McTaskInfo& t);
    friend std::ostream& prettyTaskInfoOutput(std::ostream& os, const McTaskInfo& t);
    friend std::ostream& debugTaskInfoOutput(std::ostream& os, const McTaskInfo& t);
  private:
    double indexToValue(const int i) const;

  private:
    mutable double                   *mArray;
    std::map<std::string, double>     mEntries; // name of entry -> value

    static std::vector<std::string>   mNames;   // name of entry <- unique index number
    static std::vector<EntryType>     mTypes;   // index -> type
    static std::vector<std::string>   mParameterNames;   // name of entry <- unique index number
    static ParameterInfoMap           mParameterInfos;
    static ParameterInfoMap           mFixedParameterInfos;
    static std::vector<std::string>   mLogLikeNames;   // name of entry <- unique index number
    static std::vector<std::string>   mExtraInfoNames;   // name of entry <- unique index number
};

std::ostream& prettyPrint(std::ostream& os);
std::ostream& noPrettyPrint(std::ostream& os);
std::ostream& convertNanToZero(std::ostream& os);
std::ostream& keepNan(std::ostream& os);
std::ostream& debugOn(std::ostream& os);
std::ostream& debugOff(std::ostream& os);
std::ostream& prettyTaskInfoOutput(std::ostream& os, const McTaskInfo& t);
std::ostream& debugTaskInfoOutput(std::ostream& os, const McTaskInfo& t);
std::ostream& operator<<(std::ostream& os, const McTaskInfo& t);
std::istream& operator>>(std::istream& is, McTaskInfo& t);

#endif // MCTASKINFO_H
