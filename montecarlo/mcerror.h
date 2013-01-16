#ifndef MCERROR_H
#define MCERROR_H

#include <string>

struct McError
{
  McError(std::string t): errorText(t) {}

  std::string errorText;
};

#endif // MCERROR_H
