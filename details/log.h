/*  Created on: Jan 2, 2016
 *      Author: T. Delame (tdelame@gmail.com)
 */

# ifndef KMEAN_LOG_H_
# define KMEAN_LOG_H_

# include "namespace.h"
# include <ctime>
# include <iostream>

BEGIN_PROJECT_NAMESPACE 

  enum severity_level {
    trace, debug, info, warning, error, fatal
  };
  extern const std::string severity_names[6];
  # define LOG( level, message )                                              \
  {                                                                           \
    if( severity_level::level > severity_level::warning )                     \
    {                                                                         \
      auto time = std::time( nullptr );                                       \
      char mbstr[ 30 ];                                                       \
      std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&time));       \
      std::cout<< '[' << mbstr << ']'                                         \
               << severity_names[severity_level::level]<< " "                 \
               << message                                                     \
               << " (in "<< __FILE__<<":"<< __LINE__ <<")"                    \
               << std::endl;                                                  \
    }                                                                         \
    else                                                                      \
    {                                                                         \
      auto time = std::time( nullptr );                                       \
      char mbstr[ 30 ];                                                       \
      std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&time));       \
      std::cout<< '[' << mbstr << ']'                                         \
               << severity_names[severity_level::level]<< " "                 \
               << message << std::endl;                                       \
    }                                                                         \
  }
END_PROJECT_NAMESPACE

# endif
