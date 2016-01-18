/* Created on: Dec 28, 2015
 *     Author: T.Delame (tdelame@gmail.com)
 */
# include <project.h>
# include <chrono>
# include <random>
BEGIN_PROJECT_NAMESPACE
const std::string severity_names[6] = {
    "[ trace ]",
    "[ debug ]",
    "[ info  ]",
    "[warning]",
    "[ error ]",
    "[ fatal ]"
};

real unit_random()
{
  //fixme:
  static std::default_random_engine generator( 7 );
//  static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count() );
  static std::uniform_real_distribution<real> distribution( real(0), real(1));
  return distribution(generator);
}

END_PROJECT_NAMESPACE
