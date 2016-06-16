#include <map>
#include <boost/variant.hpp>


namespace types {
  typedef std::map<std::string, boost::variant<double, int, bool, std::string> > SettingsMap;
}
