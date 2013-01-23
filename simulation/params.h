#ifndef PARAMS_H
#define PARAMS_H
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include "range.h"
#include <algorithm>
#include "pos.h"

#ifndef uint
typedef unsigned int uint;
#endif

struct ParamFileHandler;

class Params {
public:

  typedef boost::property_tree::ptree ptree;
  typedef struct ParamFileHandler ParamFileHandler;

  /// Auxilary template class in case we need to deal with container classes
  /// Some parts may need to be specialized in these cases
  template<typename T>
  struct is_container {   //SFINAE container test
    typedef char yes[1];  //sizeof() == 1
    typedef char no[2];   //sizeof() == 2
    template<typename C>
    static yes& test(typename C::iterator*);
    template<typename>
    static no& test(void*);

    enum { value = sizeof(test<T>(0)) == sizeof(yes) };
  };

  /// Placeholder for description information of parameter
  template<typename T>
  struct ParamDescriptor {
    boost::optional<Range<T> > range; /// Valid range for the parameter
    boost::optional<T> def;           /// Default value for parameter
    bool wasRead; // Whether this parameter had its value read from a parmeter file or not.
    std::string name, desc, units, path;
    ParamDescriptor(T min, T max, std::string n, boost::optional<T> defau, std::string d, std::string u, std::string p)
    : def(defau), wasRead(false), name(n), desc(d), units(u), path(p) { if(min != max) range = Range<T>(min, max); }
  };

private:

  /// Convert token string to ptree path (_ -> .)
  static const std::string make_path(const char path[]) {
    std::string s(path);
    std::replace(s.begin(), s.end(), '_', '.');
    return s;
  }
  /// Runs through all the parameters and checks they are in an appropriate range
  struct ValidatorVisitor {
    bool good;
    ValidatorVisitor() : good(true) {}
    template<typename T>
    void visit(const T& val, const Params::ParamDescriptor<T>& desc) {
      bool res = !(!!desc.range && !desc.range->contains(val));
      if(!res) std::clog<<"Invalid parameter for "<<desc.path<<'.'<<desc.name<<std::endl;
      good &= res;
    }
  };

  // Define any computed parameters including any that are computed based on other parameters.
  void computeParams();

  /* Definitions */

#define P(type, name, path, def, desc, units, min, max) type _##path##_##name; ParamDescriptor<type> _##path##_##name##Desc;
#include "params.def"

  // Needed for some parameters computed based on other parameter values and the simulation dimensions.
  // For example, source count or initial mac count based on density values.
  Pos _dim; 


public:
  Params(const Pos& dim=Pos(-1,-1));

  /// @name Accessors
  /// @{
#define P(type, name, path, def, desc, units, min, max) \
  const type& get##path##_##name() const { return _##path##_##name; } \
  const ParamDescriptor<type>& get##path##_##name##Desc() const { return _##path##_##name##Desc; } \
  type& get##path##_##name() { return _##path##_##name; } \
  void set##path##_##name(const type& v) { _##path##_##name = v; }
#include "params.def"
  /// @}

  /* Visitors */
  template<typename Visitor>
  void visitProperties(Visitor& visitor) const {
    #define P(type, name, path, def, desc, units, min, max) \
      visitor.visit(_##path##_##name, _##path##_##name##Desc);
    #include "params.def"
  }
  template<typename Visitor>
  void visitProperties(Visitor& visitor) {
    #define P(type, name, path, def, desc, units, min, max) \
      visitor.visit(_##path##_##name, _##path##_##name##Desc);
    #include "params.def"
  }
  void save(std::ostream& file, ParamFileHandler* writer, boost::property_tree::ptree& pt) const;
  void save(std::ostream& file, ParamFileHandler* writer) const {
    boost::property_tree::ptree pt;
    save(file, writer, pt);
  }
  void load(std::istream& file, ParamFileHandler* reader, boost::property_tree::ptree& pt);
  void load(std::istream& file, ParamFileHandler* reader) {
    boost::property_tree::ptree pt;
    load(file, reader, pt);
  }
  bool validate() const {
    ValidatorVisitor visitor;
    visitProperties(visitor);
    return visitor.good;
  }
};

/// Base class interface for reading file types and loading parameters
class ParamFileHandler {
protected:
  typedef boost::property_tree::ptree ptree;
private:
  ptree* pt;
  enum ERROR { NONE, RANGE, MISSING };
  ERROR error;
  const std::string root;
public:
  explicit ParamFileHandler(const std::string& _root) : root(_root) {}
  virtual ~ParamFileHandler() {}
  virtual void read(std::istream& file, ptree& pt) const = 0;
  virtual void write(std::ostream& file, const ptree& pt) const = 0;
  std::string getFullPath(const std::string& name, const std::string& path = "") const
    { return getFullPath(root, name, path); }
  virtual std::string getFullPath(const std::string& root, const std::string& name, const std::string& path) const
    {  return getPath(root, name, path); }
  static std::string getPath(const std::string& root, const std::string& name, const std::string& path)
    {  return root + '.' + (path.empty() ? std::string() : path + '.') + name; }
  bool good() { return error == NONE; }
  void visit(const Params& params, ptree& _pt) {
    error = NONE;
    pt = &_pt;
    params.visitProperties(*this);
  }
  void visit(Params& params, ptree& _pt) {
    error = NONE;
    pt = &_pt;
    params.visitProperties(*this);
  }
  //Reader visiting
  template<typename T>
  void visit(T& val, Params::ParamDescriptor<T>& desc) {

    //Read in
    using namespace boost::property_tree;
    boost::optional<T> v = pt->get_optional<T>(getFullPath(root, desc.name, desc.path));

    if(!!v)
    {
      val = *v;
      desc.wasRead = true;
    }
    else if(!!desc.def)
    {
      val = *(desc.def);
    }
    else {
      std::cerr<<"Invalid or missing parameter: "<< desc.path << '.' << desc.name << std::endl;
      error = MISSING;
      return;
    }

    if(desc.range && !desc.range->contains(val)) {
      std::cerr<<"Parameter outside of range: " << desc.path << '.' << desc.name << std::endl;
      error = RANGE;
    }
  }
  //Writer visiting
  template<typename T>
  void visit(const T& val, const Params::ParamDescriptor<T>& desc) {
    if(!desc.wasRead) return;  //Not read from a parameter file, skip it
    pt->put(getFullPath(root, desc.name, desc.path), val);
  }
};

inline void Params::save(std::ostream& file, ParamFileHandler* writer, boost::property_tree::ptree& pt) const {
  assert(writer!=NULL);
  writer->visit(*this, pt);
  writer->write(file, pt);
}

inline void Params::load(std::istream& file, ParamFileHandler* reader, boost::property_tree::ptree& pt) {
 assert(reader!=NULL);
 reader->read(file, pt);
 reader->visit(*this, pt);
 computeParams();
}

#endif // PARAMS_H
