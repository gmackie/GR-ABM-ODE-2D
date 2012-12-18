#ifndef PARAMS_H
#define PARAMS_H
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include "range.h"
#include <algorithm>

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
    std::string name, desc, units, path;
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

  /* Definitions */
#define P(type, name, path, def, desc, units, min, max) type _##path##_##name; ParamDescriptor<type> _##path##_##name##Desc;
#include "params.def"


  /* Called from constructor */
   void initialize() {
 #define P(type, n, p, defau, d, u, min, max)  \
     if(type(min) != type(max)) \
       _##p##_##n##Desc.range = Range<type>(type(min), type(max));             \
     _##p##_##n##Desc.name = std::string(#n);           \
     _##p##_##n##Desc.def = boost::optional<type>(defau); \
     _##p##_##n##Desc.desc = std::string(d);           \
     _##p##_##n##Desc.units = std::string(u);         \
     _##p##_##n##Desc.path = make_path(#p);           \
     _##p##_##n = !_##p##_##n##Desc.def ? type() : *(_##p##_##n##Desc.def);
 #include "params.def"
   }

public:

  Params() { initialize(); }

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
    if(!!v) val = *v;
    else if(!!desc.def) val = *(desc.def);
    else {
      std::clog<<"Invalid or missing parameter: "<< desc.path << '.' << desc.name << std::endl;
      error = MISSING;
      return;
    }
    if(desc.range && !desc.range->contains(val)) {
      std::clog<<"Parameter outside of range: " << desc.path << '.' << desc.name << std::endl;
      error = RANGE;
    }
  }
  //Writer visiting
  template<typename T>
  void visit(const T& val, const Params::ParamDescriptor<T>& desc) {
    if(desc.def && *(desc.def) == val) return;  //It's the same as the default, skip it
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
}

#endif // PARAMS_H
