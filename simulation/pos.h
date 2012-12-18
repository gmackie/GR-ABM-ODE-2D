#ifndef POS_H
#define POS_H

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <vector>

#define GETROW(p) (p.x)
#define GETCOL(p) (p.y)

struct Pos;
typedef std::vector<Pos> PosVector;

inline std::ostream& operator<<(std::ostream& s, const Pos& p);
inline std::istream& operator>>(std::istream& s, Pos& p);

struct Pos
{
  union
  {
    struct
    {
      int x, y;
    };
    int v[2];
  };
  Pos(int _x=0, int _y=0)
  {
    x=_x, y=_y;
  }
  const Pos operator+(int c) const
  {
    return Pos(x+c, y+c);
  }
  const Pos operator/(int c) const
  {
    return Pos(x/c, y/c);
  }
  bool operator==(const Pos& p) const
  {
    return !memcmp(v, p.v, sizeof(int)*2);
  }
  Pos(const boost::property_tree::ptree& p, const Pos& center)
  {
    boost::optional<std::string> pos = p.get_optional<std::string>("pos");
    if(pos)
      {
        if(*pos == std::string("center"))
          {
            *this = center;
          }
        else
          {
            std::istringstream s(*pos);
            s >> *this;
          }
      }
    else
      {
        boost::optional<int> _x = p.get_optional<int>("row");
        boost::optional<int> _y = p.get_optional<int>("col");
        x = _x ? *_x : 0;
        y = _y ? *_y : 0;
      }
    boost::optional<std::string> offset = p.get_optional<std::string>("offset");
    if(offset)
      {
        std::stringstream s(*offset);
        Pos p;
        s >> p;
        x += p.x;
        y += p.y;
      }
  }

  double distance(const Pos& p) const
  {
    double deltaX = x - p.x;
    double deltaY = y - p.y;
    double dist = sqrt(deltaX * deltaX + deltaY * deltaY);
    return dist;
  }

  /**
  * @copydoc GrSimulation::serialize
  */
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int /*version*/) {
    ar & BOOST_SERIALIZATION_NVP(x);
    ar & BOOST_SERIALIZATION_NVP(y);
  }

};
inline std::ostream& operator<<(std::ostream& s, const Pos& p)
{
  return s<< '(' << p.x << ", " << p.y <<')';
}
inline std::istream& operator>>(std::istream& s, Pos& p)
{
  unsigned char tmp;
  s>>tmp;
  assert(tmp == '(');
  s>>p.x;
  s>>tmp;
  assert(tmp == ',');
  s>>p.y;
  s>>tmp;
  assert(tmp == ')');
  return s;
}

#endif //POS_H

