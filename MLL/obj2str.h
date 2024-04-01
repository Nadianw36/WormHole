#ifndef OBJ2STR_H
#define OBJ2STR_H

#include <sstream>
#include <vector>
#include <map>
#include <set>

template <typename T, typename checkNum = std::enable_if_t<std::is_arithmetic_v<T>>>
std::string str(const T &obj)
{
  return std::to_string(obj);
}

template <typename T1, typename T2>
std::string str(const std::pair<T1, T2> &pair);

template <typename T1, typename T2>
std::string str(const std::map<T1, T2> &mp);

template <typename T>
std::string str(const std::set<T> &st);

template <typename T>
std::string str(const std::vector<T> &arr)
{
  std::string s = "[";
  for (size_t i = 0; i < arr.size(); ++i)
  {
    s.append(str(arr[i]));
    if (i < arr.size() - 1)
      s.append(",");
  }
  s.append("]");
  return s;
}

template <
    typename C,
    typename T = std::decay_t<decltype(*begin(std::declval<C>()))>>
auto str(C const &container)
{
  std::string s = "{";
  auto it = container.begin();
  if (it != container.end())
  {
    s.append(str(*it));
    it++;
  }
  while (it != container.end())
  {
    s.append(",");
    s.append(str(*it));
    it++;
  }
  s.append("}");
  return s;
}

template <typename T1, typename T2>
std::string str(const std::pair<T1, T2> &pair)
{
  return str(pair.first).append(":").append(str(pair.second));
}

template <template <typename...> class Map, typename T1, typename T2>
std::string str(const Map<T1, T2> &mp)
{
  std::string s = "{";
  if (!mp.empty())
  {
    auto it = mp.cbegin();
    s.append(str(*it));
    ++it;
    for (; it != mp.cend(); ++it)
    {
      s.append(",").append(str(*it));
    }
  }
  s.append("}");

  return s;
}

template <typename T>
std::string str(const std::set<T> &st)
{
  std::string s = "{";
  auto it = st.cbegin();
  s.append(str(*it));
  ++it;
  for (; it != st.cend(); ++it)
  {
    s.append(",").append(str(*it));
  }
  s.append("}");

  return s;
}

// template<>
// std::string str(const std::string &obj) {
//     std::ostringstream oss;
//     oss << obj;
//     return oss.str();
// }
#endif // OBJ2STR_H
