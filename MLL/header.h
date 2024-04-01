#ifndef LCTD_HEADERS_HEADER_H_
#define LCTD_HEADERS_HEADER_H_
#include <ctime>

#include <algorithm>
#include <cassert>
#include <ctime>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

// double get_current_time() {
//   timeval t;
//   gettimeofday(&t, 0);
//   return (double)t.tv_sec + (double)t.tv_usec/1000000;
// }

/**
 * @return current time in seconds, cpu clock
 */
[[maybe_unused]] inline double cpu_time()
{
  return clock() / static_cast<double>(CLOCKS_PER_SEC);
}

inline double wall_time()
{
  timespec ts{};
  timespec_get(&ts, TIME_UTC);
  return ts.tv_sec + ts.tv_nsec / 1.0e9;
}

#define get_current_time wall_time

void print_trace(const std::string &msg)
{ // compile with -rddynamic too see the demangled symbols
  char **strings;
  size_t i, size;
  enum Constexpr
  {
    MAX_SIZE = 1024
  };
  void *array[MAX_SIZE];
  size = backtrace(array, MAX_SIZE);
  strings = backtrace_symbols(array, size);
  for (i = 0; i < size; i++)
    printf("%s. msg:%s\n", strings[i], msg.c_str());
  puts("");
  free(strings);
}

inline std::string now()
{
  time_t t = time(nullptr);
  char tmp[64];
  strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&t));
  return std::string(tmp);
}

inline std::string getfilename(std::string path)
{
  size_t pos = path.rfind('/');
  if (pos == std::string::npos)
  {
    return path;
  }
  else
  {
    return path.substr(pos + 1, path.size() - pos);
  }
}

//=========================tools for parse command line=======================
inline char *getCmdOption(int argc, char **argv, const std::string &option)
{
  char **itr = std::find(argv, argv + argc, option);
  if (itr != (argv + argc) && ++itr != (argv + argc))
  {
    return *itr;
  }
  return nullptr;
}
inline bool cmdOptionExists(int argc, char **argv, const std::string &option)
{
  return std::find(argv, argv + argc, option) != (argv + argc);
}

template <typename T, typename T2>
inline void cmd_opt_set(int argc, char **argv, const std::string &opt, T &var, const T2 &dft)
{
  char **itr = std::find(argv, argv + argc, opt);
  static_assert(std::is_convertible<T2, T>::value, "Cannot convert from T2 to T");
  if (itr != (argv + argc) && ++itr != (argv + argc))
  {
    if constexpr (std::is_integral<T>::value)
      var = std::stoi(*itr);
    else if constexpr (std::is_floating_point<T>::value)
      var = std::stof(*itr);
    else if constexpr (std::is_same<T, std::string>::value)
      var = *itr;
  }
  else
    var = dft; // or var=T(dft);
}

#define CMTOPTSET_MACROF(opt, var, dft) cmd_opt_set(argc, argv, opt, var, dft)
#define CMTOPTFIND_MACROF(opt) cmdOptionExists(argc, argv, opt)

#define mlogn(format, ...) \
  printf("[%s %s:%d] " format, now().c_str(), getfilename(__FILE__).c_str(), __LINE__, ##__VA_ARGS__)
#define mlog(format, ...) mlogn(format "\n", ##__VA_ARGS__)

//#define flogn(format, ...)                                                                                             \
//  fmt::print("[{} {}:{}] " format, now().c_str(), getfilename(__FILE__).c_str(), __LINE__, ##__VA_ARGS__)
// #define flog(format, ...) flogn(format "\n", ##__VA_ARGS__)

//==================random generators===================

/**
 * generate type T values uniformly in [min,max]
 * @tparam T
 * @tparam URBG
 */
template <typename T = int, typename URBG = std::minstd_rand>
class Random
{
  using DT = std::conditional_t<std::is_integral<T>::value, std::uniform_int_distribution<T>,
                                std::uniform_real_distribution<T>>;
  using paramT = typename DT::param_type;

  URBG _rng;
  DT _dist;

public:
  inline Random() : Random(T(std::is_integral<T>::value ? std::numeric_limits<T>::max() : T(1))){};
  explicit Random(T max) : Random(T(0), max)
  {
  }
  inline Random(T min, T max) : _rng(0), _dist(min, max)
  {
  }

  inline void seed(unsigned int s)
  {
    _rng.seed(s);
  }
  [[maybe_unused]] inline void seedr()
  {
    seed(std::random_device()());
  }

  inline void set(T max)
  {
    set(T(0), max);
  }
  inline void set(T min, T max)
  {
    _dist = DT(min, max);
  }

  inline T min()
  {
    return _dist.min();
  }
  inline T max()
  {
    return _dist.max();
  }
  inline T operator()()
  {
    return _dist(_rng);
  }
  inline T operator()(T max)
  {
    return operator()(0, max);
  }
  inline T operator()(T min, T max)
  {
    return _dist(_rng, paramT(min, max));
  }
};

[[deprecated("Use Random class instead")]] unsigned long rand_long()
{
  static std::mt19937_64 gen(0);
  return gen();
}

template <typename Container, typename RandomGenerator>
inline typename Container::iterator select_randomly(Container &c, RandomGenerator &g)
{
  auto it = c.begin();
  std::advance(it, g(0, std::distance(c.begin(), c.end()) - 1));
  return it;
}

#endif // LCTD_HEADERS_HEADER_H_
