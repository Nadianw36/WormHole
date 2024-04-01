//
// Created by 24148 on 8/9/2021.
//

#ifndef PRUNEDTD_UTILS_FILE_IO_H_
#define PRUNEDTD_UTILS_FILE_IO_H_
#include <string>
#include <type_traits>
using std::string;

FILE *write_file(const string &path, bool exit_on_failure = true)
{
  FILE *ofile = fopen(path.c_str(), "w");
  if (ofile == nullptr)
  {
    fprintf(stderr, "open file %s failed.\n", path.c_str());
    if (exit_on_failure)
    {
      exit(EXIT_FAILURE);
    }
  }
  return ofile;
}

FILE *read_file(const string &path, bool exit_on_failure = true)
{
  FILE *ifile = fopen(path.c_str(), "r");
  if (ifile == nullptr)
  {
    fprintf(stderr, "open file %s failed.\n", path.c_str());
    if (exit_on_failure)
    {
      exit(EXIT_FAILURE);
    }
  }
  return ifile;
}

template <typename T>
void write(FILE *ofile, const T &v)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  fwrite(&v, sizeof(T), 1, ofile);
}

template <typename T>
void read(FILE *ifile, T &v)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  fread(&v, sizeof(T), 1, ifile);
}

template <typename T>
void write(FILE *ofile, const std::vector<T> &arr)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  int size = arr.size(); // size should be unsigned long
  write(ofile, size);
  fwrite(arr.data(), sizeof(T), size, ofile);
}

template <typename T>
void read(FILE *ifile, std::vector<T> &arr)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  int size; // size should be unsigned long
  read(ifile, size);
  arr.resize(size);
  fread(arr.data(), sizeof(T), size, ifile);
}

template <typename T>
void write(FILE *ofile, const std::vector<std::vector<T>> &arr)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  int size = arr.size(); // size should be unsigned long
  write(ofile, size);
  for (int i = 0; i < size; ++i)
  {
    write(ofile, arr[i]);
  }
}

template <typename T>
void read(FILE *ifile, std::vector<std::vector<T>> &arr)
{
  static_assert(std::is_trivially_copyable<T>::value || std::is_trivially_copy_constructible<T>::value);
  int size; // size should be unsigned long
  read(ifile, size);
  arr.resize(size);
  for (int i = 0; i < size; ++i)
  {
    read(ifile, arr[i]);
  }
}

void check_eof(FILE *fin)
{
  if (feof(fin))
  {
    mlog("read exceed the file size");
    fclose(fin);
  }
  else if (fgetc(fin) != EOF)
  {
    mlog("file is not fully read");
    fclose(fin);
  }
  else
  {
    fclose(fin);
    return;
  }
  exit(-1);
}
// template<typename T>
// void write(FILE *ofile, const vector<T> &arr) {
//   printf("invoke write(arr)\n");
//   write_arr(ofile, arr);
// }
//
// template<typename T>
// void read(FILE *ifile, vector<T> &arr) {
//   printf("invoke read(arr)\n");
//   read_arr(ifile, arr);
// }

#endif // PRUNEDTD_UTILS_FILE_IO_H_
