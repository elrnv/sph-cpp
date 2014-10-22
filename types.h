#ifndef TYPES_H
#define TYPES_H

// This file contains type defines used throughout the program. Most importantly,
// floating point types, integral types and some constants

// One radian measured in degrees (conversion constant)
#define RADIAN 0.017453292519943

typedef double Real;
typedef unsigned int Size;
typedef Size Index;

#if 0
// Type ID factory
// This is used to avoid virtual function calls and dynamic casts throughout the
// code_base by keeping a type-unique id in the base class
class TypeIdFactory
{
public:
  template <class T>
  static size_t GetTypeId(T const&) // argument deduction
  {
    static size_t const Id = GetTypeIdImpl();
    return Id;
  }

private:
  static size_t GetTypeIdImpl()
  {
    static size_t Id = 0;
    return ++Id;
  }
}; // class TypeIdFactory
#endif

#endif // TYPES_H
