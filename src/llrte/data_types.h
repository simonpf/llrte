#ifndef _LLRTE_DATA_TYPES_H_
#define _LLRTE_DATA_TYPES_H_

#include <string>

namespace llrte {

class DataType {
 public:
  DataType() = default;
  virtual std::string get_name() const {
    throw std::runtime_error(
        "You are trying to get the name of an unspecified datatype. "
        "Something is likely going wrong.");
  };
};

template <typename Backend>
class Integer : public DataType {
 public:
  Integer(size_t size) : size_(size) {}
  std::string get_name() const { return Backend::get_name(*this); }

 private:
  size_t size_;
};

template <typename Backend>
class String : public DataType {
 public:
  String(size_t size) {}
  std::string get_name() const { return Backend::get_name(*this); }

 private:
};

template <typename Backend>
class Float : public DataType {
 public:
  Float(size_t size) : size_(size) {}
  std::string get_name() const { return Backend::get_name(*this); }
  size_t get_size() const { return size_; }

 private:
  size_t size_;
};

template <typename NestedType, typename Backend>
class Array : public DataType {
 public:
  Array(NestedType nested_type) : nested_type_(nested_type) {}
  std::string get_name() const { return Backend::get_name(*this); }
  NestedType get_nested_type() { return nested_type_; }

 private:
  NestedType nested_type_;
};

}  // namespace llrte

#endif
