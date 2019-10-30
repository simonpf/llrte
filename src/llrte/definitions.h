#ifndef _LLRTE_DEFINITIONS_H_
#define _LLRTE_DEFINITIONS_H_

#include <memory>
#include <vector>

#include "llrte/data_types.h"

namespace llrte {

enum class Event { step, absorption, scattering, left_domain };

class ClassAttribute {
 public:
  template <typename T>
  ClassAttribute(std::string name, const T& type)
      : name_(name), type_(std::make_shared<T>(type)) {
    // Nothing to do here.
  }

  const DataType& get_type() const { return *type_; }

  std::string get_name() const { return name_; }

 private:
  std::string name_;
  std::shared_ptr<DataType> type_;
};

using TypeVector = std::vector<std::shared_ptr<DataType>>;

class ClassMethod {
 public:
  template <typename T, typename... TT>
  ClassMethod(std::string name, const T& return_type,
              std::tuple<TT...> input_types)
      : name_(name), return_type_(std::make_shared<T>(return_type)) {
    input_types_ = TypeVector{};
    for (size_t i = 0; i < std::tuple_size<std::tuple<TT...>>::value; ++i) {
      // input_types_.push_back(std::make_shared<decltype(std::get<i>(input_types))>(std::get<i>(input_types)));
    }

    // Nothing to do here.
  }

 private:
  std::string name_;
  std::shared_ptr<DataType> return_type_;
  TypeVector input_types_;
};

}  // namespace llrte
#endif
