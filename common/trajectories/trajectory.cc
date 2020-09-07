#include "drake/common/trajectories/trajectory.h"

#include "fmt/format.h"

#include "drake/common/unused.h"

namespace drake {
namespace trajectories {

template <typename T>
std::string Trajectory<T>::get_coordinate_name(int row, int col) const {
  DRAKE_DEMAND(row >=0 && row < rows());
  DRAKE_DEMAND(col >=0 && col < cols());
  const auto name = coordinate_names_.find({row, col});
  if (name != coordinate_names_.end()) {
    return name->second();
  }
  return fmt::format("{}({},{})", name_, row, col);
}

template <typename T>
std::string Trajectory<T>::get_coordinate_name(int index) const {
  DRAKE_DEMAND(index >=0 && index < rows()*cols());
  int row = index % rows();
  int col = index / cols();
  const auto name = coordinate_names_.find({row, col});
  if (name != coordinate_names_.end()) {
    return name->second();
  }
  return fmt::format("{}({},{})", name_, row, col);
}

template <typename T>
std::string Trajectory<T>::get_name(void) const {
  return name_;
}

template <typename T>
void Trajectory<T>::set_coordinate_name(int row, int col, std::string name) {
  DRAKE_DEMAND(row >=0 && row < rows());
  DRAKE_DEMAND(col >=0 && col < cols());
  coordinate_name_.insert_or_assign({row, col}) = name;
}

template <typename T>
void Trajectory<T>::set_coordinate_name(int index, std::string name) {
  DRAKE_DEMAND(index >=0 && index < rows()*cols());
  int row = index % rows();
  int col = index / cols();
  coordinate_name_.insert_or_assign({row, col}) = name;
}

template <typename T>
void Trajectory<T>::set_name(std::string name) {
  name_ = name;
}

template <typename T>
void Trajectory<T>::set_names_from(const Trajectory<T>& other) {
  name_ = other.name_;
  coordinate_names_ = other.coordinate_names_;
}

template <typename T>
MatrixX<T> Trajectory<T>::vector_values(const std::vector<T>& t) const {
  if (cols() != 1 && rows() != 1) {
    throw std::runtime_error(
        "This method only supports vector-valued trajectories.");
  }
  if (cols() == 1) {
    MatrixX<T> values(rows(), t.size());
    for (int i = 0; i < static_cast<int>(t.size()); i++) {
      values.col(i) = value(t[i]);
    }
    return values;
  }
  MatrixX<T> values(t.size(), cols());
  for (int i = 0; i < static_cast<int>(t.size()); i++) {
    values.row(i) = value(t[i]);
  }
  return values;
}

template <typename T>
bool Trajectory<T>::has_derivative() const {
  return do_has_derivative();
}

template <typename T>
bool Trajectory<T>::do_has_derivative() const {
  return false;
}

template <typename T>
MatrixX<T> Trajectory<T>::EvalDerivative(const T& t,
                                         int derivative_order) const {
  return DoEvalDerivative(t, derivative_order);
}

template <typename T>
MatrixX<T> Trajectory<T>::DoEvalDerivative(const T& t,
                                           int derivative_order) const {
  unused(t);
  unused(derivative_order);
  if (has_derivative()) {
    throw std::logic_error(
        "Trajectory classes that promise derivatives via do_has_derivative() "
        "must implement DoEvalDerivative().");
  } else {
    throw std::logic_error(
        "You asked for derivatives from a class that does not support "
        "derivatives.");
  }
}

template <typename T>
std::unique_ptr<Trajectory<T>> Trajectory<T>::MakeDerivative(
    int derivative_order) const {
  return DoMakeDerivative(derivative_order);
}

template <typename T>
std::unique_ptr<Trajectory<T>> Trajectory<T>::DoMakeDerivative(
    int derivative_order) const {
  unused(derivative_order);
  if (has_derivative()) {
    throw std::logic_error(
        "Trajectory classes that promise derivatives via do_has_derivative() "
        "must implement DoMakeDerivative().");
  } else {
    throw std::logic_error(
        "You asked for derivatives from a class that does not support "
        "derivatives.");
  }
}

}  // namespace trajectories
}  // namespace drake

DRAKE_DEFINE_CLASS_TEMPLATE_INSTANTIATIONS_ON_DEFAULT_SCALARS(
    class drake::trajectories::Trajectory)
