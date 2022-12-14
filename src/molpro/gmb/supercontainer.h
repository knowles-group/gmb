#ifndef GMB_SRC_MOLPRO_GMB_SUPERCONTAINER_H
#define GMB_SRC_MOLPRO_GMB_SUPERCONTAINER_H

#include <molpro/iostream.h>
#include "container.h"
#include "get_integral.h"

template<typename T=double>
class supercontainer {
protected:
  std::map<std::string, std::unique_ptr<container<2,T>>> m_m2;
  std::map<std::string, std::unique_ptr<container<4,T>>> m_m4;
public:
  using value_type = T;

  supercontainer() = default;
  virtual ~supercontainer() = default;

  supercontainer(const std::map<size_t,value_type>& source) {
    throw std::logic_error("supercontainer(const std::map<size_t,value_type>& source) unimplemented");
  }

  supercontainer(const supercontainer &sc) {
    for (const auto &im2 : sc.m_m2) {
      m_m2.insert(std::make_pair(im2.first, std::make_unique<container<2,T>> (*im2.second)));
    }
    for (const auto &im4 : sc.m_m4)
      m_m4.insert(std::make_pair(im4.first, std::make_unique<container<4,T>> (*im4.second)));
  }

  supercontainer& operator=(const supercontainer &sc) {
    for (const auto &im2 : sc.m_m2) {
      set(im2.first, *im2.second);
    }
    for (const auto &im4 : sc.m_m4)
      set(im4.first, *im4.second);
    return *this;
  }

  void set(std::string key, const container<2,T> &c2) {
    if (m_m2.find(key) == m_m2.end()) {
      m_m2.insert(std::make_pair(key, std::make_unique<container<2,T>> (c2)));
    }
    else {
      m_m2[key] = std::make_unique<container<2,T>> (c2);
    }
  }

  void set(std::string key, const container<4,T> &c4) {

     if (m_m4.find(key) == m_m4.end())
      m_m4.insert(std::make_pair(key, std::make_unique<container<4,T>> (c4)));
    else
      m_m4[key] = std::make_unique<container<4,T>> (c4);
   }

  const std::map<std::string, std::unique_ptr<container<2,T>>>& get_m2() const { return m_m2; }
  const std::map<std::string, std::unique_ptr<container<4,T>>>& get_m4() const { return m_m4; }

  container<2,T>& m2get(std::string key) {
    if (m_m2.find(key) == m_m2.end()) molpro::cout << key <<" not found!\n";
    return *m_m2[key];
  }

  container<4,T>& m4get(std::string key) {
    if (m_m4.find(key) == m_m4.end()) molpro::cout << key <<" not found!\n";
    return *m_m4[key]; 
  }


  void fill(T a) {
    for (auto &im2 : m_m2)
      im2.second->fill(a);
    for (auto &im4 : m_m4)
      im4.second->fill(a);
  }

  void scal(T a) {
    for (auto &im2 : m_m2)
      im2.second->scale(a);
    for (auto &im4 : m_m4)
      im4.second->scale(a);
  }

  T dot(const supercontainer &sc) const {
    double product(0.0);
    for (auto &im2 : m_m2)
      for (const auto &isc : sc.get_m2())
        product += im2.second->dot(*isc.second);
    for (auto &im4 : m_m4)
      for (const auto &isc : sc.get_m4())
        product += im4.second->dot(*isc.second);
    return product;
  }

  void axpy(T a, const supercontainer &sc) {
    for (auto &im2 : m_m2)
      for (const auto &isc : sc.get_m2())
        im2.second->axpy(a, *isc.second);
    for (auto &im4 : m_m4)
      for (const auto &isc : sc.get_m4())
        im4.second->axpy(a, *isc.second);
    }

  std::map<size_t, T> select_max_dot(size_t n, const supercontainer& sc) const {
    throw std::logic_error("supercontainer::select_max_dot(size_t n, const supercontainer& sc) unimplemented");
  }

  std::map<size_t, T> select(size_t n, bool max = false, bool ignore_sign = false) const {
    throw std::logic_error("supercontainer::select(size_t n, bool max = false, bool ignore_sign = false) unimplemented");
  }

};

#endif // GMB_SRC_MOLPRO_GMB_SUPERCONTAINER_H
