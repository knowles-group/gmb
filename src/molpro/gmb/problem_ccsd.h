#ifndef GMB_PROBLEM_CCSD_H_
#define GMB_PROBLEM_CCSD_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "expressions/ccsd/ccsd.h"
#include "expressions/eom-ccsd/precond_oovv.h"
#include "expressions/eom-ccsd/precond_ov.h"
#include "expressions/diag_oovv.h"
#include "expressions/diag_ov.h"
#include "update.h"


class problem_ccsd : public problem_gen {
public:
  problem_ccsd(const hamiltonian<> &ham)
  : problem_gen(ham) {}

  virtual ~problem_ccsd() {}

  void precondition(const VecRef<container_t> &g,
                    const std::vector<value_t> &shift) const override {
    for (int k = 0; k < g.size(); k++) {
      auto &a = g[k].get();     
      // a.m2get(t1).print();
      // auto t1_new = update_t1(m_ham.m2get(f_oo), m_ham.m2get(f_vv), a.m2get(t1));
      // auto t2_new = update_t2(m_ham.m2get(f_oo), m_ham.m2get(f_vv), a.m4get(t2));
      auto d_ov = diag_ov(m_ham.m2get(f_oo), m_ham.m2get(f_ov), m_ham.m2get(f_vv));
      auto d_oovv = diag_oovv(d_ov, m_ham.m4get(i_oovv));
      auto t1_new = precond_ov(a.m2get(t1), d_ov, 0.0);
      auto t2_new = precond_oovv(a.m4get(t2),d_oovv, 0.0);
      a.set(t1, t1_new);
      a.set(t2, t2_new);
      if (false) {
        std::cout << "This is the precondition" << std::endl;
        std::cout << "t1:" << std::endl;
        t1_new.print();
        std::cout << "t2:" << std::endl;
        t2_new.print();
      }
    }
  }


  value_t residual(const container_t &v, container_t &a) const override {
    value_t value = 0;
    auto &ccv = const_cast<container_t&> (v);        
    auto if_oo = ccsd_if_oo(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_oo),
                   m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto if_ov = ccsd_if_ov(ccv.m2get(t1), m_ham.m2get(f_ov), m_ham.m4get(i_oovv)); 
    auto if_vv = ccsd_if_vv(ccv.m2get(t1), ccv.m4get(t2), if_ov, m_ham.m2get(f_vv), 
                   m_ham.m4get(i_oovv), m_ham.m4get(i_ovvv));

    auto if2_oo = ccsd_if2_oo(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_ov),
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto tau = ccsd_tau(ccv.m2get(t1), ccv.m4get(t2));
    auto iw_oooo = ccsd_iw_oooo(ccv.m2get(t1), tau,   
                    m_ham.m4get(i_oooo), m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));
    auto iw_vvvv = ccsd_iw_vvvv(ccv.m2get(t1), 
                    m_ham.m4get(i_ovvv), m_ham.m4get(i_vvvv));
    auto iw_ovov = ccsd_iw_ovov(ccv.m2get(t1), ccv.m4get(t2), 
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
                         
    auto t1_new = ccsd_t1(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_ov), if_oo, if_ov, if_vv, 
                    m_ham.m4get(i_ooov),m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
    auto t2_new = ccsd_t2(ccv.m2get(t1), ccv.m4get(t2), tau, if2_oo, if_vv, 
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv), 
                    iw_oooo, iw_ovov, iw_vvvv);


    a.set(t1, t1_new);
    a.set(t2, t2_new); 
    
    #if 0 // projector
      auto proj1 = t1_new;
     {
        libtensor::block_tensor_wr_ctrl<2, double> ctrl(proj1);
      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
                 it != ol.end(); it++) {
          libtensor::index<2> bidx;
          ol.get_index(it, bidx);
          if (false) {
            std::cout << "bidx[0] = " << bidx[0]<< "\n";
            std::cout << "bidx[1] = " << bidx[1]<< "\n";
          }
          if ((bidx[0] != 2 && bidx[1] == 2)
            || (bidx[0] == 2 && bidx[1] != 2)) {
              if (false) std::cout << "block should be zero" << std::endl;
              ctrl.req_zero_block(bidx);
          } 
          continue;
      }
     }
    a.set(t1, proj1); 


      auto proj2 = t2_new;
      libtensor::block_tensor_wr_ctrl<4, double> ctrl(proj2);
      libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<4, double>::iterator it = ol.begin();
                 it != ol.end(); it++) {
          libtensor::index<4> bidx;
          ol.get_index(it, bidx);
          if (false) {
            std::cout << "bidx[0] = " << bidx[0]<< "\n";
            std::cout << "bidx[1] = " << bidx[1]<< "\n";
            std::cout << "bidx[2] = " << bidx[2]<< "\n";
            std::cout << "bidx[3] = " << bidx[3]<< "\n";
          }
          if ((bidx[0] != 2 && bidx[1] == 2  && bidx[2] == 2 && bidx[3] != 2)
            || (bidx[0] == 2 && bidx[1] != 2  && bidx[2] != 2 && bidx[3] == 2)
            || (bidx[0] == 2 && bidx[1] == 2  && bidx[2] != 2 && bidx[3] != 2)
            || (bidx[0] != 2 && bidx[1] != 2  && bidx[2] == 2 && bidx[3] == 2)) {
              if (false) std::cout << "block should be zero" << std::endl;
              ctrl.req_zero_block(bidx);
          } 
          continue;
      }
    a.set(t2, proj2); 
     #endif // projector 
    return value;
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {}

  void print(std::ostream& s) const {
    s << "CCSD";
  }
};

#endif //GMB_PROBLEM_CCSD_H_
