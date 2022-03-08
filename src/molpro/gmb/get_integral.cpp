#include "get_integral.h"
#include "utils.h"
#include "expressions/anti.h"
#include "expressions/add_d2.h"

#include <molpro/iostream.h>
#include <numeric>

double get_integral(const std::string &filename) {
  molpro::FCIdump dump(filename);
  int i, j, k, l;
  double integral{0.0};
  molpro::FCIdump::integralType type;
  dump.rewind();
  while ((type = dump.nextIntegral(i, j, k, l, integral)) != molpro::FCIdump::endOfFile) {
    if (type == molpro::FCIdump::I0)
      break;
  }
  return integral;
}

  container<2,double> get_integral(const std::string &fname_integrals, const std::string &fname_header, 
    const std::vector<std::shared_ptr<polariton>> &v_ppol, 
    const std::vector<std::unique_ptr<vibration>> &v_pvib, 
    const orb_type &o1, const orb_type &o2, bool add_ph) 
  {
    std::vector<orb_type> v_orb_type = {o1,o2}; // vector containing orbital types

    std::vector<spin> v_spin = {alpha, beta}; // vector containing possible spins
    for (size_t i = 0; i < v_ppol.size(); i++)
      v_spin.push_back(photon);
    for (size_t i = 0; i < v_pvib.size(); i++)
      v_spin.push_back(vib);
    std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (v_orb_type.size())); // vector containing bra and ket
    std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (v_orb_type.size())); // vector containing number of orbitals in each bra/ket
    std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (v_orb_type.size(), std::vector<int> (nsym,0))); // vector containing symmetry shift 
    std::vector<libtensor::bispace<1>> v_sp; // vector containing 1D spaces for each bra/ket
    std::vector<std::vector<bool>> v_exist(v_spin.size(), std::vector<bool> (v_orb_type.size(), true)); // vector containing if block exists or not
    bool uhf{false};
  
    read_dump(fname_header, v_ppol, v_pvib, v_exist, v_norb, v_orb_type, v_psi,  v_shift, v_sp, v_spin, uhf);
    auto integral = set_space(v_orb_type, v_sp);
    gmb::zero(integral);
    get_one_electron_part(integral, fname_integrals, v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
    
    if (v_ppol.size() > 0 && add_ph) {
      get_one_photon_part(integral, v_ppol, v_exist, v_orb_type);
      for (const auto &i_ppol : v_ppol) {
        double fact = i_ppol->omega*i_ppol->gamma*i_ppol->gamma;
        auto rnuc = get_integral(i_ppol->fname_dm);
        #if 1 // add self-energy
        if (i_ppol->self_energy) {
          container<2> sm(integral.get_space()); // second moment of charges
          gmb::zero(sm);
          get_one_electron_part(sm, i_ppol->fname_sm, v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
          integral.axpy(-fact, sm);
          container<2> dm(integral.get_space()); // dipole moment 
          gmb::zero(dm);
          get_one_electron_part(dm, i_ppol->fname_dm, v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
          integral.axpy(2.0*fact*rnuc, dm);
        }
        #endif
      }
    }
    for (const auto &i_pvib : v_pvib) {
      get_one_vibration_part(integral, v_ppol, v_pvib, v_exist, v_orb_type);
      container<2> PI(integral.get_space()); // PI integral
      gmb::zero(PI);
      get_one_electron_part(PI, i_pvib->integral_files[3], v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
      integral.axpy(0.5, PI);
    }
  
    return integral;
  }


  container<4,double> get_i(const std::string &filename, 
                            const std::vector<std::shared_ptr<polariton>> &v_ppol,
                            const std::vector<std::unique_ptr<vibration>> &v_pvib,
                            const orb_type &o1, const orb_type &o2, const orb_type &o3, const orb_type &o4, const bool &add_ph) {
  
  std::unique_ptr<container<4>> tmp_o1o2o3o4, h2_o1o3o2o4, h2_o1o4o2o3;
  
#if 1
  h2_o1o3o2o4 = std::make_unique<container<4>> (get_integral(filename, v_ppol, v_pvib, o1, o3, o2, o4)); 
  if (o3 == o4) 
    h2_o1o4o2o3 = std::make_unique<container<4>>(*h2_o1o3o2o4);
  else 
    h2_o1o4o2o3 = std::make_unique<container<4>> (get_integral(filename, v_ppol, v_pvib, o1, o4, o2, o3)); 
  
  if (o2 == o4 && o2 == o3) 
    tmp_o1o2o3o4 = std::make_unique<container<4>>(*h2_o1o4o2o3);
  else 
    tmp_o1o2o3o4 = std::make_unique<container<4>> (get_integral(filename, v_ppol, v_pvib, o1, o2, o3, o4)); 

#endif
  container<4,double> h2_o1o2o3o4(tmp_o1o2o3o4->get_space());
  gmb::zero(h2_o1o2o3o4);
#if 1

  // set symmetry
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(h2_o1o2o3o4);
  libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();
  libtensor::scalar_transf<double> tr_sym(1.0);
  libtensor::scalar_transf<double> tr_asym(-1.0);
  if (o1 == o2) {
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::se_perm<4, double> se_01(p01, tr_asym);
    sym.insert(se_01);
  }
  if (o2 == o3) {
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::se_perm<4, double> se_23(p23, tr_asym);
    sym.insert(se_23);
  }
  if (o1 == o3 && o2 == o4) {
    libtensor::permutation<4> p0213; p0213.permute(0, 2).permute(1, 3);
    libtensor::se_perm<4, double> se_0213(p0213, tr_sym);
    sym.insert(se_0213);
  }

  anti(h2_o1o2o3o4, *h2_o1o3o2o4, *h2_o1o4o2o3);

  #endif
  #if 1 // add self-energy if needed
  if (add_ph) {
    for (size_t i = 0; i < v_ppol.size(); i++) {

      if (v_ppol[i]->self_energy) {
        std::unique_ptr<container<2>> pd_o1o3, pd_o2o4, pd_o2o3, pd_o1o4;
        double fact = v_ppol[i]->omega*v_ppol[i]->gamma*v_ppol[i]->gamma;

        // add dipole integrals to two-electron part
        pd_o1o3 = std::make_unique<container<2>> (get_integral(v_ppol[i]->fname_dm, filename, v_ppol, v_pvib, o1, o3, false));
        if (o1 == o2 && o3 == o4)
          pd_o2o4 = std::make_unique<container<2>> (*pd_o1o3);
        else 
          pd_o2o4 = std::make_unique<container<2>> (get_integral(v_ppol[i]->fname_dm, filename, v_ppol, v_pvib, o2, o4, false));

        if (o1 == o2)
          pd_o2o3 = std::make_unique<container<2>> (*pd_o1o3);
        else if(o3 == o4)
          pd_o2o3 = std::make_unique<container<2>> (*pd_o2o4);
        else
          pd_o2o3 = std::make_unique<container<2>> (get_integral(v_ppol[i]->fname_dm, filename, v_ppol, v_pvib, o2, o3, false));

        if (o3 == o4)
          pd_o1o4 = std::make_unique<container<2>> (*pd_o1o3);
        else if (o1 == o2)
          pd_o1o4 = std::make_unique<container<2>> (*pd_o2o4);
        else if (o1 == o2 && o3 == o4)
          pd_o1o4 = std::make_unique<container<2>> (*pd_o2o3);
        else
          pd_o1o4 = std::make_unique<container<2>> (get_integral(v_ppol[i]->fname_dm, filename, v_ppol, v_pvib, o1, o4, false));

        add_d2(fact, *pd_o1o3, *pd_o2o4, *pd_o2o3, *pd_o1o4, h2_o1o2o3o4);
      }
    }
  }
  #endif
  return h2_o1o2o3o4;
}

  void read_dump(const std::string &filename, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::unique_ptr<vibration>> &v_pvib,
               std::vector<std::vector<bool>>& v_exist,
               std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type>& v_orb_type, 
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>> &v_sp,
               const std::vector<spin>& v_spin,
               bool &uhf) {


  // read parameters from fcidump file
  molpro::FCIdump dump{filename};
  unsigned int nb = dump.parameter("NORB")[0];
  std::vector<unsigned int> nphoton(v_ppol.size(), 1); // occupied is always 1 (vacuum orbital)
  std::vector<unsigned int> ivib(v_pvib.size(), 1); // occupied is always 1 (vacuum orbital)


  std::vector<int> orbsym = dump.parameter("ORBSYM");
  // size_t ms2 = dump.parameter("MS2")[0];
  uhf = dump.parameter("IUHF")[0];
  
  syms_t empty = {0, 0, 0, 0, 0, 0, 0, 0};
  syms_t fermi(nsym);
  for (size_t i = 0; i < nsym; i++)
    fermi[i] = static_cast<unsigned int>(dump.parameter("OCC")[i]);
  syms_t closed(nsym);
  for (size_t i = 0; i < nsym; i++)
    closed[i] = static_cast<unsigned int>(dump.parameter("CLOSED")[i]);
  
  syms_t full(nsym);
  for (const auto &os : orbsym) full[os-1] += 1;

  sym_t nalpha = std::accumulate(fermi.cbegin(),fermi.cend(),0);
  sym_t nbeta = std::accumulate(closed.cbegin(),closed.cend(),0);
  std::vector<size_t> no = {nalpha, nbeta}, nv = {nb - nalpha, nb - nbeta};

  std::vector<std::pair<syms_t, syms_t>> 
    occ = {{empty, fermi},{empty, closed}}, 
    vir = {{fermi, full},{closed, full}};

  std::pair<syms_t, syms_t> bas = {empty, full};
  

  // photon space
  if (!v_ppol.empty() ) {
    for (size_t i = 0; i < v_ppol.size(); i++) {
      syms_t fermi_ph = {nphoton[i], 0, 0, 0, 0, 0, 0, 0};
      syms_t full_ph = {nphoton[i] + v_ppol[i]->nmax, 0, 0, 0, 0, 0, 0, 0};
      no.push_back(nphoton[i]);
      nv.push_back(v_ppol[i]->nmax);
      occ.emplace_back(std::pair<syms_t,syms_t>({empty, fermi_ph}));
      vir.emplace_back(std::pair<syms_t,syms_t>({fermi_ph, full_ph}));
    }
  }

  // vib space
  if (!v_pvib.empty() ) {
    for (size_t i = 0; i < v_pvib.size(); i++) {
      syms_t fermi_vib = {ivib[i], 0, 0, 0, 0, 0, 0, 0};
      syms_t full_vib = {ivib[i] + v_pvib[i]->nmax, 0, 0, 0, 0, 0, 0, 0};
      no.push_back(ivib[i]);
      nv.push_back(v_pvib[i]->nmax);
      occ.emplace_back(std::pair<syms_t,syms_t>({empty, fermi_vib}));
      vir.emplace_back(std::pair<syms_t,syms_t>({fermi_vib, full_vib}));
    }
  }

  // fill in vectors
  for (size_t iot = 0; iot < v_orb_type.size(); iot++) {
    for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
      switch (v_orb_type[iot]) {
        case (o): { 
          v_norb[ispin][iot] = no[ispin]; 
          v_psi[ispin][iot] = occ[ispin]; 
        } break;
        case (v): { 
          v_norb[ispin][iot] = nv[ispin]; 
          v_psi[ispin][iot] = vir[ispin]; 
        } break;
        case (b): { 
          v_norb[ispin][iot] = nb; 
          v_psi[ispin][iot] = bas; 
        } break;
      }
    }
    size_t sp{0};
    for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
      sp += v_norb[ispin][iot]; 
    }
    libtensor::bispace<1> space(sp); 
    size_t cumspace{0};
    for (size_t i = 0; i < v_spin.size()-1; i++) {
      if (v_norb[i][iot] != 0) {
        cumspace += v_norb[i][iot];
        space.split(cumspace);
      }
    }
    v_sp.push_back(std::move(space));
  }

  for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
    for (size_t iot = 0; iot < v_orb_type.size(); iot++) {
      size_t count{0};
      for (sym_t isym = 0; isym < nsym; isym++) {
        if ((v_psi[ispin][iot].second[isym] - v_psi[ispin][iot].first[isym]) == 0) 
          ++count;
      }
      if (count == nsym) {
        v_exist[ispin][iot] = false;
        molpro::cout << "orb number: " << iot
                  << " ispin: " << ispin
                  << " does not exist."
                  << "\n";
      }
    }
  }

  for (sym_t isym = 0; isym < nsym; isym++) {
    for (size_t ino = 0; ino < v_orb_type.size(); ino++) 
      for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
        if (isym == 0) v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym];
        else v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym] + v_psi[ispin][ino].second[isym-1] + v_shift[ispin][ino][isym-1];
    }
  }
}

  void get_one_electron_part(container<2,double> &integral, 
               const std::string &filename, 
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift, 
               const bool &uhf) 
  {
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(integral);
    libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      libtensor::index<2> bidx;
      ol.get_index(it, bidx);
      std::vector<size_t> bidx_cp(v_orb_type.size());
      for (size_t i = 0; i < v_orb_type.size(); i++) {
        bidx_cp[i] = bidx[i];
        if (!v_exist[0][i]) ++bidx_cp[i]; // if alpha block doesn't exist
      }
      spin spin{alpha};
      auto itype = molpro::FCIdump::I1a;
      if (bidx_cp[0] !=  bidx_cp[1] ) { 
          ctrl.req_zero_block(bidx);
          continue;
      } else if (bidx_cp[0] == alpha) { 
        itype = molpro::FCIdump::I1a;
        spin = alpha; 
      } else if (bidx_cp[0] == beta) { 
        if (uhf) itype = molpro::FCIdump::I1b;
        spin = beta;
      } else if (bidx_cp[0] > beta) { 
          continue;
      }
      
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      double *ptr = tc.req_dataptr();
      size_t i, j, k, l;
      unsigned int symi, symj, symk, syml;
      double value{0.0};
      molpro::FCIdump::integralType type;
      
      molpro::FCIdump dump(filename);
      dump.rewind();
      while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
        if (type == itype) {
          if ((((i) >= v_psi[spin][0].first[symi] && (i) < v_psi[spin][0].second[symi]) 
            && ((j) >= v_psi[spin][1].first[symj] && (j) < v_psi[spin][1].second[symj]))) {
            ptr[gmb::get_offset(i+v_shift[spin][0][symi], j+v_shift[spin][1][symj], v_norb[spin][1])] 
              = value;
          }
          if ((((i) >= v_psi[spin][1].first[symi] && (i) < v_psi[spin][1].second[symi])
            && ((j) >= v_psi[spin][0].first[symj] && (j) < v_psi[spin][0].second[symj]))) {
            ptr[gmb::get_offset(j+v_shift[spin][0][symj], i+v_shift[spin][1][symi], v_norb[spin][1])] 
              = value;
          }
        }
      }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
  }

  void get_one_photon_part(container<2,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::vector<bool>>& v_exist,
               const std::vector<orb_type>& v_orb_type) 
  {
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(integral);
    libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      
      libtensor::index<2> bidx;
      ol.get_index(it, bidx);

      if (bidx[0] != bidx[1] || bidx[0] < photon || bidx[0] >= photon+v_ppol.size()) 
        continue;

      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      const libtensor::dimensions<2> &tdims = blk.get_dims();
      double *ptr = tc.req_dataptr();

      double fact{v_ppol[bidx[0]-2]->gamma*v_ppol[bidx[0]-2]->omega};
      auto nmax = v_ppol[bidx[0]-2]->nmax;
      auto omega = v_ppol[bidx[0]-2]->omega;
      auto rnuc = get_integral(v_ppol[bidx[0]-2]->fname_dm);

      if (v_orb_type[0] != v_orb_type[1]) { // ov block - only one element
        #if 1 //coupling
        if (v_ppol[bidx[0]-2]->coupling) {
          for (size_t i = 0; i < 1; i++) 
             ptr[0] = - fact*rnuc;
        }
        #endif
      } else {
        switch (v_orb_type[0]) {
          case (o): // oo block - only one diagonal element
            for (size_t i = 0; i < 1; i++) 
              ptr[gmb::get_offset(i,i,nmax)] = i*omega;          
            break;
          case (v): // vv block
            for (int p = 1; p < nmax + 1; p++) {
              auto q = p+1;
              // index in block
              auto ip = p-1; 
              auto iq = q-1;
              // diagonal elements 
              ptr[gmb::get_offset(ip,ip,nmax)] = p*omega;        
              // off-diagonal elements
              #if 1 //coupling
              if (v_ppol[bidx[0]-2]->coupling) {
                if (!(q > nmax)) {
                  ptr[gmb::get_offset(ip,iq,nmax)] = - fact*rnuc*sqrt(q);
                  ptr[gmb::get_offset(iq,ip,nmax)] = - fact*rnuc*sqrt(p+1);
                }
              }
              #endif
            }
            break;
          case (b):
            for (int p = 0; p < nmax + 1; p++) {
              auto q = p+1;
              // diagonal 
              ptr[p+p*nmax] = p*omega;        
              // off-diagonal 
              if (!(q > nmax)) {
                ptr[gmb::get_offset(p,q,nmax)] = - fact*rnuc*sqrt(q);
                ptr[gmb::get_offset(q,p,nmax)] = - fact*rnuc*sqrt(p+1);
              }
            }      
          break;  
        }
      }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
  }

  void get_one_vibration_part(container<2,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::unique_ptr<vibration>> &v_pvib,
               const std::vector<std::vector<bool>>& v_exist,
               const std::vector<orb_type>& v_orb_type) 
  {

    libtensor::block_tensor_wr_ctrl<2, double> ctrl(integral);
    libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      
      libtensor::index<2> bidx;
      ol.get_index(it, bidx);

      if (bidx[0] != bidx[1] || bidx[0] < photon+v_ppol.size()) 
        continue;

      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      const libtensor::dimensions<2> &tdims = blk.get_dims();
      double *ptr = tc.req_dataptr();

      auto nmax = v_pvib[bidx[0]-2-v_ppol.size()]->nmax;
      auto omega = v_pvib[bidx[0]-2-v_ppol.size()]->omega;
      auto c = get_integral(v_pvib[bidx[0]-2-v_ppol.size()]->integral_files[0]);
      
      double K{omega}; // K = mw/hbar
      double fact{c/sqrt(2*K)};

      if (v_orb_type[0] != v_orb_type[1]) { // ov block - only one element
        #if 1 //coupling
        if (v_pvib[bidx[0]-2-v_ppol.size()]->coupling) {
          for (size_t i = 0; i < 1; i++) 
             ptr[0] = fact;
        }
        #endif
      } else {
        switch (v_orb_type[0]) {
          case (o): // oo block - only one diagonal element
            for (size_t i = 0; i < 1; i++) 
              ptr[gmb::get_offset(i,i,nmax)] = i*omega;          
            break;
          case (v): // vv block
            for (int p = 1; p < nmax + 1; p++) {
              auto q = p+1;
              // index in block
              auto ip = p-1; 
              auto iq = q-1;
              // diagonal elements 
              ptr[gmb::get_offset(ip,ip,nmax)] = p*omega;        
              // off-diagonal elements
              #if 1 //coupling
              if (v_pvib[bidx[0]-2-v_ppol.size()]->coupling) {
                if (!(q > nmax)) {
                  ptr[gmb::get_offset(ip,iq,nmax)] = fact*sqrt(q);
                  ptr[gmb::get_offset(iq,ip,nmax)] = fact*sqrt(p+1);
                }
              }
              #endif
            }
            break;
          case (b):
            for (int p = 0; p < nmax + 1; p++) {
              auto q = p+1;
              // diagonal 
              ptr[p+p*nmax] = p*omega;        
              // off-diagonal 
              if (!(q > nmax)) {
                ptr[gmb::get_offset(p,q,nmax)] = fact*sqrt(q);
                ptr[gmb::get_offset(q,p,nmax)] = fact*sqrt(p+1);
              }
            }      
          break;  
        }
      }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
  }


 container<2,double> set_space(const std::vector<orb_type> &v_orb_type, const std::vector<libtensor::bispace<1>> &v_sp) 
  {
    // set integral space
    std::unique_ptr<libtensor::bispace<2>> pspace;
    if (v_orb_type[0] == v_orb_type[1]) {
      libtensor::bispace<2> space(v_sp[0]&v_sp[1]);
      pspace = std::make_unique<libtensor::bispace<2>>(space);
    } else {
      libtensor::bispace<2> space(v_sp[0]|v_sp[1]);
      pspace = std::make_unique<libtensor::bispace<2>>(space);
    }
    container<2,double> integral(*pspace);
    
    // set integral permutational symmetry
    if (v_orb_type[0] == v_orb_type[1]) 
      gmb::set_sym_pp(integral);

    gmb::zero(integral);

    return integral;
  
  }


  void get_two_electron_part(container<4,double> &integral, 
               const std::string &filename, 
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift, 
               const bool &uhf) 
  {
 
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(integral);
    libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<4, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      libtensor::index<4> bidx;
      ol.get_index(it, bidx);
      std::vector<size_t> bidx_cp(v_orb_type.size());
      for (size_t i = 0; i < v_orb_type.size(); i++) {
        bidx_cp[i] = bidx[i];
        if (!v_exist[0][i]) ++bidx_cp[i]; // if alpha block doesn't 
      }
      bool block1{true}, block2{true};
      auto itype = molpro::FCIdump::I2aa;
      spin spin1{alpha}, spin2{alpha}; 
      if (bidx_cp[0] == alpha && bidx_cp[1] == alpha  && bidx_cp[2] == alpha && bidx_cp[3] == alpha) { 
        itype = molpro::FCIdump::I2aa;
        spin1 = alpha; spin2 = alpha; 
      } else if ((bidx_cp[0] == alpha && bidx_cp[1] == alpha  && bidx_cp[2] == beta && bidx_cp[3] == beta)) { 
          spin1 = alpha; spin2 = beta;
          if (uhf) { 
            itype = molpro::FCIdump::I2ab;
            block2 = false;
          }
      } else if ((bidx_cp[0] == beta && bidx_cp[1] == beta  && bidx_cp[2] == alpha && bidx_cp[3] == alpha)) { 
          spin1 = beta; spin2 = alpha;
          if (uhf) {
            itype = molpro::FCIdump::I2ab;
            block1 = false;
          }
      } else if (bidx_cp[0] == beta && bidx_cp[1] == beta  && bidx_cp[2] == beta && bidx_cp[3] == beta) { 
          if (uhf) itype = molpro::FCIdump::I2bb;
          spin1 = beta;
          spin2 = beta;   
      } else if ((bidx_cp[0] <= beta && bidx_cp[1] <= beta  && bidx_cp[2] > beta && bidx_cp[3] > beta)
        || ((bidx_cp[0] > beta && bidx_cp[1] > beta  && bidx_cp[2] <= beta && bidx_cp[3] <= beta))){
          continue;
      } else {
          ctrl.req_zero_block(bidx);
          continue;
      } 
  
    libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
    libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
    const libtensor::dimensions<4> &tdims = blk.get_dims();
    double *ptr = tc.req_dataptr();
    molpro::FCIdump dump(filename);
    size_t i, j, k, l;
    unsigned int symi, symj, symk, syml;
    double value{0.0};
    molpro::FCIdump::integralType type;
    dump.rewind();     
    while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
      if (type == itype) {
        if (block1) {
          // 1 (ij|kl)
          if ( ((v_psi[spin1][0].first[symi] <= i && i < v_psi[spin1][0].second[symi]) && (v_psi[spin1][1].first[symj] <= j && j < v_psi[spin1][1].second[symj]))
            && ((v_psi[spin2][2].first[symk] <= k && k < v_psi[spin2][2].second[symk]) && (v_psi[spin2][3].first[syml] <= l && l < v_psi[spin2][3].second[syml]))) {
            ptr[gmb::get_offset(i+v_shift[spin1][0][symi], j+v_shift[spin1][1][symj], k+v_shift[spin2][2][symk], l+v_shift[spin2][3][syml],
                                v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 2 (ji|lk)
          if ( ((v_psi[spin1][0].first[symj] <= j && j < v_psi[spin1][0].second[symj]) && (v_psi[spin1][1].first[symi] <= i && i < v_psi[spin1][1].second[symi]))
            && ((v_psi[spin2][2].first[syml] <= l && l < v_psi[spin2][2].second[syml]) && (v_psi[spin2][3].first[symk] <= k && k < v_psi[spin2][3].second[symk]))) {
            ptr[gmb::get_offset(j+v_shift[spin1][0][symj], i+v_shift[spin1][1][symi], l+v_shift[spin2][2][syml], k+v_shift[spin2][3][symk],
                                v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 3 (ji|kl)
          if ( ((v_psi[spin1][0].first[symj] <= j && j < v_psi[spin1][0].second[symj]) && (v_psi[spin1][1].first[symi] <= i && i < v_psi[spin1][1].second[symi]))
            && ((v_psi[spin2][2].first[symk] <= k && k < v_psi[spin2][2].second[symk]) && (v_psi[spin2][3].first[syml] <= l && l < v_psi[spin2][3].second[syml]))) {
            ptr[gmb::get_offset(j+v_shift[spin1][0][symj], i+v_shift[spin1][1][symi], k+v_shift[spin2][2][symk], l+v_shift[spin2][3][syml],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 4 (ij|lk)
          if ( ((v_psi[spin1][0].first[symi] <= i && i < v_psi[spin1][0].second[symi]) && (v_psi[spin1][1].first[symj] <= j && j < v_psi[spin1][1].second[symj]))
            && ((v_psi[spin2][2].first[syml] <= l && l < v_psi[spin2][2].second[syml]) && (v_psi[spin2][3].first[symk] <= k && k < v_psi[spin2][3].second[symk]))) {
            ptr[gmb::get_offset(i+v_shift[spin1][0][symi], j+v_shift[spin1][1][symj], l+v_shift[spin2][2][syml], k+v_shift[spin2][3][symk],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
        }
        if (block2) {
          // 5 (kl|ij)
          if ( ((v_psi[spin1][0].first[symk] <= k && k < v_psi[spin1][0].second[symk]) && (v_psi[spin1][1].first[syml] <= l && l < v_psi[spin1][1].second[syml]))
            && ((v_psi[spin2][2].first[symi] <= i && i < v_psi[spin2][2].second[symi]) && (v_psi[spin2][3].first[symj] <= j && j < v_psi[spin2][3].second[symj]))) {
            ptr[gmb::get_offset(k+v_shift[spin1][0][symk], l+v_shift[spin1][1][syml], i+v_shift[spin2][2][symi], j+v_shift[spin2][3][symj],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 6 (lk|ji)
          if ( ((v_psi[spin1][0].first[syml] <= l && l < v_psi[spin1][0].second[syml]) && (v_psi[spin1][1].first[symk] <= k && k < v_psi[spin1][1].second[symk]))
            && ((v_psi[spin2][2].first[symj] <= j && j < v_psi[spin2][2].second[symj]) && (v_psi[spin2][3].first[symi] <= i && i < v_psi[spin2][3].second[symi])) ) {
            ptr[gmb::get_offset(l+v_shift[spin1][0][syml], k+v_shift[spin1][1][symk], j+v_shift[spin2][2][symj], i+v_shift[spin2][3][symi],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 7 (lk|ij)
          if ( ((v_psi[spin1][0].first[syml] <= l && l < v_psi[spin1][0].second[syml]) && (v_psi[spin1][1].first[symk] <= k && k < v_psi[spin1][1].second[symk]))
            && ((v_psi[spin2][2].first[symi] <= i && i < v_psi[spin2][2].second[symi]) && (v_psi[spin2][3].first[symj] <= j && j < v_psi[spin2][3].second[symj]))) {
            ptr[gmb::get_offset(l+v_shift[spin1][0][syml], k+v_shift[spin1][1][symk], i+v_shift[spin2][2][symi], j+v_shift[spin2][3][symj],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }
          // 8 (kl|ji)
          if ( ((v_psi[spin1][0].first[symk] <= k && k < v_psi[spin1][0].second[symk]) && (v_psi[spin1][1].first[syml] <= l && l < v_psi[spin1][1].second[syml]))
            && ((v_psi[spin2][2].first[symj] <= j && j < v_psi[spin2][2].second[symj]) && (v_psi[spin2][3].first[symi] <= i && i < v_psi[spin2][3].second[symi]))) {
            ptr[gmb::get_offset(k+v_shift[spin1][0][symk], l+v_shift[spin1][1][syml], j+v_shift[spin2][2][symj], i+v_shift[spin2][3][symi],
                                  v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
              = value;
          }  
        }
      }
    }
    tc.ret_dataptr(ptr);
    ctrl.ret_block(bidx);
    }
  }

  void get_electron_photon_part(container<4,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift)
  {
     libtensor::block_tensor_wr_ctrl<4, double> ctrl(integral);
    libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<4, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      libtensor::index<4> bidx;
      ol.get_index(it, bidx);
      std::vector<size_t> bidx_cp(v_orb_type.size());
      for (size_t i = 0; i < v_orb_type.size(); i++) {
        bidx_cp[i] = bidx[i];
        if (!v_exist[0][i]) ++bidx_cp[i]; // if alpha block doesn't exist
        if (!v_exist[1][i]) ++bidx_cp[i]; // if beta block doesn't exist
      }
      bool block1{true}, block2{true};
      size_t spin1{alpha}, spin2{alpha}; 
      if ((bidx_cp[0] < 2 && bidx_cp[1] < 2 && bidx_cp[2] < 2 && bidx_cp[3] < 2))
        continue; 
      else if (!((bidx_cp[0] == bidx_cp[1]) && (bidx_cp[2] == bidx_cp[3]))) {
          ctrl.req_zero_block(bidx);
          continue;
      } else if ((bidx_cp[0] == alpha) && (bidx_cp[2] >= photon) && (bidx_cp[2] < photon+v_ppol.size())) { //aapp
          spin1 = alpha;
          spin2 = bidx_cp[2];
          block2 = false;
        } else if ((bidx_cp[0] >= photon) && (bidx_cp[0] < photon+v_ppol.size()) && (bidx_cp[2] == alpha)) {// ppaa
          spin2 = bidx_cp[0];
          spin1 = alpha;
          block1 = false;
        } else if ((bidx_cp[0] == beta) && (bidx_cp[2] >= photon) && (bidx_cp[2] <photon+v_ppol.size())) {// bbpp
          spin1 = beta;
          spin2 = bidx_cp[2];
          block2 = false;
        } else if ((bidx_cp[0] >= photon) && (bidx_cp[0] < photon+v_ppol.size()) && (bidx_cp[2] == beta)) {// ppbb
          spin2 = bidx_cp[0];
          spin1 = beta;
          block1 = false;
        } else {
          ctrl.req_zero_block(bidx);
          continue;
        }   
      #if 1 // coupling
      if (v_ppol[spin2-2]->coupling) {
        libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
        libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
        const libtensor::dimensions<4> &tdims = blk.get_dims();
        double *ptr = tc.req_dataptr();

        // read dipole integrals
        std::string fname_dm{v_ppol[0]->fname_dm};

        molpro::FCIdump dump{fname_dm}; 
        size_t p, q, r, s;
        unsigned int symp, symq, symr, syms;
        double value{0.0};
        molpro::FCIdump::integralType type;
        dump.rewind();
        double fact{v_ppol[spin2-2]->gamma*v_ppol[spin2-2]->omega};
        while ((type = dump.nextIntegral(symp, p, symq, q, symr, r, syms, s, value)) != molpro::FCIdump::endOfFile) {
          if (type != molpro::FCIdump::I0)
          for (int r = 0; r < v_ppol[spin2-2]->nmax + 1; r++) {
            s = r+1;
            symr = 0;
            syms = 0;
            // 1 (pq|rs)
            if (block1) { // ppee
            if (((v_psi[spin1][0].first[symp] <= p && p < v_psi[spin1][0].second[symp]) && (v_psi[spin1][1].first[symq] <= q && q < v_psi[spin1][1].second[symq]))
              && ((v_psi[spin2][2].first[symr] <= r && r < v_psi[spin2][2].second[symr]) && (v_psi[spin2][3].first[syms] <= s && s < v_psi[spin2][3].second[syms]))) {
                ptr[gmb::get_offset(p+v_shift[spin1][0][symp], q+v_shift[spin1][1][symq], r+v_shift[spin2][2][symr], s+v_shift[spin2][3][syms],
                                    v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                  = - fact*sqrt(s)*value;
              }
            //2 (qp|rs)
            if (((v_psi[spin1][0].first[symq] <= q && q < v_psi[spin1][0].second[symq]) && (v_psi[spin1][1].first[symp] <= p && p < v_psi[spin1][1].second[symp]))
              && ((v_psi[spin2][2].first[symr] <= r && r < v_psi[spin2][2].second[symr]) && (v_psi[spin2][3].first[syms] <= s && s < v_psi[spin2][3].second[syms]))) {
                ptr[gmb::get_offset(q+v_shift[spin1][0][symq], p+v_shift[spin1][1][symp], r+v_shift[spin2][2][symr], s+v_shift[spin2][3][syms],
                                    v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                  = - fact*sqrt(s)*value;
              }
            // 3 (pq|sr)
            if (((v_psi[spin1][0].first[symp] <= p && p < v_psi[spin1][0].second[symp]) && (v_psi[spin1][1].first[symq] <= q && q < v_psi[spin1][1].second[symq]))
              && ((v_psi[spin2][3].first[symr] <= r && r < v_psi[spin2][3].second[symr]) && (v_psi[spin2][2].first[syms] <= s && s < v_psi[spin2][2].second[syms]))) {
                ptr[gmb::get_offset(p+v_shift[spin1][0][symp], q+v_shift[spin1][1][symq], s+v_shift[spin2][2][syms], r+v_shift[spin2][3][symr],
                                    v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                  = - fact*sqrt(s)*value;
              }
            // 4 (qp|sr)  
            if (((v_psi[spin1][1].first[symp] <= p && p < v_psi[spin1][1].second[symp]) && (v_psi[spin1][0].first[symq] <= q && q < v_psi[spin1][0].second[symq]))
              && ((v_psi[spin2][3].first[symr] <= r && r < v_psi[spin2][3].second[symr]) && (v_psi[spin2][2].first[syms] <= s && s < v_psi[spin2][2].second[syms]))) {
                ptr[gmb::get_offset(q+v_shift[spin1][0][symq], p+v_shift[spin1][1][symp], s+v_shift[spin2][2][syms], r+v_shift[spin2][3][symr],
                                    v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                  = - fact*sqrt(s)*value;
              }
            }
            if (block2) {
            // 5 (rs|pq)
            if (((v_psi[spin1][2].first[symp] <= p && p < v_psi[spin1][2].second[symp]) && (v_psi[spin1][3].first[symq] <= q && q < v_psi[spin1][3].second[symq]))
              && ((v_psi[spin2][0].first[symr] <= r && r < v_psi[spin2][0].second[symr]) && (v_psi[spin2][1].first[syms] <= s && s < v_psi[spin2][1].second[syms]))) {
                ptr[gmb::get_offset(r+v_shift[spin2][0][symr], s+v_shift[spin2][1][syms], p+v_shift[spin1][2][symp], q+v_shift[spin1][3][symq],
                                    v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                  = - fact*sqrt(s)*value;
              }
            // 6 (sr|pq)
            if (((v_psi[spin1][2].first[symp] <= p && p < v_psi[spin1][2].second[symp]) && (v_psi[spin1][3].first[symq] <= q && q < v_psi[spin1][3].second[symq]))
              && ((v_psi[spin2][1].first[symr] <= r && r < v_psi[spin2][1].second[symr]) && (v_psi[spin2][0].first[syms] <= s && s < v_psi[spin2][0].second[syms]))) {
                ptr[gmb::get_offset(s+v_shift[spin2][0][syms], r+v_shift[spin2][1][symr], p+v_shift[spin1][2][symp], q+v_shift[spin1][3][symq],
                                    v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                  = - fact*sqrt(s)*value;
              }
            // 7 (rs|qp)
            if (((v_psi[spin1][3].first[symp] <= p && p < v_psi[spin1][3].second[symp]) && (v_psi[spin1][2].first[symq] <= q && q < v_psi[spin1][2].second[symq]))
              && ((v_psi[spin2][0].first[symr] <= r && r < v_psi[spin2][0].second[symr]) && (v_psi[spin2][1].first[syms] <= s && s < v_psi[spin2][1].second[syms]))) {
                ptr[gmb::get_offset(r+v_shift[spin2][0][symr], s+v_shift[spin2][1][syms], q+v_shift[spin1][2][symq], p+v_shift[spin1][3][symp],
                                    v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                  = - fact*sqrt(s)*value;
              }
            // 8 (sr|qp)
            if (((v_psi[spin1][3].first[symp] <= p && p < v_psi[spin1][3].second[symp]) && (v_psi[spin1][2].first[symq] <= q && q < v_psi[spin1][2].second[symq]))
              && ((v_psi[spin2][1].first[symr] <= r && r < v_psi[spin2][1].second[symr]) && (v_psi[spin2][0].first[syms] <= s && s < v_psi[spin2][0].second[syms]))) {
                ptr[gmb::get_offset(s+v_shift[spin2][0][syms], r+v_shift[spin2][1][symr], q+v_shift[spin1][2][symq], p+v_shift[spin1][3][symp],
                                    v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                  = - fact*sqrt(s)*value;
              }
            }
          }
        }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
    #endif
    }
  }

  void get_electron_vibration_part(container<4,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::unique_ptr<vibration>> &v_pvib,
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift,
               const size_t &nfname)
  {
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(integral);
    libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());

    for (libtensor::orbit_list<4, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      libtensor::index<4> bidx;
      ol.get_index(it, bidx);
      std::vector<size_t> bidx_cp(v_orb_type.size());
      for (size_t i = 0; i < v_orb_type.size(); i++) {
        bidx_cp[i] = bidx[i];
        if (!v_exist[0][i]) ++bidx_cp[i]; // if alpha block doesn't exist
        if (!v_exist[1][i]) ++bidx_cp[i]; // if beta block doesn't exist
      }

      
      bool block1{true}, block2{true};
      size_t spin1{alpha}, spin2{alpha}; 
      if ((bidx_cp[0] < photon+v_ppol.size() && bidx_cp[1] < photon+v_ppol.size() && bidx_cp[2] < photon+v_ppol.size() && bidx_cp[3] < photon+v_ppol.size()))
        continue; 
      else if (!((bidx_cp[0] == bidx_cp[1]) && (bidx_cp[2] == bidx_cp[3]))) {
          ctrl.req_zero_block(bidx);
          continue;
      } else if ((bidx_cp[0] == alpha) && (bidx_cp[2] >= photon+v_ppol.size())) { //aavv
          spin1 = alpha;
          spin2 = bidx_cp[2];
          block2 = false;
        } else if ((bidx_cp[0] >= photon+v_ppol.size())  && (bidx_cp[2] == alpha)) {// vvaa
          spin2 = bidx_cp[0];
          spin1 = alpha;
          block1 = false;
        } else if ((bidx_cp[0] == beta) && (bidx_cp[2] >= photon+v_ppol.size())) {// bbvv
          spin1 = beta;
          spin2 = bidx_cp[2];
          block2 = false;
        } else if ((bidx_cp[0] >= photon+v_ppol.size()) && (bidx_cp[2] == beta)) {// vvbb
          spin2 = bidx_cp[0];
          spin1 = beta;
          block1 = false;
        } else {
          ctrl.req_zero_block(bidx);
          continue;
        }   

      auto ivib = spin2-2-v_ppol.size();
      auto omega = v_pvib[ivib]->omega;
      auto damping = v_pvib[ivib]->damping;
      #if 1 // zero out coupling

      if (v_pvib[ivib]->coupling) {
        libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
        libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
        const libtensor::dimensions<4> &tdims = blk.get_dims();
        double *ptr = tc.req_dataptr();

        // read coupling integrals - A & pi matrices 
        std::string fname{v_pvib[ivib]->integral_files[nfname]};

        molpro::FCIdump dump{fname}; 
        size_t p, q, r, s;
        unsigned int symp, symq, symr, syms;
        double value{0.0};
        molpro::FCIdump::integralType type;
        dump.rewind();
        
        
        double K = omega; // K = mw/hbar
        double fact{0.0};

        while ((type = dump.nextIntegral(symp, p, symq, q, symr, r, syms, s, value)) != molpro::FCIdump::endOfFile) {
          if (type != molpro::FCIdump::I0) {
            for (int r = 0; r < v_pvib[ivib]->nmax; r++) {
              for (int s = r+1; s < v_pvib[ivib]->nmax+1; s++) {
                if ( (r+s) % 2 == 0)
                  continue;
                if (nfname == 2) {
                  if (s != r+1)
                    continue;
                  else
                    fact = -1/sqrt(2.0*K)*sqrt(s);
                } else {
                  if (damping > 0 ) {
                    fact = 0;
                    for (size_t i = 0; i <= floor(r/2.0); i++) {
                      for (size_t j = 0; j <= floor(s/2.0); j++) {
                        if (s-2*j == r-2*i+1 ) {
                        fact += 1 / sqrt( pow(2.0,r+s) * gmb::factorial(r) * gmb::factorial(s))
                            * pow( K/(K+damping) , (2.0+r+s)/2.0-i-j )
                            * pow( (K/(K+damping) -1 ), i+j)
                            * gmb::factorial(r) / ( gmb::factorial(2*i)*gmb::factorial(r-2*i) )
                            * gmb::factorial(s) / (gmb::factorial(2*j)*gmb::factorial(s-2*j))
                            * gmb::factorial(2.0*i) / gmb::factorial(i)
                            * gmb::factorial(2.0*j) / gmb::factorial(j)
                            * sqrt( pow(2.0,(r-2*i+s-2*j)) * gmb::factorial(r-2.0*i)*gmb::factorial(s-2*j)) 
                            * sqrt(1/ (2.0*omega))*(sqrt(s-2*j))
                        ;
                        } else  if (s-2*j == r-2*i-1) {
                        fact += 1 / sqrt( pow(2.0,r+s) * gmb::factorial(r) * gmb::factorial(s))
                            * pow( K/(K+damping) , (2.0+r+s)/2.0-i-j )
                            * pow( (K/(K+damping) -1 ), i+j)
                            * gmb::factorial(r) / ( gmb::factorial(2*i)*gmb::factorial(r-2*i) )
                            * gmb::factorial(s) / (gmb::factorial(2*j)*gmb::factorial(s-2*j))
                            * gmb::factorial(2.0*i) / gmb::factorial(i)
                            * gmb::factorial(2.0*j) / gmb::factorial(j)
                            * sqrt( pow(2.0,(r-2*i+s-2*j)) * gmb::factorial(r-2.0*i)*gmb::factorial(s-2*j)) 
                            * sqrt(1/ (2.0*omega))*(sqrt(r-2*i))
                        ;
                        }
                      }
                    }
                } else {
                  fact = 1/sqrt(2.0*K)*sqrt(s);
                }
            }

            symr = 0;
            syms = 0;

            // std::cout << "r = " << r << ", s = " << s << "\n";
            // std::cout << "fact = " << fact << "\n";
            // std::cout << "p = " << p << ", q = " << q << "\n";
            // std::cout << "value = " << value << "\n";
            // std::cout << "fact*value = " << fact*value << "\n";
            // 1 (pq|rs)
            if (block1) { // ppee
            if (((v_psi[spin1][0].first[symp] <= p && p < v_psi[spin1][0].second[symp]) && (v_psi[spin1][1].first[symq] <= q && q < v_psi[spin1][1].second[symq]))
              && ((v_psi[spin2][2].first[symr] <= r && r < v_psi[spin2][2].second[symr]) && (v_psi[spin2][3].first[syms] <= s && s < v_psi[spin2][3].second[syms]))) {
                ptr[gmb::get_offset(p+v_shift[spin1][0][symp], q+v_shift[spin1][1][symq], r+v_shift[spin2][2][symr], s+v_shift[spin2][3][syms],
                                    v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])]
                  += fact*value;
              }
            //2 (qp|rs)
            if (p != q) {
              if (((v_psi[spin1][0].first[symq] <= q && q < v_psi[spin1][0].second[symq]) && (v_psi[spin1][1].first[symp] <= p && p < v_psi[spin1][1].second[symp]))
                && ((v_psi[spin2][2].first[symr] <= r && r < v_psi[spin2][2].second[symr]) && (v_psi[spin2][3].first[syms] <= s && s < v_psi[spin2][3].second[syms]))) {
                  ptr[gmb::get_offset(q+v_shift[spin1][0][symq], p+v_shift[spin1][1][symp], r+v_shift[spin2][2][symr], s+v_shift[spin2][3][syms],
                                      v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                    += fact*value;
                }
            }
            // 3 (pq|sr)
            if (r != s) {
              if (((v_psi[spin1][0].first[symp] <= p && p < v_psi[spin1][0].second[symp]) && (v_psi[spin1][1].first[symq] <= q && q < v_psi[spin1][1].second[symq]))
                && ((v_psi[spin2][3].first[symr] <= r && r < v_psi[spin2][3].second[symr]) && (v_psi[spin2][2].first[syms] <= s && s < v_psi[spin2][2].second[syms]))) {
                  ptr[gmb::get_offset(p+v_shift[spin1][0][symp], q+v_shift[spin1][1][symq], s+v_shift[spin2][2][syms], r+v_shift[spin2][3][symr],
                                      v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                    += fact*value;
                }
            }
            // 4 (qp|sr)  
            if (p!=q && r!=s) {
              if (((v_psi[spin1][1].first[symp] <= p && p < v_psi[spin1][1].second[symp]) && (v_psi[spin1][0].first[symq] <= q && q < v_psi[spin1][0].second[symq]))
                && ((v_psi[spin2][3].first[symr] <= r && r < v_psi[spin2][3].second[symr]) && (v_psi[spin2][2].first[syms] <= s && s < v_psi[spin2][2].second[syms]))) {
                  ptr[gmb::get_offset(q+v_shift[spin1][0][symq], p+v_shift[spin1][1][symp], s+v_shift[spin2][2][syms], r+v_shift[spin2][3][symr],
                                      v_norb[spin1][1], v_norb[spin2][2], v_norb[spin2][3])] 
                    += fact*value;
                }
            }
            }
            if (block2) {
            // 5 (rs|pq)
            if (((v_psi[spin1][2].first[symp] <= p && p < v_psi[spin1][2].second[symp]) && (v_psi[spin1][3].first[symq] <= q && q < v_psi[spin1][3].second[symq]))
              && ((v_psi[spin2][0].first[symr] <= r && r < v_psi[spin2][0].second[symr]) && (v_psi[spin2][1].first[syms] <= s && s < v_psi[spin2][1].second[syms]))) {
                ptr[gmb::get_offset(r+v_shift[spin2][0][symr], s+v_shift[spin2][1][syms], p+v_shift[spin1][2][symp], q+v_shift[spin1][3][symq],
                                    v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                  += fact*value;
              }
            // 6 (sr|pq)
            if (s != r) {
              if (((v_psi[spin1][2].first[symp] <= p && p < v_psi[spin1][2].second[symp]) && (v_psi[spin1][3].first[symq] <= q && q < v_psi[spin1][3].second[symq]))
                && ((v_psi[spin2][1].first[symr] <= r && r < v_psi[spin2][1].second[symr]) && (v_psi[spin2][0].first[syms] <= s && s < v_psi[spin2][0].second[syms]))) {
                  ptr[gmb::get_offset(s+v_shift[spin2][0][syms], r+v_shift[spin2][1][symr], p+v_shift[spin1][2][symp], q+v_shift[spin1][3][symq],
                                      v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                    += fact*value;
                }
            }
            // 7 (rs|qp)
            if (p != q) {
              if (((v_psi[spin1][3].first[symp] <= p && p < v_psi[spin1][3].second[symp]) && (v_psi[spin1][2].first[symq] <= q && q < v_psi[spin1][2].second[symq]))
                && ((v_psi[spin2][0].first[symr] <= r && r < v_psi[spin2][0].second[symr]) && (v_psi[spin2][1].first[syms] <= s && s < v_psi[spin2][1].second[syms]))) {
                  ptr[gmb::get_offset(r+v_shift[spin2][0][symr], s+v_shift[spin2][1][syms], q+v_shift[spin1][2][symq], p+v_shift[spin1][3][symp],
                                      v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                    += fact*value;
                }
            }
            // 8 (sr|qp)
            if (p!=q && r!=s) {
              if (((v_psi[spin1][3].first[symp] <= p && p < v_psi[spin1][3].second[symp]) && (v_psi[spin1][2].first[symq] <= q && q < v_psi[spin1][2].second[symq]))
                && ((v_psi[spin2][1].first[symr] <= r && r < v_psi[spin2][1].second[symr]) && (v_psi[spin2][0].first[syms] <= s && s < v_psi[spin2][0].second[syms]))) {
                  ptr[gmb::get_offset(s+v_shift[spin2][0][syms], r+v_shift[spin2][1][symr], q+v_shift[spin1][2][symq], p+v_shift[spin1][3][symp],
                                      v_norb[spin2][1], v_norb[spin1][2], v_norb[spin1][3])] 
                    += fact*value;
                }
            }
            }
          }
        }
        }
        }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
    #endif
    }
  }


  container<4,double> get_integral(const std::string &filename, 
    const std::vector<std::shared_ptr<polariton>> &v_ppol,
    const std::vector<std::unique_ptr<vibration>> &v_pvib,
    const orb_type &o1, const orb_type &o2, const orb_type &o3, const orb_type &o4) {
                                 
  std::vector<spin> v_spin = {alpha, beta}; // possible spins
  for (size_t i = 0; i < v_ppol.size(); i++)
    v_spin.push_back(photon);
  for (size_t i = 0; i < v_pvib.size(); i++)
    v_spin.push_back(vib);
  std::vector<orb_type> v_orb_type = {o1,o2,o3,o4}; // orbital types
  std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (v_orb_type.size())); // bra and ket
  std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (v_orb_type.size())); // number of orbitals in each bra/ket
  std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (v_orb_type.size(), std::vector<int> (nsym,0))); // symmetry shift 
  std::vector<libtensor::bispace<1>> v_sp; // 1D spaces for each bra/ket
  std::vector<std::vector<bool>> v_exist(v_spin.size(), std::vector<bool> (v_orb_type.size(), true)); // if block exists or not
  bool uhf{false};

  read_dump(filename, v_ppol, v_pvib, v_exist, v_norb, v_orb_type, v_psi,  v_shift, v_sp, v_spin, uhf);

  // set up integral symmetry
  std::unique_ptr<libtensor::bispace<4>> p_sp4; // pointer to 4D space

    if (v_orb_type[0] == v_orb_type[1]) {
     if (v_orb_type[1] == v_orb_type[2]) {
       if (v_orb_type[2] == v_orb_type[3]) {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]&v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       } else {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]|v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       }
     } else if (v_orb_type[2] == v_orb_type[3]) {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]&v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       } else {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]|v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       }
    } else if (v_orb_type[1] == v_orb_type[2]) {
       if (v_orb_type[2] == v_orb_type[3]) {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]&v_sp[2]&v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       } else {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]&v_sp[2]|v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       }
     } else if (v_orb_type[2] == v_orb_type[3]) {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]|v_sp[2]&v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       } else {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]|v_sp[2]|v_sp[3]);
        p_sp4 = std::make_unique<libtensor::bispace<4>>(sp4);
       }
  container<4, double> integral(*p_sp4);

  // set symmetry
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(integral);
  libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();
  libtensor::scalar_transf<double> tr(1.0);
  if (o1 == o2) {
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::se_perm<4, double> se_01(p01, tr);
    sym.insert(se_01);
  }
  if (o2 == o3) {
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::se_perm<4, double> se_23(p23, tr);
    sym.insert(se_23);
  }
  if (o1 == o3 && o2 == o4) {
    libtensor::permutation<4> p0213; p0213.permute(0, 2).permute(1, 3);
    libtensor::se_perm<4, double> se_0213(p0213, tr);
    sym.insert(se_0213);
  }
  gmb::zero(integral);
  get_two_electron_part(integral, filename, v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
  if (!v_ppol.empty() )
    get_electron_photon_part(integral, v_ppol, v_exist, v_norb, v_orb_type, v_psi, v_shift);
  if (!v_pvib.empty() ) {
    get_electron_vibration_part(integral, v_ppol, v_pvib, v_exist, v_norb, v_orb_type, v_psi, v_shift, 1); // A matrix
    get_electron_vibration_part(integral, v_ppol, v_pvib, v_exist, v_norb, v_orb_type, v_psi, v_shift, 2); // pi matrix
  }
  return integral;
}



