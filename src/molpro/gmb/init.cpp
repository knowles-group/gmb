#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"
#include "expressions/anti.h"
#include "expressions/add_d2.h"

// extern std::unique_ptr<polariton> ppol;

namespace gmb {

  void init(const std::string &filename, const std::string &method, hamiltonian<> &ham, const std::vector<std::unique_ptr<polariton>> &v_ppol) {

    // One-electron integrals
    auto h1_oo = get_hamiltonian(filename,v_ppol,o,o);
    h1_oo.print();
    auto h1_vv = get_hamiltonian(filename,v_ppol,v,v);
    h1_vv.print();

    #if 0
    // Two-electron integrals <pq||rs> 
    auto int_oooo = get_i(filename, o, o, o, o);
    auto int_oovv = get_i(filename, o, o, v, v);
    auto int_ovov = get_i(filename, o, v, o, v);
    
    // add self-energy if needed
    double fact{0};
    std::unique_ptr<container<2>> pd_oo, pd_ov, pd_vo, pd_vv;
    if (ppol != nullptr) { 

      fact = ppol->omega*ppol->gamma*ppol->gamma;

      // add second moment integrals to one-electron part
      auto sm_oo = get_integral(ppol->fname_sm,o,o,false);
      auto sm_vv = get_integral(ppol->fname_sm,v,v,false);
      h1_oo.axpy(fact, sm_oo);
      h1_vv.axpy(fact, sm_vv);

      // add dipole integrals to two-electron part
      pd_oo = std::make_unique<container<2>> (get_integral(ppol->fname_dip,o,o,false));
      pd_ov = std::make_unique<container<2>> (get_integral(ppol->fname_dip,o,v,false));
      pd_vo = std::make_unique<container<2>> (get_integral(ppol->fname_dip,v,o,false));
      pd_vv = std::make_unique<container<2>> (get_integral(ppol->fname_dip,v,v,false));
      add_d2(fact, *pd_oo, *pd_oo, *pd_oo, *pd_oo, int_oooo);
      add_d2(fact, *pd_ov, *pd_ov, *pd_ov, *pd_ov, int_oovv);
      add_d2(fact, *pd_oo, *pd_vv, *pd_vo, *pd_ov, int_ovov);
    }

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    // get fock matrix
    auto d_oo = diag_xx(h1_oo);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    if (method.find("cc") != std::string::npos) {

      auto int_vvvv = get_i(filename, v, v, v, v);
      if (ppol != nullptr) 
        add_d2(fact, *pd_vv, *pd_vv, *pd_vv, *pd_vv, int_vvvv);
      ham.set(i_vvvv, int_vvvv);

      if (method.find("ccsd") != std::string::npos) {
        auto h1_ov = get_integral(filename,o,v);

        auto int_ooov = get_i(filename, o, o, o, v);
        auto int_ovvv = get_i(filename, o, v, v, v);

        if (ppol != nullptr) {
          auto sm_ov = get_integral(ppol->fname_sm,o,v,false);
          h1_ov.axpy(fact, sm_ov);
          add_d2(fact, *pd_oo, *pd_ov, *pd_oo, *pd_ov, int_ooov);
          add_d2(fact, *pd_ov, *pd_vv, *pd_vv, *pd_ov, int_ovvv);
        }
        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        // getting fock matrix - ov block
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));

      }
    }
    #endif
  }

  container<4,double> get_i(std::string filename, 
                            orb_type o1, 
                            orb_type o2, 
                            orb_type o3, 
                            orb_type o4) {
  
  std::shared_ptr<container<4>> tmp_o1o2o3o4, h2_o1o3o2o4, h2_o1o4o2o3;
  
  h2_o1o3o2o4 = std::make_shared<container<4>> (get_integral(filename, o1, o3, o2, o4)); 

  if (o3 == o4) 
    h2_o1o4o2o3.reset(new container<4>(*h2_o1o3o2o4));
  else 
    h2_o1o4o2o3 = std::make_shared<container<4>> (get_integral(filename, o1, o4, o2, o3)); 
  
  if (o2 == o4 && o2 == o3) 
    tmp_o1o2o3o4.reset(new container<4>(*h2_o1o4o2o3));
  else 
    tmp_o1o2o3o4 = std::make_shared<container<4>> (get_integral(filename, o1, o2, o3, o4)); 

  container<4,double> h2_o1o2o3o4(tmp_o1o2o3o4->get_space());

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

  if (false) {
    std::cout << "printing integral " << o1 << o2 << o3 << o4 <<"\n";
    libtensor::bto_print<4, double>(std::cout).perform(h2_o1o2o3o4);
  }  
  
  return h2_o1o2o3o4;
}

  void read_dump(const std::string &filename, 
               const std::vector<std::unique_ptr<polariton>> &v_ppol,
               std::vector<std::vector<bool>>& v_exist,
               std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type>& v_orb_type, 
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>> &v_space,
               const std::vector<spin>& v_spin,
               bool &uhf) {

  gmb::check_file(filename);

  // read parameters from fcidump file
  molpro::FCIdump dump{filename};
  unsigned int nb = dump.parameter("NORB")[0];
  unsigned int nel = dump.parameter("NELEC")[0];
  unsigned int nbeta = nel/2;
  unsigned int nalpha = nel - nbeta;
  std::vector<unsigned int> nphoton(v_ppol.size(), 1); // occupied is always 1 (vacuum orbital)

  std::vector<size_t> no = {nalpha, nbeta}, nv = {nb - nalpha, nb - nbeta};

  std::vector<int> orbsym = dump.parameter("ORBSYM");
  // size_t ms2 = dump.parameter("MS2")[0];
  uhf = dump.parameter("IUHF")[0];
  
  syms_t empty = {0, 0, 0, 0, 0, 0, 0, 0};
  syms_t fermi(8);
  for (size_t i = 0; i < 8; i++)
    fermi[i] = (unsigned int)dump.parameter("OCC")[i];
  syms_t closed(8);
  for (size_t i = 0; i < 8; i++)
    closed[i] = (unsigned int)dump.parameter("CLOSED")[i];
  
  syms_t full(8);
  for (auto &&os : orbsym) full[os-1] += 1;

  std::vector<std::pair<syms_t, syms_t>> 
    occ = {{empty, fermi},{empty, closed}}, 
    vir = {{fermi, full},{closed, full}};

  std::pair<syms_t, syms_t> bas = {empty, full};
  

  // photon space
  if (v_ppol.size() > 0) {
    for (size_t i = 0; i < v_ppol.size(); i++) {
      syms_t fermi_ph = {nphoton[i], 0, 0, 0, 0, 0, 0, 0};
      syms_t full_ph = {nphoton[i] + v_ppol[i]->nmax, 0, 0, 0, 0, 0, 0, 0};
      no.push_back(nphoton[i]);
      nv.push_back(v_ppol[i]->nmax);
      occ.push_back({empty, fermi_ph});
      vir.push_back({fermi_ph, full_ph});
      std::pair<syms_t, syms_t> bas_ph = {empty, full_ph}; // idk if needeed, but just in case
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
        default: 
          break;
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
    v_space.push_back(std::move(space));
  }

  for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
    for (size_t iot = 0; iot < v_orb_type.size(); iot++) {
      size_t count{0};
      for (sym_t isym = 0; isym < 8; isym++)
        if ((v_psi[ispin][iot].second[isym] - v_psi[ispin][iot].first[isym]) == 0) 
          ++count;
        if (count == 8) {
          v_exist[ispin][iot] = false;
          std::cout << "orb number: " << iot 
                    << " ispin: " << ispin
                    << " does not exist."
                    << "\n";
      }
    }
  }

  for (sym_t isym = 0; isym < 8; isym++) {
    for (size_t ino = 0; ino < v_orb_type.size(); ino++) 
      for (size_t ispin = 0; ispin < v_spin.size(); ispin++) {
        if (isym == 0) v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym];
        else v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym] + v_psi[ispin][ino].second[isym-1] + v_shift[ispin][ino][isym-1];
    }
  }
}

  void get_electron_part(container<2,double> &integral, 
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
        if (!v_exist[0][i]) ++bidx_cp[i]; // if alpha block doesn't 
        if (!v_exist[1][i]) ++bidx_cp[i]; // if beta block doesn't 
      }
      spin spin{alpha};
      auto itype = molpro::FCIdump::I1a;
      bool skip{false};
      if (bidx_cp[0] == 0 && bidx_cp[1] == 0) { // aa
        itype = molpro::FCIdump::I1a;
        spin = alpha; 
      } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1) { // bb
        if (uhf) itype = molpro::FCIdump::I1b;
        spin = beta;
      } else { 
          ctrl.req_zero_block(bidx);
          continue;
      }
      
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      const libtensor::dimensions<2> &tdims = blk.get_dims();
      double *ptr = tc.req_dataptr();
      size_t i, j, k, l;
      unsigned int symi, symj, symk, syml;
      double value;
      molpro::FCIdump::integralType type;
      
      molpro::FCIdump dump(filename);
      dump.rewind();
      while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
        if (type == itype) {
          if ((((i) >= v_psi[spin][0].first[symi] & (i) < v_psi[spin][0].second[symi]) 
            && ((j) >= v_psi[spin][1].first[symj] & (j) < v_psi[spin][1].second[symj]))) {
          auto offset = gmb::get_offset(i+v_shift[spin][0][symi], j+v_shift[spin][1][symj], v_norb[spin][1]);
          ptr[offset] = value;
          }
          if ((((i) >= v_psi[spin][1].first[symi] & (i) < v_psi[spin][1].second[symi])
            && ((j) >= v_psi[spin][0].first[symj] & (j) < v_psi[spin][0].second[symj]))
            && type == itype) {
            auto offset = gmb::get_offset(j+v_shift[spin][0][symj], i+v_shift[spin][1][symi], v_norb[spin][1]);
            ptr[offset] = value;
          }
        }
      }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
  }

  void get_photon_part(container<2,double> &integral, 
               const std::vector<std::unique_ptr<polariton>> &v_ppol,
               const std::vector<std::vector<bool>>& v_exist,
               const std::vector<orb_type>& v_orb_type) 
  {
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(integral);
    libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
      libtensor::index<2> bidx;
      ol.get_index(it, bidx);
      if (bidx[0] != bidx[1] || bidx[0] < 2) 
        continue;
      if (v_orb_type[0] != v_orb_type[1]) {
        ctrl.req_zero_block(bidx);
        continue;
      }
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      const libtensor::dimensions<2> &tdims = blk.get_dims();
      double *ptr = tc.req_dataptr();
      switch (v_orb_type[0]) {
      case (o):
        for (size_t i = 0; i < 1; i++) 
            ptr[i+i*v_ppol[bidx[0]-2]->nmax] = v_ppol[bidx[0]-2]->omega*(i);          
        break;
      case (v):
        for (size_t i = 0; i < v_ppol[bidx[0]-2]->nmax; i++) 
            ptr[i+i*v_ppol[bidx[0]-2]->nmax] = v_ppol[bidx[0]-2]->omega*(1+i);        
        break;
      case (b):
        for (size_t i = 0; i < 1+v_ppol[bidx[0]-2]->nmax; i++) 
            ptr[i+i*v_ppol[bidx[0]-2]->nmax] = v_ppol[bidx[0]-2]->omega*(i);        
        break;
      }
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
    }
  }


  container<2,double> get_hamiltonian(const std::string &filename, 
    const std::vector<std::unique_ptr<polariton>> &v_ppol, const orb_type &o1, const orb_type &o2) 
  {
    std::vector<spin> v_spin = {alpha, beta}; // vector containing possible spins
    if (v_ppol.size() > 0) 
      for (size_t i = 0; i < v_ppol.size(); i++)
        v_spin.push_back(photon);
    std::vector<orb_type> v_orb_type = {o1,o2}; // vector containing orbital types
    std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (v_orb_type.size())); // vector containing bra and ket
    std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (v_orb_type.size())); // vector containing number of orbitals in each bra/ket
    std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (v_orb_type.size(), std::vector<int> (8,0))); // vector containing symmetry shift 
    std::vector<libtensor::bispace<1>> v_space; // vector containing 1D spaces for each bra/ket
    std::vector<std::vector<bool>> v_exist(v_spin.size(), std::vector<bool> (v_orb_type.size(), true)); // vector containing if block exists or not
    bool uhf{false};
  
    read_dump(filename, v_ppol, v_exist, v_norb, v_orb_type, v_psi,  v_shift, v_space, v_spin, uhf);
    
    // set integral space
    std::unique_ptr<libtensor::bispace<2>> pspace;
    if (o1 == o2) {
      libtensor::bispace<2> space(v_space[0]&v_space[1]);
      pspace = std::make_unique<libtensor::bispace<2>>(space);
    } else {
      libtensor::bispace<2> space(v_space[0]|v_space[1]);
      pspace = std::make_unique<libtensor::bispace<2>>(space);
    }
    container<2,double> integral(*pspace);
    pspace.release();
    
    // set integral symmetry
    if (o1 == o2) gmb::set_sym_pp(integral);
    gmb::zero(integral);
  
    get_electron_part(integral, filename, v_exist, v_norb, v_orb_type, v_psi, v_shift, uhf);
    if (v_ppol.size() > 0)
      get_photon_part(integral, v_ppol, v_exist, v_orb_type);
  
    return integral;
  }

} // namespace gmb

