#include "get_integral.h"
#include "utils.h"

bool noel{false}; // no electrons
bool nobeta{false}; // nobeta
bool shift{true}; // shift
bool help{false}; // print stuff
bool zerophoton{false}; // to test - zeroing photonic parts

// polaritonic parameters  
extern std::unique_ptr<polariton> ppol;

// get nuclear energy
double get_integral(std::string filename) {
  gmb::check_file(filename);
  molpro::FCIdump dump(filename);
  int i, j, k, l;
  double integral(0.0);
  molpro::FCIdump::integralType type;
  dump.rewind();
  while ((type = dump.nextIntegral(i, j, k, l, integral)) != molpro::FCIdump::endOfFile) {
    if (type == molpro::FCIdump::I0)
      if (false) std::cout << "found " <<
                "scalar integral " << integral<< "\n";
  }
  return integral;
}

void read_dump(std::string filename, 
               const std::vector<orb_type>& orb_types, 
               const std::vector<spin>& v_spin,
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<size_t>>& v_norb,
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>> &v_sp,
               std::vector<std::vector<bool>>& ssss,
               bool &uhf) {

  gmb::check_file(filename);

  // read parameters from fcidump file
  molpro::FCIdump dump{filename};
  unsigned int nb = dump.parameter("NORB")[0];
  unsigned int nel = dump.parameter("NELEC")[0];
  unsigned int nbeta = nel/2;
  unsigned int nalpha = nel - nbeta;
  unsigned int nphoton{1}; // occupied is always 1 (vacuum orbital)

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
  if (ppol != nullptr) {
    syms_t fermi_ph = {nphoton, 0, 0, 0, 0, 0, 0, 0};
    syms_t full_ph = {nphoton + ppol->nmax, 0, 0, 0, 0, 0, 0, 0};
    no.push_back(nphoton);
    nv.push_back(ppol->nmax);
    occ.push_back({empty, fermi_ph});
    vir.push_back({fermi_ph, full_ph});
    std::pair<syms_t, syms_t> bas_ph = {empty, full_ph}; // idk if needeed, but just in case
  }


  // fill in vectors
  for (size_t iot = 0; iot < orb_types.size(); iot++) {
    for (auto &&ispin : v_spin) {
      switch (orb_types[iot]) {
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
    for (auto &&ispin : v_spin) 
      sp += v_norb[ispin][iot]; //alpha+beta space 
    libtensor::bispace<1> space(sp); 
    if(v_spin.size() > 1 && v_norb[1][iot] != 0) 
      space.split(v_norb[0][iot]); // split space alpha/beta
    if (ppol != nullptr && v_norb[2][iot] != 0) 
      space.split(v_norb[0][iot]+v_norb[1][iot]); // split space electrons/photons
    v_sp.push_back(space);
  }

  for (auto &&ispin : v_spin) {
    for (size_t iot = 0; iot < orb_types.size(); iot++) {
      size_t count{0};
      for (sym_t isym = 0; isym < 8; isym++)
        if ((v_psi[ispin][iot].second[isym] - v_psi[ispin][iot].first[isym]) == 0) 
          ++count;
        if (count == 8) {
          ssss[ispin][iot] = false;
          std::cout << "orb number: " << iot 
                    << " ispin: " << ispin
                    << " does not exist."
                    << "\n";
      }
    }
  }

  for (sym_t isym = 0; isym < 8; isym++) {
    for (size_t ino = 0; ino < orb_types.size(); ino++) 
      for (auto &&ispin : v_spin) {
        if (isym == 0) v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym];
        else v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym] + v_psi[ispin][ino].second[isym-1] + v_shift[ispin][ino][isym-1];
    }
  }
}

// get one-electron integral
container<2,double> get_integral(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2,
                                 bool so_basis) {
                                 
  std::vector<spin> v_spin = {alpha}; // vector containing possible spins
  if (so_basis) 
    v_spin.push_back(beta);
  if (ppol != nullptr) 
    v_spin.push_back(photon);
  std::vector<orb_type> orb_types = {o1,o2}; // vector containing orbital types
  std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (orb_types.size())); // vector containing bra and ket
  std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (orb_types.size())); // vector containing number of orbitals in each bra/ket
  std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (orb_types.size(), std::vector<int> (8,0))); // vector containing symmetry shift 
  std::vector<libtensor::bispace<1>> v_sp; // vector containing 1D spaces for each bra/ket
  std::vector<std::vector<bool>> ssss(v_spin.size(), std::vector<bool> (orb_types.size(), true)); // vector containing if block s or not
  bool uhf{false};

  read_dump(filename, orb_types, v_spin, v_psi, v_norb, v_shift, v_sp, ssss, uhf);

  
  // set integral space
  std::unique_ptr<libtensor::bispace<2>> psp;
  if (o1 == o2) {
    libtensor::bispace<2> sp(v_sp[0]&v_sp[1]);
    psp.reset(new libtensor::bispace<2>(sp));
  } else {
    libtensor::bispace<2> sp(v_sp[0]|v_sp[1]);
    psp.reset(new libtensor::bispace<2>(sp));
  }
  container<2,double> integral(*psp);
  psp.release();
  
  // set integral symmetry
  if (o1 == o2) gmb::set_sym_pp(integral);
  gmb::zero(integral);

  libtensor::block_tensor_wr_ctrl<2, double> ctrl(integral);

  // Loop over blocks using orbit_list
  libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
  for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
         it != ol.end(); it++) {
    // Obtain the index of the current block
    libtensor::index<2> bidx;
    ol.get_index(it, bidx);
#if 1

    std::vector<size_t> bidx_cp(orb_types.size());
    for (size_t i = 0; i < orb_types.size(); i++) {
      bidx_cp[i] = bidx[i];
      if (!ssss[0][i]) ++bidx_cp[i]; // if alpha block doesn't 
    }

    spin spin = alpha;
    auto itype = molpro::FCIdump::I1a;
    bool skip(false);
    if (bidx_cp[0] == 0 && bidx_cp[1] == 0) { // aa
      itype = molpro::FCIdump::I1a;
      spin = alpha; 
    } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1) { // bb
      if (uhf) itype = molpro::FCIdump::I1b;
      spin = beta;
    } else if (bidx_cp[0] == 2 && bidx_cp[1] == 2) { // pp
      skip = true;
    } else { 
        ctrl.req_zero_block(bidx);
        continue;
    }
      
    if (!skip) {
      // Request tensor block from control object
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      // Request data pointer
       const libtensor::dimensions<2> &tdims = blk.get_dims();
      // Request data pointer
      double *ptr = tc.req_dataptr();
      // read integrals from fcidump file
      size_t i, j, k, l;
      unsigned int symi, symj, symk, syml;
      double value;
      molpro::FCIdump::integralType type;
      
      std::string h1file{filename};
      if (ppol != nullptr) {
        h1file.resize(filename.find_last_of("/")+1);
        h1file += "PERTH0MO";
        // h1file += "H0MO";
      }
      gmb::check_file(h1file);
      molpro::FCIdump dump(h1file);
      dump.rewind();
      while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
        if ((((i) >= v_psi[spin][0].first[symi] & (i) < v_psi[spin][0].second[symi]) 
          && ((j) >= v_psi[spin][1].first[symj] & (j) < v_psi[spin][1].second[symj]))
          && type == itype) {
          size_t offset = (j+v_shift[spin][1][symj])+(i+v_shift[spin][0][symi])*(v_norb[spin][1]);
          ptr[offset] = value;
          if (false) {
            std::cout << "i = " << i << " j = " << j << " k = " << k << " l = " << l<< "\n";
            std::cout << "symi = " << symi << " symj = " << symj << " symk = " << symk << " syml = " << syml<< "\n";
            std::cout << "first if\n";
            std::cout << " value = " << value<< "\n";
            std::cout << " offset = " << offset<< "\n";
          }
        }
        if ((((j) >= v_psi[spin][0].first[symj] & (j) < v_psi[spin][0].second[symj]) 
          && ((i) >= v_psi[spin][1].first[symi] & (i) < v_psi[spin][1].second[symi]))
          && type == itype) {
          size_t offset = (i+v_shift[spin][1][symi])+(j+v_shift[spin][0][symj])*(v_norb[spin][1]);
          ptr[offset] = value;
          if (false) {
            std::cout << "second if\n";
            std::cout << " value = " << value<< "\n";
            std::cout << " offset = " << offset<< "\n";
          }
        }
      }
      // Return data pointer
      tc.ret_dataptr(ptr);
      // Return the tensor block (mark as done)
      ctrl.ret_block(bidx);
    } else if (ppol != nullptr) {
      if (o1 == o2) {
        // Request tensor block from control object
        libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
        libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
        // Request data pointer
         const libtensor::dimensions<2> &tdims = blk.get_dims();
        // Request data pointer
        double *ptr = tc.req_dataptr();
        switch (o1) {
        case (o):
          for (size_t i = 0; i < 1; i++) {
            if (shift) 
              ptr[i+i*ppol->nmax] = ppol->omega*(i); // with shift          
            else
              ptr[i+i*ppol->nmax] = ppol->omega*(0.5+i);
          }
          break;
        case (v):
          for (size_t i = 0; i < ppol->nmax; i++) {
            if (shift) 
              ptr[i+i*ppol->nmax] = ppol->omega*(1+i); // with shift        
            else
              ptr[i+i*ppol->nmax] = ppol->omega*(1.5+i);
          }
          break;
        case (b):
          for (size_t i = 0; i < 1+ppol->nmax; i++) {
            if (shift) 
              ptr[i+i*ppol->nmax] = ppol->omega*(i); // with shift        
            else
              ptr[i+i*ppol->nmax] = ppol->omega*(0.5+i);
          }
          break;
        }
        // Return data pointer
        tc.ret_dataptr(ptr);
        // Return the tensor block (mark as done)
        ctrl.ret_block(bidx);
      } else { 
        ctrl.req_zero_block(bidx);
        continue;
    }
    }
#endif
  }  
  if (false) {
    std::cout << "printing integral\n";
    libtensor::bto_print<2, double>(std::cout).perform(integral);
  }

  return integral;
}

// get two-electron integral
container<4,double> get_integral(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2, 
                                 orb_type o3, 
                                 orb_type o4,
                                 bool so_basis) {
                                 
  std::vector<spin> v_spin = {alpha}; // vector containing possible spins
  if (so_basis) 
    v_spin.push_back(beta);
  if (ppol != nullptr)
    v_spin.push_back(photon);
  std::vector<orb_type> orb_types = {o1,o2,o3,o4}; // vector containing orbital types
  std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (orb_types.size())); // vector containing bra and ket
  std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (orb_types.size())); // vector containing number of orbitals in each bra/ket
  std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (orb_types.size(), std::vector<int> (8,0))); // vector containing symmetry shift 
  std::vector<libtensor::bispace<1>> v_sp; // vector containing 1D spaces for each bra/ket
  std::vector<std::vector<bool>> ssss(v_spin.size(), std::vector<bool> (orb_types.size(), true)); // vector containing if block s or not
  bool uhf{false};

  read_dump(filename, orb_types, v_spin, v_psi, v_norb, v_shift, v_sp, ssss, uhf);

  if (help) {
    std::cout << "this is orb_types: " << std::endl;                                    
    for (auto &&i : orb_types)
      std::cout << i << "  ";
    std::cout << "\n";
  }


  // set up integral symmetry
  std::unique_ptr<libtensor::bispace<4>> p_sp4; // pointer to 4D space

    if (orb_types[0] == orb_types[1]) {
     if (orb_types[1] == orb_types[2]) {
       if (orb_types[2] == orb_types[3]) {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
     } else if (orb_types[2] == orb_types[3]) {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
    } else if (orb_types[1] == orb_types[2]) {
       if (orb_types[2] == orb_types[3]) {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]&v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]&v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
     } else if (orb_types[2] == orb_types[3]) {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]|v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        libtensor::bispace<4> sp4(v_sp[0]|v_sp[1]|v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
       
  
  container<4, double> integral(*p_sp4);
  p_sp4.release();

  // set symmetry
  // Request a control object
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(integral);
  libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();
  // permutational symmetry
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
  
  if (false) {
    std::cout << "In get_integral, integral:\n";
    libtensor::bto_print<4, double>(std::cout).perform(integral);
  }

  // Loop over all blocks using orbit_list
  libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());

#if 1
  for (libtensor::orbit_list<4, double>::iterator it = ol.begin();
         it != ol.end(); it++) {

    // Obtain the index of the current block
    libtensor::index<4> bidx;
    
    ol.get_index(it, bidx);
    if (false) {
      std::cout << "bidx[0] = " << bidx[0]<< "\n";
      std::cout << "bidx[1] = " << bidx[1]<< "\n";
      std::cout << "bidx[2] = " << bidx[2]<< "\n";
      std::cout << "bidx[3] = " << bidx[3]<< "\n";
    }

    std::vector<size_t> bidx_cp(orb_types.size());
    for (size_t i = 0; i < orb_types.size(); i++) {
      bidx_cp[i] = bidx[i];
      if (!ssss[0][i]) ++bidx_cp[i]; // if alpha block doesn't 
    }

    if (false) {
      std::cout << "bidx_cp[0] = " << bidx_cp[0]<< "\n";
      std::cout << "bidx_cp[1] = " << bidx_cp[1]<< "\n";
      std::cout << "bidx_cp[2] = " << bidx_cp[2]<< "\n";
      std::cout << "bidx_cp[3] = " << bidx_cp[3]<< "\n";
    }
    
    bool block1{true}, block2{true};
    auto itype = molpro::FCIdump::I2aa;
    spin spin1{alpha}, spin2{alpha}; 
    bool skip{false};
    if (bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 0 && bidx_cp[3] == 0) { // aaaa
      itype = molpro::FCIdump::I2aa;
      spin1 = alpha; spin2 = alpha; 
    } else if ((bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 1 && bidx_cp[3] == 1)) { // aabb
        spin1 = alpha; spin2 = beta;
        if (uhf) { 
          itype = molpro::FCIdump::I2ab;
          block2 = false;
        }
        if (nobeta) skip = true;
    } else if ((bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 0 && bidx_cp[3] == 0)) { // bbaa
        spin1 = beta; spin2 = alpha;
        if (uhf) {
          itype = molpro::FCIdump::I2ab;
          block1 = false;
        }
        if (nobeta) skip = true;
    } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 1 && bidx_cp[3] == 1) { // bbbb
        if (uhf) itype = molpro::FCIdump::I2bb;
        spin1 = beta;
        spin2 = beta;
        if (nobeta) skip = true;
    } else {
        skip = true;
    } 

      if (!skip) {
      // Request tensor block from control object
      libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
      // Obtain dimensions of tensor block
      const libtensor::dimensions<4> &tdims = blk.get_dims();
      // Request data pointer
      double *ptr = tc.req_dataptr();

      molpro::FCIdump dump(filename);
      size_t i, j, k, l;
      unsigned int symi, symj, symk, syml;
      double value;
      molpro::FCIdump::integralType type;
      dump.rewind();
      
      while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
        if (false) {
        std::cout << "type: \n";
            if (type == molpro::FCIdump::I2aa) std::cout << "I2aa\n";
            if (type == molpro::FCIdump::I2ab) std::cout << "I2ab\n";
            if (type == molpro::FCIdump::I2bb) std::cout << "I2bb\n";
            std::cout << "i = " << i << " j = " << j << " k = " << k << " l = " << l<< "\n";
            std::cout << "symi = " << symi << " symj = " << symj << " symk = " << symk << " syml = " << syml<< "\n";
            std::cout << " value = " << value<< "\n";
            std::cout << "v_psi[spin1][0].first[symi]: " << v_psi[spin1][0].first[symi]<< "\n";
            std::cout << "v_psi[spin1][0].second[symi]: " << v_psi[spin1][0].second[symi]<< "\n";
        }
        if (block1) {
        // 1 (ij|kl)
        if ( (((i) >= v_psi[spin1][0].first[symi] && (i) < v_psi[spin1][0].second[symi]) 
          && ((j) >= v_psi[spin1][1].first[symj] && (j)<v_psi[spin1][1].second[symj]))
          && (((k) >= v_psi[spin2][2].first[symk] && (k) < v_psi[spin2][2].second[symk]) 
          && ((l) >= v_psi[spin2][3].first[syml] && (l)<v_psi[spin2][3].second[syml]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(i+v_shift[spin1][0][symi])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(j+v_shift[spin1][1][symj])
                        + (v_norb[spin2][3])*(k+v_shift[spin2][2][symk])
                        + (l+v_shift[spin2][3][syml]);
          ptr[offset] = value;
          if (false) std::cout << "1 offset = " << offset<< "\n";
        }
        // 2 (ji|lk)
        if ( (((j) >= v_psi[spin1][0].first[symj] && (j) < v_psi[spin1][0].second[symj]) && ((i) >= v_psi[spin1][1].first[symi] && (i)<v_psi[spin1][1].second[symi]))
          && (((l) >= v_psi[spin2][2].first[syml] && (l) < v_psi[spin2][2].second[syml]) && ((k) >= v_psi[spin2][3].first[symk] && (k)<v_psi[spin2][3].second[symk]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(j+v_shift[spin1][0][symj])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(i+v_shift[spin1][1][symi])
                        + (v_norb[spin2][3])*(l+v_shift[spin2][2][syml])
                        + (k+v_shift[spin2][3][symk]);
          ptr[offset] = value;
          if (false) std::cout << "3 offset = " << offset<< "\n";
        }
        // 3 (ji|kl)
        if ( (((j) >= v_psi[spin1][0].first[symj] && (j) < v_psi[spin1][0].second[symj]) && ((i) >= v_psi[spin1][1].first[symi] && (i)<v_psi[spin1][1].second[symi]))
          && (((k) >= v_psi[spin2][2].first[symk] && (k) < v_psi[spin2][2].second[symk]) && ((l) >= v_psi[spin2][3].first[syml] && (l)<v_psi[spin2][3].second[syml]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(j+v_shift[spin1][0][symj])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(i+v_shift[spin1][1][symi])
                        + (v_norb[spin2][3])*(k+v_shift[spin2][2][symk])
                        + (l+v_shift[spin2][3][syml]);
          ptr[offset] = value;
          if (false) std::cout << "5 offset = " << offset<< "\n";
        }
        // 4 (ij|lk)
        if ( (((i) >= v_psi[spin1][0].first[symi] && (i) < v_psi[spin1][0].second[symi]) && ((j) >= v_psi[spin1][1].first[symj] && (j)<v_psi[spin1][1].second[symj]))
          && (((l) >= v_psi[spin2][2].first[syml] && (l) < v_psi[spin2][2].second[syml]) && ((k) >= v_psi[spin2][3].first[symk] && (k)<v_psi[spin2][3].second[symk]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(i+v_shift[spin1][0][symi])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(j+v_shift[spin1][1][symj])
                        + (v_norb[spin2][3])*(l+v_shift[spin2][2][syml])
                        + (k+v_shift[spin2][3][symk]);
          ptr[offset] = value;
          if (false) std::cout << "7 offset = " << offset<< "\n";
        }
        }
        if (block2) {
        // 5 (kl|ij)
        if ( (((k) >= v_psi[spin1][0].first[symk] && (k) < v_psi[spin1][0].second[symk]) && ((l) >= v_psi[spin1][1].first[syml] && (l)<v_psi[spin1][1].second[syml]))
          && (((i) >= v_psi[spin2][2].first[symi] && (i) < v_psi[spin2][2].second[symi]) && ((j) >= v_psi[spin2][3].first[symj] && (j)<v_psi[spin2][3].second[symj]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(k+v_shift[spin1][0][symk])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(l+v_shift[spin1][1][syml])
                        + (v_norb[spin2][3])*(i+v_shift[spin2][2][symi])
                        + (j+v_shift[spin2][3][symj]);
          ptr[offset] = value;
          if (false) std::cout << "2 offset = " << offset<< "\n";
        }
        // 6 (lk|ji)
        if ( (((l) >= v_psi[spin1][0].first[syml] && (l) < v_psi[spin1][0].second[syml]) && ((k) >= v_psi[spin1][1].first[symk] && (k)<v_psi[spin1][1].second[symk]))
          && (((j) >= v_psi[spin2][2].first[symj] && (j) < v_psi[spin2][2].second[symj]) && ((i) >= v_psi[spin2][3].first[symi] && (i)<v_psi[spin2][3].second[symi])) 
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(l+v_shift[spin1][0][syml])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(k+v_shift[spin1][1][symk])
                        + (v_norb[spin2][3])*(j+v_shift[spin2][2][symj])
                        + (i+v_shift[spin2][3][symi]);
          ptr[offset] = value;
          if (false) std::cout << "4 offset = " << offset<< "\n";
        }
        // 7 (lk|ij)
        if ( (((l) >= v_psi[spin1][0].first[syml] && (l) < v_psi[spin1][0].second[syml]) && ((k) >= v_psi[spin1][1].first[symk] && (k)<v_psi[spin1][1].second[symk]))
          && (((i) >= v_psi[spin2][2].first[symi] && (i) < v_psi[spin2][2].second[symi]) && ((j) >= v_psi[spin2][3].first[symj] && (j)<v_psi[spin2][3].second[symj]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(l+v_shift[spin1][0][syml])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(k+v_shift[spin1][1][symk])
                        + (v_norb[spin2][3])*(i+v_shift[spin2][2][symi])
                        + (j+v_shift[spin2][3][symj]);
          ptr[offset] = value;
          if (false) std::cout << "6 offset = " << offset<< "\n";
        }
        // 8 (kl|ji)
        if ( (((k) >= v_psi[spin1][0].first[symk] && (k) < v_psi[spin1][0].second[symk]) && ((l) >= v_psi[spin1][1].first[syml] && (l)<v_psi[spin1][1].second[syml]))
          && (((j) >= v_psi[spin2][2].first[symj] && (j) < v_psi[spin2][2].second[symj]) && ((i) >= v_psi[spin2][3].first[symi] && (i)<v_psi[spin2][3].second[symi]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(k+v_shift[spin1][0][symk])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(l+v_shift[spin1][1][syml])
                        + (v_norb[spin2][3])*(j+v_shift[spin2][2][symj])
                        + (i+v_shift[spin2][3][symi]);
          ptr[offset] = value;
          if (false) std::cout << "8 offset = " << offset<< "\n";
        }  
      }
      }
      // Return data pointer
      tc.ret_dataptr(ptr);
      // Return the tensor block (mark as done)
      ctrl.ret_block(bidx);
    } else {
      block1 = true;
      block2 = true;
    if (ppol != nullptr) {
      if (help) {
        std::cout << "bidx_cp[0] = " << bidx_cp[0]<< "\n";
        std::cout << "bidx_cp[1] = " << bidx_cp[1]<< "\n";
        std::cout << "bidx_cp[2] = " << bidx_cp[2]<< "\n";
        std::cout << "bidx_cp[3] = " << bidx_cp[3]<< "\n";
      }
      if ((bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 2 && bidx_cp[3] == 2)) // aapp
      {
        // std::cout << "alpha" << std::endl;    
        spin1 = alpha;
        spin2 = photon;
        block2 = false;
      } else if (bidx_cp[0] == 2 && bidx_cp[1] == 2  && bidx_cp[2] == 0 && bidx_cp[3] == 0) // ppaa
      {
        spin2 = photon;
        spin1 = alpha;
        block1 = false;
      } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 2 && bidx_cp[3] == 2) // bbpp
     {
        spin1 = beta;
        spin2 = photon;
        block2 = false;
      } else if (bidx_cp[0] == 2 && bidx_cp[1] == 2  && bidx_cp[2] == 1 && bidx_cp[3] == 1) // ppbb
     {
        spin2 = photon;
        spin1 = beta;
        block1 = false;
      } else {
        ctrl.req_zero_block(bidx);
        continue;
    } 
    } else {
        ctrl.req_zero_block(bidx);
        continue;
    } 

      // Request tensor block from control object
      libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
      // Obtain dimensions of tensor block
      const libtensor::dimensions<4> &tdims = blk.get_dims();
      // Request data pointer
      double *ptr = tc.req_dataptr();

      // read dipole integrals
      std::string dipfile{filename};
      dipfile.resize(filename.find_last_of("/")+1);
      dipfile += "DMO";

      gmb::check_file(dipfile);
      molpro::FCIdump dump{dipfile}; 
      size_t p, q, r, s;
      unsigned int symp, symq, symr, syms;
      double value;
      molpro::FCIdump::integralType type;
      dump.rewind();
#if 1
      while ((type = dump.nextIntegral(symp, p, symq, q, symr, r, syms, s, value)) != molpro::FCIdump::endOfFile) {
        for (int r = 0; r < ppol->nmax + 1; r++) {
          s = r+1;
          if (help) {
            std::cout << "p = " << p  << "; q = " << q << " value = " << value << std::endl;
            std::cout << "r = " << r  << "; s = " << s <<  std::endl;
            std::cout << "symp = " << symp << " symq = " << symq << " symr = " << symr << " syms = " << syms<< "\n";
          }
          symr = 0;
          syms = 0;
          // 1
          // (pq|rs)
          if (block1) { // ppee
          if ((((p) >= v_psi[spin1][0].first[symp] && (p) < v_psi[spin1][0].second[symp]) 
            && ((q) >= v_psi[spin1][1].first[symq] && (q) < v_psi[spin1][1].second[symq]))
            && (((r) >= v_psi[spin2][2].first[symr] && (r) < v_psi[spin2][2].second[symr])
            && ((s) >= v_psi[spin2][3].first[syms] && (s) < v_psi[spin2][3].second[syms]))) {
              size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(p+v_shift[spin1][0][symp])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(q+v_shift[spin1][1][symq])
                        + (v_norb[spin2][3])*(r+v_shift[spin2][2][symr])
                        + (s+v_shift[spin2][3][syms]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "1 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          //2
          // (qp|rs)
          if ((((q) >= v_psi[spin1][0].first[symq] && (q) < v_psi[spin1][0].second[symq]) 
            && ((p) >= v_psi[spin1][1].first[symp] && (p) < v_psi[spin1][1].second[symp]))
            && (((r) >= v_psi[spin2][2].first[symr] && (r) < v_psi[spin2][2].second[symr])
            && ((s) >= v_psi[spin2][3].first[syms] && (s) < v_psi[spin2][3].second[syms]))) {
              size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(q+v_shift[spin1][0][symq])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(p+v_shift[spin1][1][symp])
                        + (v_norb[spin2][3])*(r+v_shift[spin2][2][symr])
                        + (s+v_shift[spin2][3][syms]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "2 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          #if 1
          // 3
          // (pq|sr)
          if ((((p) >= v_psi[spin1][0].first[symp] && (p) < v_psi[spin1][0].second[symp]) 
            && ((q) >= v_psi[spin1][1].first[symq] && (q) < v_psi[spin1][1].second[symq]))
            && (((r) >= v_psi[spin2][3].first[symr] && (r) < v_psi[spin2][3].second[symr])
            && ((s) >= v_psi[spin2][2].first[syms] && (s) < v_psi[spin2][2].second[syms]))) {
              size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(p+v_shift[spin1][0][symp])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(q+v_shift[spin1][1][symq])
                        + (v_norb[spin2][3])*(s+v_shift[spin2][2][syms])
                        + (r+v_shift[spin2][3][symr]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "3 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          // 4
          // (qp|sr)  
          if ((((p) >= v_psi[spin1][1].first[symp] && (p) < v_psi[spin1][1].second[symp]) 
            && ((q) >= v_psi[spin1][0].first[symq] && (q) < v_psi[spin1][0].second[symq]))
            && (((r) >= v_psi[spin2][3].first[symr] && (r) < v_psi[spin2][3].second[symr])
            && ((s) >= v_psi[spin2][2].first[syms] && (s) < v_psi[spin2][2].second[syms]))) {
              size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(q+v_shift[spin1][0][symq])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(p+v_shift[spin1][1][symp])
                        + (v_norb[spin2][3])*(s+v_shift[spin2][2][syms])
                        + (r+v_shift[spin2][3][symr]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "4 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
            #endif 
          }
          if (block2) {
          // 5
          // (rs|pq)
          if ((((p) >= v_psi[spin1][2].first[symp] && (p) < v_psi[spin1][2].second[symp]) 
            && ((q) >= v_psi[spin1][3].first[symq] && (q) < v_psi[spin1][3].second[symq]))
            && (((r) >= v_psi[spin2][0].first[symr] && (r) < v_psi[spin2][0].second[symr])
            && ((s) >= v_psi[spin2][1].first[syms] && (s) < v_psi[spin2][1].second[syms]))) {
              size_t offset = (v_norb[spin2][1]*v_norb[spin1][2]*v_norb[spin1][3])*(r+v_shift[spin2][0][symr])
                        + (v_norb[spin1][2]*v_norb[spin1][3])*(s+v_shift[spin2][1][syms])
                        + (v_norb[spin1][3])*(p+v_shift[spin1][2][symp])
                        + (q+v_shift[spin1][3][symq]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "5 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          // 6
          // (sr|pq)
          if ((((p) >= v_psi[spin1][2].first[symp] && (p) < v_psi[spin1][2].second[symp]) 
            && ((q) >= v_psi[spin1][3].first[symq] && (q) < v_psi[spin1][3].second[symq]))
            && (((r) >= v_psi[spin2][1].first[symr] && (r) < v_psi[spin2][1].second[symr])
            && ((s) >= v_psi[spin2][0].first[syms] && (s) < v_psi[spin2][0].second[syms]))) {
              size_t offset = (v_norb[spin2][1]*v_norb[spin1][2]*v_norb[spin1][3])*(s+v_shift[spin2][0][syms])
                        + (v_norb[spin1][2]*v_norb[spin1][3])*(r+v_shift[spin2][1][symr])
                        + (v_norb[spin1][3])*(p+v_shift[spin1][2][symp])
                        + (q+v_shift[spin1][3][symq]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "6 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          // 7
          // (rs|qp)
          if ((((p) >= v_psi[spin1][3].first[symp] && (p) < v_psi[spin1][3].second[symp]) 
            && ((q) >= v_psi[spin1][2].first[symq] && (q) < v_psi[spin1][2].second[symq]))
            && (((r) >= v_psi[spin2][0].first[symr] && (r) < v_psi[spin2][0].second[symr])
            && ((s) >= v_psi[spin2][1].first[syms] && (s) < v_psi[spin2][1].second[syms]))) {
              size_t offset = (v_norb[spin2][1]*v_norb[spin1][2]*v_norb[spin1][3])*(r+v_shift[spin2][0][symr])
                        + (v_norb[spin1][2]*v_norb[spin1][3])*(s+v_shift[spin2][1][syms])
                        + (v_norb[spin1][3])*(q+v_shift[spin1][2][symq])
                        + (p+v_shift[spin1][3][symp]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "7 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
          // 8
          // (sr|qp)
          if ((((p) >= v_psi[spin1][3].first[symp] && (p) < v_psi[spin1][3].second[symp]) 
            && ((q) >= v_psi[spin1][2].first[symq] && (q) < v_psi[spin1][2].second[symq]))
            && (((r) >= v_psi[spin2][1].first[symr] && (r) < v_psi[spin2][1].second[symr])
            && ((s) >= v_psi[spin2][0].first[syms] && (s) < v_psi[spin2][0].second[syms]))) {
              size_t offset = (v_norb[spin2][1]*v_norb[spin1][2]*v_norb[spin1][3])*(s+v_shift[spin2][0][syms])
                        + (v_norb[spin1][2]*v_norb[spin1][3])*(r+v_shift[spin2][1][symr])
                        + (v_norb[spin1][3])*(q+v_shift[spin1][2][symq])
                        + (p+v_shift[spin1][3][symp]);
              ptr[offset] = - ppol->gamma*ppol->omega*sqrt(s)*value;
              if (help) std::cout << "8 off set = " << offset << std::endl;
              if (help) std::cout << "ptr[offset]  = " << ptr[offset]  << std::endl;
            }
        }}
      }
    #endif
    // Return data pointer
    tc.ret_dataptr(ptr);
    // Return the tensor block (mark as done)
    ctrl.ret_block(bidx);
  }
    } 

  #endif
  if (false) {
    std::cout << "printing integral\n";
    libtensor::bto_print<4, double>(std::cout).perform(integral);
  }
  return integral;
}

  // get antisymmetrized two-electron integral <pq||rs> 
  container<4,double> get_i(std::string filename, 
  orb_type o1, orb_type o2, orb_type o3, orb_type o4) {
  
  auto h2_o1o3o2o4 = get_integral(filename, o1, o3, o2, o4); 
  auto h2_o1o4o2o3 = get_integral(filename, o1, o4, o2, o3); 
  auto tmpi = get_integral(filename, o1, o2, o3, o4); 
  container<4,double> i(tmpi.get_space());

  // set symmetry
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(i);
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
  gmb::zero(i);
  {
      libtensor::letter p,q,r,s;
      // <pq||rs> = <pq|rs> - <pq|sr> = [pr|qs] - [ps|qr]
      i(p|q|r|s) = h2_o1o3o2o4(p|r|q|s) - h2_o1o4o2o3(p|s|q|r);
  }
  if (false) {
    std::cout << "printing integral " << o1 << o2 << o3 << o4 <<"\n";
    libtensor::bto_print<4, double>(std::cout).perform(i);
  }  
  return i;
}