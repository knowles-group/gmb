#include "get_integral.h"

// get nuclear energy
double get_integral(std::string filename) {
  std::ifstream infile;
  infile.open(filename);
  //send error if output not found
  if (!infile) {
      std::cout << "Unable to open file: " << filename<< "\n";
      exit(1); // terminate with error
  }
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

// get one-electron integral
container<2,double> get_integral(std::string filename, 
                                    orb_type o1, 
                                    orb_type o2,
                                    bool so_basis) {
                                            
  std::ifstream infile;
  infile.open(filename);
  //send error if output not found
  if (!infile) {
      std::cout << "Unable to open file: " << filename<< "\n";
      exit(1); // terminate with error
  }
  

  // read parameters from fcidump file
  molpro::FCIdump dump(filename);
  size_t nb = (dump.parameter("NORB")[0]);
  size_t nel = (dump.parameter("NELEC")[0]);
  size_t nbeta = nel/2;
  size_t nalpha = nel - nbeta;

  std::vector<spin> v_spin = {alpha};
  if (so_basis) v_spin.push_back(beta);

  std::vector<size_t> no, nv;
  no.push_back(nalpha);
  no.push_back(nbeta);
  for (auto &&ispin : v_spin) nv.push_back(nb - no[ispin]);

  std::vector<int> orbsym = dump.parameter("ORBSYM");
  // size_t ms2 = dump.parameter("MS2")[0];
  bool uhf = dump.parameter("IUHF")[0];
  
  syms_t empty = {0, 0, 0, 0, 0, 0, 0, 0};
  syms_t fermi(8);
  for (size_t i = 0; i < 8; i++)
    fermi[i] = (unsigned int)dump.parameter("OCC")[i];
  syms_t closed(8);
  for (size_t i = 0; i < 8; i++)
    closed[i] = (unsigned int)dump.parameter("CLOSED")[i];
  
  syms_t full(8);
  for (auto &&os : orbsym) full[os-1] += 1;

  std::vector<std::pair<syms_t, syms_t>> occ, vir;
  occ.push_back({empty, fermi}); // alpha
  occ.push_back({empty, closed}); // beta

  vir.push_back({fermi, full}); // alpha
  vir.push_back({closed, full}); // beta
  std::pair<syms_t, syms_t> bas = {empty, full};

  std::vector<orb_type> orb_types = {o1,o2};
  
  // 2D vectors for spin and orb_type
  std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (orb_types.size())); // vector containing bra1, ket1, bra2, ket2
  std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (orb_types.size())); // vector containing number of orbitals in each bra/ket
  std::vector<libtensor::bispace<1>> v_sp; // vector containing 1D spaces for each bra/ket

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
    size_t sp(0);
    for (auto &&ispin : v_spin) 
      sp += v_norb[ispin][iot]; //space contains alpha+beta 
    libtensor::bispace<1> space(sp);
    if (so_basis) space.split(v_norb[0][iot]); // split space alpha/beta
    v_sp.push_back(space);
  }

  // vector containing if block exists or not
  std::vector<std::vector<bool>> ssss(v_spin.size(), std::vector<bool> (orb_types.size(), true));

  for (auto &&ispin : v_spin) {
    for (size_t iot = 0; iot < orb_types.size(); iot++) {
      size_t count(0);
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

  std::vector<std::vector<std::vector<int>>> v_shift(v_spin.size(), std::vector<std::vector<int>> (2, std::vector<int> (8,0))); // vector containing symmetry shift 

  for (sym_t isym = 0; isym < 8; isym++) {
    for (size_t ino = 0; ino < orb_types.size(); ino++) 
      for (auto &&ispin : v_spin) {
        if (isym == 0) v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym];
        else v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym] + v_psi[ispin][ino].second[isym-1] + v_shift[ispin][ino][isym-1];
    }
  }

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

    std::vector<size_t> bidx_cp(orb_types.size());
    for (size_t i = 0; i < orb_types.size(); i++) {
      bidx_cp[i] = bidx[i];
      if (!ssss[0][i]) ++bidx_cp[i]; // if alpha block doesn't exist
    }

    spin spin = alpha;
    auto itype = molpro::FCIdump::I1a;
    if (bidx_cp[0] == 0 && bidx_cp[1] == 0) { // aa
      itype = molpro::FCIdump::I1a;
      spin = alpha; 
    } else if ((bidx_cp[0] == 0 && bidx_cp[1] == 1) // ab 
            || (bidx_cp[0] == 1 && bidx_cp[1] == 0)) { // ba
        ctrl.req_zero_block(bidx);
        continue;
    } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1) { // bbbb
        if (uhf) itype = molpro::FCIdump::I1b;
        spin = beta;
    }
      // Request tensor block from control object
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      // Request data pointer
      double *ptr = tc.req_dataptr();
      
      // read integrals from fcidump file
      molpro::FCIdump dump(filename);
      size_t i, j, k, l;
      unsigned int symi, symj, symk, syml;
      double value;
      molpro::FCIdump::integralType type;
      dump.rewind();
      while ((type = dump.nextIntegral(symi, i, symj, j, symk, k, syml, l, value)) != molpro::FCIdump::endOfFile) {
        if ((((i) >= v_psi[spin][0].first[symi] & (i) < v_psi[spin][0].second[symi]) 
          && ((j) >= v_psi[spin][1].first[symj] & (j) < v_psi[spin][1].second[symj]))
          && type == itype) {
          size_t offset = (j+v_shift[spin][1][symj])+(i+v_shift[spin][0][symi])*(v_norb[spin][1]);
          ptr[offset] = value;
          if (false) {
            std::cout << "first if\n";
            std::cout << "i = " << i << " j = " << j << " k = " << k << " l = " << l<< "\n";
            std::cout << "symi = " << symi << " symj = " << symj << " symk = " << symk << " syml = " << syml<< "\n";
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
            std::cout << "i = " << i << " j = " << j << " k = " << k << " l = " << l<< "\n";
            std::cout << " value = " << value<< "\n";
            std::cout << " offset = " << offset<< "\n";
          }
        }
      }
      // Return data pointer
      tc.ret_dataptr(ptr);
      // Return the tensor block (mark as done)
      ctrl.ret_block(bidx);
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
                                            
  std::ifstream infile;
  infile.open(filename);
  //send error if output not found
  if (!infile) {
      std::cout << "Unable to open file: " << filename<< "\n";
      exit(1); // terminate with error
  } 
  
  // vector containing possible spins
  std::vector<spin> v_spin = {alpha}; // only alpha if mo basis
  if (so_basis) v_spin.push_back(beta); // also beta if so basis

  // read fci dump parameters
  molpro::FCIdump dump(filename);
  size_t nb = dump.parameter("NORB")[0];
  size_t nel = dump.parameter("NELEC")[0];
  size_t nbeta = nel/2;
  size_t nalpha = nel - nbeta;

  std::vector<size_t> no, nv;
  no.push_back(nalpha);
  no.push_back(nbeta);
  for (auto &&ispin : v_spin) nv.push_back(nb - no[ispin]);

  std::vector<int> orbsym = dump.parameter("ORBSYM");
  bool uhf = dump.parameter("IUHF")[0];

  // set up spaces
  syms_t empty = {0, 0, 0, 0, 0, 0, 0, 0};
  syms_t fermi(8);
  for (size_t i = 0; i < 8; i++)
    fermi[i] = (unsigned int)dump.parameter("OCC")[i];
  syms_t closed(8);
  for (size_t i = 0; i < 8; i++)
    closed[i] = (unsigned int)dump.parameter("CLOSED")[i];
  syms_t full(8);
  for (auto &&os : orbsym) full[os-1] += 1;

  std::vector<std::pair<syms_t, syms_t>> occ, vir;
  occ.push_back({empty, fermi}); // alpha
  occ.push_back({empty, closed}); // beta

  vir.push_back({fermi, full}); // alpha
  vir.push_back({closed, full}); // beta
  std::pair<syms_t, syms_t> bas = {empty, full};


  std::vector<orb_type> orb_types = {o1,o2,o3,o4};

  // 2D vectors for spin and orb_type
  std::vector<std::vector<std::pair<syms_t, syms_t>>> v_psi(v_spin.size(), std::vector<std::pair<syms_t, syms_t>> (orb_types.size())); // vector containing bra1, ket1, bra2, ket2
  std::vector<std::vector<size_t>> v_norb(v_spin.size(), std::vector<size_t> (orb_types.size())); // vector containing number of orbitals in each bra/ket
  
  std::vector<libtensor::bispace<1>> v_sp; // vector containing 1D spaces for each bra/ket


  // fill in vectors
  for (size_t iot = 0; iot < orb_types.size(); iot++) {
    for (auto &&ispin : v_spin) {
      switch (orb_types[iot]) {
        case (o): { 
          v_norb[ispin][iot] = no[ispin]; 
          v_psi[ispin][iot] = occ[ispin]; 
        } break;
        case (v): { 
          v_norb[ispin][iot] = nv[ispin] ; 
          v_psi[ispin][iot] = vir[ispin]; } break;
        case (b): { 
          v_norb[ispin][iot] = nb ; 
          v_psi[ispin][iot] = bas; } break;
        default: 
          break;
      }
  }
    // set up integral space
    size_t sp(0);
    for (auto &&ispin : v_spin) 
      sp += v_norb[ispin][iot]; //space contains alpha+beta 
    libtensor::bispace<1> space(sp);
    if (so_basis && (v_norb[0][iot] < sp)) space.split(v_norb[0][iot]); // split space alpha/beta
    v_sp.push_back(space);
  }

  // vector containing if block exists or not
  std::vector<std::vector<bool>> ssss(2, std::vector<bool> (orb_types.size(), true));


  for (auto &&ispin : v_spin) {
    for (size_t iot = 0; iot < orb_types.size(); iot++) {
      size_t count(0);
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



  std::vector<std::vector<std::vector<int>>> v_shift(2, std::vector<std::vector<int>> (4, std::vector<int> (8,0))); // vector containing symmetry shift 
  for (sym_t isym = 0; isym < 8; isym++) {
    for (size_t ino = 0; ino < orb_types.size(); ino++) 
      for (auto &&ispin : v_spin) {
        if (isym == 0) v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym];
        else v_shift[ispin][ino][isym] = - v_psi[ispin][ino].first[isym] + v_psi[ispin][ino].second[isym-1] + v_shift[ispin][ino][isym-1];
    }
  }
  

  // set up integral symmetry
  size_t sym_type(0);
  std::unique_ptr<libtensor::bispace<4>> p_sp4; // pointer to 4D space

    if (orb_types[0] == orb_types[1]) {
     if (orb_types[1] == orb_types[2]) {
       if (orb_types[2] == orb_types[3]) {
        sym_type = 1;
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        sym_type = 2;
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]&v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
     } else if (orb_types[2] == orb_types[3]) {
        sym_type = 3;
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]&v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       } else {
        sym_type = 4;
        libtensor::bispace<4> sp4(v_sp[0]&v_sp[1]|v_sp[2]|v_sp[3]);
        p_sp4.reset(new libtensor::bispace<4>(sp4));
       }
    } else if (orb_types[1] == orb_types[2]) {
       if (orb_types[2] == orb_types[3]) {
        sym_type = 5;
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
      if (!ssss[0][i]) ++bidx_cp[i]; // if alpha block doesn't exist
    }

    if (false) {
      std::cout << "bidx_cp[0] = " << bidx_cp[0]<< "\n";
      std::cout << "bidx_cp[1] = " << bidx_cp[1]<< "\n";
      std::cout << "bidx_cp[2] = " << bidx_cp[2]<< "\n";
      std::cout << "bidx_cp[3] = " << bidx_cp[3]<< "\n";
    }
    
    bool block1(true), block2(true);
    auto itype = molpro::FCIdump::I2aa;
    spin spin1(alpha), spin2(alpha); 
    int exist = 0;
    if (bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 0 && bidx_cp[3] == 0) { // aaaa
      itype = molpro::FCIdump::I2aa;
      spin1 = alpha; spin2 = alpha; 
    } else if ((bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 0 && bidx_cp[3] == 1) // aaab
            || (bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 1 && bidx_cp[3] == 0) // aaba
            || (bidx_cp[0] == 0 && bidx_cp[1] == 1  && bidx_cp[2] == 0 && bidx_cp[3] == 0) // abaa
            || (bidx_cp[0] == 0 && bidx_cp[1] == 1  && bidx_cp[2] == 0 && bidx_cp[3] == 1) // abab 
            || (bidx_cp[0] == 1 && bidx_cp[1] == 0  && bidx_cp[2] == 0 && bidx_cp[3] == 0) // baaa
            || (bidx_cp[0] == 1 && bidx_cp[1] == 0  && bidx_cp[2] == 0 && bidx_cp[3] == 1) // baab
            || (bidx_cp[0] == 0 && bidx_cp[1] == 1  && bidx_cp[2] == 1 && bidx_cp[3] == 0) // abba
            || (bidx_cp[0] == 1 && bidx_cp[1] == 0  && bidx_cp[2] == 1 && bidx_cp[3] == 0) // baba
            || (bidx_cp[0] == 0 && bidx_cp[1] == 1  && bidx_cp[2] == 1 && bidx_cp[3] == 1) // abbb
            || (bidx_cp[0] == 1 && bidx_cp[1] == 0  && bidx_cp[2] == 1 && bidx_cp[3] == 1) // babb
            || (bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 0 && bidx_cp[3] == 1) // bbab
            || (bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 1 && bidx_cp[3] == 0)) { // bbba
        ctrl.req_zero_block(bidx);
        continue;
    } else if ((bidx_cp[0] == 0 && bidx_cp[1] == 0  && bidx_cp[2] == 1 && bidx_cp[3] == 1)) { // aabb
        spin1 = alpha; spin2 = beta;
        if (uhf) { 
          itype = molpro::FCIdump::I2ab;
          block2 = false;
        }
    } else if ((bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 0 && bidx_cp[3] == 0)) { // bbaa
        spin1 = beta; spin2 = alpha;
        if (uhf) {
          itype = molpro::FCIdump::I2ab;
          block1 = false;
        }
    } else if (bidx_cp[0] == 1 && bidx_cp[1] == 1  && bidx_cp[2] == 1 && bidx_cp[3] == 1) { // bbbb
        if (uhf) itype = molpro::FCIdump::I2bb;
        spin1 = beta;
        spin2 = beta;
    }

      // Request tensor block from control object
      libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
      // Obtain dimensions of tensor block
      const libtensor::dimensions<4> &tdims = blk.get_dims();
      // Request data pointer
      double *ptr = tc.req_dataptr();
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
        if ( (((i) >= v_psi[spin1][0].first[symi] && (i) < v_psi[spin1][0].second[symi]) && ((j) >= v_psi[spin1][1].first[symj] && (j)<v_psi[spin1][1].second[symj]))
          && (((k) >= v_psi[spin2][2].first[symk] && (k) < v_psi[spin2][2].second[symk]) && ((l) >= v_psi[spin2][3].first[syml] && (l)<v_psi[spin2][3].second[syml]))
          && type == itype ) {
          size_t offset = (v_norb[spin1][1]*v_norb[spin2][2]*v_norb[spin2][3])*(i+v_shift[spin1][0][symi])
                        + (v_norb[spin2][2]*v_norb[spin2][3])*(j+v_shift[spin1][1][symj])
                        + (v_norb[spin2][3])*(k+v_shift[spin2][2][symk])
                        + (l+v_shift[spin2][3][syml]);
          ptr[offset] = value;
          if (false) std::cout << "1 offset = " << offset<< "\n";
        }
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
  }
  
  if (false) {
    std::cout << "printing integral\n";
    libtensor::bto_print<4, double>(std::cout).perform(integral);
  }
  return integral;
}

// get antisymmetrized two-electron integral <pq||rs> 
container<4,double> get_i(std::string filename, 
  orb_type o1, orb_type o2, orb_type o3, orb_type o4) {
  
  auto h2_o1o4o2o3 = get_integral(filename, o1, o4, o2, o3); 
  // std::cout << "h2_o1o4o2o3\n";
  // h2_o1o4o2o3.print();
  auto h2_o1o3o2o4 = get_integral(filename, o1, o3, o2, o4); 
  // std::cout << "h2_o1o3o2o4\n";
  // h2_o1o3o2o4.print();
  auto tmpi = get_integral(filename, o1, o2, o3, o4); 
  container<4,double> i(tmpi.get_space());
  // set symmetry
  // Request a control object
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(i);
  libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();
  // permutational symmetry
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
  return i;
}