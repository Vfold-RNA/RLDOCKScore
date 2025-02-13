#pragma once

// #include <iostream>
// #include <sstream>
// #include <istream>
// #include <ostream>
// #include <fstream>
// #include <map>
// #include <string>
// #include <iomanip>
// #include <algorithm>
#include <vector>
#include <cassert>
// #include <set>
// #include <thread>
// #include <numeric>
// #include <mutex>
// #include <cmath>
// #include <unordered_map>

#include <exception> // std::bad_alloc

#include "common.h"
#include "vec3d.h"
#include "mol.h"

namespace rlscore {

// enum DISTANCE_TYPE {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_FLEXIBLE};

inline Size_Type checked_multiply(Size_Type i, Size_Type j) {
    #ifdef MODE_DEBUG
    assert(i > 0 && j > 0);
    #endif
	if(i == 0 || j == 0) return 0;
	const Size_Type tmp = i * j;
	if(tmp < i || tmp < j || tmp / i != j)
		throw std::bad_alloc(); // can't alloc if the size makes sz wrap around
	return tmp;
}

inline Size_Type checked_multiply(Size_Type i, Size_Type j, Size_Type k) {
	return checked_multiply(checked_multiply(i, j), k);
}

class Cell {
  public:
  	const Float cell_size;
    const Vec3d coords_min;
    const Vec3d coords_max;
	std::vector<Size_Type> links;

	Cell(const const Float bf, const Float cs) : coords_min(rec.min_xyz-std::abs(buffer)), coords_max(rec.max_xyz+std::abs(buffer)), cell_size(cs) {
		const Vec3d box_size = coords_max - coords_min;
		const Size_Type c_i = static_cast<Size_Type>(box_size[0]/cell_size) + 1;
		const Size_Type c_j = static_cast<Size_Type>(box_size[1]/cell_size) + 1;
		const Size_Type c_k = static_cast<Size_Type>(box_size[2]/cell_size) + 1;
		rec_heads = Heads(c_i, c_j, c_k);
		lig_heads = Heads(c_i, c_j, c_k);
		rec_links.resize(rec.atoms.size(), k_max_size_type);
		lig_links.resize(lig.atoms.size(), k_max_size_type);
    }

	void update_rec_part(const Receptor& rec) {
		rec_heads = Heads(c_i, c_j, c_k);
	}

	void update_lig_part(const Ligand& lig) {
		lig_heads = Heads(c_i, c_j, c_k);
	}


		Size_Type c_i, c_j, c_k;
	std::vector<Size_Type> data;
  public:
	Heads() : c_i(0), c_j(0), c_k(0) {}
	Heads(const Size_Type i, const Size_Type j, const Size_Type k) : c_i(i), c_j(j), c_k(k), data(checked_multiply(i, j, k), k_max_size_type) {}
	Size_Type dim0() const { return c_i; }
	Size_Type dim1() const { return c_j; }
	Size_Type dim2() const { return c_k; }

	Size_Type& operator()(Size_Type i, Size_Type j, Size_Type k) {
        #ifdef MODE_DEBUG
        assert(i >= 0 && j >= 0 && k >= 0);
        #endif
        return data[c_k*(c_j*i + j) + k];
    }
	const Size_Type& operator()(Size_Type i, Size_Type j, Size_Type k) const {
        #ifdef MODE_DEBUG
        assert(i >= 0 && j >= 0 && k >= 0);
        #endif
        return data[c_k*(c_j*i + j) + k];
    }
	void reset() {
		std::fill(this->data.begin(), this->data.end(), k_max_size_type);
	}
	// Size_Type dim(Size_Type i) const {
	// 	switch(i) {
	// 		case 0: return m_i;
	// 		case 1: return m_j;
	// 		case 2: return m_k;
	// 		default: assert(false); return 0; // to get rid of the warning
	// 	}
	// }
	// void resize(Size_Type i, Size_Type j, Size_Type k) { // data is essentially garbled
    //     assert(i >= m_i && j >= m_j && k >= m_k);
    //     assert(i*j*k >= m_i*m_j*m_k);
	// 	m_i = i;
	// 	m_j = j;
	// 	m_k = k;
	// 	m_data.resize(checked_multiply(i, j, k));
	// }
};



// struct Grid {
//     bool flag = false;
//     std::vector<int> rigid_atoms;
//     Grid() {}
//     Grid(const bool f) : flag(f) {}
// };
// const Grid k_false_grid(false);
// struct Grids {
//     Float width;
//     Array3d<Grid> g_data;
//     int g_i, g_j, g_k;
//     int min_i, min_j, min_k;
// public:
//     //column-major
//     Grids() : g_i(0), g_j(0), g_k(0), min_i(0), min_j(0), min_k(0), width(0) {
//         // std::cout << "grids init" << std::endl;
//     }
//     Grids(const int ii, const int jj, const int kk, const int si, const int sj, const int sk, const Float w, const Grid& filler_val) : g_data(ii,jj,kk,filler_val), g_i(ii), g_j(jj), g_k(kk), min_i(si), min_j(sj), min_k(sk), width(w) {
//         assert(ii >= 0 && jj >= 0 && kk >= 0);
//         assert(w >= 0.0);
//     }
//     const Grid& at(int i, int j, int k) const {
//         if(this->out_boundary(i,j,k)) {
//             return k_false_grid;
//         }
//         else {
//             return g_data(i-this->min_i,j-this->min_j,k-this->min_k);
//         }
//     }
//     // const Grid& operator()(int i, int j, int k) const {
//     //     if(this->out_boundary(i,j,k)) {
//     //         return k_false_grid;
//     //     }
//     //     else {
//     //         #ifdef MODE_DEBUG
//     //         assert(i-this->min_i >= 0 && j-this->min_j >= 0 && k-this->min_k >= 0);
//     //         #endif
//     //         return g_data(i-this->min_i,j-this->min_j,k-this->min_k);
//     //     }
//     // }
//     // be careful this operator will conflict with the const one
//     Grid& operator()(int i, int j, int k) {
//         if(this->out_boundary(i,j,k)) {
//             std::cout << i << " " << j << " " << k << " is out of boundary, size: " << g_i << " " << g_j << " " << g_k << " min: " << min_i << " " << min_j << " " << min_k  << std::endl;
//             assert(false);
//         }
//         return g_data(i-this->min_i,j-this->min_j,k-this->min_k);
//     }
//     void assign(int i, int j, int k, const Grid& g) {
//         if(this->out_boundary(i,j,k)) {
//             std::cout << i << " " << j << " " << k << " is out of boundary, size: " << g_i << " " << g_j << " " << g_k << " min: " << min_i << " " << min_j << " " << min_k  << std::endl;
//             assert(false);
//         }
//         #ifdef MODE_DEBUG
//         assert(i-this->min_i >= 0 && j-this->min_j >= 0 && k-this->min_k >= 0);
//         #endif
//         g_data(i-this->min_i,j-this->min_j,k-this->min_k) = g;
//     }

//     bool out_boundary(int i, int j, int k) const {
//         const int& ii = i - this->min_i;
//         const int& jj = j - this->min_j;
//         const int& kk = k - this->min_k;
//         if(kk >= g_k || kk < 0 || jj >= g_j || jj < 0 || ii >= g_i || ii < 0) {
//             return true;
//         }
// 		return false;
//     }
//     int dim_1() const { return g_i; }
// 	int dim_2() const { return g_j; }
//     int dim_3() const { return g_k; }
// };



// struct Interacting_Pair {
//     const Atom_Index i;
//     const Atom_Index j;
//     Interacting_Pair(const Atom_Index ii, const Atom_Index jj) : i(ii), j(jj) {}
// };
// using Interacting_Pairs = std::vector<Interacting_Pair>;
// //precaluate interacting pairs
// Atom_Index num_of_ligand_interacing_pairs = 0;
// const Atom_Index num_lig_atom = lig.atom_num();
// const Atoms& lig_atoms = lig.get_atoms_reference();
// this->ligand_pairs.resize((num_lig_atom*(num_lig_atom-1))/2);
// for(Atom_Index i = 0; i < num_lig_atom; ++i) {
//     for(Atom_Index j = i+1; j < num_lig_atom; ++j) {
//         ++num_of_ligand_interacing_pairs;
//         const Float rij = lig_atoms[i].get_sybyl_vdw_radius() + lig_atoms[j].get_sybyl_vdw_radius();
//         const Float eij = std::sqrt(lig_atoms[i].get_sybyl_well_depth()*lig_atoms[j].get_sybyl_well_depth());
//         this->ligand_pairs[num_of_ligand_interacing_pairs-1] = std::move(Interacting_Pair(i,j,eij,rij));
//     }
// }

// class Interaction_Matrix {

//     Triangular_Matrix<Float> vdw_ACOEFF_tri_matrix;
//     Triangular_Matrix<Float> vdw_BCOEFF_tri_matrix;
//     Triangular_Matrix<Float> vdw_characteristic_dis_tri_matrix;
//     Triangular_Matrix<Float> vdw_well_depth_tri_matrix;
//     Strictly_Triangular_Matrix<BOND_TYPE> bond_

// public:
//     Interaction_Matrix(const RNA& r, const Ligand& l) : vdw_ACOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_BCOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_characteristic_dis_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_well_depth_tri_matrix(r.atom_num()+l.atom_num(),0.0) {
//         Atom_Index rna_atom_num = r.atom_num();
//         Atom_Index lig_atom_num = l.atom_num();
//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = 0; j < rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = r.get_atoms_reference()[j];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = rna_atom_num; i < lig_atom_num+rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = l.get_atoms_reference()[i-rna_atom_num];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }


//     }

// };

// // enum INTERACTION_TYPE {ELECTROSTATIC,LJ,SOLVATION,POLARIZATION};
// class Interaction_Matrix {

//     // Triangular_Matrix<Float> vdw_ACOEFF_tri_matrix;
//     // Triangular_Matrix<Float> vdw_BCOEFF_tri_matrix;
//     // Triangular_Matrix<Float> vdw_characteristic_dis_tri_matrix;
//     // Triangular_Matrix<Float> vdw_well_depth_tri_matrix;
//     Strictly_Triangular_Matrix<BOND_TYPE> bond_strict_tri_matrix;

// public:
//     Interaction_Matrix(const RNA& r, const Ligand& l) : vdw_ACOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_BCOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_characteristic_dis_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_well_depth_tri_matrix(r.atom_num()+l.atom_num(),0.0) {
//         Atom_Index rna_atom_num = r.atom_num();
//         Atom_Index lig_atom_num = l.atom_num();
//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = 0; j < rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = r.get_atoms_reference()[j];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = rna_atom_num; i < lig_atom_num+rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = l.get_atoms_reference()[i-rna_atom_num];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }


//     }

// };

}
