#pragma once

#include <iostream>
// #include <fstream>
// #include <sstream>
// #include <algorithm>
// #include <functional>
#include <vector>
#include <array>
#include <map>
#include <string>
// #include <iomanip>
// #include <chrono>
#include <cmath>
#include "common.h"
#include "parameter.h"
#include "matrix.h"
#include "mol.h"
// #include "cell.h"

namespace rlscore {

struct Complex {
  public:
    Size_Type lig_heavy_atom_num = 0;
    Size_Type rec_heavy_atom_num = 0;
    Matrix<Float> lig2rec; //[LIG][RNA]
    Matrix<Float> lig2lig;
    Matrix<Float> rec2rec;

    std::vector<std::vector<Size_Type>> sasa_neighboring_lig_indices_for_rec_atoms;
    std::vector<std::vector<Size_Type>> sasa_neighboring_rec_indices_for_lig_atoms;
    std::vector<std::vector<Size_Type>> sasa_neighboring_rec_indices_for_rec_atoms;
    std::vector<std::vector<Size_Type>> sasa_neighboring_lig_indices_for_lig_atoms;

    struct Born {
        Float alpha = 0.0;
        Float intb = 0.0;
    };
    std::vector<Born> rec_born_list_before;
    std::vector<Born> rec_born_list_after;
    std::vector<Born> lig_born_list_before;
    std::vector<Born> lig_born_list_after;

    struct Terms {
        Float U_lig_lj_before = 0.0;
        Float U_lig_lj_after = 0.0;

        Float U_lig_ele_before = 0.0;
        Float U_lig_ele_after = 0.0;

        Float lig_sasa_before = 0.0;
        Float lig_sasa_after = 0.0;
        Float delta_lig_sasa = 0.0;
        Float delta_rec_sasa = 0.0;
        Float delta_sasa = 0.0;

        Float U_rec_mutual_pol_before = 0.0;
        // Float rec_mutual_pol_after = 0.0;
        Float U_rec_self_pol_before = 0.0;
        Float U_rec_self_pol_after = 0.0;

        Float U_lig_mutual_pol_before = 0.0;
        // Float lig_mutual_pol_after = 0.0;
        Float U_lig_self_pol_before = 0.0;
        Float U_lig_self_pol_after = 0.0;
        Float U_complex_mutual_pol = 0.0;

        Float deltaU_inter_stack = 0.0;
        Float deltaU_inter_lj = 0.0;
        Float deltaU_inter_ele = 0.0;
        Float deltaU_inter_hbond = 0.0;
        Float deltaU_mutual_pol = 0.0;
        Float deltaU_rec_self_pol = 0.0;
        Float deltaU_lig_self_pol = 0.0;
        Float deltaU_lig_ele = 0.0;
        Float deltaU_lig_lj = 0.0;
        Float deltaU_sasa = 0.0;
    } terms;

    void initialize_receptor_and_ligand(const Parameter& param, const Receptor& rec, const Ligand& lig) {
        this->sasa_neighboring_rec_indices_for_rec_atoms = std::vector<std::vector<Size_Type>>(rec.atoms.size());
        this->sasa_neighboring_lig_indices_for_lig_atoms = std::vector<std::vector<Size_Type>>(lig.atoms.size());
        this->sasa_neighboring_lig_indices_for_rec_atoms = std::vector<std::vector<Size_Type>>(rec.atoms.size());
        this->sasa_neighboring_rec_indices_for_lig_atoms = std::vector<std::vector<Size_Type>>(lig.atoms.size());
        this->rec2rec.resize(rec.atoms.size(), rec.atoms.size(), 0.0);
        this->lig2lig.resize(lig.atoms.size(), lig.atoms.size(), 0.0);
        this->lig2rec.resize(lig.atoms.size(), rec.atoms.size(), 0.0);
        this->rec_heavy_atom_num = rec.atoms.size();
        this->lig_heavy_atom_num = lig.atoms.size();
        this->lig_born_list_before.resize(lig.atoms.size(), Born());
        this->lig_born_list_after.resize(lig.atoms.size(), Born());
        this->rec_born_list_before.resize(rec.atoms.size(), Born());
        this->rec_born_list_after.resize(rec.atoms.size(), Born());

        // for rec
        for(const auto& a : rec.atoms) {
            for(const auto& b : rec.atoms) {
                if(a.index < b.index) {
                    const Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                    this->rec2rec(b.index, a.index) = dis;
                    this->rec2rec(a.index, b.index) = dis;
                    if(dis < (a.radius+b.radius+2*param.r_water)) {
                        this->sasa_neighboring_rec_indices_for_rec_atoms[a.index].push_back(b.index);
                        this->sasa_neighboring_rec_indices_for_rec_atoms[b.index].push_back(a.index);
                    }
                }
            }
        }
        //for lig
        for(const auto& a : lig.atoms) {
            for(Size_Type b_index = 0; b_index < a.index; b_index++) {
                const Atom& b = lig.atoms[b_index];
                const Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                this->lig2lig(b.index, a.index) = dis;
                this->lig2lig(a.index, b.index) = dis;
                if(dis < (a.radius+b.radius+2*param.r_water)) {
                    this->sasa_neighboring_lig_indices_for_lig_atoms[a.index].push_back(b.index);
                    this->sasa_neighboring_lig_indices_for_lig_atoms[b.index].push_back(a.index);
                }
            }
        }

        // Cal born before docking
        this->cal_born_radius_before_docking<Ligand>(param, lig, this->lig_born_list_before, this->lig2lig);
        this->cal_born_radius_before_docking<Receptor>(param, rec, this->rec_born_list_before, this->rec2rec);
        this->cal_pol_before<Ligand>(param, lig, this->lig_born_list_before, this->lig2lig, this->terms.U_lig_mutual_pol_before, this->terms.U_lig_self_pol_before);
        this->cal_pol_before<Receptor>(param, rec, this->rec_born_list_before, this->rec2rec, this->terms.U_rec_mutual_pol_before, this->terms.U_rec_self_pol_before);

        //cal reference lig lj and sasa
        this->cal_lig_lj(param, lig, this->terms.U_lig_lj_before);
        this->cal_lig_ele(param, lig, this->terms.U_lig_ele_before);
        //ligand sasa before
        this->cal_mol_sasa(param, lig, this->sasa_neighboring_lig_indices_for_lig_atoms, this->terms.lig_sasa_before);
    }

    void update_ligand(const Parameter& param, const Receptor& rec, const Ligand& lig) {
        std::fill(this->sasa_neighboring_lig_indices_for_rec_atoms.begin(), this->sasa_neighboring_lig_indices_for_rec_atoms.end(), std::vector<Size_Type>());
        std::fill(this->sasa_neighboring_rec_indices_for_lig_atoms.begin(), this->sasa_neighboring_rec_indices_for_lig_atoms.end(), std::vector<Size_Type>());
        std::fill(this->sasa_neighboring_lig_indices_for_lig_atoms.begin(), this->sasa_neighboring_lig_indices_for_lig_atoms.end(), std::vector<Size_Type>());
        for(const auto& a : rec.atoms) {
            for(const auto& b : lig.atoms) {
                const Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                this->lig2rec(b.index, a.index) = dis;
                if(dis < (a.radius+b.radius+2*param.r_water)) {
                    this->sasa_neighboring_lig_indices_for_rec_atoms[a.index].push_back(b.index);
                    this->sasa_neighboring_rec_indices_for_lig_atoms[b.index].push_back(a.index);
                }
            }
        }
        for(const auto& a : lig.atoms) {
            for(Size_Type b_index = 0; b_index < a.index; b_index++) {
                const Atom& b = lig.atoms[b_index];
                const Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                this->lig2lig(b.index, a.index) = dis;
                this->lig2lig(a.index, b.index) = dis;
                if(dis < (a.radius+b.radius+2*param.r_water)) {
                    this->sasa_neighboring_lig_indices_for_lig_atoms[a.index].push_back(b.index);
                    this->sasa_neighboring_lig_indices_for_lig_atoms[b.index].push_back(a.index);
                }
            }
        }
        std::fill(this->lig_born_list_after.begin(), this->lig_born_list_after.end(), Born());
        std::fill(this->rec_born_list_after.begin(), this->rec_born_list_after.end(), Born());
        this->cal_rec_born_radius_after_docking(param, rec, lig, this->rec_born_list_before, this->lig2rec, this->rec_born_list_after);
        this->cal_lig_born_radius_after_docking(param, lig, rec, this->lig2lig, this->lig2rec, this->lig_born_list_after);
        this->cal_pol_after(param, rec, lig);
        this->terms.deltaU_mutual_pol = this->terms.U_complex_mutual_pol - this->terms.U_rec_mutual_pol_before - this->terms.U_lig_mutual_pol_before;
        this->terms.deltaU_rec_self_pol = this->terms.U_rec_self_pol_after - this->terms.U_rec_self_pol_before;
        this->terms.deltaU_lig_self_pol = this->terms.U_lig_self_pol_after - this->terms.U_lig_self_pol_before;

        // RNA delta sasa upon binding
        this->cal_delta_rec_sasa(param, rec, lig, this->sasa_neighboring_rec_indices_for_rec_atoms, this->sasa_neighboring_lig_indices_for_rec_atoms, this->terms.delta_rec_sasa);
        //cal lig lj and sasa after
        this->cal_lig_lj(param, lig, this->terms.U_lig_lj_after);
        this->cal_lig_ele(param, lig, this->terms.U_lig_ele_after);
        this->cal_mol_sasa(param, lig, rec, this->sasa_neighboring_lig_indices_for_lig_atoms, this->sasa_neighboring_rec_indices_for_lig_atoms, this->terms.lig_sasa_after);
        this->terms.deltaU_lig_lj = this->terms.U_lig_lj_after - this->terms.U_lig_lj_before;
        this->terms.delta_lig_sasa = this->terms.lig_sasa_after - this->terms.lig_sasa_before;
        this->terms.deltaU_lig_ele = this->terms.U_lig_ele_after - this->terms.U_lig_ele_before;
        this->terms.deltaU_sasa = param.gamma_SASA*(this->terms.delta_rec_sasa + this->terms.delta_lig_sasa);

        //cal intermolecular energy terms
        this->cal_inter_stack(param, rec, lig, this->terms.deltaU_inter_stack);
        this->cal_inter_lj(param, rec, lig, this->terms.deltaU_inter_lj);
        this->cal_inter_ele(param, rec, lig, this->terms.deltaU_inter_ele);
        this->cal_inter_hbond(param, rec, lig, this->terms.deltaU_inter_hbond);
    }

  private:
    void cal_inter_stack(const Parameter& param, const Receptor& rec, const Ligand& lig, Float& inter_stack)  {
        Float inter_stack_energy = 0.0;
        for(const auto& a : lig.atoms) {
            if(!a.is_aromatic) {
                continue;
            }
            for(const auto& b : rec.atoms) {
                if(!b.is_aromatic) {
                    continue;
                }
                const Float& dis = this->lig2rec(a.index, b.index);
                if(param.stack_lower_bound < dis && dis < param.stack_upper_bound) {
                    const Float tmp_dis = (dis-param.stack_lower_bound);
                    const Size_Type i_stack = tmp_dis/param.stack_step;
                    Float tmp_score = 0.0;
                    if(a.type_index == 1 && b.type_index == 1) {
                        tmp_score = (param.CC_stack_score[i_stack+1]-param.CC_stack_score[i_stack])*(tmp_dis/param.stack_step-Float(i_stack)) + param.CC_stack_score[i_stack];
                    } else if(a.type_index == 1 && b.type_index == 2) {
                        tmp_score = (param.CN_stack_score[i_stack+1]-param.CN_stack_score[i_stack])*(tmp_dis/param.stack_step-Float(i_stack)) + param.CN_stack_score[i_stack];
                    } else if(a.type_index == 2 && b.type_index == 1) {
                        tmp_score = (param.NC_stack_score[i_stack+1]-param.NC_stack_score[i_stack])*(tmp_dis/param.stack_step-Float(i_stack)) + param.NC_stack_score[i_stack];
                    } else if(a.type_index == 2 && b.type_index == 2) {
                        tmp_score = (param.NN_stack_score[i_stack+1]-param.NN_stack_score[i_stack])*(tmp_dis/param.stack_step-Float(i_stack)) + param.NN_stack_score[i_stack];
                    }
                    // std::cout << a.index << " " << b.index << " " << a.type_index << " " << b.type_index << " dis=" << dis << " stack_i=" << i_stack << " tmp_score=" << tmp_score << std::endl;
                    // if(tmp_score < 0) {
                        inter_stack_energy += tmp_score;
                    // }
                }
            }
        }
        inter_stack = param.rt_to_kcal_per_mol*inter_stack_energy*param.stack_scale; // 1 RT = 0.593 kcal/mol at T=298K
    }

    void cal_inter_lj(const Parameter& param, const Receptor& rec, const Ligand& lig, Float& inter_lj) {
        Float inter_lj_energy = 0.0;
        for(const auto& a : lig.atoms) {
            for(const auto& b : rec.atoms) {
                const Float equ_dis = a.radius + b.radius;
                if(this->lig2rec(a.index, b.index) < param.LJ_cut*equ_dis) {
                    const Size_Type iLJ_dis = this->lig2rec(a.index, b.index)/(param.LJ_cut*equ_dis)*param.LJstep+1;
                    inter_lj_energy += param.lj_table(a.type_index, b.type_index)[iLJ_dis];
                }
            }
        }
        inter_lj = inter_lj_energy;
    }

    void cal_lig_lj(const Parameter& param, const Ligand& lig, Float& lig_lj) {
        Float lig_lj_energy = 0.0;
        for(Size_Type a_id = 0; a_id < lig.atoms.size(); a_id++) {
            const auto& a = lig.atoms[a_id];
            for(Size_Type b_id = 0; b_id < a_id; b_id++) {
                const auto& b = lig.atoms[b_id];
                bool cal_flag = true;
                if(a.bonded_atom_ids.find(b.index) != a.bonded_atom_ids.end()) {
                    cal_flag = false;
                }
                const Float equ_dis = a.radius + b.radius;
                if(cal_flag && this->lig2lig(a.index, b.index) < param.LJ_cut*equ_dis) {
                    const Size_Type iLJ_dis = this->lig2lig(a.index, b.index)/(param.LJ_cut*equ_dis)*param.LJstep+1;
                    lig_lj_energy += param.lj_table(a.type_index, b.type_index)[iLJ_dis];
                }
            }
        }
        lig_lj = lig_lj_energy;
    }

    void cal_lig_ele(const Parameter& param, const Ligand& lig, Float& lig_ele) {
        Float lig_ele_energy = 0.0;
        for(Size_Type a_id = 0; a_id < lig.atoms.size(); a_id++) {
            const auto& a = lig.atoms[a_id];
            for(Size_Type b_id = 0; b_id < a_id; b_id++) {
                const auto& b = lig.atoms[b_id];
                bool cal_flag = true;
                if(a.bonded_atom_ids.find(b.index) != a.bonded_atom_ids.end()) {
                    cal_flag = false;
                }
                if(cal_flag) {
                    lig_ele_energy +=  param.lb0/param.eps_rna*a.charge*b.charge/this->lig2lig(a.index, b.index);
                }
            }
        }
        lig_ele = lig_ele_energy;
    }

    void cal_inter_ele(const Parameter& param, const Receptor& rec, const Ligand& lig, Float& inter_ele) {
        Float inter_ele_energy = 0.0;
        for(const auto& a : lig.atoms) {
            for(const auto& b : rec.atoms) {
                inter_ele_energy +=  param.lb0/param.eps_rna*a.charge*b.charge/this->lig2rec(a.index, b.index);
            }
        }
        inter_ele = inter_ele_energy;
    }

    void cal_inter_hbond(const Parameter& param, const Receptor& rec, const Ligand& lig, Float& inter_hb) {
        Float inter_hb_energy = 0.0;
        for(const auto& a : rec.atoms) {
            for(const auto& b : lig.atoms) {
                if( (a.is_don && b.is_acc) || (a.is_acc && b.is_don) ) {
                // if( a.num_h != 0 || b.num_h != 0 ) {
                    const Float equ_dis = param.Hbond_max*(a.radius+b.radius);
                    if(this->lig2rec(b.index, a.index) < equ_dis) {
                        const Size_Type iLJ_dis = (this->lig2rec(b.index, a.index)/equ_dis)*param.LJstep+1;
                        inter_hb_energy += -param.hbond_table(a.type_index, b.type_index)[iLJ_dis];
                    }
                }
            }
        }
        inter_hb = inter_hb_energy;
    }

    template<typename T>
    void cal_mol_sasa(const Parameter& param, const T& target_mol, const std::vector<std::vector<Size_Type>>& intra_near_indices_list, Float& mol_sasa) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : target_mol.atoms) {
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(intra_near_indices.size() != 0) {
                Size_Type effective_count = 0;
                for(const auto& point : param.sasa_points[a.type_index]) {
                    const Vec3d surface_point = point + a.xyz;
                    bool buried = false;
                    for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                        const Atom& intra_near_atom = target_mol.atoms[intra_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                        if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                            buried = true;
                            break;
                        }
                    }
                    if(buried == true) {
                        continue;
                    }
                    effective_count++;
                }
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.type_index].size());
            } else {
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water);
            }
        }
        mol_sasa = SumSASAtmp;
    }

    void cal_delta_rec_sasa(const Parameter& param, const Receptor& rec, const Ligand& lig, const std::vector<std::vector<Size_Type>>& intra_near_indices_list, const std::vector<std::vector<Size_Type>>& inter_near_indices_list, Float& delta_sasa) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : rec.atoms) {
            const auto& inter_near_indices = inter_near_indices_list[a.index];
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(inter_near_indices.size() == 0) {
                continue;
            }
            Size_Type effective_count = 0;
            for(const auto& point : param.sasa_points[a.type_index]) {
                const Vec3d surface_point = point + a.xyz;
                bool buried = false;
                for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                    const Atom& intra_near_atom = rec.atoms[intra_near_indices[l]];
                    const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                    if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                        buried = true;
                        break;
                    }
                }
                if(buried == true) {
                    continue;
                }
                for(Size_Type l=0; l != inter_near_indices.size(); l++) {
                    const Atom& inter_near_atom = lig.atoms[inter_near_indices[l]];
                    const Float dis_square = vec3d_distance_square(surface_point, inter_near_atom.xyz);
                    if(dis_square < ((inter_near_atom.radius+param.r_water)*(inter_near_atom.radius+param.r_water))) {
                        effective_count++;
                        break;
                    }
                }
            }
            SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.type_index].size());
        }
        delta_sasa = -SumSASAtmp;
    }

    template<typename T1, typename T2>
    void cal_mol_sasa(const Parameter& param, const T1& target_mol, const T2& near_mol, const std::vector<std::vector<Size_Type>>& intra_near_indices_list, const std::vector<std::vector<Size_Type>>& inter_near_indices_list, Float& mol_sasa) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : target_mol.atoms) {
            const auto& inter_near_indices = inter_near_indices_list[a.index];
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(inter_near_indices.size() == 0 and intra_near_indices.size() == 0) {
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water);
            } else {
                Size_Type effective_count = 0;
                for(const auto& point : param.sasa_points[a.type_index]) {
                    const Vec3d surface_point = point + a.xyz;
                    bool intra_buried = false;
                    for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                        const Atom& intra_near_atom = target_mol.atoms[intra_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                        if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                            intra_buried = true;
                            break;
                        }
                    }
                    if(intra_buried == true) {
                        continue;
                    }
                    bool inter_buried = false;
                    for(Size_Type l=0; l != inter_near_indices.size(); l++) {
                        const Atom& inter_near_atom = near_mol.atoms[inter_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, inter_near_atom.xyz);
                        if(dis_square < ((inter_near_atom.radius+param.r_water)*(inter_near_atom.radius+param.r_water))) {
                            inter_buried = true;
                            break;
                        }
                    }
                    if(inter_buried) {
                        continue;
                    }
                    effective_count++;
                }
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.type_index].size());
            }
        }
        mol_sasa = SumSASAtmp;
    }


    template<typename T>
    void cal_born_radius_before_docking(const Parameter& param, const T& self_mol, std::vector<Born>& self_born_list, const Matrix<Float>& dis_matrix) {
        for(const auto& a : self_mol.atoms) {
            Float intb = 0.0;
            for(const auto& b : self_mol.atoms) {
                if(a.index!=b.index) {
                    Float Lij, Uij;
                    const Float& dis = dis_matrix(a.index, b.index);
                    if(a.radius >= (dis+b.born_scale*b.radius)) {
                        Lij = 1.0;
                        Uij = 1.0;
                    } else {
                        if(a.radius > (dis-b.born_scale*b.radius)) {
                            Lij = a.radius;
                        } else {
                            Lij = dis-b.born_scale*b.radius;
                        }
                        Uij = dis+b.born_scale*b.radius;
                    }
                    intb += 0.5*((1.0/Lij-1.0/Uij)+(b.born_scale*b.radius*b.born_scale*b.radius/(4.0*dis)-dis/4.0)*(1.0/(Lij*Lij)-1.0/(Uij*Uij))+(1.0/(2.0*dis))*std::log(Lij/Uij));
                }
            }
            self_born_list[a.index].alpha = 1.0/(1.0/a.radius-intb);
            self_born_list[a.index].intb = intb;
        }
    }

    template<typename T>//cal pol for rigid RNA and reference ligand conformation
    void cal_pol_before(const Parameter& param, const T& self_mol, const std::vector<Born>& self_born_list, const Matrix<Float>& dis_matrix, Float& mutual_pol, Float& self_pol) {
        // const Float pol_prefix = param.lb0*(1.0/param.eps_h2o-1.0/param.eps_rna);
        const Float pol_rna_prefix = param.lb0/param.eps_rna;
        const Float pol_h2o_prefix = param.lb0/param.eps_h2o;
        Float pol = 0.0;
        for(const auto& a : self_mol.atoms) {
            for(Size_Type i = a.index+1; i != self_mol.atoms.size(); ++i) {
                const Atom& b = self_mol.atoms[i];
                const Float& a_born_r = self_born_list[a.index].alpha;
                const Float& b_born_r = self_born_list[b.index].alpha;
                const Float ab_born_rr = a_born_r*b_born_r;
                const Float& dis = dis_matrix(b.index, a.index);
                const Float dis_square = dis*dis;
                const Float f_GB = std::sqrt(dis_square+ab_born_rr*std::exp(-dis_square/(4.0*ab_born_rr)));
                pol += (pol_h2o_prefix*std::exp(-param.debye_kappa*f_GB)-pol_rna_prefix)*a.charge*b.charge/f_GB;
            }
        }
        mutual_pol = pol;

        // self part before docking
        Float self = 0.0;
        for(const auto& a : self_mol.atoms) {
            const Float& born_r = self_born_list[a.index].alpha;
            self += 0.5*(pol_h2o_prefix*std::exp(-param.debye_kappa*born_r)-pol_rna_prefix)*a.charge*a.charge*(1.0/born_r);
            if(born_r == 0) {
                std::cout << a.index << " born_before " << born_r << std::endl;
                exit(2);
            }
        }
        self_pol = self;
    }

    void cal_rec_born_radius_after_docking(const Parameter& param, const Receptor& rec, const Ligand& lig, const std::vector<Born>& self_born_list_before, const Matrix<Float>& dis_matrix, std::vector<Born>& self_born_list_after) {
        for(const auto& a : rec.atoms) {
            Float intb = self_born_list_before[a.index].intb;
            for(const auto& b : lig.atoms) {
                Float Lij, Uij;
                const Float& dis = dis_matrix(b.index, a.index);
                if(a.radius >= (dis+b.born_scale*b.radius)) {
                    Lij = 1.0;
                    Uij = 1.0;
                } else {
                    if(a.radius > (dis-b.born_scale*b.radius)) {
                        Lij = a.radius;
                    } else {
                        Lij = dis-b.born_scale*b.radius;
                    }
                    Uij = dis+b.born_scale*b.radius;
                }
                intb += 0.5*((1.0/Lij-1.0/Uij)+(b.born_scale*b.radius*b.born_scale*b.radius/(4.0*dis)-dis/4.0)*(1.0/(Lij*Lij)-1.0/(Uij*Uij))+(1.0/(2.0*dis))*std::log(Lij/Uij));
            }
            self_born_list_after[a.index].alpha = 1.0/(1.0/a.radius-intb);
            self_born_list_after[a.index].intb = intb;
        }
    }

    void cal_lig_born_radius_after_docking(const Parameter& param, const Ligand& lig, const Receptor& rec, const Matrix<Float>& lig2lig_dis, const Matrix<Float>& lig2rec_dis, std::vector<Born>& self_born_list_after) {
        for(const auto& a : lig.atoms) {
            Float intb = 0.0;
            //lig-lig
            for(const auto& b : lig.atoms) {
                if(a.index!=b.index) {
                    Float Lij, Uij;
                    const Float& dis = lig2lig_dis(a.index, b.index);
                    if(a.radius >= (dis+b.born_scale*b.radius)) {
                        Lij = 1.0;
                        Uij = 1.0;
                    } else {
                        if(a.radius > (dis-b.born_scale*b.radius)) {
                            Lij = a.radius;
                        } else {
                            Lij = dis-b.born_scale*b.radius;
                        }
                        Uij = dis+b.born_scale*b.radius;
                    }
                    intb += 0.5*((1.0/Lij-1.0/Uij)+(b.born_scale*b.radius*b.born_scale*b.radius/(4.0*dis)-dis/4.0)*(1.0/(Lij*Lij)-1.0/(Uij*Uij))+(1.0/(2.0*dis))*std::log(Lij/Uij));
                }
            }
            // lig-rec
            for(const auto& b : rec.atoms) {
                Float Lij, Uij;
                const Float& dis = lig2rec_dis(a.index, b.index);
                if(a.radius >= (dis+b.born_scale*b.radius)) {
                    Lij = 1.0;
                    Uij = 1.0;
                } else {
                    if(a.radius > (dis-b.born_scale*b.radius)) {
                        Lij = a.radius;
                    } else {
                        Lij = dis-b.born_scale*b.radius;
                    }
                    Uij = dis+b.born_scale*b.radius;
                }
                intb += 0.5*((1.0/Lij-1.0/Uij)+(b.born_scale*b.radius*b.born_scale*b.radius/(4.0*dis)-dis/4.0)*(1.0/(Lij*Lij)-1.0/(Uij*Uij))+(1.0/(2.0*dis))*std::log(Lij/Uij));
            }
            self_born_list_after[a.index].alpha = 1.0/(1.0/a.radius-intb);
            self_born_list_after[a.index].intb = intb;
        }
    }

    void cal_pol_after(const Parameter& param, const Receptor& rec, const Ligand& lig) {
        // const Float pol_prefix = param.lb0*(1.0/param.eps_h2o-1.0/param.eps_rna);
        const Float pol_rna_prefix = param.lb0/param.eps_rna;
        const Float pol_h2o_prefix = param.lb0/param.eps_h2o;
        // rec-rec mutual part after docking
        Float rec_mutual_pol = 0.0;
        for(const auto& a : rec.atoms) {
            for(Size_Type i = a.index+1; i != rec.atoms.size(); ++i) {
                const Atom& b = rec.atoms[i];
                const Float& a_born_r = this->rec_born_list_after[a.index].alpha;
                const Float& b_born_r = this->rec_born_list_after[b.index].alpha;
                const Float ab_born_rr = a_born_r*b_born_r;
                const Float& dis = this->rec2rec(b.index, a.index);
                const Float dis_square = dis*dis;
                const Float f_GB = std::sqrt(dis_square+ab_born_rr*std::exp(-dis_square/(4.0*ab_born_rr)));
                rec_mutual_pol += (pol_h2o_prefix*std::exp(-param.debye_kappa*f_GB)-pol_rna_prefix)*a.charge*b.charge/f_GB;
            }
        }
        // this->terms.rec_mutual_pol_after = rec_mutual_pol;

        // rec self part after docking
        Float rec_self_pol = 0.0;
        for(const auto& a : rec.atoms) {
            const Float& born_r = this->rec_born_list_after[a.index].alpha;
            rec_self_pol += 0.5*(pol_h2o_prefix*std::exp(-param.debye_kappa*born_r)-pol_rna_prefix)*a.charge*a.charge*(1.0/born_r);
            if(born_r == 0) {
                std::cout << a.index << " born_after " << born_r << std::endl;
                exit(2);
            }
        }
        this->terms.U_rec_self_pol_after = rec_self_pol;

        // lig-lig mutual part after docking
        Float lig_mutual_pol = 0.0;
        for(const auto& a : lig.atoms) {
            for(Size_Type i = a.index+1; i != lig.atoms.size(); ++i) {
                const Atom& b = lig.atoms[i];
                const Float& dis = this->lig2lig(b.index, a.index);
                const Float dis_square = dis*dis;
                const Float& a_born_r_after = this->lig_born_list_after[a.index].alpha;
                const Float& b_born_r_after = this->lig_born_list_after[b.index].alpha;
                const Float ab_born_rr_after = a_born_r_after*b_born_r_after;
                const Float f_GB = std::sqrt(dis_square+ab_born_rr_after*std::exp(-dis_square/(4.0*ab_born_rr_after)));
                lig_mutual_pol += (pol_h2o_prefix*std::exp(-param.debye_kappa*f_GB)-pol_rna_prefix)*a.charge*b.charge/f_GB;

                // if(std::isnan(pol0) || std::isnan(pol1)) {
                // if(a_born_r_after < 0 || b_born_r_after < 0) {
                //     std::cout << a.index << " " << b.index << " " << dis << " " << a_born_r_before << " " << b_born_r_before << " " << a_born_r_after << " " << b_born_r_after << std::endl;
                //     exit(2);
                // }
            }
        }
        // this->terms.lig_mutual_pol_after = lig_mutual_pol;

        // lig self part after docking
        Float lig_self_pol = 0.0;
        for(const auto& a : lig.atoms) {
            const Float& born_r_after = this->lig_born_list_after[a.index].alpha;
            const Float self_prefix = 0.5*(pol_h2o_prefix*std::exp(-param.debye_kappa*born_r_after)-pol_rna_prefix)*a.charge*a.charge;
            lig_self_pol += self_prefix*(1.0/born_r_after);
            if(born_r_after == 0) {
                std::cout << a.index << " " << born_r_after << std::endl;
                exit(2);
            }
        }
        this->terms.U_lig_self_pol_after = lig_self_pol;

        //lig-rec part (only mutual part)
        Float lig_rec_mutual_pol = 0.0;
        for(const auto& a : rec.atoms) {
            for(const auto& b : lig.atoms) {
                const Float& a_born_r = this->rec_born_list_after[a.index].alpha;
                const Float& b_born_r = this->lig_born_list_after[b.index].alpha;
                const Float& dis = this->lig2rec(b.index, a.index);
                const Float dis_square = dis*dis;
                const Float ab_born_rr = a_born_r*b_born_r;
                const Float f_GB = std::sqrt(dis_square+ab_born_rr*std::exp(-dis_square/(4.0*ab_born_rr)));
                lig_rec_mutual_pol += (pol_h2o_prefix*std::exp(-param.debye_kappa*f_GB)-pol_rna_prefix)*a.charge*b.charge/f_GB;
            }
        }

        this->terms.U_complex_mutual_pol = lig_rec_mutual_pol + lig_mutual_pol + rec_mutual_pol;
    }
};

}