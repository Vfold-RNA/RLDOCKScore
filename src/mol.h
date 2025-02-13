#pragma once

// #include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <map>
#include <vector>
#include <set>
// #include <regex>
// #include "atom.h"
// #include "molecule.h"
#include "common.h"
#include "parameter.h"
#include "str.h"
#include "vec3d.h"
// #include "vec2d.h"

namespace rlscore {

struct Atom {
    Size_Type index;
    int serial;
    std::string name;
    Vec3d xyz;
    std::string sybyl_type;
    int res_seq;
    std::string res_name;
    Float charge;
    Float radius;
    Float born_scale;
    Size_Type type_index;
    Size_Type num_h = 0;
    Size_Type degree = 0; //heavy atom degree
    std::set<Size_Type> bonded_atom_ids;
    // std::vector<Size_Type> hydrogen_ids;

    bool is_acc = false;
    bool is_don = false;
    bool is_aromatic = false;
};

void parse_single_pose_from_stream_each_time(std::ifstream& pose_stream, std::string& name, std::vector<Atom>& atoms);
// read atoms from mol2
void parse_single_mol2_file(const Parameter& param, const std::string file_path, std::string& name, std::vector<Atom>& atoms, const std::vector<Size_Type>& aromatic_indices, const std::vector<Size_Type>& h_don_indices, const std::vector<Size_Type>& h_acc_indices);

struct Receptor {
    std::string name;
    std::string file_path;
    std::vector<Atom> atoms;
    Receptor(const Parameter& param, const std::string fp, const std::vector<Size_Type>& aromatic_indices, const std::vector<Size_Type>& h_don_indices, const std::vector<Size_Type>& h_acc_indices) : 
    file_path(fp) {
        parse_single_mol2_file(param, this->file_path, this->name, this->atoms, aromatic_indices, h_don_indices, h_acc_indices);
        // for(auto& a : this->atoms) {
        //     // if((a.sybyl_type == "C.ar") || (a.sybyl_type == "N.ar")) {
        //     //     a.is_aromatic = true;
        //     // }
        //     // if((a.sybyl_type[0] == 'N')) {
        //     //     if(a.num_h != 0) {
        //     //         a.is_don = true;
        //     //     } else {
        //     //         if(a.degree < 3) { a.is_acc = true; }
        //     //     }
        //     // } else if(a.sybyl_type[0] == 'O') {
        //     //     if(a.num_h != 0) {
        //     //         a.is_don = true;
        //     //     } else {
        //     //         if(a.degree < 3) { a.is_acc = true; }
        //     //     }
        //     // }
        // }
    }
};

struct Ligand {
    std::string name;
    std::string file_path;
    std::vector<Atom> atoms;
    Float rot_dof = 0.0; // assigned in initialization
    Ligand(const Parameter& param, const std::string fp, const Float rd,  const std::vector<Size_Type>& aromatic_indices, const std::vector<Size_Type>& h_don_indices, const std::vector<Size_Type>& h_acc_indices) : 
    file_path(fp), rot_dof(rd) {
        parse_single_mol2_file(param, this->file_path, this->name, this->atoms, aromatic_indices, h_don_indices, h_acc_indices);
        // for(auto& a : this->atoms) {
        //     // if((a.sybyl_type == "C.ar") || (a.sybyl_type == "N.ar")) {
        //     //     a.is_aromatic = true;
        //     // }
        //     // if((a.sybyl_type[0] == 'N')) {
        //     //     if(a.num_h != 0) {
        //     //         a.is_don = true;
        //     //     } else {
        //     //         if(a.degree < 3) { a.is_acc = true; }
        //     //     }
        //     // } else if(a.sybyl_type[0] == 'O') {
        //     //     if(a.num_h != 0) {
        //     //         a.is_don = true;
        //     //     } else {
        //     //         if(a.degree < 3) { a.is_acc = true; }
        //     //     }
        //     // }
        // }
    }
};

}