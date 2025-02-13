// #include <iostream>
#include <fstream>
#include <string>
#include <cassert>
// #include <vector>
// #include <regex>

#include "mol.h"

namespace rlscore {

void parse_single_pose_from_stream_each_time(std::ifstream& pose_stream, std::string& name, std::vector<Atom>& atoms) {
    std::string sline;
    Size_Type number_atoms, number_bonds;
    std::vector<Atom> atoms_temp;
    std::string name_temp;
    while(std::getline(pose_stream, sline)) {
        if(sline.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            std::getline(pose_stream, sline);
            name_temp = trim(sline);
            std::getline(pose_stream, sline);
            const std::vector<std::string> token = string2vector(sline);
            assert(token.size() >= 2);
            number_atoms = std::stoi(token[0]);
            number_bonds = std::stoi(token[1]);
        }
        if(sline.find("@<TRIPOS>ATOM") != std::string::npos) {
            Size_Type heavy_index = 0;
            for(Size_Type index = 0; index < number_atoms; ++index) {
                std::getline(pose_stream, sline);
                const std::vector<std::string> token = string2vector(sline);
                Atom mol2_atom;
                // mol2_atom.index = index;
                mol2_atom.serial = std::stoi(token[0]);
                mol2_atom.name = token[1];
                mol2_atom.xyz = Vec3d(std::stod(token[2]), std::stod(token[3]), std::stod(token[4]));
                mol2_atom.sybyl_type = token[5];
                mol2_atom.res_seq = std::stoi(token[6]);
                mol2_atom.res_name = token[7];
                mol2_atom.charge = std::stod(token[8]);
                if(mol2_atom.sybyl_type[0] != 'H') {
                    mol2_atom.index = heavy_index;
                    atoms_temp.push_back(mol2_atom);
                    heavy_index++;
                }
            }// for atom
        }

        if(name_temp != "" && atoms_temp.size() != 0) {
            name = name_temp;
            atoms = atoms_temp;
            // assert(atoms.size() == number_atoms);
            break; // break the while for this read
        }
    }
}

// read heavy atoms from mol2
void parse_single_mol2_file(const Parameter& param, const std::string file_path, std::string& name, std::vector<Atom>& atoms, const std::vector<Size_Type>& aromatic_indices, const std::vector<Size_Type>& h_don_indices, const std::vector<Size_Type>& h_acc_indices) {
    // read all atom part
    std::ifstream inmol2(file_path, std::ios::in);
    std::string sline;
    Size_Type number_atoms, number_bonds;
    std::vector<Atom> atoms_temp;
    std::map<Size_Type,Size_Type> index_to_heavy_index;
    while(std::getline(inmol2, sline)) {
        if(sline.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            std::getline(inmol2, sline);
            name = trim(sline);
            std::getline(inmol2, sline);
            const std::vector<std::string> token = string2vector(sline);
            assert(token.size() >= 2);
            number_atoms = std::stoi(token[0]);
            number_bonds = std::stoi(token[1]);
        }
        if(sline.find("@<TRIPOS>ATOM") != std::string::npos) {
            Size_Type heavy_index = 0;
            for(Size_Type index = 0; index < number_atoms; ++index) {
                std::getline(inmol2, sline);
                const std::vector<std::string> token = string2vector(sline);
                Atom mol2_atom;
                mol2_atom.index = index;
                mol2_atom.serial = std::stoi(token[0]);
                mol2_atom.name = token[1];
                mol2_atom.xyz = Vec3d(std::stod(token[2]), std::stod(token[3]), std::stod(token[4]));
                mol2_atom.sybyl_type = token[5];
                mol2_atom.res_seq = std::stoi(token[6]);
                mol2_atom.res_name = token[7];
                mol2_atom.charge = std::stod(token[8]);
                if(mol2_atom.sybyl_type[0] != 'H') {
                    index_to_heavy_index.insert({index, heavy_index});
                    heavy_index++;
                }
                atoms_temp.push_back(mol2_atom);
            }// for atom
        }
        // read bond, update degree, update num_h and merge hydrogen charge
        if(sline.find("@<TRIPOS>BOND") != std::string::npos) {
            for(Size_Type ib = 0; ib != number_bonds; ++ib) {
                std::getline(inmol2, sline);
                const std::vector<std::string> token = string2vector(sline);
                const Size_Type left_idx = std::stoi(token[1])-1;
                const Size_Type right_idx = std::stoi(token[2])-1;

                // get bonded atoms and heavy atom degree
                if(atoms_temp[left_idx].sybyl_type[0] != 'H' && atoms_temp[right_idx].sybyl_type[0] != 'H') {
                    atoms_temp[right_idx].degree++;
                    atoms_temp[left_idx].degree++;
                    atoms_temp[right_idx].bonded_atom_ids.insert(left_idx);
                    atoms_temp[left_idx].bonded_atom_ids.insert(right_idx);
                }
                // merge h
                if(atoms_temp[left_idx].sybyl_type[0] == 'H' && atoms_temp[right_idx].sybyl_type[0] != 'H') {
                    atoms_temp[right_idx].charge += atoms_temp[left_idx].charge;
                    atoms_temp[left_idx].charge = 0.0;
                    atoms_temp[right_idx].num_h++;
                } else if(atoms_temp[left_idx].sybyl_type[0] != 'H' && atoms_temp[right_idx].sybyl_type[0] == 'H') {
                    atoms_temp[left_idx].charge += atoms_temp[right_idx].charge;
                    atoms_temp[right_idx].charge = 0.0;
                    atoms_temp[left_idx].num_h++;
                }
            }
        }
    }// end reading
    inmol2.close();

    // update aromatic info
    for(const auto& aromatic_idx : aromatic_indices) {
        atoms_temp[aromatic_idx].is_aromatic = true;
    }

    // update h bond info
    for(const auto& h_don_idx : h_don_indices) {
        atoms_temp[h_don_idx].is_don = true;
    }
    for(const auto& h_acc_idx : h_acc_indices) {
        atoms_temp[h_acc_idx].is_acc = true;
    }

    // assign type, radius and born_scale to atom
    for(auto& a : atoms_temp) {
        switch(a.sybyl_type[0]) {
            case 'H': {
                a.type_index = 0;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            case 'C': {
                a.type_index = 1;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            case 'N': {
                a.type_index = 2;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            case 'O': {
                a.type_index = 3;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            case 'P': {
                a.type_index = 4;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            case 'S': {
                a.type_index = 5;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
            default: {
                //std::cout << "new atom type : " << a.type[0] << " ---> " << a.serial << " " << a.name << " set as P" << std::endl;
                a.type_index = 6;
                a.radius = param.radii[a.type_index];
                a.born_scale = param.born_scales[a.type_index];
            }
            break;
        }
    }

    Size_Type updated_index = 0;
    // extract heavy atoms from atoms_temp, and update index, and bonded ids
    for(const auto& a : atoms_temp) {
        if(a.type_index != 0) {
            atoms.push_back(a);
            auto& new_atom = atoms[updated_index];
            new_atom.index = index_to_heavy_index[a.index];
            new_atom.bonded_atom_ids.clear();
            for(const Size_Type b_id : a.bonded_atom_ids) {
                new_atom.bonded_atom_ids.insert(index_to_heavy_index[b_id]);
            }
            updated_index++;
        }
    }
}//end read mol2

}