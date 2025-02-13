// #include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
// #include <algorithm>
// #include <functional>
#include <vector>
#include <map>
#include <string>
// #include <cstring>
// #include <iomanip>
// #include <chrono>
#include <cmath>
#include <unistd.h> // for getopt
#include <getopt.h> // for getopt_long()
#include <sys/stat.h> // for checking if path exists

#include "config.h"
#include "common.h"
#include "parameter.h"
#include "mol.h"
#include "score.h"

#ifdef MODE_PROFILE
#include <gperftools/profiler.h>
#endif

// inline bool file_exists(const std::string& name) {
//     std::ifstream f(name);
//     return f.good();
// }

inline bool path_exists(const std::string& path_name, const char mode = 'p') {
    struct stat sb;
    if(stat(path_name.c_str(), &sb) == 0) {
        switch(mode) {
            case 'p':
                return true;
                break;
            case 'f':
                return S_ISREG(sb.st_mode);
                break;
            case 'd':
                return S_ISDIR(sb.st_mode);
                break;
            default:
                std::cout << "Unkown path mode: " << mode << std::endl;
                exit(2);
                break;
        }
    }
    return false;
}

//start main//////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
#ifdef MODE_PROFILE
ProfilerStart("/tmp/profile.out");
#endif

    const std::string usage_hint = "rldock_score_v2 -r in_rec -l in_ref_lig -p in_poses -i info_path -s score_mode [-t num_torsion] [-v] [-h] [--w_xx val] [--debye_length val]";
    int c;
    int verbose = 0; // 0 for false
    int help = 0; // 0 for false

    rlscore::Float rot_dof = 0.0;

    std::string in_rec_path;
    std::string in_ref_lig_path;
    std::string in_poses_path;
    std::string in_info_path;
    std::string in_score_mode;
    rlscore::Float debye_length = 11.0;
    rlscore::Float global_score_scale = 0.0;

    rlscore::Float w_lj       = 0.0; rlscore::Float w_ele    = 0.0; rlscore::Float w_pol     = 0.0;
    rlscore::Float w_lig_spol = 0.0; rlscore::Float w_sasa   = 0.0; rlscore::Float w_hbond   = 0.0;
    rlscore::Float w_rec_spol = 0.0; rlscore::Float w_lig_lj = 0.0; rlscore::Float w_lig_ele = 0.0;
    rlscore::Float w_rot      = 0.0; rlscore::Float w_stack  = 0.0;

    // first check argument of -s if -s exists, before parsing arguments
    for(int i = 0; i < argc; ++i) {
        if(std::string(argv[i]) == "-s") {
            if((i+1) < argc) {
                in_score_mode = std::string(argv[i+1]);
            }
            if(in_score_mode != "virtual_screen" && in_score_mode != "binding_mode") {
                std::cout << "-s requries an argument, it can be either virtual_screen or binding_mode." << std::endl;
                return 2;
            }
        }
    }

    if(in_score_mode == "virtual_screen") {
        // trained for both pose and affinity
        // for debye_length = 11.0; (HIV TAR Hashimi EF2%=25.0)
        global_score_scale = 0.35115732;
        w_lj       = 2.00000000; w_ele    =  0.97765769; w_pol     = 1.04964457;
        w_lig_spol = 1.51627838; w_sasa   = 10.00000000; w_hbond   = 0.20930519;
        w_rec_spol = 3.48643219; w_lig_lj =  0.72776193; w_lig_ele = 0.83557641;
        w_rot      = 0.02400000; w_stack  =  2.56815736;

        // for debye_length = 1000000.0; (HIV TAR Hashimi EF2%=25.0)
        // [0.61099581,1.55651817,2.12423939,3.30145327,15.96859012,0.58070387,4.63355800,8.09908147,2.62822413,0.02116081,4.16999834]
        // w_lj       = 0.61099581; w_ele    =  1.55651817; w_pol     = 2.12423939;
        // w_lig_spol = 3.30145327; w_sasa   = 15.96859012; w_hbond   = 0.58070387;
        // w_rec_spol = 4.63355800; w_lig_lj =  8.09908147; w_lig_ele = 2.62822413;
        // w_rot      = 0.02116081; w_stack  =  4.16999834;
    } else if(in_score_mode == "binding_mode") {
        // for debye_length = 11.0;
        // [2.47845781,0.97765769,1.04964457,1.51627838,8.27989567,0.20930519,3.48643219,0.72776193,0.83557641,0.02358925,2.56815736]
        global_score_scale = 0.38553958;
        w_lj       = 2.47845781; w_ele    = 0.97765769; w_pol     = 1.04964457;
        w_lig_spol = 1.51627838; w_sasa   = 8.27989567; w_hbond   = 0.20930519;
        w_rec_spol = 3.48643219; w_lig_lj = 0.72776193; w_lig_ele = 0.83557641;
        w_rot      = 0.02358925; w_stack  = 2.56815736;

        // for debye_length = 13.17454041;
        // [2.70351897,0.98922035,1.07559982,1.60863637,4.99995565,0.13021738,3.85680969,0.50317140,0.75873499,0.01212298,2.15269900]
        // w_lj       = 2.70351897; w_ele    = 0.98922035; w_pol     = 1.07559982;
        // w_lig_spol = 1.60863637; w_sasa   = 4.99995565; w_hbond   = 0.13021738;
        // w_rec_spol = 3.85680969; w_lig_lj = 0.50317140; w_lig_ele = 0.75873499;
        // w_rot      = 0.01212298; w_stack  = 2.15269900;
    }

    /////////////////////////////////////////

    int option_index = 0;
    struct option longopts[] = {
        { "receptor",     required_argument,  NULL,      'r' },
        { "ligand",       required_argument,  NULL,      'l' },
        { "poses",        required_argument,  NULL,      'p' },
        { "info",         required_argument,  NULL,      'i' },
        { "score",        required_argument,  NULL,      's' },

        { "verbose",      no_argument,        & verbose,  1  },
        { "help",         no_argument,        & help,     1  },

        { "torsion",      required_argument,  NULL,      't' },
        { "w_lj",         required_argument,  0,          0  },
        { "w_ele",        required_argument,  0,          0  },
        { "w_pol",        required_argument,  0,          0  },
        { "w_lig_spol",   required_argument,  0,          0  },
        { "w_rec_spol",   required_argument,  0,          0  },
        { "w_sasa",       required_argument,  0,          0  },
        { "w_hbond",      required_argument,  0,          0  },
        { "w_lig_lj",     required_argument,  0,          0  },
        { "w_lig_ele",    required_argument,  0,          0  },
        { "w_rot",        required_argument,  0,          0  },
        { "w_stack",      required_argument,  0,          0  },
        { "debye_length", required_argument,  0,          0  },
        {0, 0, 0, 0}
    };

    bool args_correct = true;
    // opterr, optopt, optind, optarg are reserved variables for getopt
    while( (c = getopt_long(argc, argv, ":r:l:p:i:s:t:vh", longopts, &option_index)) != -1 ) {
        switch(c) {
            case 'r':
                in_rec_path = std::string(optarg);
                break;
            case 'l':
                in_ref_lig_path = std::string(optarg);
                break;
            case 'p':
                in_poses_path = std::string(optarg);
                break;
            case 'i':
                in_info_path = std::string(optarg);
                break;
            case 's':
                in_score_mode = std::string(optarg);
                assert(in_score_mode == "virtual_screen" || in_score_mode == "binding_mode");
                break;
            case 't':
                rot_dof = std::stod(std::string(optarg));
                assert(rot_dof >= 0);
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                help = true;
                break;
            case 0: {
                /* getopt_long() set a variable, just keep going */
                const std::string arg_name = std::string(longopts[option_index].name);
                // std::cout << arg_name << std::endl;
                // std::cout << longopts[option_index].flag << std::endl;
                if(arg_name == "verbose" || arg_name == "help") {}
                else if(arg_name == "w_lj")        { w_lj = std::stod(std::string(optarg)); }
                else if(arg_name == "w_ele")       { w_ele = std::stod(std::string(optarg)); }
                else if(arg_name == "w_pol")       { w_pol = std::stod(std::string(optarg)); }
                else if(arg_name == "w_rec_spol")  { w_rec_spol = std::stod(std::string(optarg)); }
                else if(arg_name == "w_lig_spol")  { w_lig_spol = std::stod(std::string(optarg)); }
                else if(arg_name == "w_sasa")      { w_sasa = std::stod(std::string(optarg)); }
                else if(arg_name == "w_hbond")     { w_hbond = std::stod(std::string(optarg)); }
                else if(arg_name == "w_lig_lj")    { w_lig_lj = std::stod(std::string(optarg)); }
                else if(arg_name == "w_lig_ele")   { w_lig_ele = std::stod(std::string(optarg)); }
                else if(arg_name == "w_rot")       { w_rot = std::stod(std::string(optarg)); }
                else if(arg_name == "w_stack")     { w_stack = std::stod(std::string(optarg)); }
                else if(arg_name == "debye_length"){ debye_length = std::stod(std::string(optarg)); }
                break;
            }
            case ':':
                if (optopt == 'r' || optopt == 'l' || optopt == 'p' || optopt == 'i' || optopt == 's' || optopt == 't') {
                    std::cout << "Option -" << char(optopt) << " requires an argument." << std::endl;
                }
                args_correct = false;
                break;
            case '?':
            default:
                if (std::isprint(optopt)) {
                    std::cout << "Invalid option `-" << char(optopt) << "'." << std::endl;
                } else {
                    std::cout << "Invalid option character `\\x" << char(optopt) << "'." << std::endl;
                }
                args_correct = false;
                break;
        }
    }


    // no options are passed, or args are not correct
    if(optind == 1 || !args_correct) {
        std::cout << "Usage: " << usage_hint << std::endl;
        return 2;
    }

    if(in_score_mode != "virtual_screen" && in_score_mode != "binding_mode") {
        std::cout << "-s must be specified, it can be set to either virtual_screen or binding_mode." << std::endl;
        return 2;
    }

    std::cout << "#################################################" << std::endl;
    std::cout << "##               rldock score v2               ##" << std::endl;
    std::cout << "#################################################" << std::endl;

    if(help) {
        std::cout << "Input:" << std::endl;
        std::cout << "  --r arg            rigid receptor (mol2)" << std::endl;
        std::cout << "  --l arg            reference ligand with bond table and partial charges (mol2)" << std::endl;
        std::cout << "  --p arg            ligand poses, only use heavy atoms' coordinates, must has same atom order as the reference ligand (mol2)" << std::endl;
        std::cout << "  --i arg            atom property info file" << std::endl;
        std::cout << "  --s arg            score mode: virtual_screen or binding_mode" << std::endl;
        std::cout << std::endl;

        std::cout << "Input (optional):" << std::endl;
        std::cout << "  --t arg            number of torsions (default=0)" << std::endl;
        std::cout << "  --w_lj arg         intermolecular LJ (default=)" << std::endl;
        std::cout << "  --w_ele arg        intermolecule Estat (default=)" << std::endl;
        std::cout << "  --w_pol arg        polarization multual (default=)" << std::endl;
        std::cout << "  --w_lig_spol arg   polarization ligand self (default=)" << std::endl;
        std::cout << "  --w_rec_spol arg   polarization receptor self (default=)" << std::endl;
        std::cout << "  --w_sasa arg       delta SASA (default=)" << std::endl;
        std::cout << "  --w_hbond arg      hydrogen bond (default=)" << std::endl;
        std::cout << "  --w_lig_lj arg     delta ligand LJ (default=)" << std::endl;
        std::cout << "  --w_lig_ele arg    delta ligand Estat (default=)" << std::endl;
        std::cout << "  --w_rot arg        torsions (default=)" << std::endl;
        std::cout << "  --w_stack arg      stacking (default=)" << std::endl;
        std::cout << "  --debye_length arg Debye screening length in Å (default=)" << std::endl;

        std::cout << std::endl;
        std::cout << "Example of the folder data structure:" << std::endl;
        std::cout << "  --|root/" << std::endl;
        std::cout << "  ----|rec.mol2 (can be symlink or in different folder)" << std::endl;
        std::cout << "  ----|ref_lig.mol2 (can be symlink or in different folder)" << std::endl;
        std::cout << "  ----|poses.mol2 (can be symlink or in different folder)" << std::endl;
        std::cout << "  ----|info.txt (can be symlink or in different folder)" << std::endl;
        std::cout << "  ----|..." << std::endl;
        return 2;
    }

    const std::string rec_file_name = in_rec_path.substr(in_rec_path.rfind("/")+1, in_rec_path.size()-in_rec_path.rfind("/")-1);
    const std::string ref_lig_file_name = in_ref_lig_path.substr(in_ref_lig_path.rfind("/")+1, in_ref_lig_path.size()-in_ref_lig_path.rfind("/")-1);
    const std::string poses_file_name = in_poses_path.substr(in_poses_path.rfind("/")+1, in_poses_path.size()-in_poses_path.rfind("/")-1);
    const std::string info_file_name = in_info_path.substr(in_info_path.rfind("/")+1, in_info_path.size()-in_info_path.rfind("/")-1);


    std::cout << "in_rec_path:      " << in_rec_path     << std::endl;
    std::cout << "in_ref_lig_path:  " << in_ref_lig_path << std::endl;
    std::cout << "in_poses_path:    " << in_poses_path   << std::endl;
    std::cout << "in_info_path:     " << in_info_path    << std::endl;
    std::cout << "in_score_mode:    " << in_score_mode   << std::endl;
    std::cout << "debye length:     " << debye_length    << std::endl;

    // std::cout << "output_path:    " << output_path << std::endl;
    std::cout << "rec_file_name:      " << rec_file_name     << std::endl;
    std::cout << "ref_lig_file_name:  " << ref_lig_file_name << std::endl;
    std::cout << "poses_file_name:    " << poses_file_name   << std::endl;
    std::cout << "info_file_name:     " << info_file_name    << std::endl;
    std::cout << "verbose:  "           << verbose           << std::endl;

    for (int index = optind; index < argc; index++) {
        std::cout << "Not used argument: " << argv[index] << std::endl;
    }

    if(verbose) {
        std::cout << "weights: " << std::endl;
        std::streamsize ss = std::cout.precision();
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "w_lj       = " << w_lj       << "  w_ele      = " << w_ele      << "  w_pol      = " << w_pol   << std::endl;
        std::cout << "w_rec_spol = " << w_rec_spol << "  w_lig_spol = " << w_lig_spol << "  w_sasa     = " << w_sasa  << std::endl;
        std::cout << "w_hbond    = " << w_hbond    << "  w_lig_lj   = " << w_lig_lj   << "  w_lig_ele  = " << w_ele   << std::endl;
        std::cout << "w_rot      = " << w_rot      << "  w_stack    = " << w_stack    << std::endl;
        std::cout << std::setprecision(6);
        std::cout.precision(ss);
    }

    std::cout << "#################################################" << std::endl;

    if(!path_exists(in_rec_path))     { std::cout << "in_rec_path: "       << in_rec_path     << " not exists!" << std::endl; exit(2); }
    if(!path_exists(in_ref_lig_path)) { std::cout << "in_ref_lig_path: "   << in_ref_lig_path << " not exists!" << std::endl; exit(2); }
    if(!path_exists(in_poses_path))   { std::cout << "in_pose_list_path: " << in_poses_path   << " not exists!" << std::endl; exit(2); }
    if(!path_exists(in_info_path))    { std::cout << "in_info_path: "      << in_info_path    << " not exists!" << std::endl; exit(2); }

    ///////////////////////////////////////////////////
    // parse info file
    bool lig_torsion_info_flag = false;
    bool lig_aromatic_info_flag = false;
    bool rec_aromatic_info_flag = false;
    bool lig_h_donor_info_flag = false;
    bool lig_h_acceptor_info_flag = false;
    bool rec_h_donor_info_flag = false;
    bool rec_h_acceptor_info_flag = false;
    std::ifstream in_info_f(in_info_path, std::ios::in);
    std::string info_sline;
    std::vector<rlscore::Size_Type> rec_aromatic_indices;
    std::vector<rlscore::Size_Type> lig_aromatic_indices;
    std::vector<rlscore::Size_Type> rec_h_donor_indices;
    std::vector<rlscore::Size_Type> rec_h_acceptor_indices;
    std::vector<rlscore::Size_Type> lig_h_donor_indices;
    std::vector<rlscore::Size_Type> lig_h_acceptor_indices;
    while(std::getline(in_info_f, info_sline)) {
        if(info_sline.find("lig_openeye_torsion:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            assert(token.size() == 2);
            rot_dof = std::stod(token[1]);
            lig_torsion_info_flag = true;
        }
        if(info_sline.find("lig_aromatic_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                lig_aromatic_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == lig_aromatic_indices.size());
            lig_aromatic_info_flag = true;
        }
        if(info_sline.find("nuc_aromatic_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                rec_aromatic_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == rec_aromatic_indices.size());
            rec_aromatic_info_flag = true;
        }
        if(info_sline.find("lig_H_donor_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                lig_h_donor_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == lig_h_donor_indices.size());
            lig_h_donor_info_flag = true;
        }
        if(info_sline.find("lig_H_acceptor_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                lig_h_acceptor_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == lig_h_acceptor_indices.size());
            lig_h_acceptor_info_flag = true;
        }
        if(info_sline.find("nuc_H_donor_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                rec_h_donor_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == rec_h_donor_indices.size());
            rec_h_donor_info_flag = true;
        }
        if(info_sline.find("nuc_H_acceptor_idx:") != std::string::npos) {
            const std::vector<std::string> token = rlscore::string2vector(info_sline);
            for(rlscore::Size_Type i = 1; i != token.size(); i++) {
                rec_h_acceptor_indices.push_back(std::stoi(token[i]));
            }
            assert((token.size()-1) == rec_h_acceptor_indices.size());
            rec_h_acceptor_info_flag = true;
        }
    }
    in_info_f.close();

    if(lig_torsion_info_flag == false || lig_aromatic_info_flag == false || rec_aromatic_info_flag == false || 
       lig_h_donor_info_flag == false || lig_h_acceptor_info_flag == false || rec_h_donor_info_flag == false || rec_h_acceptor_info_flag == false) {
        std::cout << "Incomplete info from "     << in_info_path << std::endl;
        std::cout << "lig_torsion_info_flag "    << lig_torsion_info_flag << std::endl;
        std::cout << "lig_aromatic_info_flag "   << lig_aromatic_info_flag << std::endl;
        std::cout << "rec_aromatic_info_flag "   << rec_aromatic_info_flag << std::endl;
        std::cout << "lig_h_donor_info_flag "    << lig_h_donor_info_flag << std::endl;
        std::cout << "lig_h_acceptor_info_flag " << lig_h_acceptor_info_flag << std::endl;
        std::cout << "rec_h_donor_info_flag "    << rec_h_donor_info_flag << std::endl;
        std::cout << "rec_h_acceptor_info_flag " << rec_h_acceptor_info_flag << std::endl;
        exit(2);
    }

    ///////////////////////////////////////////////////
    const rlscore::Parameter param(debye_length);

    rlscore::Receptor rec(param, in_rec_path, rec_aromatic_indices, rec_h_donor_indices, rec_h_acceptor_indices);
    std::cout << "Load receptor: "      << rec.name << " from " << in_rec_path << std::endl;
    std::cout << "--> heavy atom num: " << rec.atoms.size() << std::endl;


    rlscore::Ligand ref_lig(param, in_ref_lig_path, rot_dof, lig_aromatic_indices, lig_h_donor_indices, lig_h_acceptor_indices);
    std::cout << "Load reference ligand: " << ref_lig.name << " from " << in_ref_lig_path << std::endl;
    std::cout << "--> heavy atom num: "    << ref_lig.atoms.size() << std::endl;

    rlscore::Complex complex;
    complex.initialize_receptor_and_ligand(param, rec, ref_lig);

    //get poses and score
    std::cout << "#Scoring -----------------------------------------" << std::endl;
    std::cout << "#pose total dU(lj) dU(ele) dU(mutual_pol) dU(lig_self_pol) dU(sasa) dU(hbond) dU(rec_self_pol) dU(lig_lj) dU(lig_ele) rot_dof dU(stack)" << std::endl;
    std::vector<std::string> output_strs;
    std::ifstream in_poses_f(in_poses_path, std::ios::in);
    while(true) {
        std::string pose_name;
        std::vector<rlscore::Atom> pose_atoms;
        rlscore::parse_single_pose_from_stream_each_time(in_poses_f, pose_name, pose_atoms);

        if(pose_name == "" && pose_atoms.size() == 0 && in_poses_f.eof()) {
                break;
        }

        assert(ref_lig.atoms.size() == pose_atoms.size());
        for(rlscore::Size_Type i = 0; i < ref_lig.atoms.size(); ++i) {
            ref_lig.atoms[i].xyz = pose_atoms[i].xyz;
        }
        complex.update_ligand(param, rec, ref_lig);

        const rlscore::Float& deltaU_lj = complex.terms.deltaU_inter_lj;
        const rlscore::Float& deltaU_lig_lj = complex.terms.deltaU_lig_lj;
        const rlscore::Float& deltaU_lig_ele = complex.terms.deltaU_lig_ele;
        const rlscore::Float& deltaU_ele = complex.terms.deltaU_inter_ele;
        const rlscore::Float& deltaU_hbond = complex.terms.deltaU_inter_hbond;
        const rlscore::Float& deltaU_sasa = complex.terms.deltaU_sasa;
        const rlscore::Float& deltaU_mutual_pol = complex.terms.deltaU_mutual_pol;
        const rlscore::Float& deltaU_rec_self_pol = complex.terms.deltaU_rec_self_pol;
        const rlscore::Float& deltaU_lig_self_pol = complex.terms.deltaU_lig_self_pol;
        const rlscore::Float& deltaU_stack = complex.terms.deltaU_inter_stack;


        const rlscore::Float& U_ref_lig_lj = complex.terms.U_lig_lj_before;
        const rlscore::Float& U_bound_lig_lj = complex.terms.U_lig_lj_after;
        const rlscore::Float& U_ref_lig_ele = complex.terms.U_lig_ele_before;
        const rlscore::Float& U_bound_lig_ele = complex.terms.U_lig_ele_after;
        const rlscore::Float& delta_rec_sasa = complex.terms.delta_rec_sasa;
        const rlscore::Float& ref_lig_sasa = complex.terms.lig_sasa_before;
        const rlscore::Float& bound_lig_sasa = complex.terms.lig_sasa_after;

        // ////////trained for 1A success rate////////
        // const rlscore::Float& total_energy = 0.81558241 * deltaU_lj + 0.04005243 * deltaU_ele + 0.14544507 * deltaU_mutual_pol + 2.02957049 * deltaU_lig_self_pol + 4.18490388 * deltaU_sasa + 1.13174837 * deltaU_hbond + 4.5524668 * deltaU_rec_self_pol + 2.50539453 * deltaU_lig_lj + 4.93970779 * deltaU_lig_ele;
        // const rlscore::Float deltaU_total_estat = deltaU_mutual_pol + deltaU_rec_self_pol + deltaU_lig_self_pol + deltaU_ele + deltaU_lig_ele;
        // const rlscore::Float deltaU_total_estat_weighted = 0.14544507 * deltaU_mutual_pol + 4.5524668 * deltaU_rec_self_pol + 2.02957049 * deltaU_lig_self_pol + 0.04005243 * deltaU_ele + 4.93970779 * deltaU_lig_ele;
        // ////////

        // [0.96885464, 0.14987373, 0.26545015, 3.13873365, 3.22132689, 1.29912147, 4.70148162, 0.64799766, 0.81178695]
        ////////trained for 2A success rate////////
        // const rlscore::Float& total_energy = 0.96885464 * deltaU_lj + 0.14987373 * deltaU_ele + 0.26545015 * deltaU_mutual_pol + 3.13873365 * deltaU_lig_self_pol + 3.22132689 * deltaU_sasa + 1.29912147 * deltaU_hbond + 4.70148162 * deltaU_rec_self_pol + 0.64799766 * deltaU_lig_lj + 0.81178695 * deltaU_lig_ele;
        // const rlscore::Float deltaU_total_estat = deltaU_mutual_pol + deltaU_rec_self_pol + deltaU_lig_self_pol + deltaU_ele + deltaU_lig_ele;
        // const rlscore::Float deltaU_total_estat_weighted = 0.26545015 * deltaU_mutual_pol + 4.70148162 * deltaU_rec_self_pol + 3.13873365 * deltaU_lig_self_pol + 0.14987373 * deltaU_ele + 0.81178695 * deltaU_lig_ele;
        ////////

        //////trained for 2A success rate and affinity////////
        const rlscore::Float total_energy = global_score_scale*(w_lj*deltaU_lj + w_ele*deltaU_ele + w_pol*deltaU_mutual_pol + w_lig_spol*deltaU_lig_self_pol + w_sasa*deltaU_sasa + w_hbond*deltaU_hbond + w_rec_spol*deltaU_rec_self_pol + w_lig_lj*deltaU_lig_lj + w_lig_ele*deltaU_lig_ele + w_stack*deltaU_stack) / (w_rot*rot_dof + 1.0);

        std::ostringstream oss;
        oss << pose_name << " ";
        oss << std::fixed << std::setprecision(3) << total_energy << " " << deltaU_lj << " " << deltaU_ele << " " << deltaU_mutual_pol << " " << deltaU_lig_self_pol << " " << deltaU_sasa << " " << deltaU_hbond << " " << deltaU_rec_self_pol << " " << deltaU_lig_lj << " " << deltaU_lig_ele << " " << rot_dof << " " << deltaU_stack;
        output_strs.push_back(oss.str());

        if(verbose) {
            std::cout << "######info: "           << pose_name << std::endl;
            std::cout  << "inter_lj:            " << deltaU_lj                       << std::endl;
            std::cout  << "inter_ele:           " << deltaU_ele                      << std::endl;
            std::cout  << "lig_ele:             " << U_bound_lig_ele                 << std::endl;
            std::cout  << "lig_lj:              " << U_bound_lig_lj                  << std::endl;
            std::cout  << "lig_sasa:            " << bound_lig_sasa                  << std::endl;
            std::cout  << "delta_lig_ele:       " << U_bound_lig_ele - U_ref_lig_ele << std::endl;
            std::cout  << "delta_lig_lj:        " << U_bound_lig_lj - U_ref_lig_lj   << std::endl;
            std::cout  << "delta_lig_sasa:      " << bound_lig_sasa - ref_lig_sasa   << std::endl;
            std::cout  << "delta_rec_sasa:      " << delta_rec_sasa                  << std::endl;
            std::cout  << "deltaU_sasa:         " << deltaU_sasa                     << std::endl;
            std::cout  << "rot_dof:             " << rot_dof                         << std::endl;
            std::cout  << "deltaU_stack:        " << deltaU_stack                    << std::endl;
            std::cout  << "deltaU_hbond:        " << deltaU_hbond                    << std::endl;
            std::cout  << "deltaU_mutual_pol:   " << deltaU_mutual_pol               << std::endl;
            std::cout  << "deltaU_lig_self_pol: " << deltaU_lig_self_pol             << std::endl;
            std::cout  << "deltaU_rec_self_pol: " << deltaU_rec_self_pol             << std::endl;
            std::cout << "###########################################" << std::endl;
        }
    }
    in_poses_f.close();

    for(const auto& out_str : output_strs) {
        std::cout << out_str << std::endl;
    }
    return 0;

#ifdef MODE_PROFILE
ProfilerStop();
#endif
}//end main