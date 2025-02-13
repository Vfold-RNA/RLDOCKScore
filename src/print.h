#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include "vec3d.h"
// #include "mat3x3.h"
// #include "quaternion.h"
// #include "atom.h"
// #include "molecule.h"

namespace rlscore {

inline void print(Float x, std::ostream& out = std::cout) {
	out << x;
}
inline void print(int x, std::ostream& out = std::cout) {
	out << x;
}
inline void print(Size_Type x, std::ostream& out = std::cout) {
	out << x;
}
template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
	out << "[";
	for(Size_Type i = 0; i< v.size(); ++i) {
		if(i != 0)
			out << " ";
		out << v[i];
	}
	out << "]";
}

inline void print(const Vec3d& v, std::ostream& out = std::cout) {
	out << "(";
	for(Size_Type i = 0; i< v.size(); ++i) {
		if(i != 0)
			out << ", ";
		out << v[i];
	}
	out << ")";
}

// inline void print(const Mat3x3& m, std::ostream& out = std::cout) {
// 	out << "{";
// 	out << m(0,0) << "," << m(0,1) << "," << m(0,2) << std::endl;
// 	out << m(1,0) << "," << m(1,1) << "," << m(1,2) << std::endl;
// 	out << m(2,0) << "," << m(2,1) << "," << m(2,2) << "}";
// }

// inline void print(const Quaternion& q, std::ostream& out = std::cout) { // print as an angle
//     print(q.to_vec3d().norm(), out);
// 	print(q.to_vec3d(), out);
// }

// inline void print(const Atom& a, std::ostream& out = std::cout) {
//     // std::vector<Bond> bonds;
// 	out << std::setfill(' ') << std::setw(6) << std::left << "atom" << std::right << std::setw(5) << a.get_serial() << " "
// 		    << std::setw(6) << std::left << " "+a.get_name() << std::right
// 		    << std::setw(5) << std::right << a.get_res_name() << std::right << " " << std::setw(1) << a.get_chain_name()
// 		    << std::setw(4) << a.get_res_serial() << std::setw(3) << " "
// 		    << std::setw(8) << std::fixed << std::setprecision(3) << a.get_coord()[0] << std::setw(8) << a.get_coord()[1]
// 		    << std::setw(8) << a.get_coord()[2] << std::setw(6) << a.get_sybyl_type_name() << std::setw(6) << a.get_element_type_name()
// 		    << std::setw(10) << a.get_charge();
// 	out << std::defaultfloat << std::setprecision(6) << std::endl;
// }

// inline void print(const Conf& c, std::ostream& out = std::cout) {
// 	out << "{";
// 	print(c.origin, out); out << ",";
// 	print(c.rot_axis, out); out << ",";
// 	print(c.rot_angle, out);
// 	out << "}";
// }

// inline void print(const Node_Id& ni, std::ostream& out = std::cout) {
// 	out << "NodeId[" << ni.depth_index << "," << ni.width_index << "]";
// }

// inline void print(const Node& n, std::ostream& out = std::cout) {
// 	out << "node--"; print(n.get_node_id(), out); out << "--> conf: "; print(n.conf, out);
// 	out << " range: "; print(n.atom_ranges, out);
// 	if(n.parent==nullptr) {
// 		out << " root ";
// 	} else {
// 		out << " parentID: ";
// 		print((*(n.parent)).get_node_id(), out);
// 	}
// 	out << " children_num: " << n.children.size();
// 	out << " bonded atom index[" << n.parent_bond_atom_index << "," << n.self_bond_atom_index << "]" << std::endl;
// }

// // inline void print(const Layer& l, std::ostream& out = std::cout) {
// // 	out << "------------------Layer------------------" << std::endl;
// // 	for(const auto& seg : l.segments){
// // 		print(seg, out); out << std::endl;
// // 	}
// // }

// // inline void print(const Tree& t, std::ostream& out = std::cout) {
// // 	out << "------------------Tree------------------" << std::endl;
// // 	for(const auto& l : t.layers) {
// // 		print(l, out); out << std::endl;
// // 	}
// // }

// inline void print(const Nodes& ns, std::ostream& out = std::cout) {
// 	for(const auto& n : ns) {
// 		print(n, out);
// 	}
// }

// inline void print(const Ligand& l, std::ostream& out = std::cout) {
// 	out << "#####################Ligand#####################" << std::endl;
// 	out << "Structure--->" << std::endl;
// 	print(l.nodes, out);
// 	out << "Atom--->" << std::endl;
// 	for(const auto& a : l.atoms) {
// 		print(a, out);
// 	}
// }

// inline void print(const RNA& r, std::ostream& out = std::cout) {
// 	out << "#####################RNA#####################" << std::endl;
// 	out << "Atom--->" << std::endl;
// 	for(const auto& a : r.atoms) {
// 		print(a, out);
// 	}
// }


// void print(const vec3d& v, std::ostream& out = std::cout);//from vec3d, printnl will find it when linking
// void print(const qt& q, std::ostream& out);//from quaternion, printnl will find it when linking
template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
	print(x, out);
	out << '\n';
}

}