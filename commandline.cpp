// Lipid bending modulus calculation
// Tom Joseph <thomas.joseph@@pennmedicine.upenn.edu>
// Reimplemented from code from Itay Barel, Frank Brown, and others at UC Santa Barbara
// Algorithm from: Levine et al, JACS 2014, 136, 13582-13585. dx.doi.org/10.1021/ja507910r
// References to lines in comments below refer to old_vmd_bending_modulus.cpp (which is a modified
// version of the original source code)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include "floatarray2d.h"
#include "bending_modulus.h"

using namespace std;

// Slurps up all lines from a text file
// TODO: Return false on error
bool load_all_lines(const string &filename, std::vector<string> &out_lines) {
    std::ifstream f(filename);
    std::string line;
    out_lines.clear();
    while(std::getline(f, line)) {
        out_lines.push_back(line);
    }
    return true;
}

// Converts three vectors of strings of x, y, z coordinates respectively into one
// 2D vector of floats.
// TODO: Return false on error
bool coalesce_xyz_lines(const std::vector<string> &x_lines, const std::vector<string> &y_lines,
                        const std::vector<string> &z_lines, FloatArray2D &out_array) {
    // The lengths of the input vectors should be the same
    if(x_lines.size() != y_lines.size() || y_lines.size() != z_lines.size()) {
        cerr << "coalesce_xyz_lines: Input vector length mismatch. How's that for cryptic?" << endl;
        cerr << "coalesce_xyz_lines: Remember, I don't tolerate empty lines in input files." << endl;
        exit(EXIT_FAILURE);
    }
    // Ensure we have enough space
    int num_tuples = x_lines.size();
    out_array.reset(num_tuples, 3);

    // Do the actual parsing
    for(int i = 0; i < num_tuples; i++) {
        out_array.set(i, 0, std::stof(x_lines[i]));
        out_array.set(i, 1, std::stof(y_lines[i]));
        out_array.set(i, 2, std::stof(z_lines[i]));
    }
    return true;
}


int main(int argc, char *argv[]) {
    const int num_lipids = 512; // known from VMD
    const int grid_size = 12; // to be a user argument
    const bool should_norm_z1 = false; // to be a user argument
    const float T0IN = 17.97264862; // initial guess for average thickness?

    cout << "Loading test set..." << endl;
    std::vector<string> box_x_lines, box_y_lines, box_z_lines;
    load_all_lines("from_barel/TestSet2/boxsizeX.out", box_x_lines);
    load_all_lines("from_barel/TestSet2/boxsizeY.out", box_y_lines);
    load_all_lines("from_barel/TestSet2/boxsizeZ.out", box_z_lines);
    cout << "Lines in box files: " << box_x_lines.size() << " " << box_y_lines.size() << " " << box_z_lines.size()
         << endl;
    FloatArray2D box_size;
    coalesce_xyz_lines(box_x_lines, box_y_lines, box_z_lines, box_size);

    // Load the lipids
    std::vector<string> lipid_x_lines, lipid_y_lines, lipid_z_lines;
    load_all_lines("from_barel/TestSet2/LipidX.out", lipid_x_lines);
    load_all_lines("from_barel/TestSet2/LipidY.out", lipid_y_lines);
    load_all_lines("from_barel/TestSet2/LipidZ.out", lipid_z_lines);
    cout << "Lines in lipid files: " << lipid_x_lines.size() << " " << lipid_y_lines.size() << " "
         << lipid_z_lines.size()
         << endl;
    // This lipid variable contains both heads and tails!
    FloatArray2D lipid;
    coalesce_xyz_lines(lipid_x_lines, lipid_y_lines, lipid_z_lines, lipid);

    do_bending_modulus(lipid, box_size, num_lipids, grid_size, should_norm_z1);
}
