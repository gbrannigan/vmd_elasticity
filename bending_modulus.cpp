// Lipid bending modulus calculation
// Tom Joseph <thomas.joseph@@pennmedicine.upenn.edu>
// Reimplemented from code from Itay Barel, Frank Brown, and others at UC Santa Barbara
// Algorithm from: Levine et al, JACS 2014, 136, 13582-13585. dx.doi.org/10.1021/ja507910r
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <fftw3.h>
#include "floatarray2d.h"

using namespace std;

constexpr float PHI0IN = 0.01588405482;

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

// Separate raw lipid coordinates into head and tail, and clean it up
// TODO: Return false on error
bool separate_lipid_heads_tails(const FloatArray2D &lipid, FloatArray2D &out_heads, FloatArray2D &out_tails) {
    if(lipid.num_rows() % 2 != 0) {
        cerr << "organize_lipid: Non-even number of lipid coordinates, which implies that number of heads" << endl;
        cerr << "    does not equal number of tails. Something is corrupt." << endl;
        exit(EXIT_FAILURE);
    }
    int num_lipids = lipid.num_rows() / 2;
    out_heads.reset(num_lipids, 3);
    out_tails.reset(num_lipids, 3);
    for(int lipid_i = 0; lipid_i < num_lipids; lipid_i++) {
        out_heads.copy_row(lipid, lipid_i*2, lipid_i);
        out_tails.copy_row(lipid, lipid_i*2+1, lipid_i);
    }
    return true;
}

// Fix lipids that aren't inside the box - but only in the XY plane
void wrap_into_box_xy_plane(FloatArray2D &heads, FloatArray2D &tails, float a, float b, float avg_a, float avg_b) {
    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        float hx = heads.get(lipid_i, 0), hy = heads.get(lipid_i, 1);
        float tx = tails.get(lipid_i, 0), ty = tails.get(lipid_i, 1);
        // Wrap only in the XY plane
        if(hx >= a) {
            hx -= a;
            tx -= a;
        }
        if(hx < 0) {
            hx += a;
            tx += a;
        }
        if(hy > b) {
            hy -= b;
            ty -= b;
        }
        if(hy < 0) {
            hy += b;
            ty += b;
        }

        // Fix the tail beads carried to the other side of the box
        if(fabs(hx - tx) > (avg_a * 0.5)) {
            if(hx > tx) {
                tx += a;
            } else {
                tx -= a;
            }
        }
        if(fabs(hy - ty) > (avg_b * 0.5)) {
            if(hy > ty) {
                ty += b;
            } else {
                ty -= b;
            }
        }

        heads.set(lipid_i, 0, hx);
        heads.set(lipid_i, 1, hy);
        tails.set(lipid_i, 0, tx);
        tails.set(lipid_i, 1, ty);
    }
}

void wrap_z_and_calc_director(FloatArray2D &heads, FloatArray2D &tails, float c,
        FloatArray2D &out_director, std::vector<bool> &out_lipid_good,
        int &out_num_upper_leaflet, int &out_num_lower_leaflet, float &out_head_avg_z) {
    const float CUTOFF_ANGLE = 0.0f;
    const float z_wrap_threshold = 0.6f * c;

    std::vector<bool> lipid_good(heads.num_rows());
    lipid_good.clear();
    float upper_leaflet_avg_z = 0, lower_leaflet_avg_z = 0;
    out_num_upper_leaflet = 0;
    out_num_lower_leaflet = 0;
    out_director.reset(heads.num_rows(), heads.num_cols());

    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        float dx = tails.get(lipid_i, 0) - heads.get(lipid_i, 0);
        float dy = tails.get(lipid_i, 1) - heads.get(lipid_i, 1);
        float dz = tails.get(lipid_i, 2) - heads.get(lipid_i, 2);
        float magnitude = sqrtf(dx * dx + dy * dy + dz * dz);

        // Keep track of which lipids are "good", meaning not at some crazy angle
        out_lipid_good.push_back(fabs(dz / magnitude) > CUTOFF_ANGLE);
        // cout << dz << endl;
        if (dz < 0) {
            upper_leaflet_avg_z += heads.get(lipid_i, 2);
            out_num_upper_leaflet++;
        }
        if (dz > 0) {
            lower_leaflet_avg_z += heads.get(lipid_i, 2);
            out_num_lower_leaflet++;
        }
        out_director.set(lipid_i, 0, dx);
        out_director.set(lipid_i, 1, dy);
        out_director.set(lipid_i, 2, dz);
    }
    upper_leaflet_avg_z /= float(out_num_upper_leaflet);
    lower_leaflet_avg_z /= float(out_num_lower_leaflet);

    // Swap lipids between leaflets if necessary, according to the way they point
    int num_swapped_from_upper = 0, num_swapped_from_lower = 0;
    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        float dz = out_director.get(lipid_i, 2);
        float hz = heads.get(lipid_i, 2);
        float tz = tails.get(lipid_i, 2);

        if(dz < 0 && fabs(hz - (upper_leaflet_avg_z/float(out_num_upper_leaflet))) > z_wrap_threshold) {
            hz += c;
            tz += c;
            num_swapped_from_upper++;
        }
        if(dz > 0 && fabs(hz - (lower_leaflet_avg_z/float(out_num_lower_leaflet))) > z_wrap_threshold) {
            hz -= c;
            tz -= c;
            num_swapped_from_lower++;
        }
        heads.set(lipid_i, 2, hz);
        tails.set(lipid_i, 2, tz);
    }

    // Now that we've done all our adjustments, we can calculate the average Z-value of head
    out_head_avg_z = 0;
    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        out_head_avg_z += heads.get(lipid_i, 2);
    }
    out_head_avg_z /= float(heads.num_rows());
}

// Set up the Q matrix (wavevector), which will end up going into the FFT.
// The "x" and "y" components are stored separately, to avoid making a "Float2Array2D" or some such
void create_q_matrix(int grid_size, FloatArray2D &qx_matrix, FloatArray2D &qy_matrix,
        FloatArray2D &cos_q_matrix, FloatArray2D &sin_q_matrix) {
    qx_matrix.reset(grid_size, grid_size);
    qy_matrix.reset(grid_size, grid_size);
    sin_q_matrix.reset(grid_size, grid_size);
    cos_q_matrix.reset(grid_size, grid_size);

    // Center on zero, or something?
    for(int i = 0; i < grid_size; i++) {
        for(int j = 0; j < grid_size; j++) {
            auto q0 = float(i), q1 = float(j);
            if(i >= grid_size/2) {
                q0 -= float(grid_size);
            }
            if(j >= grid_size/2) {
                q1 -= float(grid_size);
            }

            if(i != 0 || j != 0) {
                float magnitude = 1.0/sqrt(q0*q0 + q1*q1);
                cos_q_matrix.set(i, j, q0*magnitude);
                sin_q_matrix.set(i, j, q1*magnitude);
            }
            qx_matrix.set(i, j, q0);
            qy_matrix.set(i, j, q1);
        }
    }
}

// This mod function does a second round of modulo to handle mildly negative numbers.
// That is, this function requires that |x| < n
int smarter_mod(int x, int n) {
    return ((x % n) + n) % n;
}

// Returns 1 if the interpolation considered any empty adjacent cells
int interpolate_one_cell(FloatArray2D &data, FloatArray2D &good_count, int i, int j) {
    // Identify the adjacent grid squares, wrapping due to periodic boundary conditions
    // We arbitrarily call the i-direction up/down and j-direction left/right, for ease of intuition
    int grid_size = data.num_rows();
    int i_up = smarter_mod(i - 1, grid_size);
    int i_down = smarter_mod(i + 1, grid_size);
    int j_left = smarter_mod(j - 1, grid_size);
    int j_right = smarter_mod(j + 1, grid_size);
    float up_count = good_count.get(i_up, j);
    float down_count = good_count.get(i_down, j);
    float left_count = good_count.get(i, j_left);
    float right_count = good_count.get(i, j_right);

    float good_adjacent = up_count + down_count + left_count + right_count;
    // FIXME: There may not be any adjacent lipids. In this case, set the z_height to zero?
    float data_up = data.get(i_up, j);
    float data_down = data.get(i_down, j);
    float data_left = data.get(i, j_left);
    float data_right = data.get(i, j_right);
    data.set(i, j, (up_count*data_up + down_count*data_down + left_count*data_left
                        + right_count*data_right) / good_adjacent);
    if(up_count == 0 || down_count == 0 || left_count == 0 || right_count == 0) {
        return 1;
    }
    return 0;
}

// Map lipids onto discretized grid which will eventually be FFTed
// Returns total adjacent empty count
int map_onto_grid(FloatArray2D &heads, FloatArray2D &tails, FloatArray2D &director,
        std::vector<bool> &lipid_good, int grid_size,
        float box_a, float box_b, float avg_head_z,
        FloatArray2D &z_height_upper, FloatArray2D &z_height_lower,
        FloatArray2D &lipid_good_upper, FloatArray2D &lipid_good_lower,
        FloatArray2D &lipid_bad_upper, FloatArray2D &lipid_bad_lower,
        std::vector<int> &xj_vec, std::vector<int> &yj_vec) {
    z_height_upper.reset(grid_size, grid_size);
    z_height_lower.reset(grid_size, grid_size);
    // Keep count of how many good and bad lipids are represented at each grid point
    // It's nice they are floats because we will divide by them later.
    lipid_good_upper.reset(grid_size, grid_size);
    lipid_good_lower.reset(grid_size, grid_size);
    lipid_bad_upper.reset(grid_size, grid_size);
    lipid_bad_lower.reset(grid_size, grid_size);
    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        float hx = heads.get(lipid_i, 0), hy = heads.get(lipid_i, 1);
        // Discretize on the XY grid
        int this_xj = int(floor(hx/(box_a/float(grid_size))));
        int this_yj = int(floor(hy/(box_b/float(grid_size))));

        if(this_xj >= grid_size) {
            // Check for lipid coordinate on the box edge exactly
            if(hx == box_a) {
                this_xj = grid_size - 1;
            } else {
                // TODO: Print something about this for user to read
                // This must be a bad thing, because we are on the edge of the grid
            }
        }
        if(this_yj >= grid_size) {
            if(hy == box_b) {
                this_yj = grid_size - 1;
            } else {
                // TODO: print something in the else condition
            }
        }
        // TODO: Print some more stuff
        if(this_xj < 0) {
            cout << "this_xj=" << this_xj << " for hx=" << hx << endl;
        }
        if(this_yj < 0) {
            cout << "this_yj=" << this_yj << " for hy=" << hy << endl;
        }

        // Fill out the discretized Z-height map
        float hz = heads.get(lipid_i, 2);
        float dz = director.get(lipid_i, 2);
        // Upper leaflet
        if(dz < 0) {
            if(lipid_good[lipid_i]) {
                z_height_upper.add(this_xj, this_yj, hz - avg_head_z);
                // TODO: z1sq_av_frame business
                // Increment upper leaflet grid count for this cell
                lipid_good_upper.add(this_xj, this_yj, 1);
            } else {
                // Remember how many bad lipids there are at this grid point
                lipid_bad_upper.add(this_xj, this_yj, 1);
            }
        }
        // Lower leaflet
        if(dz > 0) {
            if(lipid_good[lipid_i]) {
                z_height_lower.add(this_xj, this_yj, hz - avg_head_z);
                lipid_good_lower.add(this_xj, this_yj, 1);
            } else {
                lipid_bad_lower.add(this_xj, this_yj, 1);
            }
        }

        // Save the individual grid points such that they can be indexed by lipid_i
        xj_vec.push_back(this_xj);
        yj_vec.push_back(this_yj);
    }
    // Normalize Z-height map by number of lipids
    z_height_upper.divide(lipid_good_upper);
    z_height_lower.divide(lipid_good_lower);

    // Fill empty grid cell by interpolation from adjacent cells, for both leaflets (e)
    int num_empty_total = 0;
    for(int i = 0; i < grid_size; i++) {
        for(int j = 0; j < grid_size; j++) {
            if(lipid_good_upper.get(i, j) == 0) {
                num_empty_total += interpolate_one_cell(z_height_upper, lipid_good_upper, i, j);
            }
            if(lipid_good_lower.get(i, j) == 0) {
                num_empty_total += interpolate_one_cell(z_height_lower, lipid_good_lower, i, j);
            }
        }
    }
    return num_empty_total;
}

// OK, OK, I know the number of parameters here is ridiculous
void calculate_number_densities(FloatArray2D &heads, FloatArray2D &tails, FloatArray2D &director, std::vector<bool> &lipid_good,
        int num_upper_leaflet, int num_lower_leaflet,
        FloatArray2D &qx_matrix, FloatArray2D &qy_matrix, int grid_size, float a, float b, float avg_head_z,
        FloatArray2D &h_real, FloatArray2D &h_imag,
        FloatArray2D &psi_upper_real, FloatArray2D &psi_upper_imag, FloatArray2D &psi_lower_real,
        FloatArray2D &psi_lower_imag) {
    h_real.reset(grid_size, grid_size);
    h_imag.reset(grid_size, grid_size);
    psi_upper_real.reset(grid_size, grid_size);
    psi_upper_imag.reset(grid_size, grid_size);
    psi_lower_real.reset(grid_size, grid_size);
    psi_lower_imag.reset(grid_size, grid_size);

    // FIXME: This stuff below is flagrantly vectorizable
    for(int lipid_i = 0; lipid_i < heads.num_rows(); lipid_i++) {
        float hx = heads.get(lipid_i, 0), hy = heads.get(lipid_i, 1), hz = heads.get(lipid_i, 2);
        // Do upper half of complex plane only, and use symmetries later
        for(int grid_i = 0; grid_i < grid_size/2+1; grid_i++) {
            for(int grid_j = 0; grid_j < grid_size; grid_j++) {
                float qx = 2.0f*float(M_PI)*qx_matrix.get(grid_i, grid_j) / a;
                float qy = 2.0f*float(M_PI)*qy_matrix.get(grid_i, grid_j) / b;

                // We assume AREA_tail is zero (from the original code)
                h_real.add(grid_i, grid_j, (hz - avg_head_z) * cos(qx*hx + qy*hy));
                h_imag.add(grid_i, grid_j, -((hz - avg_head_z) * sin(qx*hx + qy*hy)));

                float dz = director.get(lipid_i, 2);
                if(lipid_good[lipid_i] && (grid_i != 0 || grid_j != 0)) {
                    // Upper leaflet
                    if(dz < 0) {
                        psi_upper_real.add(grid_i, grid_j, cos(qx*hx + qy*hy));
                        psi_upper_imag.add(grid_i, grid_j, sin(qx*hx + qy*hy));
                    }
                    if(dz > 0) {
                        psi_lower_real.add(grid_i, grid_j, cos(qx*hx + qy*hy));
                        psi_lower_imag.add(grid_i, grid_j, sin(qx*hx + qy*hy));
                    }
                }
            }
        } // for loop over upper half of FFT grid

        // Normalize by XY box size
        float box_xy_magnitude = sqrtf(a*b);
        psi_upper_real.multiply(1.0f/box_xy_magnitude);
        psi_upper_imag.multiply(1.0f/box_xy_magnitude);
        psi_lower_real.multiply(1.0f/box_xy_magnitude);
        psi_lower_imag.multiply(1.0f/box_xy_magnitude);
        // Treat the 0,0 position specially
        psi_upper_real.set(0, 0, float(num_upper_leaflet)/box_xy_magnitude - PHI0IN*box_xy_magnitude);
        psi_upper_imag.set(0, 0, 0);
        psi_lower_real.set(0, 0, float(num_lower_leaflet)/box_xy_magnitude - PHI0IN*box_xy_magnitude);
        psi_lower_imag.set(0, 0, 0);
    } // for loop over lipids
}


// Do a DFT on a grid, allocating an output buffer if necessary
fftwf_complex *do_dft2d_r2c(int grid_size, fftwf_plan plan, float *raw_data, fftwf_complex *fft_data = nullptr) {
    int num_grid_pairs = grid_size*(grid_size/2+1);
    if(fft_data == nullptr) {
        fft_data = new fftwf_complex[num_grid_pairs];
    }
    std::memset(fft_data, 0, sizeof(fftwf_complex)*num_grid_pairs);
    fftwf_execute_dft_r2c(plan, raw_data, fft_data);
    return fft_data;
}


// Do an DFT on a FloatArray2D
// This is pretty much syntactic sugar.
// Caller must delete[] the returned fftwf_complex*.
fftwf_complex *do_dft2d_r2c(fftwf_plan plan, FloatArray2D &data, fftwf_complex *fft_data = nullptr) {
    assert(data.num_rows() == data.num_cols());
    int grid_size = data.num_rows();
    auto raw_data = new float[grid_size*grid_size];
    data.copy_raw_to(raw_data);
    fft_data = do_dft2d_r2c(grid_size, plan, raw_data, fft_data);

    delete[] raw_data;
    return fft_data;
}


// Converts fftw_complex array into 2D matrices, real and imaginary, normalized,
// with some weird stuff happening in first row and column.
// Fills in the right half of the full array, using the property
// h_{-q} = h*_{q}=h_{N-q}, whatever that means.
// This modifies all the array arguments!
void fill_full_array(FloatArray2D &real, FloatArray2D &imaginary, fftwf_complex *array2,
    float box_xy_magnitude) {

    // We can avoid passing in the grid size because real/imaginary should already be the correct size
    // cerr << "grid_size: " << grid_size << endl
    //      << real.num_rows() << " " << real.num_cols() << " " << imaginary.num_rows() << " "  << imaginary.num_cols() <<  endl;
    assert(real.num_rows() == real.num_cols());
    assert(imaginary.num_rows() == imaginary.num_cols());
    assert(real.num_rows() == imaginary.num_rows());
    int grid_size = real.num_rows();

    // Surely there's a more concise way to write this normalization.
    // Or else, we could just do the multiplication inline, when we're writing to real
    // and imaginary arrays below, and save a bunch of memory writes. I don't know.
    float normalization_factor = box_xy_magnitude / float(grid_size*grid_size);
    for(int i = 0; i < (grid_size*(grid_size/2+1)); i++) {
        array2[i][0] *= normalization_factor;
        array2[i][1] *= normalization_factor;
    }

    // XXX: Can pull out the first row/col special cases to avoid all those if statements in
    // the inner loop. On the other hand, CPU branch prediction may make it a moot point
    // performance-wise and the only remaining reason to do so is clarity of presentation
    for(int row = 0; row < grid_size; row++) {
        for(int col = 0; col < grid_size; col++) {
            // Do something special for the top row
            if(row == 0) {
                if(col <= grid_size/2) {
                    // Left half of the top row
                    real.set(row, col,       array2[row*(grid_size/2+1) + col][0]);
                    imaginary.set(row, col,  array2[row*(grid_size/2+1) + col][1]);
                } else {
                    // Right half of the top row
                    real.set(row, col,       array2[row*(grid_size/2+1) + (grid_size-col)][0]);
                    imaginary.set(row, col, -array2[row*(grid_size/2+1) + (grid_size-col)][1]);
                }
            }

            // Do something special for the leftmost column, except the first row
            if(row != 0 && col == 0) {
                real.set(row, col,       array2[row*(grid_size/2+1) + col][0]);
                imaginary.set(row, col,  array2[row*(grid_size/2+1) + col][1]);
            }

            if(row > 0 && col > 0) {
                if(col <= grid_size/2) {
                    // Left half of the top row
                    real.set(row, col,       array2[row*(grid_size/2+1) + col][0]);
                    imaginary.set(row, col,  array2[row*(grid_size/2+1) + col][1]);
                } else {
                    // Right half of the top row
                    real.set(row, col,       array2[(grid_size-row)*(grid_size/2+1) + (grid_size-col)][0]);
                    imaginary.set(row, col, -array2[(grid_size-row)*(grid_size/2+1) + (grid_size-col)][1]);
                }
            }
        }
    }
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
    cout << "Lines in lipid files: " << lipid_x_lines.size() << " " << lipid_y_lines.size() << " " << lipid_z_lines.size()
         << endl;
    // This lipid variable contains both heads and tails!
    FloatArray2D lipid;
    coalesce_xyz_lines(lipid_x_lines, lipid_y_lines, lipid_z_lines, lipid);

    // Calculate average box size in X and Y dimensions - we'll need this below to wrap lipids
    float avg_a = 0.0, avg_b = 0.0;
    for(int frame_i = 0; frame_i < box_size.num_rows(); frame_i++) {
        avg_a += box_size.get(frame_i, 0);
        avg_b += box_size.get(frame_i, 1);
    }
    avg_a /= float(box_size.num_rows());
    avg_b /= float(box_size.num_rows());

    // Start setting up for the FFT
    FloatArray2D qx_matrix, qy_matrix, cos_q_matrix, sin_q_matrix;
    create_q_matrix(grid_size, qx_matrix, qy_matrix, cos_q_matrix, sin_q_matrix);
    auto z_height_raw = new float[grid_size*grid_size];
    // Number of unique grid pairs, I guess
    int num_grid_pairs = grid_size*(grid_size/2+1);

    // Preallocate this memory so we don't do it every frame
    // ...although we could, and it wouldn't be noticeably slower, probably
    auto z_height_fft = new fftwf_complex[num_grid_pairs];
    std::memset(z_height_fft, 0, sizeof(fftwf_complex) * num_grid_pairs);

    // Probably some of these we could declare within the frame loop, for clarity.
    // After all, there aren't that many frames
    // FIXME: In the original code, many of these are zeroed out for each frame
    // Original code called these dz1xqS, dz2xqS, dz1yqS, dz2yzS
    auto dz_upper_xqS = new fftwf_complex[num_grid_pairs];
    auto dz_upper_yqS = new fftwf_complex[num_grid_pairs];
    auto dz_lower_xqS = new fftwf_complex[num_grid_pairs];
    auto dz_lower_yqS = new fftwf_complex[num_grid_pairs];
    auto dz_upper_xqS_raw = new float[grid_size*grid_size];
    auto dz_upper_yqS_raw = new float[grid_size*grid_size];
    auto dz_lower_xqS_raw = new float[grid_size*grid_size];
    auto dz_lower_yqS_raw = new float[grid_size*grid_size];
    auto dz_xqS = new fftwf_complex[num_grid_pairs];
    auto dz_yqS = new fftwf_complex[num_grid_pairs];
    auto dz_xqS_copy = new fftwf_complex[num_grid_pairs];
    auto dz_yqS_copy = new fftwf_complex[num_grid_pairs];
    auto dz_xqS_raw = new float[grid_size*grid_size];
    auto dz_yqS_raw = new float[grid_size*grid_size];
    double norm_avg = 0.0, dot_cum = 0.0, director_avg = 0.0;
    auto norm_upper_x_raw = new float[grid_size*grid_size];
    auto norm_upper_y_raw = new float[grid_size*grid_size];
    auto norm_upper_z_raw = new float[grid_size*grid_size];
    auto norm_lower_x_raw = new float[grid_size*grid_size];
    auto norm_lower_y_raw = new float[grid_size*grid_size];
    auto norm_lower_z_raw = new float[grid_size*grid_size];


    const int fft_size[2] = {grid_size, grid_size};
    auto spectrum_plan = fftwf_plan_many_dft_r2c(2, fft_size, 1, z_height_raw, nullptr, 1, 0,
            z_height_fft, nullptr, 1, 0, FFTW_MEASURE);
    auto inv_plan = fftwf_plan_many_dft_c2r(2, fft_size, 1, dz_upper_xqS, NULL, 1, 0,
                                               dz_upper_xqS_raw, NULL, 1, 0, FFTW_MEASURE);

    // Go through individual frames
    // Lipid number density, or phi0 is, I guess, number of lipids per XY area divided by 2
    float lipid_number_density = 0;

    for(int frame_i = 0; frame_i < box_size.num_rows(); frame_i++) {
        cout << "== Frame " << frame_i << " ==" << endl;
        // Each frame contains 2 * num_lipids coordinates because there are heads and tails...duh
        FloatArray2D lipid_frame;
        lipid_frame.init_from(lipid, 2*num_lipids*frame_i, 2*num_lipids);

        // Deinterlace the blob of lipid coordinates into heads and tails, and wrap them into the periodic box
        FloatArray2D heads, tails;
        separate_lipid_heads_tails(lipid_frame, heads, tails);
        float box_a = box_size.get(frame_i, 0), box_b = box_size.get(frame_i, 1);
        cout << "Box XY: " << box_a << ", " << box_b << endl;
        wrap_into_box_xy_plane(heads, tails, box_a, box_b, avg_a, avg_b);

        // Calculate director and associated quantities, and wrap lipids in Z direction
        FloatArray2D director;
        // Keep track of which lipids are "good" and therefore can be considered in the calculation
        std::vector<bool> lipid_good;
        // How many lipids in upper and lower leaflets
        int num_upper_leaflet, num_lower_leaflet;
        float avg_head_z;
        float box_z = box_size.get(frame_i, 2);
        wrap_z_and_calc_director(heads, tails, box_z, director, lipid_good, num_upper_leaflet,
                num_lower_leaflet, avg_head_z);

        cout << "Upper and lower leaflet lipid counts: " << num_upper_leaflet << " " << num_lower_leaflet << endl;

        lipid_number_density += 0.5 * (num_upper_leaflet + num_lower_leaflet) / box_size.get(frame_i, 0)
                / box_size.get(frame_i, 1);

        // "CALCULATE NUMBER DENSITIES/////"
        // Not sure what these are yet
        FloatArray2D h_real, h_imag, psi_upper_real, psi_upper_imag, psi_lower_real, psi_lower_imag;
        // OMG this function does a lot of stuff
        calculate_number_densities(heads, tails, director, lipid_good, num_upper_leaflet, num_lower_leaflet,
                qx_matrix, qy_matrix, grid_size, box_a, box_b, avg_head_z,
                h_real, h_imag, psi_upper_real, psi_upper_imag, psi_lower_real, psi_lower_imag);

        // Quantize into coarse-grained grid (c)
        FloatArray2D z_height_upper, z_height_lower;
        // Keep count of how many good and bad lipids are represented at each grid point
        // It's nice they are floats because we will divide by them later.
        FloatArray2D lipid_good_upper, lipid_good_lower, lipid_bad_upper, lipid_bad_lower;
        // Index of grid points by lipid index
        std::vector<int> xj_vec, yj_vec;
        int num_empty_total = map_onto_grid(heads, tails, director, lipid_good, grid_size, box_a, box_b, avg_head_z,
                z_height_upper, z_height_lower, lipid_good_upper, lipid_good_lower, lipid_bad_upper, lipid_bad_lower,
                xj_vec, yj_vec);
        cout << "At least one empty adjacent cell encountered: " << num_empty_total << endl;

        // Calculate height and thickness (d)
        FloatArray2D z_height(grid_size, grid_size); // previously called just 'h'
        z_height.add(z_height_upper).add(z_height_lower);
        FloatArray2D z_thickness(grid_size, grid_size); // previously called just 't'
        z_thickness.add(z_height_upper).subtract(z_height_lower);

        // TODO: calculate t0_frame, tq0_frame, t0, tq0

        // Step (f)
        // "CALCULATE NORMAL VECTORS//////"
        // Copy raw floats into 1D array, suitable for FFTW, then calculate Fourier spectra
        // We copy so we don't have to expose the internal representation of FloatArray2D. Probably a bit silly.
        auto z_height_upper_fft = do_dft2d_r2c(spectrum_plan, z_height_upper);
        auto z_height_lower_fft = do_dft2d_r2c(spectrum_plan, z_height_lower);

        // set wave vector: (Line 1040)
        // Convert from vectors in Fourier space, which must be normalized by box length, to tilt angles

        // Seems like there should be a way to vectorize this
        for(int i = 0; i < grid_size; i++) {
            for(int j = 0; j <= grid_size/2; j++) {
                // Get wavevector component
                float qi = qx_matrix.get(i, j), qj = qy_matrix.get(i, j);
                // Offset into raw FFT output
                int k = i*(grid_size/2 + 1) + j;

                // Handle Nyquist element (aliased between -N/2 and N/2, set the N/2 term of the
                // derivative to zero)
                if(i == grid_size/2) {
                    qi = 0;
                }
                if(j == grid_size/2) {
                    qj = 0;
                }

                float box_xy_magnitude = sqrtf(box_a*box_b);
                float norm = box_xy_magnitude / (grid_size*grid_size);
                // Note the double subscript on the fftwf_complex values, since these are
                // really just float[2], for real and complex components
                dz_upper_xqS[k][0] = -norm*qi*z_height_upper_fft[k][1]*2*M_PI*box_a; // dz1xqS
                dz_upper_xqS[k][1] = -norm*qi*z_height_upper_fft[k][0]*2*M_PI*box_a;
                dz_upper_yqS[k][0] = -norm*qj*z_height_upper_fft[k][1]*2*M_PI*box_b; // dz1yqS
                dz_upper_yqS[k][1] = -norm*qj*z_height_upper_fft[k][0]*2*M_PI*box_b;

                dz_lower_xqS[k][0] = -norm*qi*z_height_lower_fft[k][1]*2*M_PI*box_a; // dz2xqS
                dz_lower_xqS[k][1] = -norm*qi*z_height_lower_fft[k][0]*2*M_PI*box_a;
                dz_lower_yqS[k][0] = -norm*qj*z_height_lower_fft[k][1]*2*M_PI*box_b; // dz2yqS
                dz_lower_yqS[k][1] = -norm*qj*z_height_lower_fft[k][0]*2*M_PI*box_b;
            }
        }

        delete[] z_height_upper_fft;
        delete[] z_height_lower_fft;

        // "*** Defining the Normal as the 2d gradient of the height field ***"
        // "first calculate the gradient in Fourier space"
        do_dft2d_r2c(spectrum_plan, z_height, z_height_fft);

        for(int i = 0; i < grid_size; i++) {
            for(int j = 0; j <= grid_size/2; j++) {
                float qi = qx_matrix.get(i, j), qj = qy_matrix.get(i, j);
                int k = i*(grid_size/2 + 1) + j;

                // For some reason we don't worry about the Nyquist element here

                dz_xqS[k][0] = -qi*z_height_fft[k][1]*2*M_PI*box_a;
                dz_xqS[k][1] =  qi*z_height_fft[k][0]*2*M_PI*box_a;
                dz_yqS[k][0] = -qj*z_height_fft[k][1]*2*M_PI*box_b;
                dz_yqS[k][1] =  qj*z_height_fft[k][0]*2*M_PI*box_b;
            }
        }

        // Backward transform to get derivatives in real space (Line: 1103)
        fftwf_execute_dft_c2r(inv_plan, dz_upper_xqS, dz_upper_xqS_raw);
        fftwf_execute_dft_c2r(inv_plan, dz_upper_yqS, dz_upper_yqS_raw);
        fftwf_execute_dft_c2r(inv_plan, dz_lower_xqS, dz_lower_xqS_raw);
        fftwf_execute_dft_c2r(inv_plan, dz_lower_yqS, dz_lower_yqS_raw);

        // Now multiply box_a and box_b again by dz_xqS and dz_yqS
        // (Line: 1113)
        // dhxqSd -> dz_xqS_copy; dhxqS -> dz_xqS
        for(int k = 0; k < num_grid_pairs; k++) {
            dz_xqS_copy[k][0] = box_a*dz_xqS[k][0];
            dz_xqS_copy[k][1] = box_a*dz_xqS[k][1];
            dz_yqS_copy[k][0] = box_b*dz_yqS[k][0];
            dz_yqS_copy[k][0] = box_b*dz_yqS[k][1];
        }

        fftwf_execute_dft_c2r(inv_plan, dz_xqS_copy, dz_xqS_raw);
        fftwf_execute_dft_c2r(inv_plan, dz_yqS_copy, dz_yqS_raw);

        // normalize: f = (1/L) Sum f_q exp(iq.r) (Line: 1123)
        // dz1x1D -> dz_upper_xqS_raw
        // Original code iterates over i, j for some reason even though these are never used: k = i*ngrid+j
        for(int k = 0; k < (grid_size*grid_size); k++) {
            dz_upper_xqS_raw[k] /= box_a;
            dz_upper_yqS_raw[k] /= box_b;
            dz_lower_xqS_raw[k] /= box_a;
            dz_lower_yqS_raw[k] /= box_b;

            // Normals
            dz_xqS_raw[k] *= 1.0 / ((float) (grid_size * grid_size)) / box_a;
            dz_yqS_raw[k] *= 1.0 / ((float) (grid_size * grid_size)) / box_b;

            float root_ginv_upper = 1.0, root_ginv_lower = 1.0;
            if (!should_norm_z1) {
                root_ginv_upper = 1.0 / sqrt(1.0 + dz_upper_xqS_raw[k] * dz_upper_xqS_raw[k]
                                             + dz_upper_yqS_raw[k] * dz_upper_yqS_raw[k]);
                root_ginv_lower = 1.0 / sqrt(1.0 + dz_lower_xqS_raw[k] * dz_lower_xqS_raw[k]
                                             + dz_lower_yqS_raw[k] * dz_lower_yqS_raw[k]);
            }

            norm_upper_x_raw[k] = dz_upper_xqS_raw[k] * root_ginv_upper;
            norm_upper_y_raw[k] = dz_upper_yqS_raw[k] * root_ginv_upper;
            norm_upper_z_raw[k] = -root_ginv_upper;
            // Note signs reversed for lower leaflet
            norm_lower_x_raw[k] = -dz_lower_xqS_raw[k] * root_ginv_lower;
            norm_lower_y_raw[k] = -dz_lower_yqS_raw[k] * root_ginv_lower;
            norm_lower_z_raw[k] = root_ginv_lower;

            norm_avg += sqrt(norm_lower_x_raw[k] * norm_lower_x_raw[k]
                             + norm_lower_y_raw[k] * norm_lower_y_raw[k]
                             + norm_lower_z_raw[k] * norm_lower_z_raw[k])
                        + sqrt(norm_lower_x_raw[k] * norm_lower_x_raw[k]
                               + norm_lower_y_raw[k] * norm_lower_y_raw[k]
                               + norm_lower_z_raw[k] * norm_lower_z_raw[k]);
            float inv_root = 1 / sqrt(1.0 + dz_xqS_raw[k] * dz_xqS_raw[k] + dz_yqS_raw[k] * dz_yqS_raw[k]);

            // Nnormh is never read in the original source code, so we won't bother to calculate it here
            // Nnormh[k][0] = dhx1D[k]*invRoot;
            // Nnormh[k][1] = dhy1D[k]*invRoot;

            dz_xqS_raw[k] *= inv_root; // dhx1D
            dz_yqS_raw[k] *= inv_root; // dhy1D
        }

        // Number of lipids per grid cell involved in tilt calculations
        // Ostensibly for normalization purposes later
        FloatArray2D lipids_per_cell_upper(grid_size, grid_size);
        FloatArray2D lipids_per_cell_lower(grid_size, grid_size);

        // The grid cell tilt vectors
        FloatArray2D tilt_upper_x(grid_size, grid_size); // t1
        FloatArray2D tilt_upper_y(grid_size, grid_size);
        FloatArray2D tilt_upper_z(grid_size, grid_size);
        FloatArray2D tilt_lower_x(grid_size, grid_size); // t2
        FloatArray2D tilt_lower_y(grid_size, grid_size);
        FloatArray2D tilt_lower_z(grid_size, grid_size);
        FloatArray2D director_upper_x(grid_size, grid_size); // n1: director binned field
        FloatArray2D director_upper_y(grid_size, grid_size);
        FloatArray2D director_upper_z(grid_size, grid_size);
        FloatArray2D director_lower_x(grid_size, grid_size); // n2
        FloatArray2D director_lower_y(grid_size, grid_size);
        FloatArray2D director_lower_z(grid_size, grid_size);

        // Calculate tilt vectors by iterating over lipids
        for(int lipid_i = 0; lipid_i < num_lipids; lipid_i++) {
            // Get grid coordinates for this particular lipid
            int xi = xj_vec[lipid_i], yi = yj_vec[lipid_i];
            // "raw" offset into grid
            int k = xi * grid_size + yi;

            auto dx = director.get(lipid_i, 0);
            auto dy = director.get(lipid_i, 1);
            auto dz = director.get(lipid_i, 2);

            // Look at upper leaflet
            // Tilt vector m = n/(n.N) - N
            if(dz < 0 && lipid_good[lipid_i]) {
                dot_cum += dx*norm_upper_x_raw[k] + dy*norm_upper_y_raw[k] + dz*norm_upper_z_raw[k];
                lipids_per_cell_upper.add(xi, yi, 1.0);
                // nt1++ (we can just sum num_upper_tilt, and the same for num_lower_tilt)
                tilt_upper_x.add(xi, yi, dx - norm_upper_x_raw[k]);
                tilt_upper_y.add(xi, yi, dy - norm_upper_y_raw[k]);
                tilt_upper_z.add(xi, yi, dz - norm_upper_z_raw[k]);
                director_upper_x.add(xi, yi, director.get(lipid_i, 0)); // n1
                director_upper_y.add(xi, yi, director.get(lipid_i, 1));
                director_upper_z.add(xi, yi, director.get(lipid_i, 2));

                // TODO: There's a bunch of stuff with variables u, v, tmag, ut, vt, etc
                // TODO: But I don't know what it is, and it's not repeated for the lower leaflet
                // TODO: So I'm leaving it out for now.
            }

            // Look at lower leaflet
            if(dz > 0 && lipid_good[lipid_i]) {
                dot_cum += dx*norm_lower_x_raw[k] + dy*norm_lower_y_raw[k] + dz*norm_lower_z_raw[k];
                lipids_per_cell_lower.add(xi, yi, 1.0);
                // Does anyone care about nt2?
                tilt_lower_x.add(xi, yi, dx - norm_lower_x_raw[k]);
                tilt_lower_y.add(xi, yi, dy - norm_lower_y_raw[k]);
                tilt_lower_z.add(xi, yi, dz - norm_lower_z_raw[k]);
                director_lower_x.add(xi, yi, director.get(lipid_i, 0)); // n2
                director_lower_y.add(xi, yi, director.get(lipid_i, 1));
                director_lower_z.add(xi, yi, director.get(lipid_i, 2));
            }
        } // loop over lipids

        // Normalize tilt and director values over each patch
        // The original code checked that it wasn't going to cause NaNs
        tilt_upper_x.divide_ignore_zero_denom(lipids_per_cell_upper);
        tilt_upper_y.divide_ignore_zero_denom(lipids_per_cell_upper);
        // Original code doesn't touch Z here
        // tilt_upper_z.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_upper_x.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_upper_y.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_upper_z.divide_ignore_zero_denom(lipids_per_cell_upper);

        tilt_lower_x.divide_ignore_zero_denom(lipids_per_cell_upper);
        tilt_lower_y.divide_ignore_zero_denom(lipids_per_cell_upper);
        // Original code doesn't touch Z here
        // tilt_lower_z.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_lower_x.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_lower_y.divide_ignore_zero_denom(lipids_per_cell_upper);
        director_lower_z.divide_ignore_zero_denom(lipids_per_cell_upper);

        // Fill empty grid cells by interpolation from adjacent cells, for both leaflets (e)
        for(int i = 0; i < grid_size; i++) {
            for(int j = 0; j < grid_size; j++) {
                if(lipid_good_upper.get(i, j) == 0) {
                    interpolate_one_cell(tilt_upper_x, lipids_per_cell_upper, i, j);
                    interpolate_one_cell(tilt_upper_y, lipids_per_cell_upper, i, j);
                    interpolate_one_cell(tilt_upper_z, lipids_per_cell_upper, i, j);
                    interpolate_one_cell(director_upper_x, lipids_per_cell_upper, i, j);
                    interpolate_one_cell(director_upper_y, lipids_per_cell_upper, i, j);
                    interpolate_one_cell(director_upper_z, lipids_per_cell_upper, i, j);
                }
                if(lipid_good_lower.get(i, j) == 0) {
                    interpolate_one_cell(tilt_lower_x, lipids_per_cell_lower, i, j);
                    interpolate_one_cell(tilt_lower_y, lipids_per_cell_lower, i, j);
                    interpolate_one_cell(tilt_lower_z, lipids_per_cell_lower, i, j);
                    interpolate_one_cell(director_lower_x, lipids_per_cell_lower, i, j);
                    interpolate_one_cell(director_lower_y, lipids_per_cell_lower, i, j);
                    interpolate_one_cell(director_lower_z, lipids_per_cell_lower, i, j);
                }
            }
        }

        // Renormalize director grid field ("Added by Itay")
        // 
        for(int i = 0; i < grid_size; i++) {
            for(int j = 0; j < grid_size; j++) {
                // By default normalize only with the Z component
                double norm_factor_upper = abs(director_upper_z.get(i, j));
                double norm_factor_lower = abs(director_lower_z.get(i, j));
                // ...unless we want to use the overall magnitude
                if(!should_norm_z1) { // normZ1
                    norm_factor_upper = sqrt(director_upper_x.get(i, j)*director_upper_x.get(i, j)
                            + director_upper_y.get(i, j)*director_upper_y.get(i, j)
                            + director_upper_z.get(i, j)*director_upper_z.get(i, j));
                    norm_factor_lower = sqrt(director_lower_x.get(i, j)*director_lower_x.get(i, j)
                            + director_lower_y.get(i, j)*director_lower_y.get(i, j)
                            + director_lower_z.get(i, j)*director_lower_z.get(i, j));
                }
                director_upper_x.divide(i, j, norm_factor_upper);
                director_upper_y.divide(i, j, norm_factor_upper);
                director_upper_z.divide(i, j, norm_factor_upper);
                director_lower_x.divide(i, j, norm_factor_lower);
                director_lower_y.divide(i, j, norm_factor_lower);
                director_lower_z.divide(i, j, norm_factor_lower);

                // TODO: dir_av - keep an average of the average z component of the director grid
                // XXX: For some reason we just overwrite only the xy components of the tilt vectors.
                int k = i*grid_size + j;
                tilt_upper_x.set(i, j, director_upper_x.get(i, j) - norm_upper_x_raw[k]);
                tilt_upper_y.set(i, j, director_upper_y.get(i, j) - norm_upper_y_raw[k]);
                tilt_lower_x.set(i, j, director_lower_x.get(i, j) - norm_lower_x_raw[k]);
                tilt_lower_y.set(i, j, director_lower_y.get(i, j) - norm_lower_y_raw[k]);
            }
        } // iteration over grid

        // Accumulate real space orientations of normals and tilts (Line 1358)
        // t1xR_cum -> tilt_upper_real_x_cum, etc...but is never used after being calculated
        // in the original code, so we won't include it at all here
        // Tilt sum and difference: dp and dm
        // Director sum and difference: up and um
        FloatArray2D tilt_sum_x, tilt_sum_y, director_sum_x, director_sum_y;
        FloatArray2D tilt_diff_x, tilt_diff_y, director_diff_x, director_diff_y;
        tilt_sum_x.init_from(tilt_upper_x).add(tilt_lower_x);
        tilt_sum_y.init_from(tilt_upper_y).add(tilt_lower_y);
        tilt_diff_x.init_from(tilt_upper_x).subtract(tilt_lower_x);
        tilt_diff_y.init_from(tilt_upper_y).subtract(tilt_lower_y);
        director_sum_x.init_from(director_upper_x).add(director_lower_x);
        director_sum_y.init_from(director_upper_y).add(director_lower_y);
        director_diff_x.init_from(director_upper_x).subtract(director_lower_x);
        director_diff_y.init_from(director_upper_y).subtract(director_lower_y);
        // TODO: Update umAcc matrices, but that looks messed up in the original source:
        //   --snip--
        //      umAcc[i][j][0] = umAcc[i][j][0] =  + pow(um[i][j][0],2);
        //      umAcc[i][j][1] = umAcc[i][j][1] =  + pow(um[i][j][1],2);
        //   --snip--
        // I guess the intent is to sum the squares of the director differences, over all frames

        // "ACCUMULATE SPECTRA//////" (Line 1385)

        // 
        // float **Nnormh = init_matrix<float>(ngrid*ngrid,2); // top normal vector
        // float **normh = init_matrix<float>(ngrid*ngrid,2); // normal calculate from h

        // TODO: (Line 1399)
        // /************* Normals *****************************************/
        // normhx1D[i*ngrid+j]=normh[i*ngrid + j][0];
        // normhy1D[i*ngrid+j]=normh[i*ngrid + j][1];
        
        // Nnormhx1D[i*ngrid+j] =dhx1D[i*ngrid+j];// Nnormh[i*ngrid+j][0];
        // Nnormhy1D[i*ngrid+j] = dhy1D[i*ngrid+j];//Nnormh[i*ngrid+j][1];
        // /***************************************************************/

        // The FFT of z_height is commented out in the original source for some reason,
        // but the corresponding call to fullArray() is not
        do_dft2d_r2c(spectrum_plan, z_height, z_height_fft);
        // Do FFT of Z-thickness (t1x, t1y)
        auto z_thickness_fft = do_dft2d_r2c(spectrum_plan, z_thickness);

        // hqR: real; hqI: imaginary
        // hqS: ?; Lxy: box lengths
        // We did the FFTs for z_height and z_thickness above
        float box_xy_magnitude = sqrtf(box_a*box_b);
        FloatArray2D z_height_real(grid_size, grid_size), z_height_imaginary(grid_size, grid_size);
        FloatArray2D z_thickness_real(grid_size, grid_size), z_thickness_imaginary(grid_size, grid_size);
        //      fullArray(ngrid, hqR, hqI, hqS, Lxy); // multiply by lx/N^2 factor inside
        fill_full_array(z_height_real, z_height_imaginary, z_height_fft, box_xy_magnitude);
        //      fullArray(ngrid, tqR, tqI, tqS, Lxy);
        fill_full_array(z_thickness_real, z_thickness_imaginary, z_thickness_fft, box_xy_magnitude);

        // if(AREA) (Line 1424)
        // Fill in lower half of complex plane
        // These variables came from calculate_number_densities()
        // psiRU -> psi_upper_real, psiIU -> psi_upper_imag
        // psiRD -> psi_lower_real, psiID -> psi_lower_imag
        for(int i = 1; i < grid_size/2; i++) {
            for(int j = 0; j < grid_size; j++) {
                psi_upper_real.set(grid_size - i, j,  psi_upper_real.get(i, j));
                psi_upper_imag.set(grid_size - i, j, -psi_upper_imag.get(i, j));
                psi_lower_real.set(grid_size - i, j,  psi_lower_real.get(i, j));
                psi_lower_imag.set(grid_size - i, j, -psi_lower_imag.get(i, j));
                h_real.set(grid_size - i, j,  h_real.get(i, j));
                // XXX: For some reason no negation on this one in the original source
                h_imag.set(grid_size - i, j, -h_imag.get(i, j));
            }
        }

        // Convert a bunch of stuff to raw arrays in preparation for more FFTs


        // Do FFT of tilt sum and difference vectors
        // And director sum and difference vectors
        // t1x1D, t1y1D
        // t1x1D is the raw version of t1 (top tilt vector, x component)
        // Tilt fields: dpx1D, dpy1D, dmx1D, dmy1D
        // "Symm and antisymm director fields": upx1D, upy1D, umx1D, umy1D
        // p: sum; m: difference

        
        // fftwf_execute_dft_r2c(spectrum_plan, t1x1D, t1xqS);
        // fftwf_execute_dft_r2c(spectrum_plan, t1y1D, t1yqS);
        auto tilt_upper_x_fft = do_dft2d_r2c(spectrum_plan, tilt_upper_x);
        auto tilt_upper_y_fft = do_dft2d_r2c(spectrum_plan, tilt_upper_y);

        auto tilt_sum_x_fft = do_dft2d_r2c(spectrum_plan, tilt_sum_x);
        auto tilt_sum_y_fft = do_dft2d_r2c(spectrum_plan, tilt_sum_y);
        auto tilt_diff_x_fft = do_dft2d_r2c(spectrum_plan, tilt_diff_x);
        auto tilt_diff_y_fft = do_dft2d_r2c(spectrum_plan, tilt_diff_y);
        auto director_sum_x_fft = do_dft2d_r2c(spectrum_plan, director_sum_x);
        auto director_sum_y_fft = do_dft2d_r2c(spectrum_plan, director_sum_y);
        auto director_diff_x_fft = do_dft2d_r2c(spectrum_plan, director_diff_x);
        auto director_diff_y_fft = do_dft2d_r2c(spectrum_plan, director_diff_y);

        // /////// Normals ////////////////////////////////////////////////
        //  ~Line 1451
        // dhx1D, dhy1D
        // inv_plan(dhxqSd) -> dhx1D
        // dhxqSd -> dz_xqS_copy; dhxqS -> dz_xqS
        // NnormhxqS -> norm_x_fft
        // NnormhyqS -> norm_y_fft
        // normhxR,I comes from fullArray(..., dhxqS, ...)
        //     dhxqS is "height field derivative"
        //     NdhxqS is "normalized height field derivative"
        // NnormhxR,I comes from fullArray(..., NnormhxqS, ...)
        //     NnormhxqS is commented as "tilt of top monolayer" (x-component) (Line: 297)
        //     but I don't believe that, since t1yqS is also annotated that way.
        auto norm_x_fft = do_dft2d_r2c(grid_size, spectrum_plan, dz_xqS_raw);
        auto norm_y_fft = do_dft2d_r2c(grid_size, spectrum_plan, dz_yqS_raw);

        // We use fullArray() now to, I guess, convert an array of
        // fftw_complex back into separate real and imaginary components in a matrix
        FloatArray2D norm_x_fft_real(grid_size, grid_size), norm_x_fft_imaginary(grid_size, grid_size);
        FloatArray2D norm_y_fft_real(grid_size, grid_size), norm_y_fft_imaginary(grid_size, grid_size);
        fill_full_array(norm_x_fft_real, norm_x_fft_imaginary, norm_x_fft, box_xy_magnitude);
        fill_full_array(norm_y_fft_real, norm_y_fft_imaginary, norm_y_fft, box_xy_magnitude);

        // Not sure what the difference is between this and the z_height variables
//        fullArray(ngrid,normhxR,normhxI,dhxqS,Lxy);
//        fullArray(ngrid,normhyR,normhyI,dhyqS,Lxy);
        fill_full_array(...);


//        fullArray(ngrid,t1xR,t1xI,t1xqS,Lxy);
//        fullArray(ngrid,t1yR,t1yI,t1yqS,Lxy);
        FloatArray2D tilt_upper_x_real(grid_size, grid_size), tilt_upper_x_imaginary(grid_size, grid_size);
        FloatArray2D tilt_upper_y_real(grid_size, grid_size), tilt_upper_y_imaginary(grid_size, grid_size);
        fill_full_array(tilt_upper_x_real, tilt_upper_x_imaginary, tilt_upper_x_fft, box_xy_magnitude);
        fill_full_array(tilt_upper_y_real, tilt_upper_y_imaginary, tilt_upper_y_fft, box_xy_magnitude);




        // ngrid -> grid_size

        delete[] norm_x_fft;
        delete[] norm_y_fft;
        delete[] z_thickness_fft;
        delete[] tilt_sum_y_fft;
        delete[] tilt_diff_x_fft;
        delete[] tilt_diff_y_fft;
    } // iteration over frames

    lipid_number_density /= (float) box_size.num_rows();
    cout << "Lipid number density: " << lipid_number_density << endl;

    // Clean up FFTW plans
    fftwf_destroy_plan(inv_plan);
    fftwf_destroy_plan(spectrum_plan);

    // As a plugin, we don't want to pollute the heap
    delete[] z_height_raw;
    delete[] dz_upper_xqS;
    delete[] dz_upper_yqS;
    delete[] dz_lower_xqS;
    delete[] dz_lower_yqS;
    delete[] dz_upper_xqS_raw;
    delete[] dz_upper_yqS_raw;
    delete[] dz_lower_xqS_raw;
    delete[] dz_lower_yqS_raw;
    delete[] dz_xqS;
    delete[] dz_yqS;
    delete[] dz_xqS_copy;
    delete[] dz_yqS_copy;
    delete[] dz_xqS_raw;
    delete[] dz_yqS_raw;
    delete[] norm_upper_x_raw;
    delete[] norm_upper_y_raw;
    delete[] norm_upper_z_raw;
    delete[] norm_lower_y_raw;
    delete[] norm_lower_x_raw;
    delete[] norm_lower_z_raw;
}