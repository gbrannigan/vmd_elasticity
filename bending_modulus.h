#pragma once
#include <vector>
#include <string>
#include <fftw3.h>
#include "floatarray2d.h"

using namespace std;

float *average_q_in_fft(FloatArray2D &fft_array, int nyquist = 0);
void calculate_number_densities(FloatArray2D &heads, FloatArray2D &tails, FloatArray2D &director, std::vector<bool> &lipid_good,
                                int num_upper_leaflet, int num_lower_leaflet,
                                FloatArray2D &qx_matrix, FloatArray2D &qy_matrix, int grid_size, float a, float b, float avg_head_z,
                                FloatArray2D &h_real, FloatArray2D &h_imag,
                                FloatArray2D &psi_upper_real, FloatArray2D &psi_upper_imag, FloatArray2D &psi_lower_real,
                                FloatArray2D &psi_lower_imag);
void create_q_matrix(int grid_size, FloatArray2D &qx_matrix, FloatArray2D &qy_matrix,
                     FloatArray2D &cos_q_matrix, FloatArray2D &sin_q_matrix);
fftwf_complex *do_dft2d_r2c(fftwf_plan plan, FloatArray2D &data, fftwf_complex *fft_data = nullptr);
fftwf_complex *do_dft2d_r2c(int grid_size, fftwf_plan plan, float *raw_data, fftwf_complex *fft_data = nullptr);
void fill_full_array(FloatArray2D &real, FloatArray2D &imaginary, fftwf_complex *array2,
                     float box_xy_magnitude);
int interpolate_one_cell(FloatArray2D &data, FloatArray2D &good_count, int i, int j);
int map_onto_grid(FloatArray2D &heads, FloatArray2D &tails, FloatArray2D &director,
                  std::vector<bool> &lipid_good, int grid_size,
                  float box_a, float box_b, float avg_head_z,
                  FloatArray2D &z_height_upper, FloatArray2D &z_height_lower,
                  FloatArray2D &good_count_upper, FloatArray2D &good_count_lower,
                  FloatArray2D &bad_count_upper, FloatArray2D &bad_count_lower,
                  std::vector<int> &xj_vec, std::vector<int> &yj_vec);
bool separate_lipid_heads_tails(const FloatArray2D &lipid, FloatArray2D &out_heads, FloatArray2D &out_tails);
void wrap_into_box_xy_plane(FloatArray2D &heads, FloatArray2D &tails, float a, float b, float avg_a, float avg_b);
void wrap_z_and_calc_director(FloatArray2D &heads, FloatArray2D &tails, float c,
                              FloatArray2D &out_director, std::vector<bool> &out_lipid_good,
                              int &out_num_upper_leaflet, int &out_num_lower_leaflet, float &out_head_avg_z);

void do_bending_modulus(FloatArray2D &lipid, FloatArray2D &box_size, int num_lipids, int grid_size,
                        const bool should_norm_z1 = false);

