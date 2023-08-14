// dph_abs_vis_images_calc_sinus_fit_low_level_1D_CPU.cpp : This file contains the 'main' function. Program execution begins and ends there.
// this program evaluates differential phase, absorption and vibility images from 4D phase CT data using sinu fitting

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <math.h>
#include <chrono>

# define M_PI 3.14159265358979323846

using namespace std;

int levmar_sinus(double* t_data_inp, double* y_data_inp, const int M_inp, double* x0_inp, double* x_fit_outp, int k_max_inp, double eps_1_inp, double eps_2_inp, double tau_inp) {
    // initial iteration variable
    int k = 0;
    int ni = 2;
    // initialize Jacobian 1D matrix
    double* J = new double[3 * M_inp];
    // fill Jacobian matrix
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        J[0 + index_0 * 3] = (-1.0f) * std::sin(t_data_inp[index_0] + x0_inp[1]);
        J[1 + index_0 * 3] = (-1.0f) * x0_inp[0] * std::cos(t_data_inp[index_0] + x0_inp[1]);
        J[2 + index_0 * 3] = -1.0f;
    }
    // initialize transpose of Jacobian matrix in 1D form
    double* J_transpose = new double[3 * M_inp];
    // fill transpose of Jacobian matrix
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        for (int index_1 = 0; index_1 < 3; index_1++) {
            J_transpose[index_0 + index_1 * M_inp] = J[index_1 + index_0 * 3];
        }
    }
    // calculate A matrix
    // initialize A matrix
    double* A = new double[9];
    // initialize A matrix to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        A[index_0] = 0.0f;
    }
    // calculate A matrix as J_transpose * J
    // use multiplication of 2D matrices
    for (int index_0 = 0; index_0 < 3; index_0++) {
        for (int index_1 = 0; index_1 < 3; index_1++) {
            for (int index_2 = 0; index_2 < M_inp; index_2++) {
                A[index_1 + index_0 * 3] = A[index_1 + index_0 * 3] + J_transpose[index_2 + index_0 * M_inp] * J[index_1 + index_2 * 3];
            }
        }
    }
    // calculate f function
    // initialize f function
    double* f = new double[M_inp];
    // fill f function
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        f[index_0] = y_data_inp[index_0] - x0_inp[0] * std::sin(t_data_inp[index_0] + x0_inp[1]) - x0_inp[2];
    }
    // calculate transpose of f
    // initialize transpose of f
    double* f_transpose = new double[M_inp];
    // fill transpose of f
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        f_transpose[index_0] = f[index_0];
    }
    // calculate g as J_transpose * f_transpose
    // initialize g
    double* g = new double[3];
    // initialize g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g[index_0] = 0.0f;
    }
    for (int index_0 = 0; index_0 < 3; index_0++) {
        for (int index_1 = 0; index_1 < 1; index_1++) {
            for (int index_2 = 0; index_2 < M_inp; index_2++) {
                g[index_1 + index_0 * 1] = g[index_1 + index_0 * 1] + J_transpose[index_2 + index_0 * M_inp] * f_transpose[index_1 + index_2 * 1];
            }
        }
    }
    // calculate norm of g
    // initialize norm of g
    double g_norm = 0.0f;
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g_norm = g_norm + g[index_0] * g[index_0];
    }
    g_norm = std::sqrt(g_norm);
    // boolean variable
    bool found_bool = (g_norm <= eps_1_inp);
    // initialize mi
    double mi = 0.0f;
    double A_diag_max = std::max(A[0], A[4]);
    A_diag_max = std::max(A_diag_max, A[8]);
    // calculate mi
    mi = tau_inp * A_diag_max;

    // initialize x vector
    double* x = new double[3];
    // fill x vector
    for (int index_0 = 0; index_0 < 3; index_0++) {
        x[index_0] = x0_inp[index_0];
    }
    // initialize x_new vector
    double* x_new = new double[3];
    // initialize x_new vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        x_new[index_0] = 0.0f;
    }
    // initialize transpose of g
    double* g_transpose = new double[3];
    // initialize transpose of g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g_transpose[index_0] = 0.0f;
    }
    // initialize B matrix
    double* B = new double[9];
    // initialize B matrix to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B[index_0] = 0.0f;
    }
    // initialize inversion matrix of B
    double* B_inv = new double[9];
    // initialize inversion matrix of B to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B_inv[index_0] = 0.0f;
    }
    // initialize adjoint matrix of B
    double* B_adj = new double[9];
    // initialize adjoint matrix of B to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B_adj[index_0] = 0.0f;
    }
    // initialize value of determinant of B
    double B_det = 0.0f;
    // initialize minors of matrix B
    double B_minor_11 = 0.0f;
    double B_minor_12 = 0.0f;
    double B_minor_13 = 0.0f;
    double B_minor_21 = 0.0f;
    double B_minor_22 = 0.0f;
    double B_minor_23 = 0.0f;
    double B_minor_31 = 0.0f;
    double B_minor_32 = 0.0f;
    double B_minor_33 = 0.0f;
    // initialize h_lm vector
    double* h_lm = new double[3];
    // initialize h_lm vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        h_lm[index_0] = 0.0f;
    }
    // initialize transpose of h_lm vector
    double* h_lm_transpose = new double[3];
    // initialize transpose of h_lm vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        h_lm_transpose[index_0] = 0.0f;
    }
    // initialize mi * h_lm_transpose - g
    double* mi_h_lm_trans_minus_g = new double[3];
    // initialize mi * h_lm_transpose - g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        mi_h_lm_trans_minus_g[index_0] = 0.0f;
    }
    // initialize norm of h_lm vector
    double h_lm_norm = 0.0f;
    // initilize norm of x vector
    double x_norm = 0.0f;
    // initialize F_x value
    double F_x = 0.0f;
    // initialize F_x_new value
    double F_x_new = 0.0f;
    // initialize ro_denominator
    double ro_denominator = 0.0f;
    // intialize ro value - gain ratio
    double ro = 0.0f;

    while (!found_bool && (k < k_max_inp)) {
        // increase iteration by one
        k++;
        // calculate matrix B
        // insert matrix A to matrix B
        for (int index_0 = 0; index_0 < 9; index_0++) {
            B[index_0] = A[index_0];
        }
        B[0] = B[0] + mi;
        B[4] = B[4] + mi;
        B[8] = B[8] + mi;
        // calculate transpose of g
        for (int index_0 = 0; index_0 < 3; index_0++) {
            g_transpose[index_0] = g[index_0];
        }
        // calculate inversion of B
        // calculate minor values
        B_minor_11 = (+1.0f) * (B[4] * B[8] - B[5] * B[7]);
        B_minor_12 = (-1.0f) * (B[3] * B[8] - B[5] * B[6]);
        B_minor_13 = (+1.0f) * (B[3] * B[7] - B[4] * B[6]);
        B_minor_21 = (-1.0f) * (B[1] * B[8] - B[2] * B[7]);
        B_minor_22 = (+1.0f) * (B[0] * B[8] - B[2] * B[6]);
        B_minor_23 = (-1.0f) * (B[0] * B[7] - B[1] * B[6]);
        B_minor_31 = (+1.0f) * (B[1] * B[5] - B[2] * B[4]);
        B_minor_32 = (-1.0f) * (B[0] * B[5] - B[2] * B[3]);
        B_minor_33 = (+1.0f) * (B[0] * B[4] - B[1] * B[3]);
        // calculate adjoint matrix of B
        B_adj[0] = B_minor_11;
        B_adj[1] = B_minor_21;
        B_adj[2] = B_minor_31;
        B_adj[3] = B_minor_12;
        B_adj[4] = B_minor_22;
        B_adj[5] = B_minor_32;
        B_adj[6] = B_minor_13;
        B_adj[7] = B_minor_23;
        B_adj[8] = B_minor_33;
        // calculate determinant value of matrix B
        B_det = B[0] * B_minor_11 + B[1] * B_minor_12 + B[2] * B_minor_13;
        // calculate inversion of matrix B
        for (int index_0 = 0; index_0 < 9; index_0++) {
            B_inv[index_0] = B_adj[index_0] / B_det;
        }
        // calculate h_lm vector
        // initialize h_lm vector to zero values
        for (int index_0 = 0; index_0 < 3; index_0++) {
            h_lm[index_0] = 0.0f;
        }
        for (int index_0 = 0; index_0 < 1; index_0++) {
            for (int index_1 = 0; index_1 < 3; index_1++) {
                for (int index_2 = 0; index_2 < 3; index_2++) {
                    h_lm[index_1 + index_0 * 3] = h_lm[index_1 + index_0 * 3] + (-1.0f) * g_transpose[index_2 + index_0 * 3] * B_inv[index_1 + index_2 * 3];
                }
            }
        }
        // calculate norm of h_lm vector
        // fill norm of h_lm vector with zero value
        h_lm_norm = 0.0f;
        for (int index_0 = 0; index_0 < 3; index_0++) {
            h_lm_norm = h_lm_norm + h_lm[index_0] * h_lm[index_0];
        }
        h_lm_norm = std::sqrt(h_lm_norm);
        // calculate norm of x vector
        // fill norm of x vector with zero value
        x_norm = 0.0f;
        for (int index_0 = 0; index_0 < 3; index_0++) {
            x_norm = x_norm + x[index_0] * x[index_0];
        }
        x_norm = std::sqrt(x_norm);
        // main condition
        if (h_lm_norm <= eps_2_inp * (x_norm + eps_2_inp)) {
            found_bool = true;
        }
        else {
            // calculate vector x_new
            for (int index_0 = 0; index_0 < 3; index_0++) {
                x_new[index_0] = x[index_0] + h_lm[index_0];
            }
            // print iteration result
            //std::cout << k << " " << double(x_new[0]) << " " << double(x_new[1]) << " " << double(x_new[2]) << std::endl;
            // calculate F(x)
            // caclulate function f
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                f[index_0] = y_data_inp[index_0] - x[0] * std::sin(t_data_inp[index_0] + x[1]) - x[2];
            }
            // calculate F_x value
            // initialize F_x to zero value 
            F_x = 0.0f;
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                F_x = F_x + f[index_0] * f[index_0];
            }
            F_x = 0.5f * F_x;

            // calculate F(x_new)
            // calculate function f
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                f[index_0] = y_data_inp[index_0] - x_new[0] * std::sin(t_data_inp[index_0] + x_new[1]) - x_new[2];
            }
            // calculate F_x_new value
            // initialize F_x_new to zero value
            F_x_new = 0.0f;
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                F_x_new = F_x_new + f[index_0] * f[index_0];
            }
            F_x_new = 0.5f * F_x_new;

            // calculate ro_denominator part
            // initialize transpose of h_lm vector
            for (int index_0 = 0; index_0 < 3; index_0++) {
                h_lm_transpose[index_0] = h_lm[index_0];
            }
            // calculate mi * h_lm_transpose - g
            // initialize mi * h_lm_transpose - g to zero values
            for (int index_0 = 0; index_0 < 3; index_0++) {
                mi_h_lm_trans_minus_g[index_0] = 0.0f;
            }
            for (int index_0 = 0; index_0 < 3; index_0++) {
                mi_h_lm_trans_minus_g[index_0] = mi * h_lm_transpose[index_0] - g[index_0];
            }
            // calculate ro_denominator
            ro_denominator = 0.0f;
            for (int index_0 = 0; index_0 < 1; index_0++) {
                for (int index_1 = 0; index_1 < 1; index_1++) {
                    for (int index_2 = 0; index_2 < 3; index_2++) {
                        ro_denominator = ro_denominator + h_lm[index_2 + index_0 * 3] * mi_h_lm_trans_minus_g[index_1 + index_2 * 1];
                    }
                }
            }
            ro_denominator = 0.5f * ro_denominator;
            // calculate ro value - gain ratio
            ro = (F_x - F_x_new) / ro_denominator;
            if (ro > 0.0f) {
                // insert vector x_new into the vector x
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    x[index_0] = x_new[index_0];
                }
                // fill Jacobian matrix
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    J[0 + index_0 * 3] = (-1.0f) * std::sin(t_data_inp[index_0] + x[1]);
                    J[1 + index_0 * 3] = (-1.0f) * x[0] * std::cos(t_data_inp[index_0] + x[1]);
                    J[2 + index_0 * 3] = -1.0f;
                }
                // fill transpose of Jacobian matrix
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    for (int index_1 = 0; index_1 < 3; index_1++) {
                        J_transpose[index_0 + index_1 * M_inp] = J[index_1 + index_0 * 3];
                    }
                }
                // calculate A matrix
                // initialize A matrix to zero values
                for (int index_0 = 0; index_0 < 9; index_0++) {
                    A[index_0] = 0.0f;
                }
                // calculate A matrix as J_transpose * J
                // use multiplication of 2D matrices
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    for (int index_1 = 0; index_1 < 3; index_1++) {
                        for (int index_2 = 0; index_2 < M_inp; index_2++) {
                            A[index_1 + index_0 * 3] = A[index_1 + index_0 * 3] + J_transpose[index_2 + index_0 * M_inp] * J[index_1 + index_2 * 3];
                        }
                    }
                }
                // calculate f function
                // fill f function
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    f[index_0] = y_data_inp[index_0] - x[0] * std::sin(t_data_inp[index_0] + x[1]) - x[2];
                }
                // calculate transpose of f
                // fill transpose of f
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    f_transpose[index_0] = f[index_0];
                }
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    g[index_0] = 0.0f;
                }
                // calculate g as J_transpose * f_transpose
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    for (int index_1 = 0; index_1 < 1; index_1++) {
                        for (int index_2 = 0; index_2 < M_inp; index_2++) {
                            g[index_1 + index_0 * 1] = g[index_1 + index_0 * 1] + J_transpose[index_2 + index_0 * M_inp] * f_transpose[index_1 + index_2 * 1];
                        }
                    }
                }
                // calculate norm of g
                // initialize norm of g
                g_norm = 0.0f;
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    g_norm = g_norm + g[index_0] * g[index_0];
                }
                g_norm = std::sqrt(g_norm);
                // calculate boolean variable
                found_bool = (g_norm <= eps_1_inp);
                // calculate mi
                mi = mi * std::max(double(0.33333333333333333333f), double(1.0f - std::pow((2.0f * ro - 1), 3.0f)));
                // define ni
                ni = 2;
            }
            else {
                // calculate mi
                mi = mi * double(ni);
                // calculate ni
                ni = 2 * ni;
            }
        }
    }

    ////std::cout << x_new << std::endl;
    //// convert phase shift from fitting to interval (0, 2*pi)
    //if (x_new[1] > 0.0f && x_new[0] > 0.0f) {
    //    x_new[1] = x_new[1] - (2.0f * M_PI) * int(x_new[1] / (2 * M_PI));
    //}
    //else if (x_new[1] < 0.0f && x_new[0] > 0.0f) {
    //    x_new[1] = x_new[1] - (2.0f * M_PI) * (int(x_new[1] / (2 * M_PI)) - 1);
    //}
    //else if (x_new[1] > 0.0f && x_new[0] < 0.0f) {
    //    x_new[0] = (-1.0f) * x_new[0];
    //    x_new[1] = x_new[1] + 1.0f * M_PI;
    //    x_new[1] = x_new[1] - (2.0f * M_PI) * int(x_new[1] / (2 * M_PI));
    //}
    //else if (x_new[1] < 0.0f && x_new[0] < 0.0f) {
    //    x_new[0] = (-1.0f) * x_new[0];
    //    x_new[1] = x_new[1] - 1.0f * M_PI;
    //    x_new[1] = x_new[1] - (2.0f * M_PI) * (int(x_new[1] / (2 * M_PI)) - 1);
    //}
    //else {
    //    x_new[1] = 0.0f;
    //}
    ////std::cout << x_new << std::endl;

    // store fitting results to output 1D double array
    x_fit_outp[0] = x_new[0];
    x_fit_outp[1] = x_new[1];
    x_fit_outp[2] = x_new[2];

    delete[] J;
    delete[] J_transpose;
    delete[] A;
    delete[] f;
    delete[] f_transpose;
    delete[] g;
    delete[] g_transpose;
    delete[] x;
    delete[] x_new;
    delete[] B;
    delete[] B_inv;
    delete[] B_adj;
    delete[] h_lm;
    delete[] h_lm_transpose;
    delete[] mi_h_lm_trans_minus_g;

    return 0;
}

int main()
{
    // define path to folder with all foreground data or all subfolders
    string path_to_fg_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/pp/fg/");
    // define path to folder with all background data or all subfolders
    string path_to_bg_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/bg/");
    // print path to folder with all foreground folders
    std::cout << path_to_fg_folder << "\n";
    // print path to folder with all background folders
    std::cout << path_to_bg_folder << "\n";

    // define path to output folder with output differential phase (dph) images
    string path_to_output_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/pp_dph_abs_vis_Momose_sinus_fitting_low_level_1D_CPU/");

    // output image name root - differential phase image or dph image
    string image_output_dph_name_root = "dph";
    // output image name root - absorption image or abs image
    string image_output_abs_name_root = "abs";
    // output image name root - visibility image or vis image
    string image_output_vis_name_root = "vis";
    // final ouput image name - differential phase image or dph image
    string image_output_dph_name;
    // final ouput image name - absorption image or abs image
    string image_output_abs_name;
    // final ouput image name - visibility image or vis image
    string image_output_vis_name;
    // extension of output image
    string image_output_extension = ".raw";

    // define size of the raw unsigned 16 bit images
    const unsigned int no_cols = 1536; // in pixels, in horizontal direction
    const unsigned int no_rows = 512; // in pixels, in vertical direction
    // total number of pixels in single image
    const unsigned int no_pixels = no_cols * no_rows;
    // total number of bytes in single image, we consider 16 bit values per pixel = 2 bytes
    const unsigned int no_bytes = 2 * no_pixels;

    // define number of initial and final subfolder for foreground
    unsigned int no_subfolder_fg_initial = 1;
    unsigned int no_subfolder_fg_final = 11200;

    // define number of initial and final folder for background
    unsigned int no_subfolder_bg_initial = 1;
    unsigned int no_subfolder_bg_final = 11200;

    // number of digits in subfolder name for foreground
    string::size_type no_subfolder_digits_fg = 6;
    // number of digits in subfolder name for background
    string::size_type no_subfolder_digits_bg = 6;

    // number of steps in fringe scanning technique
    const unsigned int M = 5;

    // calculate differential phase image for foreground
    // fringe scanning defined from initial value
    const unsigned int M_fg_initial = 1;
    // fringe scanning defined to final value
    const unsigned int M_fg_final = M;
    // number of steps in fringe scanning for foreground
    const unsigned int M_fg = M_fg_final - M_fg_initial + 1;

    // calculate differential phase image for background
    // fringe scanning defined from initial value
    const unsigned int M_bg_initial = 1;
    // fringe scanning defined to final value
    const unsigned int M_bg_final = M;
    // number of steps in fringe scanning for background
    const unsigned int M_bg = M_bg_final - M_bg_initial + 1;

    // define root name of images for foreground
    string root_image_name_fg("a");
    // define root name of images for background
    string root_image_name_bg("a");

    // number of digits in image name for foreground
    string::size_type no_image_digits_fg = 6;
    // number of digits in image name for background
    string::size_type no_image_digits_bg = 6;

    // define image extensions
    // image extension for foreground
    string image_extension_fg = ".raw";
    // image extension for background
    string image_extension_bg = ".raw";

    // allocate image buffer for foreground
    auto image_buffer_fg = new unsigned short int[no_pixels][M_fg];
    // allocate image buffer for background
    auto image_buffer_bg = new unsigned short int[no_pixels][M_fg];

    // allocate phase buffer for foreground
    double* phase_buffer_fg = new double[no_pixels];
    // allocate phase buffer for background
    double* phase_buffer_bg = new double[no_pixels];
    // allocate amplitude buffer for foreground
    double* amp_buffer_fg = new double[no_pixels];
    // allocate amplitude buffer for background
    double* amp_buffer_bg = new double[no_pixels];
    // allocate offset buffer for foreground
    double* offset_buffer_fg = new double[no_pixels];
    // allocate offset buffer for background
    double* offset_buffer_bg = new double[no_pixels];

    // allocate memory for differential phase image
    double* dph_image = new double[no_pixels];
    // allocate memory for absorption image
    double* abs_image = new double[no_pixels];
    // allocate memory for visibility image
    double* vis_image = new double[no_pixels];

    // define phase for foreground
    double phase_step_fg = (2 * M_PI) / M_fg;
    // define phase_step for background
    double phase_step_bg = (2 * M_PI) / M_bg;

    // auxiliary variables for iteration through subfolder name for foreground
    string subfolder_name(no_subfolder_digits_fg, '0');
    string subfolder_number = "";
    string::size_type counter_digits = 0;
    string::size_type difference = 0;
    string::size_type counter = 0;
    string path_to_fg_subfolder = "";

    // auxiliary variables for iteration through M_fg images
    int counter_image = 0;
    string image_name = root_image_name_fg;
    string image_name_number(no_image_digits_fg, '0');
    string image_number = "";
    // counter_digits, difference and counter variables are taken from iterations through subfolders
    string path_to_fg_image = "";

    // auxiliary variables for iteration through subfolder name for background
    string path_to_bg_subfolder = "";

    // auxiliary variables for iteration through M_bg images
    image_name = root_image_name_bg;
    image_name_number = string(no_image_digits_bg, '0');
    image_number = "";
    // counter_digits, difference and counter variables are taken from iterations through subfolders
    string path_to_bg_image = "";

    // declare auxiliary variable for output image
    string path_to_output_image = "";

    // declare auxiliary variable for output dph image
    string path_to_output_dph_image = "";
    // declare auxiliary variable for output abs image
    string path_to_output_abs_image = "";
    // declare auxiliary variable for output vis image
    string path_to_output_vis_image = "";

    // define intial parameters of sinusoidal function for fitting
    // intial amplitude
    double x01 = 0.0f;
    // initial phase shift
    double x02 = M_PI / 2.0f;
    // initial offset
    double x03 = 0.0f;
    // define initial parameter vector
    double* x0 = new double[3];
    // fill initial parameters vector
    x0[0] = x01;
    x0[1] = x02;
    x0[2] = x03;
    // initialize variables for calculation initial parameters x0
    // for foreground
    double y_data_fg_max = 0.0f;
    double y_data_fg_min = 0.0f;
    // for background
    double y_data_bg_max = 0.0f;
    double y_data_bg_min = 0.0f;
    // define input t_data 1D array for foreground, on the t axis 
    double* t_data_fg = new double[M_fg];
    for (unsigned int index_0 = 0; index_0 < M_fg; index_0++) {
        t_data_fg[index_0] = double(index_0) * ((2.0f * M_PI) / double(M_fg));
    }
    // define input t_data 1D array for background, on the t axis
    double* t_data_bg = new double[M_bg];
    for (unsigned int index_0 = 0; index_0 < M_bg; index_0++) {
        t_data_bg[index_0] = double(index_0) * ((2.0f * M_PI) / double(M_bg));
    }
    // define input y_data 1D array for foreground
    double* y_data_fg = new double[M_fg];
    // define input y_data 1D array for background
    double* y_data_bg = new double[M_bg];
    // maximum number of iterations in fitting
    unsigned int k_max = 1000;
    // auxiliar fitting variable epsilon 1
    double eps_1 = 1.0E-8;
    // auxiliar fitting variable epsilon 2
    double eps_2 = 1.0E-8;
    // auxiliar fitting variable tau
    double tau = 1.0E-3;
    // create output 1D array where fitting results will be stored
    double* x_fit = new double[3];

    // go through all foreground subfolders and foreground images
    for (unsigned int index_0 = no_subfolder_fg_initial; index_0 <= no_subfolder_fg_final; index_0++) {
        // start to measure elapsed time at the beginning
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // initialize subfolder name to "000000"
        subfolder_name = string(no_subfolder_digits_fg, '0');
        // typcast integer value to string, convert integer value to string
        subfolder_number = std::to_string(index_0);
        // initialize digits counter
        counter_digits = subfolder_number.size();
        // initialize difference
        difference = no_subfolder_digits_fg - counter_digits;
        // initialize counter
        counter = 0;
        // generate subfolder name
        for (string::size_type index_1 = difference; index_1 < no_subfolder_digits_fg; index_1++) {
            subfolder_name[index_1] = subfolder_number[counter];
            counter++;
        }
        // generate path to foreground subfolder 
        path_to_fg_subfolder = path_to_fg_folder + subfolder_name;
        // print final path to foreground subfolder
        std::cout << path_to_fg_subfolder << "\n";
        // initilize counter for one of M images
        counter_image = 0;
        // open images in foreground subfolder
        for (unsigned int index_2 = M_fg_initial; index_2 <= M_fg_final; index_2++) {
            // initialize image name to root value
            image_name = root_image_name_fg;
            // initialize image number to "000000"
            image_name_number = string(no_image_digits_fg, '0');
            // typcast integer value to string, convert integer value to string
            image_number = std::to_string(index_2);
            // initialize digits counter
            counter_digits = image_number.size();
            // initialize difference
            difference = no_image_digits_fg - counter_digits;
            // initialize counter
            counter = 0;
            // generate image name number
            for (string::size_type index_3 = difference; index_3 < no_image_digits_fg; index_3++) {
                image_name_number[index_3] = image_number[counter];
                counter++;
            }
            // concatenate root image value and image number
            image_name += image_name_number;
            // generate path to foreground image 
            path_to_fg_image = path_to_fg_subfolder + "/" + image_name + image_extension_fg;
            // print path to image for foregrouund
            //std::cout << path_to_fg_image << "\n";
            //************************************************************************************
            // read binary images
            //************************************************************************************
            ifstream raw_image(path_to_fg_image, ios::out | ios::binary);
            /*streampos begin, end;
            begin = raw_image.tellg();
            raw_image.seekg(0, ios::end);
            end = raw_image.tellg();*/
            //std::cout << "Size of the raw image is: " << (end - begin) << " bytes.\n";
            if (raw_image.is_open())
            {
                //unsigned int counter_pixel = 0;
                //raw_image.seekg(0, ios::beg);
                //while (raw_image.read(reinterpret_cast<char*>(&image_buffer_fg[counter_pixel][counter_image]), sizeof(uint16_t))) { // Read 16-bit integer values from file
                //    counter_pixel++;
                //}
                //raw_image.close();
                for (unsigned int counter_pixel = 0; counter_pixel < no_pixels; counter_pixel++) {
                    raw_image.read((char*)&image_buffer_fg[counter_pixel][counter_image], sizeof(unsigned short int));
                }
                raw_image.close();
            }
            else {
                std::cout << "Warning: Unable to open raw image file!!!" << "\n";
            }
            //************************************************************************************
            // end of reading of binary images
            //************************************************************************************
            // increase image counter by one
            counter_image++;
        }

        // go through background subfolder and background images
        // initialize subfolder name to "000000"
        subfolder_name = string(no_subfolder_digits_bg, '0');
        // typcast integer value to string, convert integer value to string
        subfolder_number = std::to_string(index_0);
        // initialize digits counter
        counter_digits = subfolder_number.size();
        // initialize difference
        difference = no_subfolder_digits_bg - counter_digits;
        // initialize counter
        counter = 0;
        // generate subfolder name
        for (string::size_type index_1 = difference; index_1 < no_subfolder_digits_bg; index_1++) {
            subfolder_name[index_1] = subfolder_number[counter];
            counter++;
        }
        // generate path to background subfolder 
        path_to_bg_subfolder = path_to_bg_folder + subfolder_name;
        // print final path to background subfolder
        std::cout << path_to_bg_subfolder << "\n";
        // initilize counter for one of M images
        counter_image = 0;
        // open images in background subfolder
        for (unsigned int index_2 = M_bg_initial; index_2 <= M_bg_final; index_2++) {
            // initialize image name to root value
            image_name = root_image_name_bg;
            // initialize image number to "000000"
            image_name_number = string(no_image_digits_bg, '0');
            // typcast integer value to string, convert integer value to string
            image_number = std::to_string(index_2);
            // initialize digits counter
            counter_digits = image_number.size();
            // initialize difference
            difference = no_image_digits_bg - counter_digits;
            // initialize counter
            counter = 0;
            // generate image name number
            for (string::size_type index_3 = difference; index_3 < no_image_digits_bg; index_3++) {
                image_name_number[index_3] = image_number[counter];
                counter++;
            }
            // concatenate root image value and image number
            image_name += image_name_number;
            // generate path to background image 
            path_to_bg_image = path_to_bg_subfolder + "/" + image_name + image_extension_bg;
            // print path to image for background
            //std::cout << path_to_bg_image << "\n";
            //************************************************************************************
            // read binary images
            //************************************************************************************
            ifstream raw_image(path_to_bg_image, ios::out | ios::binary);
            /*streampos begin, end;
            begin = raw_image.tellg();
            raw_image.seekg(0, ios::end);
            end = raw_image.tellg();*/
            //std::cout << "Size of the raw image is: " << (end - begin) << " bytes.\n";
            if (raw_image.is_open())
            {
                //unsigned int counter_pixel = 0;
                //raw_image.seekg(0, ios::beg);
                //while (raw_image.read(reinterpret_cast<char*>(&image_buffer_bg[counter_pixel][counter_image]), sizeof(uint16_t))) { // Read 16-bit integer values from file
                //    counter_pixel++;
                //}
                //raw_image.close();
                raw_image.seekg(0, ios::beg);
                for (unsigned int counter_pixel = 0; counter_pixel < no_pixels; counter_pixel++) {
                    raw_image.read((char*)&image_buffer_bg[counter_pixel][counter_image], sizeof(unsigned short int));
                }
                raw_image.close();
            }
            else {
                std::cout << "Warning: Unable to open raw image file!!!" << "\n";
            }
            //************************************************************************************
            // end of reading of binary images
            //************************************************************************************
            // increase image counter by one
            counter_image++;
        }

        // calculate phase image for foreground
        for (unsigned int index_8 = 0; index_8 < no_pixels; index_8++) {
            x0[2] = 0.0f;
            for (unsigned int index_9 = 0; index_9 < M_fg; index_9++) {
                y_data_fg[index_9] = double(image_buffer_fg[index_8][index_9]);
                if (index_9 == 0) {
                    y_data_fg_max = y_data_fg[0];
                    y_data_fg_min = y_data_fg[0];
                }
                else {
                    y_data_fg_max = std::max(y_data_fg_max, y_data_fg[index_9]);
                    y_data_fg_min = std::min(y_data_fg_min, y_data_fg[index_9]);
                }
                x0[2] += y_data_fg[index_9] / M_fg;
            }
            x0[0] = (y_data_fg_max - y_data_fg_min) / 2;
            levmar_sinus(t_data_fg, y_data_fg, M_fg, x0, x_fit, k_max, eps_1, eps_2, tau);
            //std::cout << "foreground: " << index_8 << std::endl;
            phase_buffer_fg[index_8] = x_fit[1];
            amp_buffer_fg[index_8] = x_fit[0];
            offset_buffer_fg[index_8] = x_fit[2];
        }
        // calculate phase image for background
        for (unsigned int index_8 = 0; index_8 < no_pixels; index_8++) {
            x0[2] = 0.0f;
            for (unsigned int index_9 = 0; index_9 < M_bg; index_9++) {
                y_data_bg[index_9] = double(image_buffer_bg[index_8][index_9]);
                if (index_9 == 0) {
                    y_data_bg_max = y_data_bg[0];
                    y_data_bg_min = y_data_bg[0];
                }
                else {
                    y_data_bg_max = std::max(y_data_bg_max, y_data_bg[index_9]);
                    y_data_bg_min = std::min(y_data_bg_min, y_data_bg[index_9]);
                }
                x0[2] += y_data_bg[index_9] / M_bg;
            }
            x0[0] = (y_data_bg_max - y_data_bg_min) / 2;
            levmar_sinus(t_data_bg, y_data_bg, M_bg, x0, x_fit, k_max, eps_1, eps_2, tau);
            //std::cout << "background: " << index_8 << std::endl;
            phase_buffer_bg[index_8] = x_fit[1];
            amp_buffer_bg[index_8] = x_fit[0];
            offset_buffer_bg[index_8] = x_fit[2];
        }

        // calculate differential phase image or dph image
        for (unsigned int index_10 = 0; index_10 < no_pixels; index_10++) {
            dph_image[index_10] = phase_buffer_fg[index_10] - phase_buffer_bg[index_10];
            abs_image[index_10] = offset_buffer_fg[index_10] / offset_buffer_bg[index_10];
            vis_image[index_10] = (amp_buffer_fg[index_10] / offset_buffer_fg[index_10]) / (amp_buffer_bg[index_10] / offset_buffer_bg[index_10]);
        }

        // define name for output dph image for current subfolder
        image_output_dph_name = image_output_dph_name_root + "_" + subfolder_name + image_output_extension;
        // define name for output abs image for current subfolder
        image_output_abs_name = image_output_abs_name_root + "_" + subfolder_name + image_output_extension;
        // define name for output vis image for current subfolder
        image_output_vis_name = image_output_vis_name_root + "_" + subfolder_name + image_output_extension;
        // define path to the output dph image
        path_to_output_dph_image = path_to_output_folder + image_output_dph_name;
        // define path to the output abs image
        path_to_output_abs_image = path_to_output_folder + image_output_abs_name;
        // define path to the output vis image
        path_to_output_vis_image = path_to_output_folder + image_output_vis_name;
        // write differential phase (dph) image
        // set for output, binary data, trunc
        fstream output_dph_image(path_to_output_dph_image, ios::out | ios::binary | ios::trunc);
        if (output_dph_image.is_open())
        {
            // set pointer to the beginning of the image
            output_dph_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_dph_image.write((char*)&dph_image[index_11], sizeof(double));
            }
            output_dph_image.close();
        }
        else {
            std::cout << "Warning: Unable to open dph image file!!!" << "\n";
        }
        // write absorption (abs) image
        // set for output, binary data, trunc
        fstream output_abs_image(path_to_output_abs_image, ios::out | ios::binary | ios::trunc);
        if (output_abs_image.is_open())
        {
            // set pointer to the beginning of the image
            output_abs_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_abs_image.write((char*)&abs_image[index_11], sizeof(double));
            }
            output_abs_image.close();
        }
        else {
            std::cout << "Warning: Unable to open abs image file!!!" << "\n";
        }
        // write visibility (vis) image
        // set for output, binary data, trunc
        fstream output_vis_image(path_to_output_vis_image, ios::out | ios::binary | ios::trunc);
        if (output_vis_image.is_open())
        {
            // set pointer to the beginning of the image
            output_vis_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_vis_image.write((char*)&vis_image[index_11], sizeof(double));
            }
            output_vis_image.close();
        }
        else {
            std::cout << "Warning: Unable to open vis image file!!!" << "\n";
        }

        // stop to measure elapsed time at the end
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // print elapsed time in milliseconds, microseconds and nanoseconds
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[seconds]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[millisec]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microsec]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[nanosec]" << std::endl;
    }

    // delete image buffer for foreground
    delete[] image_buffer_fg;
    // delete image buffer for background
    delete[] image_buffer_bg;
    // delete phase buffer for foreground
    delete[] phase_buffer_fg;
    // delete phase buffer for background
    delete[] phase_buffer_bg;
    // delete amplitude buffer for foreground
    delete[] amp_buffer_fg;
    // delete amplitude buffer for background
    delete[] amp_buffer_bg;
    // delete offset buffer for foreground
    delete[] offset_buffer_fg;
    // delete offset buffer for background
    delete[] offset_buffer_bg;
    // delete buffer for differential phase image
    delete[] dph_image;
    // delete buffer for absorption image
    delete[] abs_image;
    // delete buffer for visibility image
    delete[] vis_image;

    // delete buffers created for fitting
    delete[] x0;
    delete[] t_data_fg;
    delete[] t_data_bg;
    delete[] y_data_fg;
    delete[] y_data_bg;
    delete[] x_fit;

    return 0;
}
