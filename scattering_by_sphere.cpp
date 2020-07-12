//
// Copyright (c) 2020 Muneyoshi Suzuki. All Rights Reserved. 
//
// 2020-7-9
//
// scattering_by_sphere.cpp: 
//

//
// Numerical analysis of frequency transfer function between a plane wave and scattering by a spherical obstacle.
//

#include <numbers>
#include <cmath>
#include <complex>
#include <array>
#include <iostream>
#include <fstream>

using namespace std;

// max value of n
constexpr unsigned int TBLSIZE {400};

// step of angle [deg]
constexpr unsigned int step_angle {3};

// number of angles
constexpr unsigned int n_angle {180 / step_angle + 1};

// change rate in convergence test
constexpr double epsilon {0.0001};

// number of consecutive convergence in truncation
constexpr unsigned int n_conv {6};

// j: imaginary unit
constexpr complex<double> j {1.0i};

// pi: circle ratio
constexpr double pi {numbers::pi};

// c: acoustic velocity [m/s]
constexpr double c {343.7};

// a: radius of spherical obstacle [m]
constexpr double a {0.0715};

// f: frequency [Hz]
constexpr double f[] {100, 200, 300, 400, 500, 600, 700, 800, 900,
                    1'000, 1'200, 1'400, 1'600, 1'800, 2'000,
                    3'000, 4'000, 5'000, 6'000, 7'000, 8'000, 9'000,
                    10'000, 11'000, 12'000, 13'000, 14'000, 15'000,
                    16'000, 17'000, 18'000, 19'000, 20'000};

constexpr size_t fsize {sizeof(f) / sizeof(*f)};

// ka = (k * a); k: wave number, k = (2 * pi * f / c)  
array<double, fsize> ka_tbl;

// spherical Bessel functions table 
array<array<double, TBLSIZE + 1>, fsize> j_tbl;

// spherical Neumann functions table 
array<array<double, TBLSIZE + 1>, fsize> n_tbl;

// spherical Hankel functions (2) table
array<array<complex<double>, TBLSIZE + 1>, fsize> h2_tbl;

// derivative of spherical Bessel functions table 
array<array<double, TBLSIZE>, fsize> Dj_tbl;

// derivative of spherical Neumann functions table 
array<array<double, TBLSIZE>, fsize> Dn_tbl;

// derivative of spherical Hankel functions (2) table
array<array<complex<double>, TBLSIZE>, fsize> Dh2_tbl;

// for test
// array<array<complex<double>, TBLSIZE>, fsize> test1;
// array<array<complex<double>, TBLSIZE>, fsize> test2;

// Legendre polynomials table
array<array<double, TBLSIZE>, n_angle> P_tbl;

// frequency transfer function table
array<array<double, fsize>, n_angle> G_tbl;

int main()
{
    // initialize ka_tbl
    for (int x = 0; x < fsize; x++)
        ka_tbl.at(x) = f[x] * 2.0 * pi * a / c;

    // initialize j_tbl, n_tbl, and h2_tbl
    for (int x = 0; x < fsize; x++)
        for (int n = 0; n < TBLSIZE + 1; n++) {
            j_tbl.at(x).at(n) = sph_bessel(n, ka_tbl.at(x));
            n_tbl.at(x).at(n) = sph_neumann(n, ka_tbl.at(x));
            h2_tbl.at(x).at(n) = j_tbl.at(x).at(n) - j * n_tbl.at(x).at(n);
        }

    // initialize Dj_tbl, Dn_tbl, and Dh2_tbl
    for (int x = 0; x < fsize; x++)
        for (int n = 0; n < TBLSIZE; n++) {
            Dj_tbl.at(x).at(n) = j_tbl.at(x).at(n) * n / ka_tbl.at(x) - j_tbl.at(x).at(n + 1);
            Dn_tbl.at(x).at(n) = n_tbl.at(x).at(n) * n / ka_tbl.at(x) - n_tbl.at(x).at(n + 1);
            Dh2_tbl.at(x).at(n) = Dj_tbl.at(x).at(n) - j * Dn_tbl.at(x).at(n);
        }

    // test
    // for (int x = 0; x < fsize; x++)
    //    for (int n = 0; n < TBLSIZE; n++) {
    //        test1.at(x).at(n) = j_tbl.at(x).at(n) - Dj_tbl.at(x).at(n) * h2_tbl.at(x).at(n) / Dh2_tbl.at(x).at(n);
    //        test2.at(x).at(n) = -1.0 * j / pow (ka_tbl.at(x), 2.0)/ Dh2_tbl.at(x).at(n);
    //    }

    // initialize P_tbl
    for (int y = 0; y < n_angle; y++)
        for (int n = 0; n < TBLSIZE; n++)
            P_tbl.at(y).at(n) = legendre(n, cos(2.0 * pi * step_angle * y / 360.0));


    for (int y = 0; y < n_angle; y++) {
        for (int x = 0; x < fsize; x++) {
            complex<double> nume {0.0 + j * 0.0};   // numerator
            complex<double> deno {0.0 + j * 0.0};   // denominator
            unsigned int conv {0};                  // number of consecutive convergence
            for (int n = 0; n < TBLSIZE; n++) {
                nume += pow(j, n + 1.0) * (2.0 * n + 1) * P_tbl.at(y).at(n) / Dh2_tbl.at(x).at(n);
                deno += pow(j, n + 0.0) * (2.0 * n + 1) * P_tbl.at(y).at(n) * j_tbl.at(x).at(n);
                double G = abs(nume / deno);
                double delta = G - G_tbl.at(y).at(x);
                G_tbl.at(y).at(x) = G;
                if (abs(delta / G) < epsilon) {     // convergence test
                    if (++conv > n_conv) {
                        cout << y << "\t" << x << "\t" << n << endl;
                        break;
                    }
                } else {
                    conv = 0;
                }
            }
            G_tbl.at(y).at(x) /= pow(ka_tbl.at(x), 2.0);
        }
    }

    ofstream fd("sphere.txt", ios_base::out | ios_base::trunc);
    for (int x = 0; x < fsize; x++)
        fd << "\t" << f[x];
    fd << endl;
    for (int y = 0; y < n_angle; y++) {
        fd << y * step_angle;
        for (int x = 0; x < fsize; x++)
            fd << "\t" << G_tbl.at(y).at(x);
        fd << endl;
    }
    for (int y = n_angle; y < 2 * n_angle - 1; y++) {
        fd << y * step_angle;
        for (int x = 0; x < fsize; x++)
            fd << "\t" << G_tbl.at(2 * n_angle - y - 2).at(x);
        fd << endl;
    }
    fd.close();
}