#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>

#include <hdf5.h>

#include "include/Declare_functions.hh"
#include "include/Check_parameters.hh"

#include "Parameters.hh"

using namespace std;


int main() {
    // Print some info
    #if (INITIAL_DATA == GAUSSIAN)
    cout << "Gaussian initial data"     << endl
         << "Initial center: " << X0    << endl
         << "Initial sigma:  " << SIGMA << endl
         << endl;
    #else
    #error "Invalid initial data"
    #endif


    // Build the initial data
    constexpr double dx = L/static_cast<double>(NX);
    vector<double> x(NX), u(NX);

    for (auto i = decltype(NX){0}; i < NX; ++i) {
        x[i] = i*dx;  // NOTE: only used for output

        #if (INITIAL_DATA == GAUSSIAN)
        u[i] = gaussian(x[i]);
        #else
        #error "Invalid initial data"
        #endif

        if (u[i] < 0.0) {
            cout << "ERROR: u = " << u[i] << " < 0 at x = " << x[i]
                 << ". The porous medium equation only admits non-negative solutions; please provide non-negative initial data." << endl;
            return 1;
        }
    }


    /* ------------------------
     * Initialize the HDF5 file
     * ------------------------ */
    const auto file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    constexpr hsize_t fdims = NX;
    const auto fspace_id    = H5Screate_simple(1, &fdims, nullptr);
    assert(fspace_id >= 0);


    // ***** Output frequency ****
    constexpr hsize_t one   = 1;
    const auto space_one_id = H5Screate_simple(1, &one, nullptr);
    assert(space_one_id >= 0);

    const auto out_freq_dset_id = H5Dcreate(file_id, "Output frequency", H5T_NATIVE_INT, space_one_id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(out_freq_dset_id >= 0);
    assert(H5Dwrite(out_freq_dset_id, H5T_NATIVE_INT, space_one_id, space_one_id, H5P_DEFAULT, &OUT_EVERY) >= 0);
    assert(H5Dclose(out_freq_dset_id) >= 0);
    assert(H5Sclose(space_one_id) >= 0);


    // ***** Position *****
    const auto x_dset_id = H5Dcreate(file_id, "Position",  H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(x_dset_id >= 0);
    assert(H5Dwrite(x_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, x.data()) >= 0);
    assert(H5Dclose(x_dset_id) >= 0);


    // ***** Solution *****
    const auto u_group_id = H5Gcreate(file_id, "Solution",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u_group_id >= 0); 

    const auto u0_dset_id = H5Dcreate(u_group_id, "Iteration 0",
                                      H5T_NATIVE_DOUBLE, fspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u0_dset_id >= 0);
    assert(H5Dwrite(u0_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, u.data()) >= 0);
    assert(H5Dclose(u0_dset_id) >= 0);



    /* ====================
     * Begin time evolution
     * ==================== */
    constexpr double nxm  = NX - 1;
    constexpr double cnst = ALPHA*(GAMMA + 1.0)*DT/(2.0*dx*dx);

    vector<double> neg_res(NX);                                // Negative residual
    vector<double> dri_dui(NX), dri_duim(nxm), dri_duip(nxm);  // Residual's Jacobian dR^i / du^j (tridiagonal)
    vector<double> uold(NX), du(NX);                           // Iterative correction to solution

    constexpr hsize_t tdims = NT/OUT_EVERY;  // NOTE: integer division
    vector<double> time(tdims);

    // n == 0 has been written to file already
    for (auto n = decltype(NT){1}; n < NT; ++n) {
        size_t nrcount = 0;

        uold = u;

        /* -------------------
         * Newton-Raphson loop
         * ------------------- */
        do {
            if (nrcount > NRMAXIT) {
                cout << "ERROR (iteration " << n << ", time " << n*DT
                     << "): number of Newton-Raphson iterations exceeded maximum limit "
                     << NRMAXIT << endl;
                return 1;
            }

            // Build the negative residual and the residuals's Jacobian
            for (auto i = decltype(NX){0}; i < NX; ++i) {
                const auto ui  = u[i];
                const auto uim = (i > 0)   ? u[i-1] : U0;
                const auto uip = (i < nxm) ? u[i+1] : UL;

                if (ui < 0.0 or uim < 0.0 or uip < 0.0) {
                    cout << "ERROR (time iteration " << n << ", time " << n*DT << ", Newton-Raphon iteration " << nrcount
                         << "): u[i] = " << ui << ", u[i-1] = " << uim << ", u[i+1] = " << uip << " at x = " << x[i]
                         << ". The porous medium equation only admits non-negative solutions."
                         << " The solution might have gone negative after Newton-Raphson overshooting due to an overly large time step."
                         << " Try reducing the time step parameter DT (currently DT = " << DT << ")."
                         << endl;
                    return 1;
                }

                const auto uig  = pow(ui,  GAMMA);
                const auto uimg = pow(uim, GAMMA);
                const auto uipg = pow(uip, GAMMA);

                const auto uimg_plus_uig = uimg + uig;
                const auto uipg_plus_uig = uipg + uig;

                const auto uim_minus_ui = uim - ui;
                const auto uip_minus_ui = uip - ui;

                neg_res[i] = uold[i] - ui + cnst*(uipg_plus_uig*uip_minus_ui + uimg_plus_uig*uim_minus_ui);

                const auto uigm = uig/ui;  // pow(u[i], GAMMA - 1.0)
                const auto uigp = uig*ui;  // pow(u[i], GAMMA + 1.0)

                const auto uimgm = uimg/ui;  // pow(u[i-1], GAMMA - 1.0)
                const auto uipgm = uipg/ui;  // pow(u[i+1], GAMMA - 1.0)

                dri_dui[i] = 1.0 - cnst*(GAMMA*uigm*(uim - 2.0*ui + uip) - (uimg_plus_uig + uipg));

                if (i < nxm) {  // dr^i/du^(i-1) undefined at i==NX-1
                    dri_duim[i] = -cnst*(GAMMA*uimgm*uim_minus_ui + uimg_plus_uig);
                }

                if (i > 0) {  // dr^i/du^(i+1) undefined at i==0
                    dri_duip[i-1] = -cnst*(GAMMA*uipgm*uip_minus_ui + uipg_plus_uig);
                }
            }

            // Solve J*du = -res where J = dr^i/du^j (linearization of res = 0)
            Thomas_tridiagonal_solver(dri_duim, dri_dui, dri_duip, neg_res, du);

            // Update the solution
            for (auto i = decltype(NX){0}; i < NX; ++i) {
                u[i] += du[i];
            }

            ++nrcount;
        } while (
            // Stop when |max_i(res[i])| <= NRTOL
            abs(*(max_element(neg_res.begin(), neg_res.end(),
                [](const double &a, const double &b) {
                    return (abs(a) < abs(b));
                }))) > NRTOL
        );


        /* ------------------
         * Write data to file
         * ------------------ */
        static_assert(sizeof(std::size_t) <= sizeof(long long),
                      "size_t too big for long long");
        const auto [quot, rem] = lldiv(static_cast<long long>(n), static_cast<long long>(OUT_EVERY));
        if (rem == 0) {
            time[quot] = n*DT;

            ostringstream iteration_ss;
            iteration_ss << "Iteration " << n;

            const auto un_dset_id = H5Dcreate(u_group_id, iteration_ss.str().c_str(),
                                              H5T_NATIVE_DOUBLE, fspace_id,
                                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert(un_dset_id >= 0);
            assert(H5Dwrite(un_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, u.data()) >= 0);
            assert(H5Dclose(un_dset_id) >= 0);

            if constexpr (VERBOSE) {
                cout << "Iteration " << n << ": data written to file" << endl;
            }
        }

        if constexpr (VERBOSE) {
            cout << "Iteration " << n << ": done" << endl;
        }
    }


    // Write time to file 
    const auto tspace_id = H5Screate_simple(1, &tdims, nullptr);
    assert(tspace_id >= 0);

    const auto t_dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, tspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(t_dset_id >= 0);

    assert(H5Dwrite(t_dset_id, H5T_NATIVE_DOUBLE, tspace_id, tspace_id, H5P_DEFAULT, time.data()) >= 0);

    assert(H5Dclose(t_dset_id) >= 0);
    assert(H5Sclose(tspace_id) >= 0);
 

    // Close up
    assert(H5Sclose(fspace_id)  >= 0);
    assert(H5Gclose(u_group_id) >= 0);
    assert(H5Fclose(file_id)    >= 0);

    cout << "Data written to '" << FILENAME << "'" << endl;

    return 0;
}
