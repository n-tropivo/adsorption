#include <iostream>
#include <type_traits>
#include "functions.hpp"


// This is an implementation of the FastIAS method
// with respect to the adsorbed molar fractions x[i]
//
// The respective Jacobian and vector look like this:
// ( 1           1            1            ...    1           | Σx[i] - 1   )
// ( C[0]/x[0]   -C[1]/x[1]   0            ...    0           | z[0] - z[1] )
// ( C[0]/x[0]   0            -C[2]/x[2]   ...    0           | z[0] - z[2] )
// ( .           .            .                   .           | .           )
// ( .           .            .                   .           | .           )
// ( .           .            .                   .           | .           )
// ( C[0]/x[0]   0            0            ...   -C[n]/x[n]   | z[0] - z[n] )
//
// In terms of the variables defined in the following function:
// ( 1      1       1       ...    1       | Σx[i] - 1   )
// ( j0     1/j[1]  0       ...    0       | z[0] - z[1] )
// ( j0     0       1/j[2]  ...    0       | z[0] - z[2] )
// ( .      .       .              .       | .           )
// ( .      .       .              .       | .           )
// ( .      .       .              .       | .           )
// ( j0     0       0       ...    1/j[n]  | z[0] - z[n] )

template<class Conc>
void fastIAS(dimless* xout, const int array_size, AdsorptionFunction<Conc>** v, const pressure& total_pressure, const dimless* y, int max_iter = 1000, double eps = 1e-7){
    // Get inverse concentration type
    typedef typename conditional<is_same_v<Conc, conc_molpkg>, inv_conc_molpkg, inv_conc_molpm3>::type inv_conc;

    // Set initial approximation for x.
    // Each x[i] is proportional to H[i] * y[i]. To keep x as a double, divide each summand by H[0] * y[0].
    double x[array_size];
    x[0] = 1.0;
    double s = 1.0;
    pressure p = total_pressure * y[0];
    Conc c = v[0]->henry_law_concentration(p);
    for(int i = 1; i < array_size; i++){
        p = total_pressure * y[i];
        x[i] = v[i]->henry_law_concentration(p) / c;
        s += x[i];
    }
    // Normalize the x[i] so that they actually sum to zero.
    for(int i = 0; i < array_size; i++) x[i] /= s;

    int iter = 1;
    double g[array_size];           // Will contain the correction to x[n]
    Conc z0;                        // rsp for faster computation
    Conc j0;                        // jacobian[i][0] element
    inv_conc j[array_size], jsum;   // *inverse* jacobian [i][i] elements and their sum.
                                    // Unused j[0] takes up memory, but recalulating indices slows computation.
    
    for(int i = 0; i < array_size; i++) if(x[i] <= 0) x[i] = eps;  // In case Henry's law has coefficient zero

    // Start iterating
    while(true){
        iter++;
        // Diagonalize the jacobian
        p = total_pressure * y[0] / x[0];
        z0 = v[0]->reduced_spreading_pressure(p);
        j0 = v[0]->surface_concentration(p) / x[0];
        jsum = 0;
        g[0] = x[0] - 1;
        // Gauss forward propagation
        for(int i = 1; i < array_size; i++){
            p = total_pressure * y[i] / x[i];
            j[i] = x[i] / v[i]->surface_concentration(p);
            jsum += j[i];

            g[i] = (z0 - v[i]->reduced_spreading_pressure(p)) * j[i];
            g[0] += x[i] - g[i];
        }

        g[0] /= 1 + (double) (j0 * jsum);
        j0 *= g[0];

        // Gauss back propagation
        for(int i = 1; i < array_size; i++) g[i] += j0 * j[i];

        // Correct x[i], find relative error
        s = 0;
        for(int i = 0; i < array_size; i++){
            x[i] -= g[i];
            if(x[i] < 0) x[i] = eps;
            s += pow(x[i]/g[i], 2);
        }

        if(sqrt(s) < eps) break;
        if(iter > max_iter) break;
    }

    for(int i = 0; i < array_size; i++) xout[i] = x[i];
}

template<class Conc>
void test(AdsorptionFunction<Conc>** v, const dimless* y){
    typename Conc::unit_type conc_unit;
    typedef typename conditional<is_same_v<Conc, conc_molpkg>, inv_conc_molpkg, inv_conc_molpm3>::type inv_conc;
    Conc a = 10.0 * conc_unit;
    inv_conc f = 10.0 / conc_unit;
    pressure p = (pressure) (1 * atm);
    // cout << v[0]->surface_concentration(p) << " " << y[0] << endl;
    cout << a * f << endl;
}

int main(int argc, char* argv[]){

    pressure total_pressure = static_cast<pressure>(7 * pascal);

    dimless y[4] = {0.37, 0.13, 0.01, 0.49};
    conc_molpkg cs[4] = {32.47 * molpkg, 82.34 * molpkg, 584.43 * molpkg, 30.15 * molpkg};
    inv_pressure b[4] = {0.255 / (mega * pascal), 2.7682 / (mega * pascal), 97.7962 / (mega * pascal), 23.703 / (mega * pascal)};
    dimless t[4] = {0.777, 0.323, 0.134, 0.600};

    AdsorptionFunction<conc_molpkg>** v = new AdsorptionFunction<conc_molpkg>* [4];

    v[0] = new Langmuir<conc_molpkg>(cs[0], b[0]);
    v[1] = new Toth<conc_molpkg>(cs[1], b[1], t[1]);
    v[2] = new Sips<conc_molpkg>(cs[2], b[2], t[2]);
    v[3] = new Unilan<conc_molpkg>(cs[3], b[3], t[3]);

    dimless x[4]; 
    fastIAS<conc_molpkg>(x, 4, v, total_pressure, y);

    for(int i = 0; i < 4; i++){
        cout << x[i] << " ";
    }
    cout << endl;

    // pressure p = (pressure)(1e1 * atm);

    // conc_molpkg z = 10 * molpkg;
    
    // for(int i = 0; i < v.size(); i++){
    //     cout << v[i]->surface_concentration(p) << " "
    //         << v[i]->henry_law_concentration(p) 
    //         << " " << v[i]->reduced_spreading_pressure(p) 
    //         << " " << v[i]->pure_component_pressure(z)
    //         << endl << endl;
    // }

    // test<conc_molpkg>(v, y);



    // cout << name_format << engineering_prefix << static_cast<double>(langmuir(cs[0], p, b[0]) / cs[0]) << endl;
}