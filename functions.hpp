#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include <boost/units/pow.hpp>

#include <boost/variant.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "units.hpp"

using namespace std;
using namespace boost::units;


template<class Conc>
class AdsorptionFunction {
public:
    virtual Conc surface_concentration(const pressure& p){
        return 0;
    };

    virtual Conc henry_law_concentration(const pressure& p){
        return 0;
    };

    virtual Conc reduced_spreading_pressure(const pressure& p){
        return 0;
    };

    virtual pressure pure_component_pressure(const Conc& z){
        return 0;
    };
};

template<class Conc>
using AFunique_ptr = unique_ptr<AdsorptionFunction<Conc>>;

template<class Conc>
using AFvector = vector<AFunique_ptr<Conc>>;

template<class Conc>
class Langmuir : public AdsorptionFunction<Conc> {
public:
    Conc cs;
    inv_pressure b;

    Langmuir(Conc& _cs, inv_pressure& _b){
        cs = _cs;
        b = _b;
    }

    Conc surface_concentration(const pressure& p){
        return cs * (b*p) / (1 + b*p);
    }

    Conc henry_law_concentration(const pressure& p){
        return cs * b * p;
    }

    Conc reduced_spreading_pressure(const pressure& p){
        return cs * log(1 + b * p);
    }

    pressure pure_component_pressure(const Conc& z){
        return (exp(z / cs) - 1) / b;
    }
};

template<class Conc>
class Toth : public AdsorptionFunction<Conc> {
public:

    Conc cs;
    inv_pressure b;
    dimless t;

    Toth(const Conc& _cs, const inv_pressure& _b, const dimless& _t){
        cs = _cs;
        b = _b;
        t = _t;
    }

    Conc surface_concentration(const pressure& p){
        return cs * (b*p) * pow(1 + pow(b*p, t), -1/t);
    }

    Conc henry_law_concentration(const pressure& p){
        return cs * b * p;
    }

    Conc reduced_spreading_pressure(const pressure& p){
        double td = (double) t;
        return cs * reduced_spreading_pressure_diml((double)(b * p), &td);
    }

    pressure pure_component_pressure(const Conc& z){
        pcp_params params = {(double) t, (double) (z/cs)};
        return pure_component_pressure_diml(&params) / b;
    }

private:
    
    struct pcp_params {
        double t;
        double zcs;
    }; 

    static double rsp_deriv_diml(double x, void* t){
        double ts = *(double*) t;
        return pow(1 + pow(x, ts), -1/ts);
    }

    static double reduced_spreading_pressure_diml(double x, void* t){
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        double result, error;
        
        gsl_function F;
        F.function = &rsp_deriv_diml;
        F.params = t;

        gsl_integration_qags(&F, 0.0, x, 0, 1e-7, 1000, w, &result, &error);

        gsl_integration_workspace_free(w);
        return result;
    }
    
    static double f (double x, void* params){
        struct pcp_params *p = (struct pcp_params *) params;
        return reduced_spreading_pressure_diml(x, &(p->t)) - p->zcs;
    }
    
    static void fdf (double x, void* params, double* y, double* dy){
        struct pcp_params *p = (struct pcp_params *) params;
        *y = f(x, params);
        *dy = rsp_deriv_diml(x, &(p->t));
    }

    static double pure_component_pressure_diml(void* params){
        const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
        gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
        gsl_function_fdf FDF;

        // struct pcp_params *p = ;

        FDF.f = &f;
        FDF.df = &rsp_deriv_diml;
        FDF.fdf = &fdf;
        FDF.params = params;

        double x0, x = exp(((struct pcp_params *) params) -> zcs) - 1;
        gsl_root_fdfsolver_set(s, &FDF, x);

        int status, iter = 0, max_iter = 100;
        do{
            iter++;
            status = gsl_root_fdfsolver_iterate(s);
            x0 = x;
            x = gsl_root_fdfsolver_root(s);

            status = gsl_root_test_delta(x, x0, 0, 1e-7);            
        }while(status == GSL_CONTINUE && iter < max_iter);

        gsl_root_fdfsolver_free(s);

        return x;
    }
};

template<class Conc>
class Sips : public AdsorptionFunction<Conc> {
public:

    Conc cs;
    inv_pressure b;
    dimless n;

    Sips(const Conc& _cs, const inv_pressure& _b, const dimless& _n){
        cs = _cs;
        b = _b;
        n = _n;
    }

    Conc surface_concentration(const pressure& p){
        double temp = pow(b * p, 1/n);
        return cs * temp / (1 + temp);
    }

    Conc henry_law_concentration(const pressure& p){
        return 0;
    }

    Conc reduced_spreading_pressure(const pressure& p){
        return cs * n * log(1 + pow(b*p, 1/n));
    }

    pressure pure_component_pressure(const Conc& z){
        return pow(exp(z/(n * cs)) - 1, n) / b;
    }
};

template<class Conc>
class Unilan : public AdsorptionFunction<Conc>{
public:

    Conc cs;
    inv_pressure b;
    dimless s;

    Unilan(const Conc& _cs, const inv_pressure& _b, const dimless& _s){
        cs = _cs;
        b = _b;
        s = _s;
        exps = exp(s);
        expmins = exp(-s);
    }

    Conc surface_concentration(const pressure& p){
        return 0.5 * cs * log((1 + b*p*exps) / (1 + b*p*expmins)) / s;
    }

    Conc henry_law_concentration(const pressure& p){
        return cs * b * p * sinh(s) / s;
    }

    Conc reduced_spreading_pressure(const pressure& p){
        double sd = (double) s;
        return cs * reduced_spreading_pressure_diml((double)(b * p), &sd);
    }

    pressure pure_component_pressure(const Conc& z){
        pcp_params params = {(double) s, (double) exps, (double) expmins, (double) (z/cs)};
        return pure_component_pressure_diml(&params) / b;
    }

private:

    dimless exps, expmins;

    struct pcp_params {
        double s;
        double exps;
        double expmins;
        double zcs;
    }; 

    static double rsp_helper(double s, void* x){
        return log(1 + exp(s) * (*(double*)x));
    }
    
    static double reduced_spreading_pressure_diml(double x, void* params){
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        double result, error;

        double s = *(double*) params;
        
        gsl_function F;
        F.function = &rsp_helper;
        F.params = &x;

        gsl_integration_qags(&F, -s, s, 0, 1e-7, 1000, w, &result, &error);

        gsl_integration_workspace_free(w);
        return 0.5 * result / s;
    }
    
    static double f (double x, void* params){
        struct pcp_params *p = (struct pcp_params *) params;
        return reduced_spreading_pressure_diml(x, &(p->s)) - p->zcs;
    }

    static double df (double x, void* params){
        struct pcp_params *p = (struct pcp_params *) params;
        return log((1 + x * p->exps) / (1 + x * p->expmins)) / (2 * x * p->s);
    }
    
    static void fdf (double x, void* params, double* y, double* dy){
        *y = f(x, params);
        *dy = df(x, params);
    }

    static double pure_component_pressure_diml(void* params){
        const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
        gsl_root_fdfsolver *solv = gsl_root_fdfsolver_alloc(T);
        gsl_function_fdf FDF;

        FDF.f = &f;
        FDF.df = &df;
        FDF.fdf = &fdf;
        FDF.params = params;

        double x0, x = exp(((struct pcp_params *) params) -> zcs) - 1;
        gsl_root_fdfsolver_set(solv, &FDF, x);

        int status, iter = 0, max_iter = 100;
        do{
            iter++;
            status = gsl_root_fdfsolver_iterate(solv);
            x0 = x;
            x = gsl_root_fdfsolver_root(solv);

            status = gsl_root_test_delta(x, x0, 0, 1e-7);            
        }while(status == GSL_CONTINUE && iter < max_iter);

        gsl_root_fdfsolver_free(solv);

        return x;
    }
};