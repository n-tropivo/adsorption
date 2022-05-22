import numpy as np
from pint import UnitRegistry
from functions import *

class IAS():

    def __init__(self, funcs):
        mod = __import__("functions")
        # for f in funcs:
        #     print(f["kw"])
        # self.funcs = getattr(mod, np.array([f["fname"] for f in funcs]))(np.array([**(f["kw"]) for f in funcs]))

    def simple_IAS(self, total_pressure, y):
        z = newton(lambda z: total_pressure * sum([y[i] / self.funcs[i].pure_comp_pressure(z) for i in range(len(y))]) - 1 ,\
            sum([y[i] * self.funcs[i].red_spr_pressure(total_pressure) for i in range(len(y))]), \
            fprime = lambda z : - total_pressure * sum([y[i] / ((p0 := self.funcs[i].pure_comp_pressure(z)) * self.funcs[i].surface_conc(p0)) for i in range(len(y))])
        )

        return [(total_pressure * y[i] / self.funcs[i].pure_comp_pressure(z)).to('').magnitude for i in range(len(y))]

    def fast_IAS(self, total_pressure, y, max_iter = 100, eps = 1e-7):
        print(getattr(self.funcs, "henry_law_const"))
        s = np.sum(self.funcs.henry_law_const * y)
        x = self.funcs.henry_law_const * y / s
        print(s)
        print(x)

        it = 0
        while(True):
            p0 = [total_pressure * y[i] / x[i] for i in range(len(y))]
            z0 = self.funcs[0].red_spr_pressure(p0[0])
            g = [sum(x) - 1] + \
                [z0 - self.funcs[i].red_spr_pressure(p0[i]) for i in range(1, len(y))]

            j = [self.funcs[0].surface_conc(p0[0]) / x[0]] + \
                [x[i] / self.funcs[i].surface_conc(p0[i]) for i in range(1, len(y))]

            s = 0
            for i in range(1, len(y)):
                s += j[i]
                g[0] -= j[i] * g[i]
            
            g[0] /= 1 + j[0] * s
            g[1:] = [(g[i] + g[0] * j[0]) * j[i] for i in range(1, len(y))]

            x = [k if (k := x[i] - g[i]).magnitude > 0 else eps for i in range(len(y))]

            e = sqrt(sum([((g[i] / x[i]).to("").magnitude)**2 for i in range(len(y))]))
            it += 1
            if (e < eps):
                return x
            if (it > max_iter):
                raise TimeoutError("Number of iterations exceeded the maximum.")
            

if __name__ == "__main__":

    ureg = UnitRegistry()
    total_pressure = 7 * ureg.Pa
    y = [0.37, 0.13, 0.01, 0.49] * ureg.dimensionless
    cmus = [32.47, 82.34, 584.43, 30.15] * ureg.mmol / ureg.gram
    b = [0.255, 2.7682, 97.7962, 23.703] * (1 / ureg.MPa)
    t = [0.777, 0.323, 0.134, 0.660] * ureg.dimensionless

    to = Toth(cmus[0], b[0], t[0])
    la = Langmuir(cmus[0], b[0])
    si = Sips(cmus[0], b[0], t[0])
    un = Unilan(cmus[0], b[0], t[0])

    p = 10 * ureg.atm
    f = un
    print(f.surface_conc(p).to("mol/kg"))
    print(f.henry_law(p).to("mol/kg"))
    print(f.red_spr_pressure(p))
    print(f.pure_comp_pressure(10 * ureg.mol/ureg.kg).to("pascal"))

    # dat = [{"fname" : "Toth", "kw" : {"sat_surface_conc" : cmus[i], "b" : b[i], "t" : t[i]}} for i in range(len(y))]

    # ias = IAS(dat)

    # print([p.to('').magnitude for p in ias.fast_IAS(total_pressure, y)])



