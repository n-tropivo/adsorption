from pint import Quantity
from numpy import log, exp, sinh
from scipy.integrate import quad
from scipy.optimize import newton

class AdsorptionFunction():
    surface_conc_func = None
    henry_law_const = None
    red_spr_pressure_func = None
    pure_comp_pressure_func = None

    def surface_conc(self, p : Quantity):
        if(p.check('[pressure]')):
            return self.surface_conc_func(p)
        else:
            raise ValueError("p must be in units of pressure")

    def henry_law(self, p : Quantity):
        if(p.check('[pressure]')):
            return self.henry_law_const * p
        else:
            raise ValueError("p must be in units of pressure")

    def red_spr_pressure(self, p0 : Quantity):
        if(p0.check('[pressure]')):
            return self.red_spr_pressure_func(p0)
        else:
            raise ValueError("p0 must be in units of pressure")

    def pure_comp_pressure(self, z : Quantity):
        if(z.check('[substance]/[mass]')):
            return self.pure_comp_pressure_func(z)
        else:
            raise ValueError("z must be in units of molar concentration")

class Langmuir(AdsorptionFunction):
    sat_surface_conc : Quantity
    b : Quantity

    def __init__(self, sat_surface_conc, b : Quantity):
        if(sat_surface_conc.check('[substance]/[mass]')):
            self.sat_surface_conc = sat_surface_conc
        else:
            raise ValueError(f"{self.__class__.__name__} function a parameter must be in units of substance over mass")

        if(b.check('[pressure]^(-1)')):
            self.b = b
        else:
            raise ValueError(f"{self.__class__.__name__} function b parameter must be in units of inverse pressure")

        self.surface_conc_func = lambda p : (self.sat_surface_conc * ((self.b * p) / (1 + self.b * p)).to('dimensionless'))
        self.henry_law_const = self.sat_surface_conc * self.b
        self.red_spr_pressure_func = lambda p0 : self.sat_surface_conc * log(1 + self.b * p0)
        self.pure_comp_pressure_func = lambda z : (exp(z/self.sat_surface_conc) - 1) / self.b

    def __str__(self):
        return f"<{self.__class__.__name__}(a = {str(self.sat_surface_conc)}, b = {str(self.b)})>"

class Unilan(AdsorptionFunction):
    sat_surface_conc : Quantity
    b : Quantity
    s : Quantity

    def __init__(self, sat_surface_conc, b, s : Quantity):
        if(sat_surface_conc.check('[substance]/[mass]')):
            self.sat_surface_conc = sat_surface_conc
        else:
            raise ValueError(f"{self.__class__.__name__} function sat_surf_conc parameter must be in units of substance over mass")

        if(b.check('[pressure]^(-1)')):
            self.b = b
        else:
            raise ValueError(f"{self.__class__.__name__} function b parameter must be in units of inverse pressure")

        if(s.check('[]')):
            self.s = s
            self.exps = exp(s)
            self.expmins = exp(-s)
        else:
            raise ValueError(f"{self.__class__.__name__} function s parameter must be dimensionless")

        self.surface_conc_func = lambda p : 0.5 * self.sat_surface_conc * log((1 + self.b * p * self.exps) / (1 + self.b * p * self.expmins)) / self.s
        self.henry_law_const = self.sat_surface_conc * self.b * sinh(self.s) / self.s
        self.red_spr_pressure_func = lambda p0 : (0.5 * self.sat_surface_conc / self.s) * quad(lambda x: log(1 + self.b * p0 * exp(x)), -self.s, self.s)[0]
        self.pure_comp_pressure_func = lambda z : newton(\
            lambda pdiml: quad(lambda x: log(1 + pdiml * exp(x)) / pdiml, -self.s, self.s)[0] - 2 * z * self.s / self.sat_surface_conc, \
            (self.s / sinh(self.s)) * (exp(z / self.sat_surface_conc) - 1),
            fprime = lambda pdiml: log((1 + pdiml * self.exps) / (1 + pdiml * self.expmins)) / pdiml 
        ) / self.b

    def __str__(self):
        return f"<{self.__class__.__name__}(sat_surf_conc = {str(self.sat_surface_conc)}, b = {str(self.b)}, s = {str(self.s)})>"

class Sips(AdsorptionFunction):
    sat_surface_conc : Quantity
    b : Quantity
    n : Quantity

    def __init__(self, sat_surface_conc, b, n : Quantity):
        if(sat_surface_conc.check('[substance]/[mass]')):
            self.sat_surface_conc = sat_surface_conc
        else:
            raise ValueError(f"{self.__class__.__name__} function sat_surf_conc parameter must be in units of substance over mass")

        if(b.check('[pressure]^(-1)')):
            self.b = b
        else:
            raise ValueError(f"{self.__class__.__name__} function b parameter must be in units of inverse pressure")

        if(n.check('[]')):
            self.n = n
            self.ninv = 1/n
        else:
            raise ValueError(f"{self.__class__.__name__} function n parameter must be dimensionless")

        self.surface_conc_func = lambda p : self.sat_surface_conc * (q := (self.b * p).to("") ** self.ninv) / (1 + q)
        self.henry_law_const = 0 
        self.red_spr_pressure_func = lambda p0 : self.n * self.sat_surface_conc * log(1 + (self.b * p0).to("") ** self.ninv)
        self.pure_comp_pressure_func = lambda z : (exp(z / (self.sat_surface_conc * self.n)) - 1)**self.n

    def __str__(self):
        return f"<{self.__class__.__name__}(sat_surf_conc = {str(self.sat_surface_conc)}, b = {str(self.b)}, n = {str(self.n)})>"

class Toth(AdsorptionFunction):
    sat_surface_conc = None
    b = None
    t = None

    def __init__(self, sat_surface_conc, b, t : Quantity):
        if(sat_surface_conc.check('[substance]/[mass]')):
            self.sat_surface_conc = sat_surface_conc
        else:
            raise ValueError(f"{self.__class__.__name__} function sat_surf_conc parameter must be in units of substance over mass")

        if(b.check('[pressure]^(-1)')):
            self.b = b
        else:
            raise ValueError(f"{self.__class__.__name__} function b parameter must be in units of inverse pressure")

        if(t.check('[]')):
            self.t = t
        else:
            raise ValueError(f"{self.__class__.__name__} function t parameter must be dimensionless")

        fdiml = lambda x: (1 + x**self.t) ** (-1/self.t)
        self.surface_conc_func = lambda p : self.sat_surface_conc * (self.b * p) / (1 + (self.b * p)**self.t)**(1/self.t)
        self.henry_law_const = self.sat_surface_conc * self.b
        self.red_spr_pressure_func = lambda p0 : self.sat_surface_conc * quad(fdiml, 0, (self.b * p0).to(""))[0]
        self.pure_comp_pressure_func = lambda z : newton(\
            lambda pdiml : quad(fdiml, 0, pdiml)[0] - z / self.sat_surface_conc, \
            exp(z / self.sat_surface_conc) - 1,
            fprime = fdiml ,\
        ) / self.b

    def __str__(self):
        return f"<{self.__class__.__name__}(sat_surf_conc = {self.sat_surface_conc.__str__()}, b = {self.b.__str__()})>"
