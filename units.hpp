#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/io.hpp>

#define kilo 1e3
#define mega 1e6
#define giga 1e9


using namespace boost::units;

typedef quantity<si::dimensionless> dimless;

#include <boost/units/base_units/metric/atmosphere.hpp>
BOOST_UNITS_STATIC_CONSTANT(pascal, si::pressure);
typedef quantity<si::pressure> pressure;

BOOST_UNITS_STATIC_CONSTANT(atm, metric::atmosphere_base_unit::unit_type);
// typedef quantity<metric::atmosphere_base_unit::unit_type, double> pressure_atm;

#include <boost/units/base_units/metric/bar.hpp>
BOOST_UNITS_STATIC_CONSTANT(bar, metric::bar_base_unit::unit_type);
// typedef quantity<metric::bar_base_unit::unit_type, double> pressure_bar;


// BOOST_UNITS_DEFINE_BASE_UNIT_WITH_CONVERSIONS(metric, inv_pascal, "pascal^-1", "atm", 1.01325e5, si::pressure, 33);
typedef unit<derived_dimension<time_base_dimension, 2, length_base_dimension, 1, mass_base_dimension, -1>::type, si::system> inv_pressure_unit;
typedef quantity<inv_pressure_unit, double> inv_pressure;
BOOST_UNITS_STATIC_CONSTANT(inv_pascal, inv_pressure_unit);

#include <boost/units/base_units/si/mole.hpp>
BOOST_UNITS_STATIC_CONSTANT(mol, si::mole_base_unit::unit_type);
typedef quantity<si::mole_base_unit::unit_type, double> quantity_mol;

#include <boost/units/base_units/cgs/gram.hpp>
BOOST_UNITS_STATIC_CONSTANT(gram, cgs::gram_base_unit::unit_type);
typedef quantity<cgs::gram_base_unit::unit_type, double> mass_gram;

typedef unit<derived_dimension<amount_base_dimension, 1, mass_base_dimension, -1>::type, si::system>::unit_type conc_molpkg_unit;
typedef quantity<conc_molpkg_unit, double> conc_molpkg;
BOOST_UNITS_STATIC_CONSTANT(molpkg, conc_molpkg_unit);

typedef unit<derived_dimension<amount_base_dimension, -1, mass_base_dimension, 1>::type, si::system>::unit_type inv_conc_molpkg_unit;
typedef quantity<inv_conc_molpkg_unit, double> inv_conc_molpkg;
BOOST_UNITS_STATIC_CONSTANT(kgpmol, inv_conc_molpkg_unit);

typedef unit<derived_dimension<amount_base_dimension, 1, length_base_dimension, -3>::type, si::system>::unit_type conc_molpm3_unit;
typedef quantity<conc_molpm3_unit, double> conc_molpm3;
BOOST_UNITS_STATIC_CONSTANT(molpm3, conc_molpm3_unit);

typedef unit<derived_dimension<amount_base_dimension, -1, length_base_dimension, 3>::type, si::system>::unit_type inv_conc_molpm3_unit;
typedef quantity<inv_conc_molpm3_unit, double> inv_conc_molpm3;
BOOST_UNITS_STATIC_CONSTANT(m3pmol, inv_conc_molpm3_unit);


// typedef make_scaled_unit<pressure, scale<10, static_rational<3> > >::type kilopascal;