#ifndef AVFORM_ATMOSPHERE_HPP
#define AVFORM_ATMOSPHERE_HPP

#include "core.hpp"
#include "units.hpp"

namespace avform::atmosphere {

// ============================================================================
// 1976 International Standard Atmosphere (ISA) Constants (#Altimetry)
// ============================================================================

namespace isa_constants {

template <typename T = double>
inline constexpr T t0_celsius = static_cast<T>(15.0);

template <typename T = double>
inline constexpr T t0_kelvin = static_cast<T>(288.15);

template <typename T = double>
inline constexpr T t0_fahrenheit = static_cast<T>(59.0);

template <typename T = double>
inline constexpr T p0_inhg = static_cast<T>(29.92126);

template <typename T = double>
inline constexpr T p0_mbar = static_cast<T>(1013.250);

template <typename T = double>
inline constexpr T p0_psf = static_cast<T>(2116.2166);

template <typename T = double>
inline constexpr T p0_mmhg = static_cast<T>(760.0);

template <typename T = double>
inline constexpr T p0_pa = static_cast<T>(101325.0);

template <typename T = double>
inline constexpr T p0_psi = static_cast<T>(14.69595);

template <typename T = double>
inline constexpr T p0_atm = static_cast<T>(1.0);

template <typename T = double>
inline constexpr T rho0_kg_m3 = static_cast<T>(1.2250);

template <typename T = double>
inline constexpr T rho0_slugs_ft3 = static_cast<T>(0.002376892);

template <typename T = double>
inline constexpr T lapse_rate_c_per_m = static_cast<T>(0.0065);

template <typename T = double>
inline constexpr T lapse_rate_c_per_ft = static_cast<T>(0.0019812);

template <typename T = double>
inline constexpr T tropopause_km = static_cast<T>(11.0);

template <typename T = double>
inline constexpr T tropopause_ft = static_cast<T>(36089.24);

template <typename T = double>
inline constexpr T t_tropopause_celsius = static_cast<T>(-56.5);

template <typename T = double>
inline constexpr T t_tropopause_kelvin = static_cast<T>(216.65);

template <typename T = double>
inline constexpr T cs0_knots = static_cast<T>(661.4786);

} // namespace isa_constants

// ============================================================================
// Standard Atmosphere Profile (#Altimetry)
// ============================================================================

/**
 * @brief Standard temperature (in Celsius) at altitude `h_ft` (feet), valid up to 20 km (~65,617 ft).
 * Formula:
 *   T_s = 15 - 0.0019812 * h (h < 36089.24 ft)
 *   T_s = -56.5              (h >= 36089.24 ft)
 */
template <typename T = double>
[[nodiscard]] constexpr T isa_temperature_c(T h_ft) noexcept {
    if (h_ft < isa_constants::tropopause_ft<T>) {
        return isa_constants::t0_celsius<T> - isa_constants::lapse_rate_c_per_ft<T> * h_ft;
    }
    return isa_constants::t_tropopause_celsius<T>;
}

/**
 * @brief Standard temperature (in Kelvin) at altitude `h_ft` (feet).
 */
template <typename T = double>
[[nodiscard]] constexpr T isa_temperature_k(T h_ft) noexcept {
    return units::celsius_to_kelvin(isa_temperature_c(h_ft));
}

/**
 * @brief Standard pressure at altitude `h_ft` (feet), in the same pressure unit as `p0` (default: inHg).
 * Formulas:
 *   p = P_0 * (1 - 6.8755856e-6 * h)^5.2558797       (h < 36,089.24 ft)
 *   p_Tr = 0.2233609 * P_0
 *   p = p_Tr * exp(-4.806346e-5 * (h - 36089.24))    (h >= 36,089.24 ft)
 */
template <typename T = double>
[[nodiscard]] inline T isa_pressure(T h_ft, T p0 = isa_constants::p0_inhg<T>) noexcept {
    if (h_ft <= isa_constants::tropopause_ft<T>) {
        return p0 * std::pow(static_cast<T>(1) - static_cast<T>(6.8755856e-6) * h_ft,
                             static_cast<T>(5.2558797));
    }
    const T p_tr = static_cast<T>(0.2233609) * p0;
    return p_tr * std::exp(static_cast<T>(-4.806346e-5) * (h_ft - isa_constants::tropopause_ft<T>));
}

/**
 * @brief Standard air density at altitude `h_ft` (feet), in the same density unit as `rho0` (default: kg/m^3).
 * Formulas:
 *   rho = rho_0 * (1 - 6.8755856e-6 * h)^4.2558797   (h < 36,089.24 ft)
 *   rho_Tr = 0.2970756 * rho_0
 *   rho = rho_Tr * exp(-4.806346e-5 * (h - 36089.24)) (h >= 36,089.24 ft)
 */
template <typename T = double>
[[nodiscard]] inline T isa_density(T h_ft, T rho0 = isa_constants::rho0_kg_m3<T>) noexcept {
    if (h_ft <= isa_constants::tropopause_ft<T>) {
        return rho0 * std::pow(static_cast<T>(1) - static_cast<T>(6.8755856e-6) * h_ft,
                               static_cast<T>(4.2558797));
    }
    const T rho_tr = static_cast<T>(0.2970756) * rho0;
    return rho_tr * std::exp(static_cast<T>(-4.806346e-5) * (h_ft - isa_constants::tropopause_ft<T>));
}

// ============================================================================
// Altimetry (#Altimetry)
// ============================================================================

/**
 * @brief Exact pressure altitude correction (in feet) given altimeter setting `alt_set_inhg` (inHg).
 * Formula: P_alt_corr = 145442.2 * (1 - (alt_set / 29.92126)^0.190261)
 */
template <typename T = double>
[[nodiscard]] inline T pressure_altitude_correction_ft(T alt_set_inhg) noexcept {
    return static_cast<T>(145442.2) *
           (static_cast<T>(1) - std::pow(alt_set_inhg / isa_constants::p0_inhg<T>, static_cast<T>(0.190261)));
}

/**
 * @brief Simple approximation for pressure altitude correction (in feet):
 * Formula: P_alt_corr = (29.92 - alt_set) * 1000
 */
template <typename T = double>
[[nodiscard]] constexpr T pressure_altitude_correction_approx_ft(T alt_set_inhg) noexcept {
    return (static_cast<T>(29.92) - alt_set_inhg) * static_cast<T>(1000);
}

/**
 * @brief Pressure altitude (in feet) from indicated altitude `ind_alt_ft` and altimeter setting `alt_set_inhg`.
 */
template <typename T = double>
[[nodiscard]] inline T pressure_altitude_ft(T ind_alt_ft, T alt_set_inhg, bool approximate = false) noexcept {
    if (approximate) {
        return ind_alt_ft + pressure_altitude_correction_approx_ft(alt_set_inhg);
    }
    return ind_alt_ft + pressure_altitude_correction_ft(alt_set_inhg);
}

/**
 * @brief Exact density altitude (in feet) given pressure altitude `p_alt_ft` (feet)
 * and actual outside air temperature `oat_c` (Celsius).
 * Formula: D_Alt = P_alt + (T_s / T_r) * (1 - (T_s / T)^0.2349690)
 */
template <typename T = double>
[[nodiscard]] inline T density_altitude_ft(T p_alt_ft, T oat_c) noexcept {
    const T ts_k = isa_temperature_k(p_alt_ft);
    const T t_k = units::celsius_to_kelvin(oat_c);
    return p_alt_ft + (ts_k / isa_constants::lapse_rate_c_per_ft<T>) *
                      (static_cast<T>(1) - std::pow(ts_k / t_k, static_cast<T>(0.2349690)));
}

/**
 * @brief Approximate density altitude (in feet):
 * Formula: D_Alt = P_Alt + 118.6 * (T - T_s)
 */
template <typename T = double>
[[nodiscard]] constexpr T density_altitude_approx_ft(T p_alt_ft, T oat_c) noexcept {
    const T ts_c = isa_temperature_c(p_alt_ft);
    return p_alt_ft + static_cast<T>(118.6) * (oat_c - ts_c);
}

/**
 * @brief True altitude above sea-level (in feet) from calibrated altitude `ca_ft`,
 * station field elevation `fe_ft`, average ISA temperature deviation `isadev_c` (Celsius),
 * and outside air temperature at altitude `oat_c` (Celsius).
 * Formula: TA = CA + (CA - FE) * (ISADEV) / (273 + OAT)
 */
template <typename T = double>
[[nodiscard]] constexpr T true_altitude_ft(T ca_ft, T fe_ft, T isadev_c, T oat_c) noexcept {
    return ca_ft + (ca_ft - fe_ft) * isadev_c / (static_cast<T>(273.0) + oat_c);
}

// ============================================================================
// Mach Numbers, True vs Calibrated Airspeeds (#Mach)
// ============================================================================

/**
 * @brief Speed of sound (in knots) given static Outside Air Temperature `oat_c` (Celsius).
 * Formula: CS = 38.967854 * sqrt(T + 273.15)
 */
template <typename T = double>
[[nodiscard]] inline T speed_of_sound_kt(T oat_c) noexcept {
    return static_cast<T>(38.967854) * std::sqrt(oat_c + static_cast<T>(273.15));
}

/**
 * @brief Mach number M = TAS / CS given TAS (knots) and OAT (Celsius).
 */
template <typename T = double>
[[nodiscard]] inline T mach_number(T tas_kt, T oat_c) noexcept {
    return tas_kt / speed_of_sound_kt(oat_c);
}

/**
 * @brief Approximate Indicated Air Temperature (IAT, Celsius) from OAT (Celsius), TAS (knots),
 * and probe recovery factor `k` (typically 0.95 to 1.0, down to 0.7).
 * Formula: IAT = OAT + K * TAS^2 / 7592
 */
template <typename T = double>
[[nodiscard]] constexpr T iat_from_oat_tas(T oat_c, T tas_kt, T k = static_cast<T>(1.0)) noexcept {
    return oat_c + k * tas_kt * tas_kt / static_cast<T>(7592.0);
}

/**
 * @brief Approximate OAT (Celsius) from IAT (Celsius), TAS (knots), and recovery factor `k`.
 */
template <typename T = double>
[[nodiscard]] constexpr T oat_from_iat_tas(T iat_c, T tas_kt, T k = static_cast<T>(1.0)) noexcept {
    return iat_c - k * tas_kt * tas_kt / static_cast<T>(7592.0);
}

/**
 * @brief True OAT (Celsius) from Indicated Air Temperature `iat_c` (Celsius), Mach number `mach`,
 * and probe recovery factor `k`.
 * Formula: OAT = (IAT + 273.15) / (1 + 0.2 * K * M^2) - 273.15
 */
template <typename T = double>
[[nodiscard]] constexpr T oat_from_iat_mach(T iat_c, T mach, T k = static_cast<T>(1.0)) noexcept {
    return (iat_c + static_cast<T>(273.15)) /
           (static_cast<T>(1.0) + static_cast<T>(0.2) * k * mach * mach) -
           static_cast<T>(273.15);
}

/**
 * @brief Low-speed (incompressible, M < 0.3) True Airspeed from Calibrated Airspeed `cas_kt`
 * and density altitude `da_ft` (< 36,089.24 ft).
 * Formula: TAS = CAS / (1 - 6.8755856e-6 * DA)^2.127940
 */
template <typename T = double>
[[nodiscard]] inline T tas_low_speed_kt(T cas_kt, T da_ft) noexcept {
    return cas_kt / std::pow(static_cast<T>(1) - static_cast<T>(6.8755856e-6) * da_ft,
                             static_cast<T>(2.127940));
}

template <typename T = double>
struct CompressibleAirspeedResult {
    T dp_inhg{}; ///< Differential pitot-static pressure (inHg)
    T p_inhg{};  ///< Static pressure at pressure altitude (inHg)
    T mach{};    ///< True Mach number (subsonic or Rayleigh supersonic)
    T oat_c{};   ///< Static outside air temperature (Celsius)
    T cs_kt{};   ///< Speed of sound at altitude (knots)
    T tas_kt{};  ///< True airspeed (knots)
};

/**
 * @brief Computes compressible Mach number, OAT, speed of sound, and True Airspeed (TAS)
 * given CAS/IAS (knots), pressure altitude `pa_ft` (feet), indicated air temperature `iat_c` (Celsius),
 * and recovery factor `k`.
 * Automatically applies Rayleigh's Supersonic Pitot iterative solver when M > 1.
 */
template <typename T = double>
[[nodiscard]] inline CompressibleAirspeedResult<T>
tas_and_mach_from_cas(T cas_kt, T pa_ft, T iat_c, T k = static_cast<T>(1.0),
                      T supersonic_tol = static_cast<T>(1e-9), int max_iter = 100) noexcept {
    CompressibleAirspeedResult<T> res{};
    constexpr T p0 = isa_constants::p0_inhg<T>;
    constexpr T cs0 = isa_constants::cs0_knots<T>;

    const T cas_ratio = cas_kt / cs0;
    res.dp_inhg = p0 * (std::pow(static_cast<T>(1) + static_cast<T>(0.2) * cas_ratio * cas_ratio,
                                 static_cast<T>(3.5)) - static_cast<T>(1));
    res.p_inhg = isa_pressure(pa_ft, p0);

    const T dp_over_p_plus_1 = res.dp_inhg / res.p_inhg + static_cast<T>(1);
    T m = std::sqrt(static_cast<T>(5) *
                    (std::pow(dp_over_p_plus_1, static_cast<T>(2.0 / 7.0)) - static_cast<T>(1)));

    // If M > 1 (supersonic flight), iterate Rayleigh's Supersonic Pitot equation:
    // M = 0.881285 * sqrt((DP/P + 1) * (1 - 1/(7*M^2))^(5/2))
    if (m > static_cast<T>(1)) {
        for (int i = 0; i < max_iter; ++i) {
            const T term = static_cast<T>(1) - static_cast<T>(1) / (static_cast<T>(7) * m * m);
            const T m_next = static_cast<T>(0.881285) *
                             std::sqrt(dp_over_p_plus_1 * std::pow(term, static_cast<T>(2.5)));
            if (std::abs(m_next - m) < supersonic_tol) {
                m = m_next;
                break;
            }
            m = m_next;
        }
    }

    res.mach = m;
    res.oat_c = oat_from_iat_mach(iat_c, res.mach, k);
    res.cs_kt = speed_of_sound_kt(res.oat_c);
    res.tas_kt = res.mach * res.cs_kt;
    return res;
}

/**
 * @brief Computes IAS/CAS (in knots) given Mach number `mach` (for M <= 1) and pressure altitude `pa_ft` (feet).
 * Formulas:
 *   x = (1 - 6.8755856e-6 * PA)^5.2558797
 *   ias = 661.4786 * (5 * ((1 + x * ((1 + M^2/5)^3.5 - 1))^(2/7) - 1))^0.5
 */
template <typename T = double>
[[nodiscard]] inline T ias_from_mach(T mach, T pa_ft) noexcept {
    const T x = isa_pressure(pa_ft, static_cast<T>(1));
    const T inner = static_cast<T>(1) +
                    x * (std::pow(static_cast<T>(1) + (mach * mach) / static_cast<T>(5),
                                  static_cast<T>(3.5)) - static_cast<T>(1));
    return isa_constants::cs0_knots<T> *
           std::sqrt(static_cast<T>(5) *
                     (std::pow(inner, static_cast<T>(2.0 / 7.0)) - static_cast<T>(1)));
}

// ============================================================================
// Relative Humidity, Dewpoint, Frostpoint (#RH)
// ============================================================================

/**
 * @brief Relative humidity `f` (as a fraction in [0, 1]) from temperature `t_c` (Celsius)
 * and dewpoint `td_c` (Celsius).
 * Formula: f = exp(17.27 * (Td / (Td + 237.3) - T / (T + 237.3)))
 */
template <typename T = double>
[[nodiscard]] inline T relative_humidity_from_dewpoint(T t_c, T td_c) noexcept {
    return std::exp(static_cast<T>(17.27) *
                    (td_c / (td_c + static_cast<T>(237.3)) - t_c / (t_c + static_cast<T>(237.3))));
}

/**
 * @brief Relative humidity `f` (as a fraction in [0, 1]) from temperature `t_c` (Celsius)
 * and frostpoint `tf_c` (Celsius).
 * Formula: f = exp(21.87 * (Tf / (Tf + 265.5) - T / (T + 265.5)))
 */
template <typename T = double>
[[nodiscard]] inline T relative_humidity_from_frostpoint(T t_c, T tf_c) noexcept {
    return std::exp(static_cast<T>(21.87) *
                    (tf_c / (tf_c + static_cast<T>(265.5)) - t_c / (t_c + static_cast<T>(265.5))));
}

/**
 * @brief Tetens (1930) empirical fit for saturation vapor pressure `e_s` (in mbar) at temperature `t_c` (Celsius).
 * Formula: e_s = 6.11 * exp(b * T / (T + a))
 *   over water: a = 237.3, b = 17.27
 *   over ice:   a = 265.5, b = 21.87
 */
template <typename T = double>
[[nodiscard]] inline T saturation_vapor_pressure_tetens(T t_c, bool over_ice = false) noexcept {
    const T a = over_ice ? static_cast<T>(265.5) : static_cast<T>(237.3);
    const T b = over_ice ? static_cast<T>(21.87) : static_cast<T>(17.27);
    return static_cast<T>(6.11) * std::exp(b * t_c / (t_c + a));
}

/**
 * @brief Lowe (1977) polynomial fit for saturation vapor pressure `e_s` (in mbar) over water at temperature `t_c` (Celsius).
 */
template <typename T = double>
[[nodiscard]] constexpr T saturation_vapor_pressure_lowe(T t_c) noexcept {
    return static_cast<T>(6.10779) +
           t_c * (static_cast<T>(4.43652e-1) +
           t_c * (static_cast<T>(1.42894e-2) +
           t_c * (static_cast<T>(2.65064e-4) +
           t_c * (static_cast<T>(3.03124e-6) +
           t_c * (static_cast<T>(2.03408e-8) +
                  static_cast<T>(6.13682e-11) * t_c)))));
}

/**
 * @brief Dewpoint `Td` (Celsius) from temperature `t_c` (Celsius) and relative humidity fraction `f` (in (0, 1]).
 * Formula: Td = 237.3 / (1 / (ln(f)/17.27 + T/(T + 237.3)) - 1)
 */
template <typename T = double>
[[nodiscard]] inline T dewpoint_from_rh(T t_c, T f) noexcept {
    const T term = std::log(f) / static_cast<T>(17.27) + t_c / (t_c + static_cast<T>(237.3));
    return static_cast<T>(237.3) / (static_cast<T>(1) / term - static_cast<T>(1));
}

/**
 * @brief Frostpoint `Tf` (Celsius) from temperature `t_c` (Celsius) and relative humidity fraction `f` (in (0, 1]).
 * Formula: Tf = 265.5 / (1 / (ln(f)/21.87 + T/(T + 265.5)) - 1)
 */
template <typename T = double>
[[nodiscard]] inline T frostpoint_from_rh(T t_c, T f) noexcept {
    const T term = std::log(f) / static_cast<T>(21.87) + t_c / (t_c + static_cast<T>(265.5));
    return static_cast<T>(265.5) / (static_cast<T>(1) / term - static_cast<T>(1));
}

template <typename T = double>
struct PsychrometricResult {
    T ed{}; ///< Saturation vapor pressure at dry-bulb temp (mbar)
    T ew{}; ///< Saturation vapor pressure at wet-bulb temp (mbar)
    T wd{}; ///< Saturation mixing ratio at T
    T ww{}; ///< Saturation mixing ratio at Tw
    T w{};  ///< Mixing ratio
    T f{};  ///< Relative humidity as a fraction
    T e{};  ///< Vapor pressure (mbar)
    T td{}; ///< Dewpoint (Celsius)
};

/**
 * @brief Computes relative humidity, vapor pressure, mixing ratio, and dewpoint from
 * dry-bulb temperature `t_c` (°C), wet-bulb temperature `tw_c` (°C), and pressure `p_mbar` (mbar).
 */
template <typename T = double>
[[nodiscard]] inline PsychrometricResult<T> psychrometric_wet_bulb(T t_c, T tw_c, T p_mbar) noexcept {
    PsychrometricResult<T> res{};
    res.ed = static_cast<T>(6.11) * std::exp(static_cast<T>(17.27) * t_c / (t_c + static_cast<T>(237.3)));
    res.ew = static_cast<T>(6.11) * std::exp(static_cast<T>(17.27) * tw_c / (tw_c + static_cast<T>(237.3)));
    res.wd = static_cast<T>(0.62197) * res.ed / (p_mbar - res.ed);
    res.ww = static_cast<T>(0.62197) * res.ew / (p_mbar - res.ew);
    res.w = (static_cast<T>(2500.0) * res.ww - static_cast<T>(1.0046) * (t_c - tw_c)) /
            (static_cast<T>(2500.0) + static_cast<T>(1.81) * (t_c - tw_c));
    res.f = res.w / res.wd;
    res.e = p_mbar * res.w / (static_cast<T>(0.62197) + res.w);
    const T log10_e = std::log10(res.e);
    res.td = (static_cast<T>(237.3) * log10_e - static_cast<T>(186.527)) /
             (static_cast<T>(8.286) - log10_e);
    return res;
}

/**
 * @brief Increase in effective density altitude (in feet) due to humidity, given relative humidity fraction `rh_fraction`,
 * temperature `t_c` (Celsius), and pressure altitude `h_ft` (feet).
 * Formula: Increase(ft) = 0.267 * RH * (T + 273) * exp(17.3 * T / (T + 237)) * (1 - 0.00000688 * H)^(-5.26)
 */
template <typename T = double>
[[nodiscard]] inline T density_altitude_increase_from_rh(T rh_fraction, T t_c, T h_ft) noexcept {
    return static_cast<T>(0.267) * rh_fraction * (t_c + static_cast<T>(273.0)) *
           std::exp(static_cast<T>(17.3) * t_c / (t_c + static_cast<T>(237.0))) *
           std::pow(static_cast<T>(1.0) - static_cast<T>(0.00000688) * h_ft, static_cast<T>(-5.26));
}

/**
 * @brief Increase in effective density altitude (in feet) due to humidity, given dewpoint `td_c` (Celsius),
 * temperature `t_c` (Celsius), and pressure altitude `h_ft` (feet).
 * Formula: Increase(ft) = 0.267 * (T + 273) * exp(17.3 * Td / (Td + 237)) * (1 - 0.00000688 * H)^(-5.26)
 */
template <typename T = double>
[[nodiscard]] inline T density_altitude_increase_from_dewpoint(T td_c, T t_c, T h_ft) noexcept {
    return static_cast<T>(0.267) * (t_c + static_cast<T>(273.0)) *
           std::exp(static_cast<T>(17.3) * td_c / (td_c + static_cast<T>(237.0))) *
           std::pow(static_cast<T>(1.0) - static_cast<T>(0.00000688) * h_ft, static_cast<T>(-5.26));
}

} // namespace avform::atmosphere

#endif // AVFORM_ATMOSPHERE_HPP

