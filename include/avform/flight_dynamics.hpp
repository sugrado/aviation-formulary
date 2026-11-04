#ifndef AVFORM_FLIGHT_DYNAMICS_HPP
#define AVFORM_FLIGHT_DYNAMICS_HPP

#include "core.hpp"
#include "units.hpp"

namespace avform::flight_dynamics {

// ============================================================================
// Turns and Pivotal Altitude (#Turns)
// ============================================================================

/**
 * @brief Radius of a steady, no-wind turn in consistent physical units.
 * Formula: R = v^2 / (g * tan(b))
 */
template <typename T = double>
[[nodiscard]] inline T turn_radius_physical(T v, T bank_rad, T g) noexcept {
    return (v * v) / (g * std::tan(bank_rad));
}

/**
 * @brief Rate of a steady, no-wind turn (in rad/s) in consistent physical units.
 * Formula: w = v / R
 */
template <typename T = double>
[[nodiscard]] constexpr T turn_rate_physical(T v, T radius) noexcept {
    return v / radius;
}

/**
 * @brief Pivotal altitude in consistent physical units.
 * Formula: h_p = v^2 / g
 */
template <typename T = double>
[[nodiscard]] constexpr T pivotal_altitude_physical(T v, T g) noexcept {
    return (v * v) / g;
}

/**
 * @brief Radius of turn (in feet) given airspeed `v_kt` (knots) and bank angle `bank_deg` (degrees).
 * Formula: R = v^2 / (11.23 * tan(0.01745 * b))
 * (If `use_exact_deg_to_rad` is true, uses pi/180 instead of the rounded 0.01745 constant).
 */
template <typename T = double>
[[nodiscard]] inline T turn_radius_ft(T v_kt, T bank_deg, bool use_exact_deg_to_rad = false) noexcept {
    const T bank_rad = use_exact_deg_to_rad ? units::deg_to_rad(bank_deg)
                                            : static_cast<T>(0.01745) * bank_deg;
    return (v_kt * v_kt) / (static_cast<T>(11.23) * std::tan(bank_rad));
}

/**
 * @brief Rate of turn `w` (in degrees/sec) given airspeed `v_kt` (knots) and turn radius `radius_ft` (feet).
 * Formula: w = 96.7 * v / R
 */
template <typename T = double>
[[nodiscard]] constexpr T turn_rate_deg_per_sec(T v_kt, T radius_ft) noexcept {
    return static_cast<T>(96.7) * v_kt / radius_ft;
}

/**
 * @brief Bank angle `b_s` (in degrees) for a standard-rate turn (3 deg/sec) given airspeed `v_kt` (knots).
 * Formula: b_s = 57.3 * atan(v / 362.1)
 */
template <typename T = double>
[[nodiscard]] inline T standard_rate_bank_angle_deg(T v_kt) noexcept {
    return static_cast<T>(57.3) * std::atan(v_kt / static_cast<T>(362.1));
}

/**
 * @brief Rule-of-thumb bank angle (in degrees) for a standard-rate turn (accurate to ~1 deg for v < 250 kt):
 * Formula: b_s = v / 7
 */
template <typename T = double>
[[nodiscard]] constexpr T standard_rate_bank_angle_rule_of_thumb_deg(T v_kt) noexcept {
    return v_kt / static_cast<T>(7.0);
}

/**
 * @brief Pivotal altitude `h_p` (in feet) given groundspeed `v_kt` (knots).
 * Formula: h_p = v^2 / 11.23
 */
template <typename T = double>
[[nodiscard]] constexpr T pivotal_altitude_ft(T v_kt) noexcept {
    return (v_kt * v_kt) / static_cast<T>(11.23);
}

// ============================================================================
// Distance to Horizon (#Horizon)
// ============================================================================

/**
 * @brief Distance to the horizon in consistent length units given height `h`, Earth radius `R`,
 * and atmospheric refraction factor `b` (default 0.8279 for standard lapse rate).
 * Formula: d = sqrt(2 * R * h / b)
 */
template <typename T = double>
[[nodiscard]] inline T distance_to_horizon_physical(T h, T radius,
                                                    T refraction_b = static_cast<T>(0.8279)) noexcept {
    return std::sqrt(static_cast<T>(2) * radius * h / refraction_b);
}

/**
 * @brief Practical distance to the horizon (in nautical miles) given height `h_ft` (in feet).
 * Formula: d = 1.17 * sqrt(h)
 */
template <typename T = double>
[[nodiscard]] inline T distance_to_horizon_nm(T h_ft) noexcept {
    return static_cast<T>(1.17) * std::sqrt(h_ft);
}

// ============================================================================
// Approximate Magnetic Variation Formulae (#Var)
// ============================================================================

/**
 * @brief Least-squares polynomial fit for magnetic variation (in degrees, West positive)
 * for the Continental US (24 < lat_n_deg < 50, 66 < lon_w_deg < 125).
 */
template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_conus_deg(T lat_n_deg, T lon_w_deg) noexcept {
    const T x = lat_n_deg;
    const T y = lon_w_deg;
    const T x2 = x * x;
    const T x3 = x2 * x;
    const T y2 = y * y;
    const T y3 = y2 * y;
    return static_cast<T>(-65.6811) +
           static_cast<T>(0.99) * x +
           static_cast<T>(0.0128899) * x2 -
           static_cast<T>(0.0000905928) * x3 +
           static_cast<T>(2.87622) * y -
           static_cast<T>(0.0116268) * x * y -
           static_cast<T>(0.00000603925) * x2 * y -
           static_cast<T>(0.0389806) * y2 -
           static_cast<T>(0.0000403488) * x * y2 +
           static_cast<T>(0.000168556) * y3;
}

template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_conus_deg(const LatLon<T>& p) noexcept {
    return magnetic_variation_conus_deg(p.lat_deg(), p.lon_west_deg());
}

/**
 * @brief Least-squares polynomial fit for magnetic variation (in degrees, West positive)
 * for Alaska (lat_n_deg > 54, 130 < lon_w_deg < 172).
 */
template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_alaska_deg(T lat_n_deg, T lon_w_deg) noexcept {
    const T x = lat_n_deg;
    const T y = lon_w_deg;
    const T x2 = x * x;
    const T x3 = x2 * x;
    const T y2 = y * y;
    const T y3 = y2 * y;
    return static_cast<T>(618.854) +
           static_cast<T>(2.76049) * x -
           static_cast<T>(0.556206) * x2 +
           static_cast<T>(0.00251582) * x3 -
           static_cast<T>(12.7974) * y +
           static_cast<T>(0.408161) * x * y +
           static_cast<T>(0.000434097) * x2 * y -
           static_cast<T>(0.00602173) * y2 -
           static_cast<T>(0.00144712) * x * y2 +
           static_cast<T>(0.000222521) * y3;
}

template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_alaska_deg(const LatLon<T>& p) noexcept {
    return magnetic_variation_alaska_deg(p.lat_deg(), p.lon_west_deg());
}

/**
 * @brief Least-squares polynomial fit for magnetic variation (in degrees) for Western Europe
 * (36 < lat_n_deg < 68, -10 < lon_e_deg < 28, where *East* longitude is positive!).
 */
template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_western_europe_deg(T lat_n_deg, T lon_e_deg) noexcept {
    const T lat = lat_n_deg;
    const T lon = lon_e_deg;
    const T lat2 = lat * lat;
    const T lat3 = lat2 * lat;
    const T lon2 = lon * lon;
    const T lon3 = lon2 * lon;
    return static_cast<T>(10.4768771667158) -
           static_cast<T>(0.507385322418858) * lon +
           static_cast<T>(0.00753170031703826) * lon2 -
           static_cast<T>(1.40596203924748e-05) * lon3 -
           static_cast<T>(0.535560699962353) * lat +
           static_cast<T>(0.0154348808069955) * lat * lon -
           static_cast<T>(8.07756425110592e-05) * lat * lon2 +
           static_cast<T>(0.00976887198864442) * lat2 -
           static_cast<T>(0.000259163929798334) * lat2 * lon -
           static_cast<T>(3.69056939266123e-05) * lat3;
}

template <typename T = double>
[[nodiscard]] constexpr T magnetic_variation_western_europe_deg(const LatLon<T>& p) noexcept {
    return magnetic_variation_western_europe_deg(p.lat_deg(), p.lon_east_deg());
}

} // namespace avform::flight_dynamics

#endif // AVFORM_FLIGHT_DYNAMICS_HPP

