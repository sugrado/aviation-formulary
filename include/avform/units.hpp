#ifndef AVFORM_UNITS_HPP
#define AVFORM_UNITS_HPP

#include "core.hpp"

namespace avform::units {

/**
 * @brief Degrees, Minutes, Seconds representation.
 */
template <typename T = double>
struct DMS {
    int degrees{};
    int minutes{};
    T seconds{};
};

// ============================================================================
// Angle Conversions
// ============================================================================

template <typename T = double>
[[nodiscard]] constexpr T deg_to_rad(T angle_degrees) noexcept {
    return (constants::pi<T> / static_cast<T>(180)) * angle_degrees;
}

template <typename T = double>
[[nodiscard]] constexpr T rad_to_deg(T angle_radians) noexcept {
    return (static_cast<T>(180) / constants::pi<T>) * angle_radians;
}

/**
 * @brief Converts degrees, minutes, and seconds to decimal degrees.
 * If `degrees` is negative, `minutes` and `seconds` are subtracted (or can be passed with same sign).
 */
template <typename T = double>
[[nodiscard]] constexpr T dms_to_deg(T degrees, T minutes, T seconds = static_cast<T>(0)) noexcept {
    T sign = (degrees < static_cast<T>(0)) ? static_cast<T>(-1) : static_cast<T>(1);
    T abs_deg = (degrees < static_cast<T>(0)) ? -degrees : degrees;
    T abs_min = (minutes < static_cast<T>(0)) ? -minutes : minutes;
    T abs_sec = (seconds < static_cast<T>(0)) ? -seconds : seconds;
    return sign * (abs_deg + (abs_min / static_cast<T>(60)) + (abs_sec / static_cast<T>(3600)));
}

template <typename T = double>
[[nodiscard]] constexpr T dms_to_rad(T degrees, T minutes, T seconds = static_cast<T>(0)) noexcept {
    return deg_to_rad(dms_to_deg(degrees, minutes, seconds));
}

/**
 * @brief Converts decimal degrees to degrees, minutes, and seconds.
 * Reference: "Wind Triangles":
 *   degrees = int(angle_degrees)
 *   minutes = int(60 * (angle_degrees - degrees))
 *   seconds = 60 * (60 * (angle_degrees - degrees) - minutes)
 */
template <typename T = double>
[[nodiscard]] constexpr DMS<T> deg_to_dms(T angle_degrees) noexcept {
    int deg = static_cast<int>(angle_degrees);
    T rem_deg = angle_degrees - static_cast<T>(deg);
    int min = static_cast<int>(static_cast<T>(60) * rem_deg);
    T sec = static_cast<T>(60) * (static_cast<T>(60) * rem_deg - static_cast<T>(min));
    return DMS<T>{deg, min, sec};
}

// ============================================================================
// Great Circle Distance Conversions (Radians <-> Nautical Miles / Kilometers)
// ============================================================================

/**
 * @brief Converts distance in nautical miles to radians.
 * Formula: distance_radians = (pi / (180 * 60)) * distance_nm
 */
template <typename T = double>
[[nodiscard]] constexpr T nm_to_rad(T distance_nm) noexcept {
    return (constants::pi<T> / static_cast<T>(180 * 60)) * distance_nm;
}

/**
 * @brief Converts angular distance in radians to nautical miles.
 * Formula: distance_nm = ((180 * 60) / pi) * distance_radians
 */
template <typename T = double>
[[nodiscard]] constexpr T rad_to_nm(T distance_radians) noexcept {
    return (static_cast<T>(180 * 60) / constants::pi<T>) * distance_radians;
}

/**
 * @brief Converts angular distance in radians to kilometers using a specified spherical radius.
 * Formula: distance_km = radius_km * distance_radians
 */
template <typename T = double>
[[nodiscard]] constexpr T rad_to_km(T distance_radians, T radius_km = constants::earth_radius_km_nm<T>) noexcept {
    return radius_km * distance_radians;
}

/**
 * @brief Converts distance in kilometers to angular distance in radians using a specified spherical radius.
 * Formula: distance_radians = distance_km / radius_km
 */
template <typename T = double>
[[nodiscard]] constexpr T km_to_rad(T distance_km, T radius_km = constants::earth_radius_km_nm<T>) noexcept {
    return distance_km / radius_km;
}

// ============================================================================
// Linear Distance Conversions
// ============================================================================

template <typename T = double>
[[nodiscard]] constexpr T nm_to_km(T nm) noexcept {
    return nm * static_cast<T>(1.852);
}

template <typename T = double>
[[nodiscard]] constexpr T km_to_nm(T km) noexcept {
    return km / static_cast<T>(1.852);
}

template <typename T = double>
[[nodiscard]] constexpr T ft_to_m(T ft) noexcept {
    return ft * static_cast<T>(0.3048);
}

template <typename T = double>
[[nodiscard]] constexpr T m_to_ft(T m) noexcept {
    return m / static_cast<T>(0.3048);
}

// ============================================================================
// Speed Conversions ("Unit conversions, etc.")
// ============================================================================

// Exact conversion factors from the Formulary:
// 1 knot = 1.852000 km/hr*
// 1 knot = 185200/109728 ft/sec* = 1.687810 ft/sec
// 1 knot = 1852000/1609344 mph* = 1.150779 mph
// 1 mph  = 1.609344 km/hr*

template <typename T = double>
[[nodiscard]] constexpr T kt_to_kmh(T kt) noexcept {
    return kt * static_cast<T>(1.852);
}

template <typename T = double>
[[nodiscard]] constexpr T kmh_to_kt(T kmh) noexcept {
    return kmh / static_cast<T>(1.852);
}

template <typename T = double>
[[nodiscard]] constexpr T kt_to_fps(T kt) noexcept {
    return kt * (static_cast<T>(185200.0L) / static_cast<T>(109728.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T fps_to_kt(T fps) noexcept {
    return fps * (static_cast<T>(109728.0L) / static_cast<T>(185200.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T kt_to_mph(T kt) noexcept {
    return kt * (static_cast<T>(1852000.0L) / static_cast<T>(1609344.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T mph_to_kt(T mph) noexcept {
    return mph * (static_cast<T>(1609344.0L) / static_cast<T>(1852000.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T mph_to_kmh(T mph) noexcept {
    return mph * static_cast<T>(1.609344);
}

template <typename T = double>
[[nodiscard]] constexpr T kmh_to_mph(T kmh) noexcept {
    return kmh / static_cast<T>(1.609344);
}

template <typename T = double>
[[nodiscard]] constexpr T mph_to_fps(T mph) noexcept {
    return mph * (static_cast<T>(5280.0L) / static_cast<T>(3600.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T fps_to_mph(T fps) noexcept {
    return fps * (static_cast<T>(3600.0L) / static_cast<T>(5280.0L));
}

template <typename T = double>
[[nodiscard]] constexpr T kmh_to_fps(T kmh) noexcept {
    return kmh * (static_cast<T>(1000.0L) / static_cast<T>(3600.0L * 0.3048L));
}

template <typename T = double>
[[nodiscard]] constexpr T fps_to_kmh(T fps) noexcept {
    return fps * (static_cast<T>(3600.0L * 0.3048L) / static_cast<T>(1000.0L));
}

// ============================================================================
// Temperature & Pressure Conversions
// ============================================================================

template <typename T = double>
[[nodiscard]] constexpr T celsius_to_kelvin(T c) noexcept {
    return c + static_cast<T>(273.15);
}

template <typename T = double>
[[nodiscard]] constexpr T kelvin_to_celsius(T k) noexcept {
    return k - static_cast<T>(273.15);
}

template <typename T = double>
[[nodiscard]] constexpr T fahrenheit_to_celsius(T f) noexcept {
    return (f - static_cast<T>(32)) * (static_cast<T>(5) / static_cast<T>(9));
}

template <typename T = double>
[[nodiscard]] constexpr T celsius_to_fahrenheit(T c) noexcept {
    return c * (static_cast<T>(9) / static_cast<T>(5)) + static_cast<T>(32);
}

template <typename T = double>
[[nodiscard]] constexpr T inhg_to_mb(T inhg) noexcept {
    return inhg * (static_cast<T>(1013.25) / static_cast<T>(29.92126));
}

template <typename T = double>
[[nodiscard]] constexpr T mb_to_inhg(T mb) noexcept {
    return mb * (static_cast<T>(29.92126) / static_cast<T>(1013.25));
}

} // namespace avform::units

#endif // AVFORM_UNITS_HPP

