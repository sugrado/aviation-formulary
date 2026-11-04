#ifndef AVFORM_CORE_HPP
#define AVFORM_CORE_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

namespace avform {

namespace constants {

template <typename T = double>
inline constexpr T pi = static_cast<T>(3.141592653589793238462643383279502884L);

template <typename T = double>
inline constexpr T two_pi = static_cast<T>(2.0L * 3.141592653589793238462643383279502884L);

template <typename T = double>
inline constexpr T half_pi = static_cast<T>(0.5L * 3.141592653589793238462643383279502884L);

template <typename T = double>
inline constexpr T quarter_pi = static_cast<T>(0.25L * 3.141592653589793238462643383279502884L);

template <typename T = double>
inline constexpr T default_eps = static_cast<T>(1e-15);

template <typename T = double>
inline constexpr T default_tol = static_cast<T>(1e-15);

// Earth radius implied by 1 nm = 1 arcminute: (180 * 60) / pi nautical miles
template <typename T = double>
inline constexpr T earth_radius_nm = static_cast<T>((180.0L * 60.0L) / 3.141592653589793238462643383279502884L);

// Earth radius in km implied by 1 nm = 1.852 km (~6366.7070195 km)
template <typename T = double>
inline constexpr T earth_radius_km_nm = static_cast<T>(1.852L * (180.0L * 60.0L) / 3.141592653589793238462643383279502884L);

// Former FAI standard sphere radius (6371.0 km)
template <typename T = double>
inline constexpr T earth_radius_km_fai = static_cast<T>(6371.0);

// WGS84 equatorial and polar radii in km
template <typename T = double>
inline constexpr T earth_radius_km_wgs84_eq = static_cast<T>(6378.137);

template <typename T = double>
inline constexpr T earth_radius_km_wgs84_polar = static_cast<T>(6356.752);

} // namespace constants

/**
 * @brief Safe arc-sine function clamped to [-1, 1] to prevent NaN from rounding errors.
 * Reference: "Implementation notes - Notes on mathematical functions"
 */
template <typename T>
[[nodiscard]] inline T asin_safe(T x) noexcept {
    return std::asin(std::clamp(x, static_cast<T>(-1), static_cast<T>(1)));
}

/**
 * @brief Safe arc-cosine function clamped to [-1, 1] to prevent NaN from rounding errors.
 * Reference: "Implementation notes - Notes on mathematical functions"
 */
template <typename T>
[[nodiscard]] inline T acos_safe(T x) noexcept {
    return std::acos(std::clamp(x, static_cast<T>(-1), static_cast<T>(1)));
}

/**
 * @brief Euclidean modulo function where the result always lies in [0, x) for x > 0.
 * Reference: "Implementation notes - Note on the mod function":
 *   mod(y, x) = y - x * floor(y / x)
 *   For instance: mod(2.3, 2.) = 0.3 and mod(-2.3, 2.) = 1.7
 */
template <typename T>
[[nodiscard]] inline T mod(T y, T x) noexcept {
    T m = y - x * std::floor(y / x);
    if (m < static_cast<T>(0)) {
        m += x;
    }
    return m;
}

/**
 * @brief Normalizes an angle/course in radians to the range [0, 2*pi).
 */
template <typename T>
[[nodiscard]] inline T normalize_course(T angle_rad) noexcept {
    return mod(angle_rad, constants::two_pi<T>);
}

/**
 * @brief Normalizes a longitude in radians to the range [-pi, pi).
 * Formula: mod(lon + pi, 2*pi) - pi
 */
template <typename T>
[[nodiscard]] inline T normalize_lon(T lon_rad) noexcept {
    return mod(lon_rad + constants::pi<T>, constants::two_pi<T>) - constants::pi<T>;
}

/**
 * @brief Fallback mathematical implementations using only atan / atan2,
 * as given in "Implementation notes".
 */
namespace fallback_math {

template <typename T>
[[nodiscard]] inline T acos_via_atan2(T x) noexcept {
    x = std::clamp(x, static_cast<T>(-1), static_cast<T>(1));
    return std::atan2(std::sqrt(static_cast<T>(1) - x * x), x);
}

template <typename T>
[[nodiscard]] inline T asin_via_atan2(T x) noexcept {
    x = std::clamp(x, static_cast<T>(-1), static_cast<T>(1));
    return std::atan2(x, std::sqrt(static_cast<T>(1) - x * x));
}

template <typename T>
[[nodiscard]] inline T asin_via_atan(T x) noexcept {
    x = std::clamp(x, static_cast<T>(-1), static_cast<T>(1));
    return static_cast<T>(2) * std::atan(x / (static_cast<T>(1) + std::sqrt(static_cast<T>(1) - x * x)));
}

template <typename T>
[[nodiscard]] inline T acos_via_atan(T x) noexcept {
    x = std::clamp(x, static_cast<T>(-1), static_cast<T>(1));
    if (x <= static_cast<T>(-1)) {
        return constants::pi<T>;
    }
    if (x >= static_cast<T>(0)) {
        return static_cast<T>(2) * std::atan(std::sqrt((static_cast<T>(1) - x) / (static_cast<T>(1) + x)));
    }
    return constants::pi<T> - static_cast<T>(2) * std::atan(std::sqrt((static_cast<T>(1) + x) / (static_cast<T>(1) - x)));
}

template <typename T>
[[nodiscard]] inline std::optional<T> atan2_via_atan(T y, T x) noexcept {
    if (x > static_cast<T>(0)) {
        return std::atan(y / x);
    }
    if (x < static_cast<T>(0) && y >= static_cast<T>(0)) {
        return std::atan(y / x) + constants::pi<T>;
    }
    if (x == static_cast<T>(0) && y > static_cast<T>(0)) {
        return constants::half_pi<T>;
    }
    if (x < static_cast<T>(0) && y < static_cast<T>(0)) {
        return std::atan(y / x) - constants::pi<T>;
    }
    if (x == static_cast<T>(0) && y < static_cast<T>(0)) {
        return -constants::half_pi<T>;
    }
    return std::nullopt; // atan2(0, 0) is undefined
}

} // namespace fallback_math

/**
 * @brief Geographic coordinate on a sphere/ellipsoid in radians.
 *
 * Sign convention (per Ed Williams' Aviation Formulary):
 * - Latitude (`lat`): North is positive (+), South is negative (-).
 * - Longitude (`lon`): West is positive (+), East is negative (-).
 *
 * Helper factory methods are provided for both West-positive and East-positive inputs.
 */
template <typename T = double>
struct LatLon {
    T lat{}; ///< Latitude in radians (North positive)
    T lon{}; ///< Longitude in radians (West positive per Formulary convention)

    constexpr LatLon() noexcept = default;
    constexpr LatLon(T lat_rad, T lon_west_rad) noexcept : lat(lat_rad), lon(lon_west_rad) {}

    [[nodiscard]] static constexpr LatLon from_rad(T lat_rad, T lon_west_rad) noexcept {
        return LatLon(lat_rad, lon_west_rad);
    }

    [[nodiscard]] static constexpr LatLon from_rad_east_positive(T lat_rad, T lon_east_rad) noexcept {
        return LatLon(lat_rad, -lon_east_rad);
    }

    [[nodiscard]] static constexpr LatLon from_deg(T lat_north_deg, T lon_west_deg) noexcept {
        constexpr T deg2rad = constants::pi<T> / static_cast<T>(180);
        return LatLon(lat_north_deg * deg2rad, lon_west_deg * deg2rad);
    }

    [[nodiscard]] static constexpr LatLon from_deg_east_positive(T lat_north_deg, T lon_east_deg) noexcept {
        constexpr T deg2rad = constants::pi<T> / static_cast<T>(180);
        return LatLon(lat_north_deg * deg2rad, -lon_east_deg * deg2rad);
    }

    [[nodiscard]] constexpr T lat_deg() const noexcept {
        return lat * (static_cast<T>(180) / constants::pi<T>);
    }

    [[nodiscard]] constexpr T lon_west_deg() const noexcept {
        return lon * (static_cast<T>(180) / constants::pi<T>);
    }

    [[nodiscard]] constexpr T lon_east_deg() const noexcept {
        return -lon * (static_cast<T>(180) / constants::pi<T>);
    }

    [[nodiscard]] constexpr T lon_east_rad() const noexcept {
        return -lon;
    }
};

} // namespace avform

#endif // AVFORM_CORE_HPP

