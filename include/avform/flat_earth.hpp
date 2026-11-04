#ifndef AVFORM_FLAT_EARTH_HPP
#define AVFORM_FLAT_EARTH_HPP

#include "core.hpp"
#include "units.hpp"

namespace avform::flat_earth {

/**
 * @brief Reference ellipsoid parameters from "Ellipsoidal parameters" section.
 */
template <typename T = double>
struct Ellipsoid {
    T a{}; ///< Semi-major axis (equatorial radius) in km
    T f{}; ///< Flattening

    [[nodiscard]] constexpr T b() const noexcept {
        return a * (static_cast<T>(1) - f);
    }

    [[nodiscard]] constexpr T e_squared() const noexcept {
        return f * (static_cast<T>(2) - f);
    }

    static constexpr Ellipsoid wgs84() noexcept {
        return Ellipsoid{static_cast<T>(6378.13700), static_cast<T>(1.0L / 298.257223563L)};
    }

    static constexpr Ellipsoid grs80_nad83() noexcept {
        return Ellipsoid{static_cast<T>(6378.13700), static_cast<T>(1.0L / 298.257222101L)};
    }

    static constexpr Ellipsoid wgs66() noexcept {
        return Ellipsoid{static_cast<T>(6378.145), static_cast<T>(1.0L / 298.25L)};
    }

    static constexpr Ellipsoid grs67_iau68() noexcept {
        return Ellipsoid{static_cast<T>(6378.16000), static_cast<T>(1.0L / 298.2472L)};
    }

    static constexpr Ellipsoid wgs72() noexcept {
        return Ellipsoid{static_cast<T>(6378.135), static_cast<T>(1.0L / 298.26L)};
    }

    static constexpr Ellipsoid krasovsky() noexcept {
        return Ellipsoid{static_cast<T>(6378.245), static_cast<T>(1.0L / 298.3L)};
    }

    static constexpr Ellipsoid clarke66_nad27() noexcept {
        return Ellipsoid{static_cast<T>(6378.2064), static_cast<T>(1.0L / 294.9786982138L)};
    }
};

template <typename T = double>
struct RadiiOfCurvature {
    T r1{}; ///< Meridional radius of curvature (same unit as ellipsoid `a`)
    T r2{}; ///< Radius of curvature in the prime vertical (same unit as ellipsoid `a`)
};

template <typename T = double>
struct FlatEarthPosition {
    T distance_north{}; ///< North offset
    T distance_east{};  ///< East offset
    T distance{};       ///< Total planar distance
    T bearing{};        ///< True bearing from origin to target in radians [0, 2*pi)
};

/**
 * @brief Computes meridional radius of curvature (R1) and radius of curvature in the prime vertical (R2)
 * at latitude `lat0` (radians).
 * Formulas:
 *   R1 = a*(1-e^2) / (1 - e^2*(sin(lat0))^2)^(3/2)
 *   R2 = a / sqrt(1 - e^2*(sin(lat0))^2)
 */
template <typename T = double>
[[nodiscard]] inline RadiiOfCurvature<T> radii_of_curvature(T lat0,
                                                            const Ellipsoid<T>& ell = Ellipsoid<T>::wgs84()) noexcept {
    const T e2 = ell.e_squared();
    const T sin_lat = std::sin(lat0);
    const T denom_sq = static_cast<T>(1) - e2 * sin_lat * sin_lat;
    const T denom_sqrt = std::sqrt(denom_sq);
    const T r1 = ell.a * (static_cast<T>(1) - e2) / (denom_sq * denom_sqrt);
    const T r2 = ell.a / denom_sqrt;
    return RadiiOfCurvature<T>{r1, r2};
}

/**
 * @brief Local flat-earth approximation given changes `dlat = lat - lat0` (North positive)
 * and `dlon_east = lon_east - lon0_east` (East positive), and radii `r1`, `r2`.
 */
template <typename T = double>
[[nodiscard]] inline FlatEarthPosition<T> flat_earth_from_deltas(T lat0, T dlat, T dlon_east,
                                                                 T r1, T r2) noexcept {
    FlatEarthPosition<T> res{};
    res.distance_north = r1 * dlat;
    res.distance_east = r2 * std::cos(lat0) * dlon_east;
    res.distance = std::hypot(res.distance_north, res.distance_east);
    res.bearing = mod(std::atan2(res.distance_east, res.distance_north), constants::two_pi<T>);
    return res;
}

/**
 * @brief Local flat-earth approximation on an ellipsoid (distances in km by default, matching `ell.a`).
 * Note: `LatLon` uses the Formulary convention (West positive), so `dlon_east = p0.lon - p.lon`.
 */
template <typename T = double>
[[nodiscard]] inline FlatEarthPosition<T> flat_earth_ellipsoid(const LatLon<T>& p0, const LatLon<T>& p,
                                                               const Ellipsoid<T>& ell = Ellipsoid<T>::wgs84()) noexcept {
    const auto radii = radii_of_curvature(p0.lat, ell);
    const T dlat = p.lat - p0.lat;
    const T dlon_east = p0.lon - p.lon; // Convert West-positive to East-positive delta
    return flat_earth_from_deltas(p0.lat, dlat, dlon_east, radii.r1, radii.r2);
}

/**
 * @brief Local flat-earth approximation on a sphere of radius `R` (use R = 1 for radians,
 * or R = constants::earth_radius_nm<T> for nautical miles).
 */
template <typename T = double>
[[nodiscard]] inline FlatEarthPosition<T> flat_earth_spherical(const LatLon<T>& p0, const LatLon<T>& p,
                                                               T radius = static_cast<T>(1)) noexcept {
    const T dlat = p.lat - p0.lat;
    const T dlon_east = p0.lon - p.lon; // Convert West-positive to East-positive delta
    return flat_earth_from_deltas(p0.lat, dlat, dlon_east, radius, radius);
}

// ============================================================================
// Geodetic <-> Geocentric Coordinate Conversion (#Conv)
// ============================================================================

template <typename T = double>
struct GeocentricCoord {
    T u{}; ///< Geocentric latitude in radians
    T r{}; ///< Geocentric radius (in same length units as ellipsoid `a` and `h`)
};

/**
 * @brief Converts geodetic coordinates (geodetic latitude `v` in radians, height `h` above ellipsoid)
 * to geocentric coordinates (geocentric latitude `u` in radians, geocentric radius `r`).
 * Note: `h` must be in the same units as `ell.a` (kilometers for standard Ellipsoid presets).
 * Formulas:
 *   tan(u) = tan(v) * (h*sqrt((a*cos(v))^2 + (b*sin(v))^2) + b^2) /
 *                     (h*sqrt((a*cos(v))^2 + (b*sin(v))^2) + a^2)
 *   r^2 = h^2 + 2*h*sqrt((a*cos(v))^2 + (b*sin(v))^2) +
 *         (a^4 - (a^4 - b^4)*(sin(v))^2) / (a^2 - (a^2 - b^2)*(sin(v))^2)
 */
template <typename T = double>
[[nodiscard]] inline GeocentricCoord<T> geodetic_to_geocentric(T v, T h_km,
                                                               const Ellipsoid<T>& ell = Ellipsoid<T>::wgs84()) noexcept {
    const T a = ell.a;
    const T b = ell.b();
    const T a2 = a * a;
    const T b2 = b * b;
    const T a4 = a2 * a2;
    const T b4 = b2 * b2;

    const T cos_v = std::cos(v);
    const T sin_v = std::sin(v);
    const T sin_v2 = sin_v * sin_v;

    const T s = std::hypot(a * cos_v, b * sin_v);
    const T tan_u = std::tan(v) * (h_km * s + b2) / (h_km * s + a2);
    const T u = std::atan(tan_u);

    const T r_sq = h_km * h_km + static_cast<T>(2) * h_km * s +
                   (a4 - (a4 - b4) * sin_v2) / (a2 - (a2 - b2) * sin_v2);
    const T r = std::sqrt(std::max(static_cast<T>(0), r_sq));
    return GeocentricCoord<T>{u, r};
}

} // namespace avform::flat_earth

#endif // AVFORM_FLAT_EARTH_HPP

