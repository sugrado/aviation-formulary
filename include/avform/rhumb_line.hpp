#ifndef AVFORM_RHUMB_LINE_HPP
#define AVFORM_RHUMB_LINE_HPP

#include "core.hpp"
#include "units.hpp"

namespace avform::rhumb_line {

template <typename T = double>
struct RhumbRoute {
    T tc{};     ///< True course in radians [0, 2*pi)
    T d{};      ///< Angular distance in radians
    T dlon_w{}; ///< Westerly longitude difference in radians
    T dlon_e{}; ///< Easterly longitude difference in radians
    T dphi{};   ///< Mercator latitude difference
    T q{};      ///< Scale factor (lat2 - lat1) / dphi
};

/**
 * @brief Mercator latitude difference dphi = log(tan(lat2/2 + pi/4) / tan(lat1/2 + pi/4)).
 */
template <typename T = double>
[[nodiscard]] inline T mercator_dphi(T lat1, T lat2) noexcept {
    return std::log(std::tan(lat2 / static_cast<T>(2) + constants::quarter_pi<T>) /
                    std::tan(lat1 / static_cast<T>(2) + constants::quarter_pi<T>));
}

/**
 * @brief Total angular length (in radians) of a non-East-West rhumb line from pole to pole:
 * length = pi / abs(cos(tc)).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<T> pole_to_pole_length_rad(T tc,
                                                              T eps = constants::default_eps<T>) noexcept {
    const T abs_cos = std::abs(std::cos(tc));
    if (abs_cos < eps) {
        return std::nullopt; // E-W rhumb line follows a parallel and does not reach the poles
    }
    return constants::pi<T> / abs_cos;
}

/**
 * @brief Computes true course (in radians) between p1 and p2 without checking 180 E/W meridian wrap.
 * Formula: tc = mod(atan2(lon1-lon2, log(tan(lat2/2+pi/4)/tan(lat1/2+pi/4))), 2*pi)
 */
template <typename T = double>
[[nodiscard]] inline T course_direct(const LatLon<T>& p1, const LatLon<T>& p2) noexcept {
    return mod(std::atan2(p1.lon - p2.lon, mercator_dphi(p1.lat, p2.lat)), constants::two_pi<T>);
}

/**
 * @brief Computes angular distance (in radians) between p1 and p2 without checking 180 E/W meridian wrap.
 */
template <typename T = double>
[[nodiscard]] inline T distance_direct(const LatLon<T>& p1, const LatLon<T>& p2,
                                       T tol = constants::default_tol<T>) noexcept {
    const T dlat = p2.lat - p1.lat;
    T q{};
    if (std::abs(dlat) < std::sqrt(tol)) {
        q = std::cos(p1.lat);
    } else {
        q = dlat / mercator_dphi(p1.lat, p2.lat);
    }
    const T dlon = p2.lon - p1.lon;
    return std::sqrt(dlat * dlat + q * q * dlon * dlon);
}

/**
 * @brief Computes the true course (`tc`) and angular distance (`d`) for the shortest rhumb line
 * connecting `p1` and `p2`, allowing crossing of the 180 E/W meridian.
 */
template <typename T = double>
[[nodiscard]] inline RhumbRoute<T> shortest_route(const LatLon<T>& p1, const LatLon<T>& p2,
                                                  T tol = constants::default_tol<T>) noexcept {
    RhumbRoute<T> res{};
    res.dlon_w = mod(p2.lon - p1.lon, constants::two_pi<T>);
    res.dlon_e = mod(p1.lon - p2.lon, constants::two_pi<T>);
    res.dphi = mercator_dphi(p1.lat, p2.lat);

    const T dlat = p2.lat - p1.lat;
    if (std::abs(dlat) < std::sqrt(tol)) {
        res.q = std::cos(p1.lat);
    } else {
        res.q = dlat / res.dphi;
    }

    if (res.dlon_w < res.dlon_e) {
        // Westerly rhumb line is the shortest
        res.tc = mod(std::atan2(-res.dlon_w, res.dphi), constants::two_pi<T>);
        res.d = std::sqrt(res.q * res.q * res.dlon_w * res.dlon_w + dlat * dlat);
    } else {
        // Easterly rhumb line is the shortest
        res.tc = mod(std::atan2(res.dlon_e, res.dphi), constants::two_pi<T>);
        res.d = std::sqrt(res.q * res.q * res.dlon_e * res.dlon_e + dlat * dlat);
    }
    return res;
}

/**
 * @brief True course (in radians) along the shortest rhumb line from p1 to p2.
 */
template <typename T = double>
[[nodiscard]] inline T course(const LatLon<T>& p1, const LatLon<T>& p2,
                              T tol = constants::default_tol<T>) noexcept {
    return shortest_route(p1, p2, tol).tc;
}

/**
 * @brief True course (in degrees) along the shortest rhumb line from p1 to p2.
 */
template <typename T = double>
[[nodiscard]] inline T course_deg(const LatLon<T>& p1, const LatLon<T>& p2,
                                  T tol = constants::default_tol<T>) noexcept {
    return units::rad_to_deg(course(p1, p2, tol));
}

/**
 * @brief Angular distance (in radians) along the shortest rhumb line from p1 to p2.
 */
template <typename T = double>
[[nodiscard]] inline T distance(const LatLon<T>& p1, const LatLon<T>& p2,
                                T tol = constants::default_tol<T>) noexcept {
    return shortest_route(p1, p2, tol).d;
}

/**
 * @brief Distance (in nautical miles) along the shortest rhumb line from p1 to p2.
 */
template <typename T = double>
[[nodiscard]] inline T distance_nm(const LatLon<T>& p1, const LatLon<T>& p2,
                                   T tol = constants::default_tol<T>) noexcept {
    return units::rad_to_nm(distance(p1, p2, tol));
}

/**
 * @brief Finds the destination LatLon at true course `tc` (radians) and angular distance `d` (radians)
 * from `p1` along a rhumb line.
 * Returns std::nullopt if `abs(lat) > pi/2` ("d too large. You can't go this far along this rhumb line!")
 * or if the initial point is a pole.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<LatLon<T>> point_from_course_distance(const LatLon<T>& p1, T tc, T d,
                                                                         T tol = constants::default_tol<T>) noexcept {
    if (std::cos(p1.lat) < tol) {
        return std::nullopt; // Initial point cannot be a pole
    }
    const T lat = p1.lat + d * std::cos(tc);
    if (std::abs(lat) > constants::half_pi<T>) {
        return std::nullopt; // d too large
    }

    T q{};
    if (std::abs(lat - p1.lat) < std::sqrt(tol)) {
        q = std::cos(p1.lat);
    } else {
        const T dphi = mercator_dphi(p1.lat, lat);
        q = (lat - p1.lat) / dphi;
    }

    if (std::abs(q) < tol) {
        return std::nullopt;
    }

    const T dlon = -d * std::sin(tc) / q;
    const T lon = normalize_lon(p1.lon + dlon);
    return LatLon<T>{lat, lon};
}

/**
 * @brief Convenience overload taking true course in degrees and distance in nautical miles.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<LatLon<T>> point_from_course_distance_nm(const LatLon<T>& p1, T tc_deg, T d_nm,
                                                                            T tol = constants::default_tol<T>) noexcept {
    return point_from_course_distance(p1, units::deg_to_rad(tc_deg), units::nm_to_rad(d_nm), tol);
}

} // namespace avform::rhumb_line

#endif // AVFORM_RHUMB_LINE_HPP

