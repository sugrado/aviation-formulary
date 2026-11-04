#ifndef AVFORM_GREAT_CIRCLE_HPP
#define AVFORM_GREAT_CIRCLE_HPP

#include "core.hpp"
#include "units.hpp"
#include <utility>

namespace avform::great_circle {

// ============================================================================
// Distance between points (#Dist)
// ============================================================================

/**
 * @brief Great circle distance in radians using the spherical law of cosines.
 * Formula: d = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1-lon2))
 */
template <typename T = double>
[[nodiscard]] inline T distance_cos(const LatLon<T>& p1, const LatLon<T>& p2) noexcept {
    return acos_safe(std::sin(p1.lat) * std::sin(p2.lat) +
                     std::cos(p1.lat) * std::cos(p2.lat) * std::cos(p1.lon - p2.lon));
}

/**
 * @brief Great circle distance in radians using the Haversine formula (less subject to rounding error for short distances).
 * Formula: d = 2*asin(sqrt((sin((lat1-lat2)/2))^2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))
 */
template <typename T = double>
[[nodiscard]] inline T distance(const LatLon<T>& p1, const LatLon<T>& p2) noexcept {
    const T sin_dlat = std::sin((p1.lat - p2.lat) / static_cast<T>(2));
    const T sin_dlon = std::sin((p1.lon - p2.lon) / static_cast<T>(2));
    const T a = sin_dlat * sin_dlat + std::cos(p1.lat) * std::cos(p2.lat) * sin_dlon * sin_dlon;
    return static_cast<T>(2) * asin_safe(std::sqrt(std::max(static_cast<T>(0), a)));
}

/**
 * @brief Great circle distance in nautical miles.
 */
template <typename T = double>
[[nodiscard]] inline T distance_nm(const LatLon<T>& p1, const LatLon<T>& p2) noexcept {
    return units::rad_to_nm(distance(p1, p2));
}

/**
 * @brief Great circle distance in kilometers.
 */
template <typename T = double>
[[nodiscard]] inline T distance_km(const LatLon<T>& p1, const LatLon<T>& p2,
                                   T radius_km = constants::earth_radius_km_nm<T>) noexcept {
    return units::rad_to_km(distance(p1, p2), radius_km);
}

// ============================================================================
// Course between points (#Crs)
// ============================================================================

/**
 * @brief Initial true course (in radians) from p1 to p2 given pre-computed distance d.
 * Includes special handling when starting from a pole.
 */
template <typename T = double>
[[nodiscard]] inline T initial_course_acos(const LatLon<T>& p1, const LatLon<T>& p2, T d,
                                           T eps = constants::default_eps<T>) noexcept {
    if (std::cos(p1.lat) < eps) {
        return (p1.lat > static_cast<T>(0)) ? constants::pi<T> : constants::two_pi<T>;
    }
    if (std::sin(d) < eps) {
        return static_cast<T>(0);
    }
    const T arg = (std::sin(p2.lat) - std::sin(p1.lat) * std::cos(d)) / (std::sin(d) * std::cos(p1.lat));
    if (std::sin(p2.lon - p1.lon) < static_cast<T>(0)) {
        return acos_safe(arg);
    }
    return constants::two_pi<T> - acos_safe(arg);
}

/**
 * @brief Initial true course (in radians) from p1 to p2 without requiring pre-computed distance.
 * Formula: tc1 = mod(atan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2)), 2*pi)
 */
template <typename T = double>
[[nodiscard]] inline T initial_course(const LatLon<T>& p1, const LatLon<T>& p2,
                                      T eps = constants::default_eps<T>) noexcept {
    if (std::cos(p1.lat) < eps) {
        return (p1.lat > static_cast<T>(0)) ? constants::pi<T> : constants::two_pi<T>;
    }
    const T y = std::sin(p1.lon - p2.lon) * std::cos(p2.lat);
    const T x = std::cos(p1.lat) * std::sin(p2.lat) -
                std::sin(p1.lat) * std::cos(p2.lat) * std::cos(p1.lon - p2.lon);
    return mod(std::atan2(y, x), constants::two_pi<T>);
}

/**
 * @brief Initial true course in degrees [0, 360).
 */
template <typename T = double>
[[nodiscard]] inline T initial_course_deg(const LatLon<T>& p1, const LatLon<T>& p2,
                                          T eps = constants::default_eps<T>) noexcept {
    return units::rad_to_deg(initial_course(p1, p2, eps));
}

// ============================================================================
// Latitude of point on GC (#Int)
// ============================================================================

/**
 * @brief Finds the latitude (in radians) where the great circle connecting p1 and p2 crosses meridian `lon`.
 * Returns std::nullopt if p1 and p2 lie on the same meridian (or opposite meridians where sin(lon1 - lon2) == 0).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<T> latitude_at_longitude(const LatLon<T>& p1, const LatLon<T>& p2, T lon,
                                                            T eps = constants::default_eps<T>) noexcept {
    const T denom = std::cos(p1.lat) * std::cos(p2.lat) * std::sin(p1.lon - p2.lon);
    if (std::abs(denom) < eps) {
        return std::nullopt;
    }
    const T num = std::sin(p1.lat) * std::cos(p2.lat) * std::sin(lon - p2.lon) -
                  std::sin(p2.lat) * std::cos(p1.lat) * std::sin(lon - p1.lon);
    return std::atan(num / denom);
}

// ============================================================================
// Lat/lon given radial and distance (#LL)
// ============================================================================

/**
 * @brief Short-range formula (valid when |dlon| < pi/2) for destination point given initial point,
 * true course `tc` (radians), and angular distance `d` (radians).
 */
template <typename T = double>
[[nodiscard]] inline LatLon<T> point_from_radial_distance_short(const LatLon<T>& p1, T tc, T d,
                                                                T eps = constants::default_eps<T>) noexcept {
    const T lat = asin_safe(std::sin(p1.lat) * std::cos(d) + std::cos(p1.lat) * std::sin(d) * std::cos(tc));
    T lon{};
    if (std::abs(std::cos(lat)) < eps) {
        lon = p1.lon; // endpoint is a pole
    } else {
        lon = normalize_lon(p1.lon - asin_safe(std::sin(tc) * std::sin(d) / std::cos(lat)));
    }
    return LatLon<T>{lat, lon};
}

/**
 * @brief General formula (valid for any distance) for destination point given initial point,
 * true course `tc` (radians), and angular distance `d` (radians).
 */
template <typename T = double>
[[nodiscard]] inline LatLon<T> point_from_radial_distance(const LatLon<T>& p1, T tc, T d) noexcept {
    const T lat = asin_safe(std::sin(p1.lat) * std::cos(d) + std::cos(p1.lat) * std::sin(d) * std::cos(tc));
    const T dlon = std::atan2(std::sin(tc) * std::sin(d) * std::cos(p1.lat),
                              std::cos(d) - std::sin(p1.lat) * std::sin(lat));
    const T lon = normalize_lon(p1.lon - dlon);
    return LatLon<T>{lat, lon};
}

/**
 * @brief Convenience overload taking course in degrees and distance in nautical miles.
 */
template <typename T = double>
[[nodiscard]] inline LatLon<T> point_from_radial_distance_nm(const LatLon<T>& p1, T tc_deg, T d_nm) noexcept {
    return point_from_radial_distance(p1, units::deg_to_rad(tc_deg), units::nm_to_rad(d_nm));
}

// ============================================================================
// Intersection of two radials (#Intersection)
// ============================================================================

enum class IntersectionStatus {
    Unique,
    Infinite,
    Ambiguous
};

template <typename T = double>
struct RadialIntersectionResult {
    IntersectionStatus status{IntersectionStatus::Ambiguous};
    std::optional<LatLon<T>> point{};
    T dst12{};
    T crs12{};
    T crs21{};
    T ang1{};
    T ang2{};
    T ang3{};
    T dst13{};
};

/**
 * @brief Computes the intersection of two great circle radials: `crs13` from `p1` and `crs23` from `p2`.
 */
template <typename T = double>
[[nodiscard]] inline RadialIntersectionResult<T> intersect_radials(const LatLon<T>& p1, T crs13,
                                                                   const LatLon<T>& p2, T crs23,
                                                                   T eps = constants::default_eps<T>) noexcept {
    RadialIntersectionResult<T> res{};
    res.dst12 = distance(p1, p2);
    if (res.dst12 < eps) {
        res.status = IntersectionStatus::Infinite;
        res.point = p1;
        return res;
    }

    const T arg12 = (std::sin(p2.lat) - std::sin(p1.lat) * std::cos(res.dst12)) /
                    (std::sin(res.dst12) * std::cos(p1.lat));
    const T arg21 = (std::sin(p1.lat) - std::sin(p2.lat) * std::cos(res.dst12)) /
                    (std::sin(res.dst12) * std::cos(p2.lat));

    if (std::sin(p2.lon - p1.lon) < static_cast<T>(0)) {
        res.crs12 = acos_safe(arg12);
        res.crs21 = constants::two_pi<T> - acos_safe(arg21);
    } else {
        res.crs12 = constants::two_pi<T> - acos_safe(arg12);
        res.crs21 = acos_safe(arg21);
    }

    res.ang1 = normalize_lon(crs13 - res.crs12);
    res.ang2 = normalize_lon(res.crs21 - crs23);

    if (std::abs(std::sin(res.ang1)) < eps && std::abs(std::sin(res.ang2)) < eps) {
        res.status = IntersectionStatus::Infinite;
        return res;
    }
    if (std::sin(res.ang1) * std::sin(res.ang2) < static_cast<T>(0)) {
        res.status = IntersectionStatus::Ambiguous;
        return res;
    }

    res.ang1 = std::abs(res.ang1);
    res.ang2 = std::abs(res.ang2);
    res.ang3 = acos_safe(-std::cos(res.ang1) * std::cos(res.ang2) +
                         std::sin(res.ang1) * std::sin(res.ang2) * std::cos(res.dst12));
    res.dst13 = std::atan2(std::sin(res.dst12) * std::sin(res.ang1) * std::sin(res.ang2),
                           std::cos(res.ang2) + std::cos(res.ang1) * std::cos(res.ang3));

    res.point = point_from_radial_distance(p1, crs13, res.dst13);
    res.status = IntersectionStatus::Unique;
    return res;
}

// ============================================================================
// Clairaut's formula (#Clairaut)
// ============================================================================

/**
 * @brief Clairaut's invariant along any great circle: sin(tc) * cos(lat) = constant.
 */
template <typename T = double>
[[nodiscard]] inline T clairaut_constant(T lat, T tc) noexcept {
    return std::sin(tc) * std::cos(lat);
}

/**
 * @brief Maximum latitude (in radians) reached at the apex of a great circle passing through `lat` with course `tc`.
 * Formula: latmx = acos(abs(sin(tc)*cos(lat)))
 */
template <typename T = double>
[[nodiscard]] inline T clairaut_max_latitude(T lat, T tc) noexcept {
    return acos_safe(std::abs(std::sin(tc) * std::cos(lat)));
}

// ============================================================================
// GC crossing parallel (#Par)
// ============================================================================

/**
 * @brief Finds the two longitudes (in radians, {lon3_1, lon3_2}) where the great circle
 * through p1 and p2 crosses the parallel `lat3`.
 * Returns std::nullopt if the great circle does not reach `lat3`.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<std::pair<T, T>> crossing_parallels(const LatLon<T>& p1, const LatLon<T>& p2,
                                                                       T lat3,
                                                                       T eps = constants::default_eps<T>) noexcept {
    const T l12 = p1.lon - p2.lon;
    const T A = std::sin(p1.lat) * std::cos(p2.lat) * std::cos(lat3) * std::sin(l12);
    const T B = std::sin(p1.lat) * std::cos(p2.lat) * std::cos(lat3) * std::cos(l12) -
                std::cos(p1.lat) * std::sin(p2.lat) * std::cos(lat3);
    const T C = std::cos(p1.lat) * std::cos(p2.lat) * std::sin(lat3) * std::sin(l12);
    const T lon = std::atan2(B, A);
    const T norm_AB = std::hypot(A, B);

    if (norm_AB < eps || std::abs(C) > norm_AB) {
        return std::nullopt;
    }

    const T dlon = acos_safe(C / norm_AB);
    const T lon3_1 = normalize_lon(p1.lon + dlon + lon);
    const T lon3_2 = normalize_lon(p1.lon - dlon + lon);
    return std::make_pair(lon3_1, lon3_2);
}

// ============================================================================
// Intermediate points on a great circle (#Intermediate)
// ============================================================================

/**
 * @brief Finds the intermediate point a fraction `f` (f=0 is p1, f=1 is p2) along the great circle route from p1 to p2.
 * Returns std::nullopt if p1 and p2 are antipodal.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<LatLon<T>> intermediate_point(const LatLon<T>& p1, const LatLon<T>& p2, T f,
                                                                 T eps = constants::default_eps<T>) noexcept {
    const T d = distance(p1, p2);
    if (std::abs(d - constants::pi<T>) < eps) {
        return std::nullopt; // Antipodal points: route is undefined
    }
    if (d < eps) {
        return p1;
    }

    const T sin_d = std::sin(d);
    const T A = std::sin((static_cast<T>(1) - f) * d) / sin_d;
    const T B = std::sin(f * d) / sin_d;

    const T x = A * std::cos(p1.lat) * std::cos(p1.lon) + B * std::cos(p2.lat) * std::cos(p2.lon);
    const T y = A * std::cos(p1.lat) * std::sin(p1.lon) + B * std::cos(p2.lat) * std::sin(p2.lon);
    const T z = A * std::sin(p1.lat) + B * std::sin(p2.lat);

    const T lat = std::atan2(z, std::hypot(x, y));
    const T lon = std::atan2(y, x);
    return LatLon<T>{lat, lon};
}

// ============================================================================
// Cross track error and along track distance (#XTE)
// ============================================================================

template <typename T = double>
struct TrackErrorResult {
    T xtd{};     ///< Cross-track error in radians (positive = right of course, negative = left)
    T atd{};     ///< Along-track distance in radians
    T dist_ad{}; ///< Great circle distance from A to D in radians
    T crs_ab{};  ///< Initial course from A to B in radians
    T crs_ad{};  ///< Initial course from A to D in radians
};

/**
 * @brief Computes cross-track error (XTD) in radians given distance AD, course AD, and course AB.
 * Positive XTD means right of course; negative means left.
 */
template <typename T = double>
[[nodiscard]] inline T cross_track_error(T dist_ad, T crs_ad, T crs_ab) noexcept {
    return asin_safe(std::sin(dist_ad) * std::sin(crs_ad - crs_ab));
}

/**
 * @brief Computes cross-track error (XTD) in radians for point D relative to great circle A -> B.
 * Handles the special case where A is the North or South Pole.
 */
template <typename T = double>
[[nodiscard]] inline T cross_track_error(const LatLon<T>& pA, const LatLon<T>& pB, const LatLon<T>& pD,
                                         T eps = constants::default_eps<T>) noexcept {
    const T dist_ad = distance(pA, pD);
    if (std::cos(pA.lat) < eps) {
        const T angle_diff = (pA.lat > static_cast<T>(0)) ? (pD.lon - pB.lon) : (pB.lon - pD.lon);
        return asin_safe(std::sin(dist_ad) * std::sin(angle_diff));
    }
    const T crs_ab = initial_course(pA, pB, eps);
    const T crs_ad = initial_course(pA, pD, eps);
    return cross_track_error(dist_ad, crs_ad, crs_ab);
}

/**
 * @brief Along-track distance (ATD) in radians using the standard acos formula:
 * ATD = acos(cos(dist_AD) / cos(XTD))
 */
template <typename T = double>
[[nodiscard]] inline T along_track_distance_cos(T dist_ad, T xtd) noexcept {
    return acos_safe(std::cos(dist_ad) / std::cos(xtd));
}

/**
 * @brief Along-track distance (ATD) in radians using the short-distance formula less susceptible to rounding error:
 * ATD = asin(sqrt(sin(dist_AD)^2 - sin(XTD)^2) / cos(XTD))
 */
template <typename T = double>
[[nodiscard]] inline T along_track_distance(T dist_ad, T xtd) noexcept {
    const T sin_ad = std::sin(dist_ad);
    const T sin_xtd = std::sin(xtd);
    const T diff = std::max(static_cast<T>(0), sin_ad * sin_ad - sin_xtd * sin_xtd);
    return asin_safe(std::sqrt(diff) / std::cos(xtd));
}

/**
 * @brief Computes both cross-track error (XTD) and along-track distance (ATD) for point D relative to A -> B.
 */
template <typename T = double>
[[nodiscard]] inline TrackErrorResult<T> cross_and_along_track(const LatLon<T>& pA, const LatLon<T>& pB,
                                                               const LatLon<T>& pD,
                                                               T eps = constants::default_eps<T>) noexcept {
    TrackErrorResult<T> res{};
    res.dist_ad = distance(pA, pD);
    res.crs_ab = initial_course(pA, pB, eps);
    res.crs_ad = initial_course(pA, pD, eps);
    res.xtd = cross_track_error( pA, pB, pD, eps );
    res.atd = along_track_distance(res.dist_ad, res.xtd);
    return res;
}

// ============================================================================
// Point(s) known distance from a great circle (#POINTDME)
// ============================================================================

template <typename T = double>
struct PointsAtDistanceResult {
    LatLon<T> point1{};
    LatLon<T> point2{};
    T dp1{}; ///< Distance of point1 from A along AB (radians)
    T dp2{}; ///< Distance of point2 from A along AB (radians)
};

/**
 * @brief Finds the points on the great circle through A and B that lie an angular distance `d` from point D, if they exist.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<PointsAtDistanceResult<T>>
points_at_distance_from_gc(const LatLon<T>& pA, const LatLon<T>& pB, const LatLon<T>& pD, T d,
                           T eps = constants::default_eps<T>) noexcept {
    const T crs_ab = initial_course(pA, pB, eps);
    const T crs_ad = initial_course(pA, pD, eps);
    const T A = crs_ad - crs_ab;
    const T b = distance(pA, pD);

    const T cos_b = std::cos(b);
    const T sin_b_cos_A = std::sin(b) * std::cos(A);
    const T r_sq = cos_b * cos_b + sin_b_cos_A * sin_b_cos_A;
    const T cos_d = std::cos(d);

    if (cos_d * cos_d > r_sq || r_sq < eps) {
        return std::nullopt;
    }

    const T r = std::sqrt(r_sq);
    const T p = std::atan2(sin_b_cos_A, cos_b);
    const T delta_p = acos_safe(cos_d / r);

    PointsAtDistanceResult<T> res{};
    res.dp1 = p + delta_p;
    res.dp2 = p - delta_p;
    res.point1 = point_from_radial_distance(pA, crs_ab, res.dp1);
    res.point2 = point_from_radial_distance(pA, crs_ab, res.dp2);
    return res;
}

} // namespace avform::great_circle

#endif // AVFORM_GREAT_CIRCLE_HPP

