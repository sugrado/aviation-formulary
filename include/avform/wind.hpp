#ifndef AVFORM_WIND_HPP
#define AVFORM_WIND_HPP

#include "core.hpp"
#include "units.hpp"

namespace avform::wind {

// ============================================================================
// Wind Triangles (#Wind)
// ============================================================================

template <typename T = double>
struct WindVector {
    T ws{}; ///< Wind speed (same speed unit as TAS and GS)
    T wd{}; ///< Wind direction (from) in radians [0, 2*pi)

    [[nodiscard]] constexpr T wd_deg() const noexcept {
        return units::rad_to_deg(wd);
    }
};

template <typename T = double>
struct HeadingGroundspeed {
    T hd{};  ///< Required true heading in radians [0, 2*pi)
    T gs{};  ///< Resulting groundspeed
    T swc{}; ///< Sine of wind correction angle

    [[nodiscard]] constexpr T hd_deg() const noexcept {
        return units::rad_to_deg(hd);
    }
};

template <typename T = double>
struct CourseGroundspeed {
    T crs{}; ///< Resulting true course in radians [0, 2*pi)
    T gs{};  ///< Resulting groundspeed
    T wca{}; ///< Wind correction angle in radians

    [[nodiscard]] constexpr T crs_deg() const noexcept {
        return units::rad_to_deg(crs);
    }

    [[nodiscard]] constexpr T wca_deg() const noexcept {
        return units::rad_to_deg(wca);
    }
};

/**
 * @brief (1) Unknown Wind: Computes wind speed (WS) and wind direction from (WD, in radians)
 * given course `crs` (rad), heading `hd` (rad), true airspeed `tas`, and groundspeed `gs`.
 */
template <typename T = double>
[[nodiscard]] inline WindVector<T> unknown_wind(T crs, T hd, T tas, T gs) noexcept {
    const T diff = hd - crs;
    const T sin_half = std::sin(diff / static_cast<T>(2));
    const T ws = std::sqrt((tas - gs) * (tas - gs) +
                           static_cast<T>(4) * tas * gs * sin_half * sin_half);
    const T wd = mod(crs + std::atan2(tas * std::sin(diff), tas * std::cos(diff) - gs),
                     constants::two_pi<T>);
    return WindVector<T>{ws, wd};
}

/**
 * @brief Convenience overload taking course and heading in degrees.
 */
template <typename T = double>
[[nodiscard]] inline WindVector<T> unknown_wind_deg(T crs_deg, T hd_deg, T tas, T gs) noexcept {
    return unknown_wind(units::deg_to_rad(crs_deg), units::deg_to_rad(hd_deg), tas, gs);
}

/**
 * @brief (2) Find HD, GS: Computes required heading (HD, in radians) and groundspeed (GS)
 * given desired course `crs` (rad), true airspeed `tas`, wind speed `ws`, and wind direction from `wd` (rad).
 * Returns std::nullopt if the course cannot be flown ("wind too strong", |SWC| > 1 or GS < 0).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<HeadingGroundspeed<T>>
heading_and_groundspeed(T crs, T tas, T ws, T wd) noexcept {
    if (tas <= static_cast<T>(0)) {
        return std::nullopt;
    }
    const T swc = (ws / tas) * std::sin(wd - crs);
    if (std::abs(swc) > static_cast<T>(1)) {
        return std::nullopt; // Course cannot be flown -- wind too strong
    }
    const T hd = mod(crs + std::asin(swc), constants::two_pi<T>);
    const T gs = tas * std::sqrt(static_cast<T>(1) - swc * swc) - ws * std::cos(wd - crs);
    if (gs < static_cast<T>(0)) {
        return std::nullopt; // Course cannot be flown -- wind too strong
    }
    return HeadingGroundspeed<T>{hd, gs, swc};
}

/**
 * @brief Convenience overload taking course and wind direction in degrees.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<HeadingGroundspeed<T>>
heading_and_groundspeed_deg(T crs_deg, T tas, T ws, T wd_deg) noexcept {
    return heading_and_groundspeed(units::deg_to_rad(crs_deg), tas, ws, units::deg_to_rad(wd_deg));
}

/**
 * @brief (3) Find CRS, GS: Computes resulting course (CRS, in radians), groundspeed (GS),
 * and wind correction angle (WCA, in radians) given heading `hd` (rad), true airspeed `tas`,
 * wind speed `ws`, and wind direction from `wd` (rad).
 */
template <typename T = double>
[[nodiscard]] inline CourseGroundspeed<T> course_and_groundspeed(T hd, T tas, T ws, T wd) noexcept {
    const T diff = hd - wd;
    const T gs = std::sqrt(ws * ws + tas * tas - static_cast<T>(2) * ws * tas * std::cos(diff));
    const T wca = std::atan2(ws * std::sin(diff), tas - ws * std::cos(diff));
    const T crs = mod(hd + wca, constants::two_pi<T>);
    return CourseGroundspeed<T>{crs, gs, wca};
}

/**
 * @brief Convenience overload taking heading and wind direction in degrees.
 */
template <typename T = double>
[[nodiscard]] inline CourseGroundspeed<T> course_and_groundspeed_deg(T hd_deg, T tas, T ws, T wd_deg) noexcept {
    return course_and_groundspeed(units::deg_to_rad(hd_deg), tas, ws, units::deg_to_rad(wd_deg));
}

// ============================================================================
// Head- and cross-wind components (#Xwind)
// ============================================================================

template <typename T = double>
struct WindComponents {
    T headwind{};  ///< Headwind component (tailwind is negative)
    T crosswind{}; ///< Crosswind component (positive = wind from right, negative = from left)
};

/**
 * @brief Computes headwind and crosswind components given wind speed `ws`,
 * wind direction from `wd` (radians), and runway direction `rd` (radians).
 * Formulas:
 *   HW = WS * cos(WD - RD)   (tailwind negative)
 *   XW = WS * sin(WD - RD)   (positive = wind from right)
 */
template <typename T = double>
[[nodiscard]] inline WindComponents<T> head_and_cross_wind(T ws, T wd, T rd) noexcept {
    return WindComponents<T>{
        ws * std::cos(wd - rd),
        ws * std::sin(wd - rd)
    };
}

/**
 * @brief Convenience overload taking wind direction and runway direction in degrees.
 */
template <typename T = double>
[[nodiscard]] inline WindComponents<T> head_and_cross_wind_deg(T ws, T wd_deg, T rd_deg) noexcept {
    return head_and_cross_wind(ws, units::deg_to_rad(wd_deg), units::deg_to_rad(rd_deg));
}

// ============================================================================
// TAS and windspeed from three (GPS) groundspeeds (#3gs)
// ============================================================================

template <typename T = double>
struct ThreeGroundspeedsResult {
    T tas{}; ///< True airspeed (when TAS > WS; exchanged if WS > TAS)
    T ws{};  ///< Wind speed
    T vms{};
    T mu{};
    T bp{};
    T bm{};
};

/**
 * @brief Computes True Airspeed (TAS) and windspeed (WS) from three groundspeeds (v1, v2, v3)
 * measured on headings differing by 120 degrees (e.g., 40, 160, and 280 degrees).
 * Provided that TAS exceeds windspeed, `tas = sqrt(vms * bp)` and `ws = sqrt(vms * bm)`.
 * If `tas_exceeds_wind` is false, the two roots are exchanged.
 */
template <typename T = double>
[[nodiscard]] inline std::optional<ThreeGroundspeedsResult<T>>
tas_and_wind_from_three_gs(T v1, T v2, T v3, bool tas_exceeds_wind = true) noexcept {
    const T v1_sq = v1 * v1;
    const T v2_sq = v2 * v2;
    const T v3_sq = v3 * v3;
    const T vms = (v1_sq + v2_sq + v3_sq) / static_cast<T>(3);
    if (vms <= static_cast<T>(0)) {
        return std::nullopt;
    }

    const T a1 = v1_sq / vms - static_cast<T>(1);
    const T a2 = v2_sq / vms - static_cast<T>(1);
    const T a3 = v3_sq / vms - static_cast<T>(1);
    const T mu = (a1 * a1 + a2 * a2 + a3 * a3) / static_cast<T>(6);

    const T disc = static_cast<T>(0.25) - mu;
    if (disc < static_cast<T>(0)) {
        return std::nullopt; // Physically inconsistent groundspeeds for a constant wind
    }

    const T bp = static_cast<T>(0.5) + std::sqrt(disc);
    const T bm = mu / bp;

    T root_high = std::sqrt(vms * bp);
    T root_low = std::sqrt(std::max(static_cast<T>(0), vms * bm));

    ThreeGroundspeedsResult<T> res{};
    res.vms = vms;
    res.mu = mu;
    res.bp = bp;
    res.bm = bm;
    res.tas = tas_exceeds_wind ? root_high : root_low;
    res.ws = tas_exceeds_wind ? root_low : root_high;
    return res;
}

// ============================================================================
// Bellamy's formula (#Bellamy)
// ============================================================================

template <typename T = double>
struct BellamyResult {
    T drift_nm{}; ///< Wind drift in nautical miles (positive when p2 > p1 -> drift to the left, WCA to the right)
    T wca_deg{};  ///< Single wind correction angle in degrees
};

/**
 * @brief Bellamy's formula for wind drift (in nm) using pressure difference `p2 - p1` in inches of Hg.
 * `lat` is average route latitude in radians, `tas_kt` is true airspeed in knots.
 * Formula: Drift (nm) = 21500 * (p2 - p1) / (sin(latitude) * TAS)
 */
template <typename T = double>
[[nodiscard]] inline T bellamy_drift_inhg(T p1_inhg, T p2_inhg, T lat_rad, T tas_kt) noexcept {
    return static_cast<T>(21500) * (p2_inhg - p1_inhg) / (std::sin(lat_rad) * tas_kt);
}

/**
 * @brief Bellamy's formula for wind drift (in nm) using pressure difference `p2 - p1` in millibars (hPa).
 * Formula: Drift (nm) = 635 * (p2 - p1) / (sin(latitude) * TAS)
 */
template <typename T = double>
[[nodiscard]] inline T bellamy_drift_mbar(T p1_mbar, T p2_mbar, T lat_rad, T tas_kt) noexcept {
    return static_cast<T>(635) * (p2_mbar - p1_mbar) / (std::sin(lat_rad) * tas_kt);
}

/**
 * @brief Bellamy's formula for wind correction angle (in degrees) using pressure in inches of Hg.
 * Formula: WCA = 1230000 * (p2 - p1) / (sin(latitude) * TAS * Dist)
 */
template <typename T = double>
[[nodiscard]] inline T bellamy_wca_inhg(T p1_inhg, T p2_inhg, T lat_rad, T tas_kt, T dist_nm) noexcept {
    return static_cast<T>(1230000) * (p2_inhg - p1_inhg) / (std::sin(lat_rad) * tas_kt * dist_nm);
}

/**
 * @brief Bellamy's formula for wind correction angle (in degrees) using pressure in millibars (hPa).
 * Formula: WCA = 36300 * (p2 - p1) / (sin(latitude) * TAS * Dist)
 */
template <typename T = double>
[[nodiscard]] inline T bellamy_wca_mbar(T p1_mbar, T p2_mbar, T lat_rad, T tas_kt, T dist_nm) noexcept {
    return static_cast<T>(36300) * (p2_mbar - p1_mbar) / (std::sin(lat_rad) * tas_kt * dist_nm);
}

/**
 * @brief Computes both Bellamy drift (nm) and WCA (deg) given latitude in degrees and pressures in inHg.
 */
template <typename T = double>
[[nodiscard]] inline BellamyResult<T> bellamy_inhg_deg(T p1_inhg, T p2_inhg, T lat_deg, T tas_kt, T dist_nm) noexcept {
    const T lat_rad = units::deg_to_rad(lat_deg);
    return BellamyResult<T>{
        bellamy_drift_inhg(p1_inhg, p2_inhg, lat_rad, tas_kt),
        bellamy_wca_inhg(p1_inhg, p2_inhg, lat_rad, tas_kt, dist_nm)
    };
}

/**
 * @brief Computes both Bellamy drift (nm) and WCA (deg) given latitude in degrees and pressures in mbar.
 */
template <typename T = double>
[[nodiscard]] inline BellamyResult<T> bellamy_mbar_deg(T p1_mbar, T p2_mbar, T lat_deg, T tas_kt, T dist_nm) noexcept {
    const T lat_rad = units::deg_to_rad(lat_deg);
    return BellamyResult<T>{
        bellamy_drift_mbar(p1_mbar, p2_mbar, lat_rad, tas_kt),
        bellamy_wca_mbar(p1_mbar, p2_mbar, lat_rad, tas_kt, dist_nm)
    };
}

} // namespace avform::wind

#endif // AVFORM_WIND_HPP

