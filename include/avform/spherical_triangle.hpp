#ifndef AVFORM_SPHERICAL_TRIANGLE_HPP
#define AVFORM_SPHERICAL_TRIANGLE_HPP

#include "core.hpp"
#include <array>

namespace avform::spherical_triangle {

/**
 * @brief Represents a spherical triangle with side lengths {a, b, c} (in radians)
 * and opposite interior angles {A, B, C} (in radians).
 */
template <typename T = double>
struct SphericalTriangle {
    T a{};
    T b{};
    T c{};
    T A{};
    T B{};
    T C{};
};

template <typename T = double>
struct AmbiguousTriangleSolutions {
    std::size_t count{0};
    std::array<SphericalTriangle<T>, 2> triangles{};
};

// ============================================================================
// General Spherical Triangle Solvers (#Triangle)
// ============================================================================

/**
 * @brief Given {A, b, c}: Two sides and included angle (SAS).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<SphericalTriangle<T>> solve_sas(T A, T b, T c,
                                                                   T eps = constants::default_eps<T>) noexcept {
    const T a = acos_safe(std::cos(b) * std::cos(c) + std::sin(b) * std::sin(c) * std::cos(A));
    if (std::sin(a) < eps || std::sin(b) < eps || std::sin(c) < eps) {
        return std::nullopt;
    }
    const T B = acos_safe((std::cos(b) - std::cos(c) * std::cos(a)) / (std::sin(c) * std::sin(a)));
    const T C = acos_safe((std::cos(c) - std::cos(a) * std::cos(b)) / (std::sin(a) * std::sin(b)));
    return SphericalTriangle<T>{a, b, c, A, B, C};
}

/**
 * @brief Given {a, B, C}: Two angles and included side (ASA).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<SphericalTriangle<T>> solve_asa(T a, T B, T C) noexcept {
    const T A = acos_safe(-std::cos(B) * std::cos(C) + std::sin(B) * std::sin(C) * std::cos(a));
    const T num = std::sin(a) * std::sin(B) * std::sin(C);
    const T b = std::atan2(num, std::cos(B) + std::cos(C) * std::cos(A));
    const T c = std::atan2(num, std::cos(C) + std::cos(A) * std::cos(B));
    if (b <= static_cast<T>(0) || c <= static_cast<T>(0)) {
        return std::nullopt;
    }
    return SphericalTriangle<T>{a, b, c, A, B, C};
}

/**
 * @brief Given {a, b, c}: Three sides (SSS).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<SphericalTriangle<T>> solve_sss(T a, T b, T c,
                                                                   T eps = constants::default_eps<T>) noexcept {
    if (std::sin(a) < eps || std::sin(b) < eps || std::sin(c) < eps) {
        return std::nullopt;
    }
    // Check spherical triangle inequality
    if (a + b <= c || a + c <= b || b + c <= a || (a + b + c) >= constants::two_pi<T>) {
        return std::nullopt;
    }
    const T A = acos_safe((std::cos(a) - std::cos(b) * std::cos(c)) / (std::sin(b) * std::sin(c)));
    const T B = acos_safe((std::cos(b) - std::cos(c) * std::cos(a)) / (std::sin(c) * std::sin(a)));
    const T C = acos_safe((std::cos(c) - std::cos(a) * std::cos(b)) / (std::sin(a) * std::sin(b)));
    return SphericalTriangle<T>{a, b, c, A, B, C};
}

/**
 * @brief Given {A, B, C}: Three angles (AAA).
 */
template <typename T = double>
[[nodiscard]] inline std::optional<SphericalTriangle<T>> solve_aaa(T A, T B, T C,
                                                                   T eps = constants::default_eps<T>) noexcept {
    const T delta = (A + B + C - constants::pi<T>) / static_cast<T>(2);
    if (delta <= static_cast<T>(0) || std::sin(A) < eps || std::sin(B) < eps || std::sin(C) < eps) {
        return std::nullopt;
    }
    const T term_a = std::sin(delta) * std::sin(A - delta) / (std::sin(B) * std::sin(C));
    const T term_b = std::sin(delta) * std::sin(B - delta) / (std::sin(C) * std::sin(A));
    const T term_c = std::sin(delta) * std::sin(C - delta) / (std::sin(A) * std::sin(B));
    if (term_a < static_cast<T>(0) || term_b < static_cast<T>(0) || term_c < static_cast<T>(0)) {
        return std::nullopt;
    }
    const T a = static_cast<T>(2) * asin_safe(std::sqrt(term_a));
    const T b = static_cast<T>(2) * asin_safe(std::sqrt(term_b));
    const T c = static_cast<T>(2) * asin_safe(std::sqrt(term_c));
    return SphericalTriangle<T>{a, b, c, A, B, C};
}

namespace detail {

template <typename T>
[[nodiscard]] inline SphericalTriangle<T> complete_from_a_b_A_B(T a, T b, T A, T B) noexcept {
    const T half_AB_sum = (A + B) / static_cast<T>(2);
    const T half_AB_diff = (A - B) / static_cast<T>(2);
    const T half_ab_sum = (a + b) / static_cast<T>(2);
    const T half_ab_diff = (a - b) / static_cast<T>(2);

    const T c = mod(static_cast<T>(2) * std::atan2(std::cos(half_AB_sum) * std::sin(half_ab_sum),
                                                   std::cos(half_AB_diff) * std::cos(half_ab_sum)),
                    constants::two_pi<T>);
    const T C = mod(static_cast<T>(2) * std::atan2(std::cos(half_ab_diff) * std::cos(half_AB_sum),
                                                   std::cos(half_ab_sum) * std::sin(half_AB_sum)),
                    constants::two_pi<T>);
    return SphericalTriangle<T>{a, b, c, A, B, C};
}

} // namespace detail

/**
 * @brief Given {A, a, b}: Two sides and non-included angle (SSA).
 * May yield 0, 1, or 2 valid spherical triangles in (0, pi).
 */
template <typename T = double>
[[nodiscard]] inline AmbiguousTriangleSolutions<T> solve_ssa(T A, T a, T b,
                                                             T eps = constants::default_eps<T>) noexcept {
    AmbiguousTriangleSolutions<T> res{};
    if (std::sin(a) < eps) {
        return res;
    }
    const T x = std::sin(A) * std::sin(b) / std::sin(a);
    if (x > static_cast<T>(1) + eps) {
        return res;
    }

    auto try_add = [&](T B_candidate) {
        const auto tri = detail::complete_from_a_b_A_B(a, b, A, B_candidate);
        if (tri.c > eps && tri.c < constants::pi<T> - eps &&
            tri.C > eps && tri.C < constants::pi<T> - eps) {
            res.triangles[res.count++] = tri;
        }
    };

    if (std::abs(x - static_cast<T>(1)) <= eps) {
        try_add(constants::half_pi<T>);
    } else if (x < static_cast<T>(1)) {
        const T B1 = asin_safe(x);
        const T B2 = constants::pi<T> - B1;
        try_add(B1);
        if (std::abs(B2 - B1) > eps) {
            try_add(B2);
        }
    }
    return res;
}

/**
 * @brief Given {a, A, B}: Two angles and non-included side (AAS).
 * May yield 0, 1, or 2 valid spherical triangles in (0, pi).
 */
template <typename T = double>
[[nodiscard]] inline AmbiguousTriangleSolutions<T> solve_aas(T a, T A, T B,
                                                             T eps = constants::default_eps<T>) noexcept {
    AmbiguousTriangleSolutions<T> res{};
    if (std::sin(A) < eps) {
        return res;
    }
    const T x = std::sin(a) * std::sin(B) / std::sin(A);
    if (x > static_cast<T>(1) + eps) {
        return res;
    }

    auto try_add = [&](T b_candidate) {
        const auto tri = detail::complete_from_a_b_A_B(a, b_candidate, A, B);
        if (tri.c > eps && tri.c < constants::pi<T> - eps &&
            tri.C > eps && tri.C < constants::pi<T> - eps) {
            res.triangles[res.count++] = tri;
        }
    };

    if (std::abs(x - static_cast<T>(1)) <= eps) {
        try_add(constants::half_pi<T>);
    } else if (x < static_cast<T>(1)) {
        const T b1 = asin_safe(x);
        const T b2 = constants::pi<T> - b1;
        try_add(b1);
        if (std::abs(b2 - b1) > eps) {
            try_add(b2);
        }
    }
    return res;
}

// ============================================================================
// Spherical Excess and Surface Area
// ============================================================================

/**
 * @brief Spherical excess E = A + B + C - pi (in radians / steradians).
 */
template <typename T = double>
[[nodiscard]] constexpr T spherical_excess_angles(T A, T B, T C) noexcept {
    return A + B + C - constants::pi<T>;
}

/**
 * @brief L'Huilier's formula for spherical excess E in terms of the three sides {a, b, c}:
 * E = 4 * atan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)))
 * where s = (a + b + c) / 2.
 */
template <typename T = double>
[[nodiscard]] inline T spherical_excess_sides(T a, T b, T c) noexcept {
    const T s = (a + b + c) / static_cast<T>(2);
    const T prod = std::tan(s / static_cast<T>(2)) *
                   std::tan((s - a) / static_cast<T>(2)) *
                   std::tan((s - b) / static_cast<T>(2)) *
                   std::tan((s - c) / static_cast<T>(2));
    return static_cast<T>(4) * std::atan(std::sqrt(std::max(static_cast<T>(0), prod)));
}

/**
 * @brief Surface area enclosed by a spherical triangle given spherical excess E and sphere radius R:
 * Area = E * R^2
 */
template <typename T = double>
[[nodiscard]] constexpr T spherical_area(T excess_rad, T radius = static_cast<T>(1)) noexcept {
    return excess_rad * radius * radius;
}

// ============================================================================
// Right Spherical Triangles (Angle B = pi / 2) — Napier's Rules
// ============================================================================

namespace right_triangle {

/**
 * @brief Right spherical triangle with B = pi/2, Given legs {a, c}.
 * Formulas:
 *   cos(b) = cos(a) * cos(c)
 *   tan(A) = tan(a) / sin(c)
 *   tan(C) = tan(c) / sin(a)
 */
template <typename T = double>
[[nodiscard]] inline SphericalTriangle<T> solve_from_a_c(T a, T c) noexcept {
    const T b = acos_safe(std::cos(a) * std::cos(c));
    const T A = std::atan2(std::tan(a), std::sin(c));
    const T C = std::atan2(std::tan(c), std::sin(a));
    return SphericalTriangle<T>{a, b, c, A, constants::half_pi<T>, C};
}

/**
 * @brief Right spherical triangle with B = pi/2, Given leg `a` and hypotenuse `b`.
 * Formulas:
 *   cos(c) = cos(b) / cos(a)
 *   sin(A) = sin(a) / sin(b)
 *   cos(C) = tan(a) / tan(b)
 */
template <typename T = double>
[[nodiscard]] inline SphericalTriangle<T> solve_from_a_b(T a, T b) noexcept {
    const T c = acos_safe(std::cos(b) / std::cos(a));
    const T A = asin_safe(std::sin(a) / std::sin(b));
    const T C = acos_safe(std::tan(a) / std::tan(b));
    return SphericalTriangle<T>{a, b, c, A, constants::half_pi<T>, C};
}

/**
 * @brief Right spherical triangle with B = pi/2, Given hypotenuse `b` and leg `c`.
 */
template <typename T = double>
[[nodiscard]] inline SphericalTriangle<T> solve_from_b_c(T b, T c) noexcept {
    const T a = acos_safe(std::cos(b) / std::cos(c));
    const T C = asin_safe(std::sin(c) / std::sin(b));
    const T A = acos_safe(std::tan(c) / std::tan(b));
    return SphericalTriangle<T>{a, b, c, A, constants::half_pi<T>, C};
}

/**
 * @brief Right spherical triangle with B = pi/2, Given hypotenuse `b` and angle `A`.
 * Formulas:
 *   sin(a) = sin(A) * sin(b)
 *   tan(c) = tan(b) * cos(A)
 *   tan(C) = 1 / (cos(b) * tan(A))
 */
template <typename T = double>
[[nodiscard]] inline SphericalTriangle<T> solve_from_b_A(T b, T A) noexcept {
    const T a = asin_safe(std::sin(A) * std::sin(b));
    const T c = std::atan(std::tan(b) * std::cos(A));
    const T C = std::atan2(static_cast<T>(1), std::cos(b) * std::tan(A));
    return SphericalTriangle<T>{a, b, (c < 0 ? c + constants::pi<T> : c), A, constants::half_pi<T>, C};
}

/**
 * @brief Right spherical triangle with B = pi/2, Given angles `A` and `C`.
 * Formulas:
 *   cos(a) = cos(A) / sin(C)
 *   cos(c) = cos(C) / sin(A)
 *   cos(b) = 1 / (tan(A) * tan(C))
 */
template <typename T = double>
[[nodiscard]] inline SphericalTriangle<T> solve_from_A_C(T A, T C) noexcept {
    const T a = acos_safe(std::cos(A) / std::sin(C));
    const T c = acos_safe(std::cos(C) / std::sin(A));
    const T b = acos_safe(static_cast<T>(1) / (std::tan(A) * std::tan(C)));
    return SphericalTriangle<T>{a, b, c, A, constants::half_pi<T>, C};
}

} // namespace right_triangle

} // namespace avform::spherical_triangle

#endif // AVFORM_SPHERICAL_TRIANGLE_HPP

