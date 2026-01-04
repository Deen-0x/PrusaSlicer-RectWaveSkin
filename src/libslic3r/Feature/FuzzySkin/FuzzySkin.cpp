// src/libslic3r/Feature/FuzzySkin/FuzzySkin.cpp

#include <algorithm>
#include <cassert>
#include <cmath>

#include "libslic3r/Algorithm/LineSegmentation/LineSegmentation.hpp"
#include "libslic3r/Arachne/utils/ExtrusionJunction.hpp"
#include "libslic3r/Arachne/utils/ExtrusionLine.hpp"
#include "libslic3r/PerimeterGenerator.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/PrintConfig.hpp"

#include "FuzzySkin.hpp"

using namespace Slic3r;

namespace Slic3r::Feature::FuzzySkin {

// --------------------------
// Rect-wave “fuzz” parameters
// --------------------------

// Smooth square wave in [-1, +1] using tanh(sin()).
// sharpness ~ 6..20 (higher = more "rect", but harsher transitions).
static inline double rect_wave(double phase_radians, double sharpness = 10.0) {
    return std::tanh(sharpness * std::sin(phase_radians));
}

static inline double two_pi() { return 6.2831853071795864769; }

// Wavelength is derived from fuzzy_skin_point_distance to avoid adding new UI knobs.
// Increase to make longer waves; decrease for tighter waves.
static constexpr double k_wavelength_factor = 4.0;

// --------------------------
// Per-layer phase handling
// --------------------------

// Thread-local so parallel slicing won't cross-contaminate.
static thread_local double tl_phase_offset_radians = 0.0;

struct PhaseScope
{
    double prev = 0.0;
    explicit PhaseScope(size_t layer_idx) : prev(tl_phase_offset_radians) {
        // Half-cycle shift every other layer: 0, π, 0, π, ...
        tl_phase_offset_radians = (layer_idx & 1) ? (two_pi() * 0.5) : 0.0;
    }
    ~PhaseScope() { tl_phase_offset_radians = prev; }
};

// --------------------------
// Core implementation
// --------------------------

void fuzzy_polyline(
    Points &poly,
    const bool closed,
    const double fuzzy_skin_thickness,
    const double fuzzy_skin_point_distance
) {
    if (poly.size() < 2 || fuzzy_skin_thickness <= 0.0 || fuzzy_skin_point_distance <= 0.0)
        return;

    // Sample points along the polyline at a fixed step, but start at half-step to avoid hitting vertices.
    const double step = fuzzy_skin_point_distance;
    const double start_dist = step * 0.5;

    const double wavelength = std::max(step * k_wavelength_factor, 1e-9);

    double dist_left_over = start_dist;
    double total_dist = 0.0; // global arclength from the start of the polyline

    Points out;
    out.reserve(poly.size());

    // Skip the first point for open polyline.
    Point *p0 = closed ? &poly.back() : &poly.front();
    for (auto it_pt1 = closed ? poly.begin() : std::next(poly.begin()); it_pt1 != poly.end();
         ++it_pt1) {
        Point &p1 = *it_pt1;

        Vec2d p0p1 = (p1 - *p0).cast<double>();
        double p0p1_size = p0p1.norm();
        if (p0p1_size <= 0.0) {
            p0 = &p1;
            continue;
        }

        // Unit normal (perpendicular) for offset.
        Vec2d n = perp(p0p1).cast<double>().normalized();

        double p0pa_dist = dist_left_over;
        for (; p0pa_dist < p0p1_size; p0pa_dist += step) {
            const double s_global = total_dist + p0pa_dist;
            const double phase = two_pi() * (s_global / wavelength) + tl_phase_offset_radians;

            // Offset in [-thickness, +thickness]
            const double r = fuzzy_skin_thickness * rect_wave(phase, /*sharpness=*/10.0);

            out.emplace_back(*p0 + (p0p1 * (p0pa_dist / p0p1_size) + n * r).cast<coord_t>());
        }

        dist_left_over = p0pa_dist - p0p1_size;
        total_dist += p0p1_size;
        p0 = &p1;
    }

    while (out.size() < 3 && poly.size() >= 2) {
        size_t point_idx = poly.size() - 2;
        out.emplace_back(poly[point_idx]);
        if (point_idx == 0)
            break;
        --point_idx;
    }

    if (out.size() >= 3) {
        poly = std::move(out);
    }
}

void fuzzy_polygon(Polygon &polygon, double fuzzy_skin_thickness, double fuzzy_skin_point_distance) {
    fuzzy_polyline(polygon.points, true, fuzzy_skin_thickness, fuzzy_skin_point_distance);
}

void fuzzy_extrusion_line(
    Arachne::ExtrusionLine &ext_lines,
    const double fuzzy_skin_thickness,
    const double fuzzy_skin_point_distance
) {
    if (ext_lines.empty() || fuzzy_skin_thickness <= 0.0 || fuzzy_skin_point_distance <= 0.0)
        return;

    const double step = fuzzy_skin_point_distance;
    const double start_dist = step * 0.5;
    const double wavelength = std::max(step * k_wavelength_factor, 1e-9);

    double dist_left_over = start_dist;
    double total_dist = 0.0;

    Arachne::ExtrusionJunction *p0 = &ext_lines.front();
    Arachne::ExtrusionJunctions out;
    out.reserve(ext_lines.size());

    for (auto &p1 : ext_lines) {
        if (p0->p == p1.p) {
            // Copy the first point.
            out.emplace_back(p1.p, p1.w, p1.perimeter_index);
            continue;
        }

        Vec2d p0p1 = (p1.p - p0->p).cast<double>();
        double p0p1_size = p0p1.norm();
        if (p0p1_size <= 0.0) {
            p0 = &p1;
            continue;
        }

        Vec2d n = perp(p0p1).cast<double>().normalized();

        double p0pa_dist = dist_left_over;
        for (; p0pa_dist < p0p1_size; p0pa_dist += step) {
            const double s_global = total_dist + p0pa_dist;
            const double phase = two_pi() * (s_global / wavelength) + tl_phase_offset_radians;
            const double r = fuzzy_skin_thickness * rect_wave(phase, /*sharpness=*/10.0);

            out.emplace_back(
                p0->p + (p0p1 * (p0pa_dist / p0p1_size) + n * r).cast<coord_t>(), p1.w,
                p1.perimeter_index
            );
        }

        dist_left_over = p0pa_dist - p0p1_size;
        total_dist += p0p1_size;
        p0 = &p1;
    }

    while (out.size() < 3 && ext_lines.size() >= 2) {
        size_t point_idx = ext_lines.size() - 2;
        out.emplace_back(
            ext_lines[point_idx].p, ext_lines[point_idx].w, ext_lines[point_idx].perimeter_index
        );
        if (point_idx == 0)
            break;
        --point_idx;
    }

    if (ext_lines.back().p == ext_lines.front().p) {
        // Connect endpoints.
        out.front().p = out.back().p;
    }

    if (out.size() >= 3) {
        ext_lines.junctions = std::move(out);
    }
}

bool should_fuzzify(
    const PrintRegionConfig &config,
    const size_t layer_idx,
    const size_t perimeter_idx,
    const bool is_contour
) {
    const FuzzySkinType fuzzy_skin_type = config.fuzzy_skin.value;

    if (fuzzy_skin_type == FuzzySkinType::None || layer_idx <= 0) {
        return false;
    }

    const bool fuzzify_contours = perimeter_idx == 0;
    const bool fuzzify_holes = fuzzify_contours && fuzzy_skin_type == FuzzySkinType::All;

    return is_contour ? fuzzify_contours : fuzzify_holes;
}

Polygon apply_fuzzy_skin(
    const Polygon &polygon,
    const PrintRegionConfig &base_config,
    const PerimeterRegions &perimeter_regions,
    const size_t layer_idx,
    const size_t perimeter_idx,
    const bool is_contour
) {
    using namespace Slic3r::Algorithm::LineSegmentation;

    PhaseScope phase_scope(layer_idx);

    auto apply_fuzzy_skin_on_polygon =
        [&layer_idx, &perimeter_idx,
         &is_contour](const Polygon &polygon, const PrintRegionConfig &config) -> Polygon {
        if (should_fuzzify(config, layer_idx, perimeter_idx, is_contour)) {
            Polygon fuzzified_polygon = polygon;
            fuzzy_polygon(
                fuzzified_polygon, scaled<double>(config.fuzzy_skin_thickness.value),
                scaled<double>(config.fuzzy_skin_point_dist.value)
            );
            return fuzzified_polygon;
        } else {
            return polygon;
        }
    };

    if (perimeter_regions.empty()) {
        return apply_fuzzy_skin_on_polygon(polygon, base_config);
    }

    PolylineRegionSegments segments = polygon_segmentation(polygon, base_config, perimeter_regions);
    if (segments.size() == 1) {
        const PrintRegionConfig &config = segments.front().config;
        return apply_fuzzy_skin_on_polygon(polygon, config);
    }

    Polygon fuzzified_polygon;
    for (PolylineRegionSegment &segment : segments) {
        const PrintRegionConfig &config = segment.config;
        if (should_fuzzify(config, layer_idx, perimeter_idx, is_contour)) {
            fuzzy_polyline(
                segment.polyline.points, false, scaled<double>(config.fuzzy_skin_thickness.value),
                scaled<double>(config.fuzzy_skin_point_dist.value)
            );
        }

        assert(!segment.polyline.empty());
        if (segment.polyline.empty()) {
            continue;
        } else if (!fuzzified_polygon.empty() &&
                   fuzzified_polygon.back() == segment.polyline.front()) {
            // Remove the last point to avoid duplicate points.
            fuzzified_polygon.points.pop_back();
        }

        Slic3r::append(fuzzified_polygon.points, std::move(segment.polyline.points));
    }

    assert(!fuzzified_polygon.empty());
    if (fuzzified_polygon.front() == fuzzified_polygon.back()) {
        // Remove the last point to avoid duplicity between the first and the last point.
        fuzzified_polygon.points.pop_back();
    }

    return fuzzified_polygon;
}

Arachne::ExtrusionLine apply_fuzzy_skin(
    const Arachne::ExtrusionLine &extrusion,
    const PrintRegionConfig &base_config,
    const PerimeterRegions &perimeter_regions,
    const size_t layer_idx,
    const size_t perimeter_idx,
    const bool is_contour
) {
    using namespace Slic3r::Algorithm::LineSegmentation;
    using namespace Slic3r::Arachne;

    PhaseScope phase_scope(layer_idx);

    if (perimeter_regions.empty()) {
        if (should_fuzzify(base_config, layer_idx, perimeter_idx, is_contour)) {
            ExtrusionLine fuzzified_extrusion = extrusion;
            fuzzy_extrusion_line(
                fuzzified_extrusion, scaled<double>(base_config.fuzzy_skin_thickness.value),
                scaled<double>(base_config.fuzzy_skin_point_dist.value)
            );
            return fuzzified_extrusion;
        } else {
            return extrusion;
        }
    }

    ExtrusionRegionSegments segments =
        extrusion_segmentation(extrusion, base_config, perimeter_regions);
    ExtrusionLine fuzzified_extrusion(extrusion.inset_idx, extrusion.is_odd, extrusion.is_closed);

    for (ExtrusionRegionSegment &segment : segments) {
        const PrintRegionConfig &config = segment.config;
        if (should_fuzzify(config, layer_idx, perimeter_idx, is_contour)) {
            fuzzy_extrusion_line(
                segment.extrusion, scaled<double>(config.fuzzy_skin_thickness.value),
                scaled<double>(config.fuzzy_skin_point_dist.value)
            );
        }

        assert(!segment.extrusion.empty());
        if (segment.extrusion.empty()) {
            continue;
        } else if (!fuzzified_extrusion.empty() &&
                   fuzzified_extrusion.back().p == segment.extrusion.front().p) {
            // Remove the last point to avoid duplicate points (We don't care if the width of both
            // points is different.).
            fuzzified_extrusion.junctions.pop_back();
        }

        Slic3r::append(fuzzified_extrusion.junctions, std::move(segment.extrusion.junctions));
    }

    assert(!fuzzified_extrusion.empty());

    return fuzzified_extrusion;
}

} // namespace Slic3r::Feature::FuzzySkin
