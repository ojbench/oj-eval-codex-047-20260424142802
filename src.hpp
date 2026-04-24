#ifndef PPCA_SRC_HPP
#define PPCA_SRC_HPP

#include "math.h"

class Monitor; // forward declaration for pointer member

class Controller {

public:
    Controller(const Vec &_pos_tar, double _v_max, double _r, int _id, Monitor *_monitor) {
        pos_tar = _pos_tar;
        v_max = _v_max;
        r = _r;
        id = _id;
        monitor = _monitor;
    }

    void set_pos_cur(const Vec &_pos_cur) {
        pos_cur = _pos_cur;
    }

    void set_v_cur(const Vec &_v_cur) {
        v_cur = _v_cur;
    }

private:
    int id;
    Vec pos_tar;
    Vec pos_cur;
    Vec v_cur;
    double v_max, r;
    Monitor *monitor;

    // Check whether a candidate velocity will be safe w.r.t. other robots,
    // assuming others keep their last known velocity from monitor.
    bool is_safe(const Vec &v_candidate) const;

    // Clamp velocity magnitude to v_max and avoid NaNs/infs.
    Vec clamp_velocity(const Vec &v) const {
        Vec res = v;
        // sanitize
        if (!(res.x == res.x)) res.x = 0; // NaN guard
        if (!(res.y == res.y)) res.y = 0;
        double sp2 = res.norm_sqr();
        double vmax2 = v_max * v_max;
        if (sp2 > vmax2) {
            double sp = std::sqrt(std::max(sp2, 1e-12));
            if (sp > 0) res = res * (v_max / sp);
        }
        return res;
    }

public:

    Vec get_v_next() {
        // Desired direction towards target
        Vec to_tar = pos_tar - pos_cur;
        double dist = to_tar.norm();

        // If already at target (within epsilon), stop
        if (dist <= EPSILON) {
            return Vec();
        }

        // Ideal velocity to reach target in one step if possible
        Vec v_desired;
        double max_move = v_max * TIME_INTERVAL;
        if (dist <= max_move) {
            // Move exactly onto the target
            v_desired = to_tar / TIME_INTERVAL;
        } else {
            v_desired = to_tar.normalize() * v_max;
        }

        v_desired = clamp_velocity(v_desired);

        // If desired is safe, take it
        if (is_safe(v_desired)) {
            return v_desired;
        }

        // Try reducing speed along the same direction
        Vec dir = (dist > 0 ? to_tar.normalize() : Vec());
        const double scales[] = {0.8, 0.6, 0.4, 0.2, 0.0};
        for (double s : scales) {
            Vec cand = clamp_velocity(dir * (v_max * s));
            if (is_safe(cand)) {
                return cand;
            }
        }

        // Try sidestep by rotating direction
        const double angles[] = {PI / 6, -PI / 6, PI / 3, -PI / 3, PI / 2, -PI / 2};
        for (double ang : angles) {
            Vec side = dir.rotate(ang) * (0.5 * v_max);
            side = clamp_velocity(side);
            if (is_safe(side)) {
                return side;
            }
        }

        // Fallback: full stop
        return Vec();
    }
};


inline bool Controller::is_safe(const Vec &v_candidate) const {
    // Mirror simulator.h collision logic as closely as possible
    int n = monitor ? monitor->get_robot_number() : 0;
    for (int j = 0; j < n; ++j) {
        if (j == id) continue;
        Vec pj = monitor->get_pos_cur(j);
        Vec vj = monitor->get_v_cur(j);
        double rj = monitor->get_r(j);

        Vec delta_pos = pos_cur - pj;
        Vec delta_v = v_candidate - vj;
        double project = delta_pos.dot(delta_v);
        if (project >= 0) {
            continue; // moving apart or stationary in separating direction
        }
        double dv_norm = delta_v.norm();
        // project is negative here; divide by -|dv|
        project /= -dv_norm;
        double min_dis_sqr;
        double delta_r = r + rj;
        if (project < dv_norm * TIME_INTERVAL) {
            min_dis_sqr = delta_pos.norm_sqr() - project * project;
        } else {
            min_dis_sqr = (delta_pos + delta_v * TIME_INTERVAL).norm_sqr();
        }
        if (min_dis_sqr <= delta_r * delta_r - EPSILON) {
            return false;
        }
    }
    // speed guard and NaN guard
    if (v_candidate.norm_sqr() >= v_max * v_max + EPSILON) return false;
    if (!(v_candidate.x == v_candidate.x) || !(v_candidate.y == v_candidate.y)) return false;
    return true;
}


#endif // PPCA_SRC_HPP
