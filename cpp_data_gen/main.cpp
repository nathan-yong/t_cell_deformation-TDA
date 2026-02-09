#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
#include <cfloat>   // at top (or <limits> already included)


struct Params {
    // constantes_tarea.m
    double L = 86.0;
    double f_rep = 20.0;
    double f_adh = 0.75;
    double d_cut = (7.0/6.0)*2.0 - 2.0; // 2*(7/6-1)
    int n_drt = 1;                      // number of dendritic cells (APC)
    double amp_rad = 0.15;
    double w = 6.0;

    // Main_tarea.m
    double dt = 0.005;
    double tfinal = 3e5 * 0.005; // 1.5
    int writeEvery = 20000;        // every 200 steps
};

static inline double hypot2(double x, double y) {
    return std::sqrt(x*x + y*y);
}

static inline void apply_pbc_and_update_ids(std::vector<double>& x,
                                           std::vector<double>& y,
                                           std::vector<int64_t>& id_par,
                                           double L) {
    // Equivalent to PBC_fun_id.m:
    // wrap positions to [0,L] and if a particle crosses any boundary,
    // assign it a fresh new id = max(id_par)+1.
    int64_t current_max = 0;
    for (auto v : id_par) current_max = std::max(current_max, v);

    auto bump = [&current_max]() {
        current_max += 1;
        return current_max;
    };

    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        bool changed = false;
        if (x[i] < 0.0) { x[i] += L; changed = true; }
        else if (x[i] > L) { x[i] -= L; changed = true; }

        if (y[i] < 0.0) { y[i] += L; changed = true; }
        else if (y[i] > L) { y[i] -= L; changed = true; }

        if (changed) id_par[i] = bump();
    }
}

static inline void normalize(double& vx, double& vy) {
    double n = std::sqrt(vx*vx + vy*vy);
    if (n <= 0.0) { vx = 1.0; vy = 0.0; return; }
    vx /= n; vy /= n;
}

static void build_triangular_grid(double L, double d,
                                  std::vector<double>& x,
                                  std::vector<double>& y) {
    // TriangularGrid2.m (no plotting)
    const int nx = static_cast<int>(std::floor(L / d)) + 1;
    const int ny = static_cast<int>(std::floor(L / d)) + 1;
    x.clear(); y.clear();
    x.reserve(static_cast<size_t>(nx*ny));
    y.reserve(static_cast<size_t>(nx*ny));
    for (int iy = 0; iy < ny; ++iy) {
        double yy = iy * d;
        for (int ix = 0; ix < nx; ++ix) {
            double xx = ix * d;
            if ((iy % 2) == 1) xx += d/2.0;
            x.push_back(xx);
            y.push_back(yy);
        }
    }
}

static void init_particles(int np, int64_t seed, const Params& p,
                           std::vector<double>& x,
                           std::vector<double>& y,
                           std::vector<double>& vx0,
                           std::vector<double>& vy0,
                           std::vector<double>& r_cer,
                           std::vector<double>& phi_rad,     // phase theta_i (kept)
                           std::vector<double>& r_part,      // initial radii (all r0)
                           std::vector<double>& t_start,     // NEW: start time for oscillation (per particle)
                           std::vector<int64_t>& id_par) {

    std::mt19937_64 rng(static_cast<uint64_t>(seed));
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    const double r0   = 1.0;
    const double vdes = 1.0;

    // -------------------------
    // 1) Build triangular grid and randomly assign np particles to np grid sites
    // -------------------------
    std::vector<double> gx, gy;
    const double dgrid = 1.8; // same spacing as before; adjust if you want
    build_triangular_grid(p.L, dgrid, gx, gy);

    if (static_cast<int>(gx.size()) < np) {
        throw std::runtime_error("Triangular grid has fewer points than np.");
    }

    std::vector<size_t> idx(gx.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::shuffle(idx.begin(), idx.end(), rng);

    x.resize(np);
    y.resize(np);
    for (int i = 0; i < np; ++i) {
        x[i] = gx[idx[i]];
        y[i] = gy[idx[i]];
    }

    // -------------------------
    // 2) Choose the APC as the particle closest to the center (among the chosen grid sites)
    //    Swap it into index 0 (no need to force exactly at L/2, L/2).
    // -------------------------
    int apc = 0;
    double best = std::numeric_limits<double>::max();
    const double cx = p.L * 0.5;
    const double cy = p.L * 0.5;

    for (int i = 0; i < np; ++i) {
        const double dx = x[i] - cx;
        const double dy = y[i] - cy;
        const double d2 = dx*dx + dy*dy;
        if (d2 < best) { best = d2; apc = i; }
    }

    if (apc != 0) {
        std::swap(x[0], x[apc]);
        std::swap(y[0], y[apc]);
    }

    // -------------------------
    // 3) IDs
    // -------------------------
    id_par.resize(np);
    for (int i = 0; i < np; ++i) id_par[i] = static_cast<int64_t>(i + 1);

    // -------------------------
    // 4) Initial radii: all r0, no oscillation yet
    // -------------------------
    r_cer.assign(np, r0);
    r_part.assign(np, r0);

    // -------------------------
    // 5) Random phases, and per-particle oscillation start times:
    //    they begin oscillating at the first time t_start_i >= 0 such that
    //    sin(w*t_start_i + theta_i) = 0  (=> radius equals r0 at the start).
    // -------------------------
    phi_rad.resize(np);
    t_start.resize(np);

    for (int i = 0; i < np; ++i) {
        const double theta = U01(rng) * 2.0 * M_PI;  // individual phase
        phi_rad[i] = theta;

        // Find smallest t >= 0 such that w*t + theta = m*pi
        // m = ceil(theta/pi)
        const double m = std::ceil(theta / M_PI);
        t_start[i] = (m * M_PI - theta) / p.w;       // in [0, pi/w]
    }

    // APC: keep it non-oscillatory from the start (optional but consistent with "APC")
    t_start[0] = std::numeric_limits<double>::max();

    // -------------------------
    // 6) Desired velocities: random directions, magnitude vdes
    // -------------------------
    vx0.resize(np);
    vy0.resize(np);
    for (int i = 0; i < np; ++i) {
        double ax = 2.0*U01(rng) - 1.0;
        double ay = 2.0*U01(rng) - 1.0;
        normalize(ax, ay);
        vx0[i] = ax * vdes;
        vy0[i] = ay * vdes;
    }
}


static inline double wrap_min_image(double d, double L) {
    d -= L * std::round(d / L);
    return d;
}


struct CellList {
    int nx = 1, ny = 1;
    double cellSize = 1.0;
    std::vector<int> head; // size nx*ny
    std::vector<int> next; // size n

    inline int cidx(int ix, int iy) const { return ix + nx * iy; }

    inline void build(const std::vector<double>& x,
                      const std::vector<double>& y,
                      double L,
                      double cellSize_in) {
        cellSize = cellSize_in;
        nx = std::max(1, (int)std::floor(L / cellSize));
        ny = std::max(1, (int)std::floor(L / cellSize));

        head.assign(nx * ny, -1);
        next.assign((int)x.size(), -1);

        for (int i = 0; i < (int)x.size(); ++i) {
            int ix = (int)std::floor(x[i] / cellSize);
            int iy = (int)std::floor(y[i] / cellSize);

            // clamp (in case x==L or y==L)
            if (ix < 0) ix = 0;
            if (iy < 0) iy = 0;
            if (ix >= nx) ix = nx - 1;
            if (iy >= ny) iy = ny - 1;

            const int c = cidx(ix, iy);
            next[i] = head[c];
            head[c] = i;
        }
    }
};


static void compute_forces(const std::vector<double>& x,
                           const std::vector<double>& y,
                           const std::vector<double>& r,
                           const Params& p,
                           std::vector<double>& fx,
                           std::vector<double>& fy,
                           CellList& cl) {
    const int n = (int)x.size();
    std::fill(fx.begin(), fx.end(), 0.0);
    std::fill(fy.begin(), fy.end(), 0.0);

    // Conservative cutoff for neighbor search:
    // repulsion if dist < (r_i+r_j)
    // adhesion if dist < (r_i+r_j + d_cut) AND not wrapped
    const double rmax = 1.0 * (1.0 + p.amp_rad);   // safe upper bound in your model
    const double rcut = 2.0 * rmax + p.d_cut;
    const double rcut2 = rcut * rcut;

    // cell size = rcut => only need 3x3 neighbor cells
    cl.build(x, y, p.L, rcut);

    for (int iy = 0; iy < cl.ny; ++iy) {
        for (int ix = 0; ix < cl.nx; ++ix) {
            const int c = cl.cidx(ix, iy);

            for (int i = cl.head[c]; i != -1; i = cl.next[i]) {

                for (int ddy = -1; ddy <= 1; ++ddy) {
                    int jy = iy + ddy;
                    if (jy < 0) jy += cl.ny;
                    else if (jy >= cl.ny) jy -= cl.ny;

                    for (int ddx = -1; ddx <= 1; ++ddx) {
                        int jx = ix + ddx;
                        if (jx < 0) jx += cl.nx;
                        else if (jx >= cl.nx) jx -= cl.nx;

                        const int cn = cl.cidx(jx, jy);

                        for (int j = cl.head[cn]; j != -1; j = cl.next[j]) {
                            if (j <= i) continue; // each pair once

                            const double dx0 = x[i] - x[j];
                            const double dy0 = y[i] - y[j];

                            const double dx = wrap_min_image(dx0, p.L);
                            const double dy = wrap_min_image(dy0, p.L);

                            const bool wrapped = (dx != dx0) || (dy != dy0);

                            const double dist2 = dx*dx + dy*dy;
                            if (dist2 > rcut2) continue;

                            const double dist = std::sqrt(dist2);
                            if (dist <= 1e-12) continue;

                            const double gap = dist - (r[i] + r[j]);
                            const double nxij = dx / dist;
                            const double nyij = dy / dist;

                            double F = 0.0;
                            if (gap < 0.0) {
                                F = (-gap) * p.f_rep;
                            } else if (!wrapped && gap < p.d_cut) {
                                F = (-gap) * p.f_adh;
                            } else {
                                continue;
                            }

                            const double fxij = nxij * F;
                            const double fyij = nyij * F;

                            fx[i] += fxij;
                            fy[i] += fyij;
                            fx[j] -= fxij;
                            fy[j] -= fyij;
                        }
                    }
                }
            }
        }
    }
}


static void run_simulation(int np, int64_t seed, const Params& p) {
    std::vector<double> x, y, vx0, vy0, r_cer, phi_rad, r_part;
    std::vector<double> t_start;
    std::vector<int64_t> id_par;
    init_particles(np, seed, p, x, y, vx0, vy0, r_cer, phi_rad, r_part, t_start, id_par);

    std::vector<double> fx(np, 0.0), fy(np, 0.0);
    std::unordered_set<int64_t> visited; // unique particle IDs that ever contacted APC
    std::vector<double> t3;              // times when a NEW particle contacts APC for the first time
    t3.reserve(static_cast<size_t>(np));
    CellList cl;

    const std::string filename =
        "Results.txt";
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open output file: " + filename);

    const int nSteps = static_cast<int>(std::llround(p.tfinal / p.dt));
    int step = 0;

    for (int s = 0; s <= nSteps; ++s) {
        const double t = s * p.dt;

        // --- Contact detection with APC (particle 0) ---
        // If a particle j touches APC and j is NEW (not in visited), store time t in t3.
        for (int j = 1; j < np; ++j) {
            double dx = x[0] - x[j];
            double dy = y[0] - y[j];
            dx = wrap_min_image(dx, p.L);
            dy = wrap_min_image(dy, p.L);

            double dist = std::sqrt(dx*dx + dy*dy);
            double gap  = dist - (r_part[0] + r_part[j]);

            if (gap < 0.0) {
                auto ins = visited.insert(id_par[j]);
                if (ins.second) {
                    t3.push_back(t); // <-- time of NEW (first-time) visit/contact
                }
            }
        }

        // Forces

        compute_forces(x, y, r_part, p, fx, fy, cl);

        // Euler overdamped
        std::vector<double> vx(np), vy(np);
        for (int i = 0; i < np; ++i) {
            vx[i] = vx0[i] + fx[i];
            vy[i] = vy0[i] + fy[i];
            x[i] += vx[i] * p.dt;
            y[i] += vy[i] * p.dt;
        }

        // Update radii oscillation
        for (int i = 0; i < np; ++i) {
            r_part[i] = r_cer[i] * (1.0 + p.amp_rad * std::sin(p.w * t + phi_rad[i]));
        }

        // Fix dendritic cells (APC)
        for (int i = 0; i < p.n_drt; ++i) {
            x[i] = p.L / 2.0;
            y[i] = p.L / 2.0;
            r_part[i] = r_cer[i];
        }

        // PBC + ID changes
        apply_pbc_and_update_ids(x, y, id_par, p.L);

        // Write state every writeEvery steps (unchanged)
       if ((step % p.writeEvery) == 0) {
    	out.setf(std::ios::fixed);
    	out.precision(3);

    	// write time
    	out << "t " << t << "\n";

    	// write one particle per line
    	for (int i = 0; i < np; ++i) {
        	out << id_par[i] << " "
            << x[i] << " "
            << y[i] << " "
            << vx[i] << " "
            << vy[i] << " "
            << r_part[i] << "\n";
    	}

    	// blank line to separate snapshots
    	//out << "\n";
    }

        step++;
    }

    // Summary (unchanged)
   // std::ofstream summary("Summary_Seed_" + std::to_string(seed) + "_np_" + std::to_string(np) + ".txt");
   // summary << "np " << np << "\n";
   // summary << "seed " << seed << "\n";
   // summary << "L " << p.L << "\n";
   // summary << "f_rep " << p.f_rep << "\n";
   // summary << "f_adh " << p.f_adh << "\n";
   // summary << "d_cut " << p.d_cut << "\n";
   // summary << "dt " << p.dt << "\n";
   // summary << "tfinal " << p.tfinal << "\n";
   // summary << "unique_contacts_APC " << visited.size() << "\n";

    // --- NEW OUTPUT FILE: times of first-time contacts only ---
    const std::string t3file = "First_contact_times.txt";
    std::ofstream outT3(t3file);
    if (!outT3) throw std::runtime_error("Failed to open t3 output file: " + t3file);

    outT3.setf(std::ios::fixed);
    outT3.precision(3);
    // one time per line (only the array of times, nothing else)
    for (double tt : t3) outT3 << tt << "\n";
}

static void print_usage() {
    std::cerr << "Usage: apc_sim --np <int> --seed <int> [--dt <double>] [--tfinal <double>]\n";
}

int main(int argc, char** argv) {
    Params p;
    int np = -1;
    int64_t seed = -1;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string& key) {
            if (i+1 >= argc) {
                std::cerr << "Missing value for " << key << "\n";
                std::exit(2);
            }
        };
        if (a == "--np") { need(a); np = std::stoi(argv[++i]); }
        else if (a == "--seed") { need(a); seed = std::stoll(argv[++i]); }
        else if (a == "--dt") { need(a); p.dt = std::stod(argv[++i]); }
        else if (a == "--tfinal") { need(a); p.tfinal = std::stod(argv[++i]); }
        else if (a == "--L") { need(a); p.L = std::stod(argv[++i]); }
        else if (a == "--f_rep") { need(a); p.f_rep = std::stod(argv[++i]); }
        else if (a == "--f_adh") { need(a); p.f_adh = std::stod(argv[++i]); }
        else if (a == "--amp_rad") { need(a); p.amp_rad = std::stod(argv[++i]); }
        else if (a == "--w") { need(a); p.w = std::stod(argv[++i]); }
        else if (a == "--writeEvery") { need(a); p.writeEvery = std::stoi(argv[++i]); }
        else { print_usage(); return 2; }
    }

    if (np <= 0 || seed < 0) { print_usage(); return 2; }

    try {
        run_simulation(np, seed, p);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
