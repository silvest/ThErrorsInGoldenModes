#ifndef HISTO_H
#define HISTO_H

#include <string>
#include <map>
#include "TH1D.h"
#include "TH2D.h"

struct Histos {
    std::map<std::string, TH1D*> h1d;
    std::map<std::string, TH2D*> h2d;
    std::map<std::string, double>& obs;  // Reference to observables map

    Histos(std::map<std::string, double>& obs_ref) : obs(obs_ref) {}

    ~Histos() {
        for (auto& p : h1d) {
            if (p.second) delete p.second;
        }
        for (auto& p : h2d) {
            if (p.second) delete p.second;
        }
    }

    void createH1D(const std::string& name, int nbins, double xmin, double xmax) {
        if (h1d.count(name) > 0) {
            // Histogram already exists
            return;
        }
        h1d[name] = new TH1D(name.c_str(), name.c_str(), nbins, xmin, xmax);
    }

    void createH2D(const std::string& name_x, const std::string& name_y, 
                   int nbinsx, double xmin, double xmax,
                   int nbinsy, double ymin, double ymax) {
        std::string name = name_x + "_vs_" + name_y;
        if (h2d.count(name) > 0) {
            // Histogram already exists
            return;
        }
        h2d[name] = new TH2D(name.c_str(), name.c_str(), 
                             nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    }

    void fill() {
        for (auto& p : h1d) {
            if (obs.find(p.first) != obs.end()) {
                p.second->Fill(obs[p.first]);
            }
        }
    }

    void fillh1d() {
        for (auto& p : h1d) {
            if (obs.find(p.first) != obs.end()) {
                p.second->Fill(obs[p.first]);
            }
            else {
                std::cout << "Observable " << p.first << " not found in obs for filling histogram!" << std::endl;
            }
        }
    }

    void fillh2d() {
        // Fill 2D histograms from pairs of observables
        for (auto& p : h2d) {
            // Extract variable names from histogram name (format: "varX_vs_varY")
            std::string name = p.first;
            size_t pos = name.find("_vs_");
            if (pos != std::string::npos) {
                std::string var_x = name.substr(0, pos);
                std::string var_y = name.substr(pos + 4);
                
                if (obs.find(var_x) != obs.end() && obs.find(var_y) != obs.end()) {
                    p.second->Fill(obs[var_x], obs[var_y]);
                }
            }
        }
    }
};

#endif // HISTO_H
