#include "diagnostics/GoldenModesPTRegions.h"
#include <BAT/BCModel.h>
#include <mpi.h>
#include <iostream>
#include <limits>
#include <stdexcept>
class Regions : public BCModel {
public:
    Regions() : BCModel("regions") {}
    double LogLikelihood(const std::vector<double>&) override { return 0; }
};
int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    struct Finalize { ~Finalize() { MPI_Finalize(); } } finalize;
    try {
        if (argc > 1) {
            std::cout << "parameter,transform,threshold,true_label,false_label\n";
            for (const auto& c : GoldenModesPT::cuts())
                std::cout << c.parameter << ',' << c.transform << ',' << c.threshold << ',' << c.yes << ',' << c.no << '\n';
            return 0;
        }
        Regions model;
        model.AddParameter("unrelated", -10, 10);
        std::vector<std::pair<unsigned,unsigned>> indices;
        if (!GoldenModesPT::definitions(model,indices).empty()) throw std::runtime_error("absent channels included");
        for (const auto& c : GoldenModesPT::cuts()) model.AddParameter(c.parameter,-10,10);
        auto defs=GoldenModesPT::definitions(model,indices);
        if (defs.size()!=14 || indices.front().first!=1) throw std::runtime_error("parameter indexing mismatch");
        std::vector<double> values(model.GetNParameters(),0.);
        std::vector<unsigned> labels;
        GoldenModesPT::classify(indices,values,labels);
        for (unsigned i=0;i<labels.size();++i) if(labels[i]) throw std::runtime_error("zero classification mismatch");
        for (unsigned i=0;i<indices.size();++i) {
            const auto& c=GoldenModesPT::cuts()[i];
            values[indices[i].first] = c.transform=='i' ? c.threshold : c.transform=='c' ? 3. : -1.;
        }
        GoldenModesPT::classify(indices,values,labels);
        for (unsigned x : labels) if(x!=1) throw std::runtime_error("opposite region/boundary mismatch");
        auto saved=values;
        GoldenModesPT::classify(indices,values,labels);
        if(values!=saved) throw std::runtime_error("classifier mutated input");
        values[1]=std::numeric_limits<double>::quiet_NaN(); bool rejected=false;
        try { GoldenModesPT::classify(indices,values,labels); } catch(const std::runtime_error&) { rejected=true; }
        if (!rejected) throw std::runtime_error("nonfinite value accepted");
        std::cout<<"PASS: missing channels, indexed parameters, boundaries, pure classification\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
