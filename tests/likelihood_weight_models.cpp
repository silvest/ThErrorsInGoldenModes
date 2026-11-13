#include "goldenmodesB.h"
#include "goldenmodesB_indSU3.h"
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

template<class Model> class Inspect : public Model {
public:
    using Model::Model;
    void verify(bool constraint) {
        std::vector<double> p;
        for(unsigned i=0;i<this->GetNParameters();++i) {
            const auto& par=this->GetParameter(i);
            p.push_back(par.Fixed()?par.GetFixedValue():par.GetLowerLimit()+(.41+.031*(i%5))*par.GetRangeWidth());
        }
        const double original=this->LogLikelihood(p);
        this->SetWriteMCMCWeights();
        const double captured=this->EvaluateOutputLikelihood(p);
        if(!std::isfinite(original) || captured!=original)throw std::runtime_error("weight capture changed model likelihood");
        double sum=0.;unsigned groups=0,singles=0,members=0,constraints=0;
        for(const auto& item:this->likelihoodWeights.Entries()) {
            const auto& e=item.second;
            if(!e.active || !std::isfinite(e.value))throw std::runtime_error("invalid captured term");
            if(e.additive)sum+=e.value;
            groups+=e.kind=="correlated_block";singles+=e.kind=="uncorrelated";
            members+=e.kind=="correlated_member";constraints+=e.kind=="constraint";
        }
        if(std::abs(sum-original)>1.e-10*std::max(1.,std::abs(original)) || !groups || !singles || members<groups || constraints!=unsigned(constraint))
            throw std::runtime_error("missing or incorrectly classified model likelihood terms");
        std::cout<<"PASS: "<<groups<<" blocks, "<<singles<<" singles, "<<members<<" member contrasts, "<<constraints<<" constraints\n";
    }
};
int main(int argc,char**argv) {
    MPI_Init(&argc,&argv);int status=0,rank=0,size=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    // Pure likelihoods also run on worker ranks, which skip dense plotting histograms.
    if(rank==(size>1?1:0)) {
    try {
        double ewp=.1,dsu3=.2;
        {Inspect<goldenmodesB> m(dsu3,ewp,true,size>1,size>1);m.verify(false);}
        for(const std::string mode: {"abs","complex","reim"}) {
            Inspect<goldenmodesB_indSU3> m(ewp,true,true,true,.3,true,false,false,false,false,false);
            m.SetSU3Weight(mode);m.verify(true);
        }
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';status=1;}
    }
    MPI_Bcast(&status,1,MPI_INT,size>1?1:0,MPI_COMM_WORLD);
    MPI_Finalize();return status;
}
