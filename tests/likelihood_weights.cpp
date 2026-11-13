#include "LikelihoodWeights.h"
#include <BAT/BCLog.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

static void check(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
static CorrelatedGaussianObservables measurement(double rho = .5) {
    TMatrixDSym correlation(2); correlation.UnitMatrix();correlation(0,1)=correlation(1,0)=rho;
    return CorrelatedGaussianObservables({dato(0.,1.),dato(0.,1.)},correlation);
}
class Toy : public GoldenModesLikelihoodModel {
public:
    unsigned evaluations=0;
    CorrelatedGaussianObservables correlated=measurement();
    Toy() {
        SetName("weights_test");AddParameter("x",-8.,8.);AddParameter("y",-8.,8.);
        SetPriorGauss("x",0.,2.);SetPriorGauss("y",0.,3.);
    }
    double LogLikelihood(const std::vector<double>& p) override {
        ++evaluations;TVectorD values(2);values(0)=p[0];values(1)=p[1];
        double ll=likelihoodWeights.Gaussian("single_observable",p[0]-1.,1.5);
        ll+=CorrelatedWeight("correlated_block",correlated,values,{"x_measurement","y_measurement"});
        const double constraint=-.05*p[1]*p[1];ll+=constraint;
        likelihoodWeights.Record("constraint_test",constraint,"constraint");
        return ll;
    }
    void MCMCUserIterationInterface() override {
        for(unsigned c=0;c<GetNColdChains();++c) {
            const double ll=EvaluateOutputLikelihood(fMCMCStates[c].parameters);
            WriteLikelihoodWeights(c,ll);
        }
    }
};
static void configure(Toy& model,bool pt,bool weights,const std::string& path) {
    model.SetNChains(2);if(pt)model.SetParallelTempering({1.,2.,4.},2,2);
    model.SetWriteMCMCWeights(weights);model.SetRandomSeed(82034);
    model.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,13);
    model.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth);
    model.SetProposeMultivariate(false);model.SetNIterationsPreRunMax(120);
    model.SetNIterationsPreRunCheck(20);model.SetNIterationsRun(101);model.SetNLag(3);
    // Gaussian-prior sampling in this BAT version uses ROOT's global RNG.
    // Explicit starts isolate the trajectory comparison from that external state.
    std::vector<std::vector<double>> starts(model.GetNChains());
    for(unsigned c=0;c<model.GetNChains();++c)starts[c]={-.6+.1*c,.3-.1*c};
    model.SetInitialPositions(starts);
    model.WriteMarkovChain(path,"RECREATE",true,false);
}
static void math_and_schema_tests() {
    TVectorD r(2);r(0)=1.;r(1)=2.;auto correlated=measurement();
    const auto w=correlated.computeLeaveOneOutWeights(r);
    check(std::abs(w[0]+2.)<1.e-12,"incorrect correlated-block weight");
    check(std::abs(w[1]-2./3.)<1.e-12 && std::abs(w[2]+4./3.)<1.e-12,"HEPfit member contrast differs");
    check(w[1]>0. && std::abs(w[1]+w[2]-w[0])>1.,"correlated contrasts incorrectly assumed additive");
    for(int i=0;i<2;++i) {
        TVectorD zeroed(r);zeroed(i)=0.;
        check(std::abs(w[i+1]-(w[0]-correlated.logweight(zeroed)))<1.e-12,"residual-removal definition differs");
    }
    const auto independent=measurement(0.).computeLeaveOneOutWeights(r);
    check(std::abs(independent[1]+independent[2]-independent[0])<1.e-12,"independent limit incorrect");
    check(LikelihoodWeightOutput::SafeName("a-b")!=LikelihoodWeightOutput::SafeName("a_x2d_b"),"branch encoding collides");
    {
        TFile file("weights_schema.root","RECREATE");LikelihoodWeightOutput output;output.Enable(true);output.Initialize(&file,"schema");
        {LikelihoodWeightOutput::Capture capture(output);output.Record("a-b",-1.);}
        output.Fill(0,1,1,0,-1.,-1.,0.,-1.);
        {LikelihoodWeightOutput::Capture capture(output);output.Record("a-b",-2.);output.Record("late",-3.);}
        output.Fill(0,2,2,0,-5.,-5.,0.,-5.);
        {LikelihoodWeightOutput::Capture capture(output);output.Record("a-b",-7.);}
        output.Fill(0,3,3,0,-7.,-7.,0.,-7.);
        bool rejected=false;try{output.Fill(0,4,4,0,-8.,-8.,0.,-8.);}catch(const std::runtime_error&){rejected=true;}
        check(rejected,"mismatched likelihood accepted");file.Write();
    }
    TFile file("weights_schema.root","READ");auto tree=dynamic_cast<TTree*>(file.Get("schema_Observables"));
    double value;bool active;tree->SetBranchAddress("late_weight",&value);tree->SetBranchAddress("late_weight_active",&active);
    check(tree->GetEntries()==3,"invalid row was saved");
    for(unsigned i=0;i<3;++i){tree->GetEntry(i);check(value==(i==1?-3.:0.) && active==(i==1),"conditional term was stale or not backfilled");}
}
static void inspect(const std::string& prefix) {
    TFile file((prefix+"on.root").c_str(),"READ"), baseline((prefix+"off.root").c_str(),"READ");
    auto mcmc=dynamic_cast<TTree*>(file.Get("weights_test_mcmc"));
    auto off=dynamic_cast<TTree*>(baseline.Get("weights_test_mcmc"));
    auto output=dynamic_cast<TTree*>(file.Get("weights_test_Observables"));
    auto schema=dynamic_cast<TTree*>(file.Get("weights_test_LikelihoodWeightInfo"));
    check(mcmc && output && schema && off,"output tree missing");
    check(!baseline.Get("weights_test_Observables"),"disabled output created weight tree");
    check(output->GetEntries()==66 && mcmc->GetEntries()==66 && off->GetEntries()==66,"lag/warmup/hot sample alignment incorrect");
    double x,y,ll,ol,lp,op,ws,single,block,memberX,memberY,xt,yt,ox,oy;
    unsigned c,oc,it,oit,r,orr;int phase;
    mcmc->SetBranchAddress("x",&x);mcmc->SetBranchAddress("y",&y);mcmc->SetBranchAddress("LogLikelihood",&ll);
    mcmc->SetBranchAddress("Chain",&c);mcmc->SetBranchAddress("Iteration",&it);mcmc->SetBranchAddress("Replica",&r);
    off->SetBranchAddress("x",&xt);off->SetBranchAddress("y",&yt);
    output->SetBranchAddress("LogLikelihood",&ol);output->SetBranchAddress("LogPrior",&lp);output->SetBranchAddress("LogProbability",&op);
    output->SetBranchAddress("x",&ox);output->SetBranchAddress("y",&oy);
    output->SetBranchAddress("WeightSum",&ws);output->SetBranchAddress("Phase",&phase);
    output->SetBranchAddress("Chain",&oc);output->SetBranchAddress("Iteration",&oit);output->SetBranchAddress("Replica",&orr);
    output->SetBranchAddress("single_observable_weight",&single);output->SetBranchAddress("correlated_block_weight",&block);
    output->SetBranchAddress("correlated_block__x_measurement_weight",&memberX);output->SetBranchAddress("correlated_block__y_measurement_weight",&memberY);
    bool exchanged=false;
    for(Long64_t i=0;i<66;++i) {
        mcmc->GetEntry(i);output->GetEntry(i);off->GetEntry(i);
        check(x==xt && y==yt,"weight output changed the sampler trajectory");
        check(x==ox && y==oy,"observable-tree parameter columns are misaligned");
        check(c==oc && it==oit && r==orr && c<2 && phase==1,"chain/iteration/replica metadata misaligned");
        check(ll==ol && std::abs(ws-ll)<1.e-10 && std::abs(op-lp-ll)<1.e-10,"total log densities differ");
        check(std::abs(single+.5*std::pow((x-1.)/1.5,2))<1.e-10,"single-observable weight is not at the current state");
        const double expected=-.5*(x*x-x*y+y*y)/.75;
        check(std::abs(block-expected)<1.e-10,"correlated block incorrect at current state");
        check(std::abs(memberX+.5*x*(x-y)/.75)<1.e-10 && std::abs(memberY+.5*y*(y-x)/.75)<1.e-10,"correlated member contrast incorrect");
        exchanged|=r>=2;
    }
    if(prefix=="weights_pt_")check(exchanged,"PT test did not exercise hot-to-cold transport");
    check(schema->GetEntries()==5,"weight metadata has wrong number of terms");
}
static void inspect_file(const std::string& path, unsigned expected) {
    TFile file(path.c_str(),"READ");check(!file.IsZombie(),"cannot read ROOT output");
    std::string model;TIter keys(file.GetListOfKeys());
    while(auto key=dynamic_cast<TKey*>(keys())) {
        const std::string name=key->GetName();
        if(name.size()>5 && name.substr(name.size()-5)=="_mcmc")model=name.substr(0,name.size()-5);
    }
    auto mcmc=dynamic_cast<TTree*>(file.Get((model+"_mcmc").c_str()));
    auto values=dynamic_cast<TTree*>(file.Get((model+"_Observables").c_str()));
    auto info=dynamic_cast<TTree*>(file.Get((model+"_LikelihoodWeightInfo").c_str()));
    check(mcmc && values && info,"missing model weight trees");
    check(values->GetEntries()==expected && mcmc->GetEntries()==expected,"weight rows differ from saved chains");
    std::string* branch=nullptr;bool additive=false;
    info->SetBranchAddress("branch",&branch);info->SetBranchAddress("additive",&additive);
    std::map<std::string,double> weights;
    for(Long64_t i=0;i<info->GetEntries();++i) {
        info->GetEntry(i);
        if(additive)weights.emplace(*branch,0.);
    }
    for(auto& item:weights)values->SetBranchAddress(item.first.c_str(),&item.second);
    unsigned chain=0,otherChain=0,iteration=0,otherIteration=0,replica=0,otherReplica=0;int phase=0;
    double ll=0.,otherLL=0.,sum=0.;
    values->SetBranchAddress("Chain",&chain);values->SetBranchAddress("Iteration",&iteration);values->SetBranchAddress("Replica",&replica);
    values->SetBranchAddress("Phase",&phase);values->SetBranchAddress("LogLikelihood",&ll);values->SetBranchAddress("WeightSum",&sum);
    mcmc->SetBranchAddress("Chain",&otherChain);mcmc->SetBranchAddress("Iteration",&otherIteration);mcmc->SetBranchAddress("Replica",&otherReplica);
    mcmc->SetBranchAddress("LogLikelihood",&otherLL);
    for(unsigned i=0;i<expected;++i) {
        values->GetEntry(i);mcmc->GetEntry(i);double total=0.;
        for(const auto& item:weights) {check(std::isfinite(item.second),"nonfinite saved weight");total+=item.second;}
        check(phase==1 && chain<2 && chain==otherChain && iteration==otherIteration && replica==otherReplica,"misaligned saved metadata");
        check(std::isfinite(ll) && ll==otherLL && std::abs(total-ll)<1.e-10*std::max(1.,std::abs(ll)) && std::abs(sum-ll)<1.e-10*std::max(1.,std::abs(ll)),"weights do not sum to model likelihood");
    }
    std::cout<<"PASS: "<<expected<<" aligned saved states, "<<weights.size()<<" additive terms\n";
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        if(argc==4 && std::string(argv[1])=="check") {
            inspect_file(argv[2],std::stoul(argv[3]));MPI_Finalize();return 0;
        }
        BCLog::SetLogLevelScreen(BCLog::nothing);Toy worker;
        if(rank) {
            double recv[3],value;
            for(;;) {
                MPI_Scatter(nullptr,3,MPI_DOUBLE,recv,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(recv[0]==-1.)break;
                value=recv[0]==1.?worker.LogEval(std::vector<double>(recv+1,recv+3)):-INFINITY;
                MPI_Gather(&value,1,MPI_DOUBLE,nullptr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        } else {
            const bool pt=argc>1 && std::string(argv[1])=="pt";
            const std::string prefix=pt?"weights_pt_":"weights_ordinary_";
            math_and_schema_tests();
            Toy off;configure(off,pt,false,prefix+"off.root");off.MarginalizeAll(BCIntegrate::kMargMetropolis);off.CloseOutputFile();
            Toy on;configure(on,pt,true,prefix+"on.root");on.MarginalizeAll(BCIntegrate::kMargMetropolis);on.CloseOutputFile();
            check(on.evaluations==off.evaluations,"weight output added likelihood evaluations");
            inspect(prefix);
            std::vector<double> stop(size*3,0.);double recv[3];for(int i=1;i<size;++i)stop[i*3]=-1.;
            if(size>1)MPI_Scatter(stop.data(),3,MPI_DOUBLE,recv,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
            std::cout<<"PASS: analytic correlated contrasts, additive sums, ROOT metadata, lag, rejected states, swaps and unchanged sampling\n";
        }
    }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
