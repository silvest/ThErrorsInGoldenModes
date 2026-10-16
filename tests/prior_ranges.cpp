// Check that PT full priors join every partition and preserve other parameters.
#include "goldenmodesB_indSU3.h"
#include <BAT/BCPrior.h>
#include <BAT/BCGaussianPrior.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

using Bounds = std::pair<double,double>;
using Ranges = std::map<std::string,Bounds>;
static void check(bool ok, const std::string& reason) {
    if (!ok) throw std::runtime_error(reason);
}
static Ranges snapshotPriorRanges(bool full, bool positive, bool small) {
    double ewp = 0.1;
    goldenmodesB_indSU3 model(ewp, true, true, true, 0.3, true,
        positive, positive, positive, small, false, positive, positive,
        small, positive, positive, full);
    Ranges result;
    for (unsigned i=0; i<model.GetNParameters(); ++i) {
        const auto& p=model.GetParameter(i);
        result[p.GetName()]={p.GetLowerLimit(),p.GetUpperLimit()};
    }
    check(dynamic_cast<const BCGaussianPrior*>(model.GetParameter("CKM_Vud").GetPrior().get()), "CKM Gaussian prior changed");
    check(dynamic_cast<const BCGaussianPrior*>(model.GetParameter("theta_P").GetPrior().get()), "theta_P Gaussian prior changed");
    check(dynamic_cast<const BCGaussianPrior*>(model.GetParameter("dP4EW_ucd_BPJPSI_abs").GetPrior().get()), "EWP Gaussian prior changed");
    return result;
}
int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int status=0;
    try {
        const Ranges off=snapshotPriorRanges(false,false,false);
        const Ranges small=snapshotPriorRanges(false,true,true);
        const Ranges large=snapshotPriorRanges(false,true,false);
        const Ranges fullOff=snapshotPriorRanges(true,false,false);
        const Ranges fullOn=snapshotPriorRanges(true,true,true);
        const Ranges fullLarge=snapshotPriorRanges(true,true,false);
        check(fullOff==fullOn && fullOff==fullLarge,"region flags affect full priors");
        const Ranges fullExpected={
            {"G2t_scd_BJPSIP_arg",{-M_PI,M_PI}},
            {"G2t_dcd_BJPSIP_arg",{-M_PI,M_PI}},
            {"E2t_ccdd_BJPSIP_abs",{0.,4.}},
            {"EA1_sdcd_BPJPSI_abs",{0.,4.}},
            {"EA1_sdcd_BPJPSI_arg",{-M_PI,M_PI}},
            {"EA1_ddcd_BPJPSI_abs",{0.,4.}},
            {"EA1_ddcd_BPJPSI_arg",{-M_PI,M_PI}},
            {"EA2t_ccds_BJPSIP_arg",{-M_PI,M_PI}},
            {"E2t_ccds_BJPSIV_paral_arg",{-M_PI,M_PI}},
            {"G2t_dcs_BJPSIV_0_arg",{-M_PI,M_PI}},
            {"G2t_dcs_BJPSIV_paral_arg",{-M_PI,M_PI}},
            {"G2t_dcs_BJPSIV_perp_arg",{-M_PI,M_PI}},
            {"G1t_scd_BDDb_arg",{-M_PI,M_PI}},
            {"A1_scdd_BDDb_arg",{-M_PI,M_PI}}
        };
        for (const auto& entry:fullExpected)
            check(fullOff.at(entry.first)==entry.second,"incorrect full interval: "+entry.first);
        for (const auto& entry:off)
            if (!fullExpected.count(entry.first))
                check(fullOff.at(entry.first)==entry.second,"unrelated interval changed: "+entry.first);
        check(off.at("G2t_scd_BJPSIP_arg")==Bounds(M_PI/2.,3.*M_PI/2.),"legacy negative region changed");
        check(small.at("G2t_scd_BJPSIP_arg")==Bounds(-M_PI/2.,M_PI/2.),"legacy positive region changed");
        check(off.at("E2t_ccdd_BJPSIP_abs")==Bounds(0.75,4.),"legacy large E2 changed");
        check(small.at("E2t_ccdd_BJPSIP_abs")==Bounds(0.,0.75),"legacy small E2 changed");
        check(small.at("EA1_sdcd_BPJPSI_abs")==Bounds(0.,1.5),"legacy small EA1 changed");
        check(large.at("EA1_sdcd_BPJPSI_abs")==Bounds(1.5,4.),"legacy large EA1 changed");
        std::cout<<"PASS: full domains for "<<fullExpected.size()<<" parameters, precedence over region flags, legacy intervals and Gaussian priors"<<std::endl;
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<std::endl;
        status=1;
    }
    MPI_Finalize();
    return status;
}
