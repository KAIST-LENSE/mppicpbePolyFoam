/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "KinematicParcel.H"
#include "forceSuSp.H"
#include "integrationScheme.H"
#include "meshTools.H"
#include "fluxLimiterMinmod.H"
#include <math.h>
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::KinematicParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    tetIndices tetIs = this->currentTetIndices();

    td.rhoc() = td.rhoInterp().interpolate(this->coordinates(), tetIs);

    if (td.rhoc() < cloud.constProps().rhoMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed density in cell " << this->cell()
                << " to " << cloud.constProps().rhoMin() <<  nl << endl;
        }

        td.rhoc() = cloud.constProps().rhoMin();
    }

    td.Uc() = td.UInterp().interpolate(this->coordinates(), tetIs);

    td.muc() = td.muInterp().interpolate(this->coordinates(), tetIs);

    td.alphaParticlec() = td.alphaParticleInterp().interpolate(this->coordinates(), tetIs);

    td.epsc() = td.epsInterp().interpolate(this->coordinates(), tetIs);


}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::calcDispersion
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() = cloud.dispersion().update
    (
        dt,
        this->cell(),
        U_,
        td.Uc(),
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() += cloud.UTrans()[this->cell()]/massCell(td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(td);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Su, dUTrans, Spu);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;

        // Update momentum transfer coefficient
        cloud.UCoeff()[this->cell()] += np0*Spu;
    }
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::PBE
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    int node = noNode() ;
    scalar deltar = dr();
    scalar phid = td.alphaParticlec() ;
    scalar particleRho = rho() ;
    scalar eps = td.epsc() ;

    scalarList r(node+2) ;
    scalarList fj(node+4) ;
    scalarList fw(node) ;
    scalarList fr0(node) ;
    scalarList fr1(node) ;
    scalarList Birth_b(node) ;
    scalarList Death_b(node) ;
    scalarList Birth_c(node) ;
    scalarList Death_c(node) ;

    forAll(Birth_b, i)
    {
        Birth_b[i] = 0.0 ;
        Death_b[i] = 0.0 ;
        Birth_c[i] = 0.0 ;
        Death_c[i] = 0.0 ;
    }

    forAll(r, L)
    {
        r[L] = L * deltar ;
    }

    fj[0] = 0.0 ;
    fj[1] = 0.0 ;
    forAll(fw, L)
    {
        fj[L+2] = F(L) ;
    }
    fj[node+2] = fj[node+1] ;
    fj[node+3] = fj[node+1] ;

    forAll(fw, L)
    {
        fr0[L] = minmod(fj[L], fj[L+1], fj[L+2], deltar, 1.5) ;
        fr1[L] = minmod(fj[L+1], fj[L+2], fj[L+3], deltar, 1.5) ;
    }

    int j = 0;

    scalarList r_v(2) ;
    scalarList f_v(2) ;
    scalarList betaKernel(2) ;
    scalarList bKernel(2) ;
    scalarList cKernel(2) ;

    forAll(betaKernel, i)
    {
        betaKernel[i] = 0.0;
        bKernel[i] = 0.0;
        cKernel[i] = 0.0;
    }

    forAll(fw, L)
    {
        for(j=L ; j<node ; j++)
        {
            betaKernel[0] = 0.5*4*3.14*sqr(0.5*r[L+1])*4.6/stabilise(pow3(r[j+1]),SMALL)*exp(-4.5*sqr(2.0*pow3(r[L+1])-pow3(r[j+1]))/stabilise(pow6(r[j+1]),SMALL));
            betaKernel[1] = 0.5*4*3.14*sqr(0.5*r[L+1])*4.6/stabilise(pow3(r[j+2]),SMALL)*exp(-4.5*sqr(2.0*pow3(r[L+1])-pow3(r[j+2]))/stabilise(pow6(r[j+2]),SMALL));
            bKernel[0] = 0.01031*pow(eps,0.33333)/stabilise(pow(r[j+1],0.66666),SMALL)/(1.0+phid)*exp(-0.06354*0.013*sqr(1.0+phid)/particleRho/stabilise(pow(r[j+1],1.66666),SMALL)/pow(eps,0.66666));
            bKernel[1] = 0.01031*pow(eps,0.33333)/stabilise(pow(r[j+2],0.66666),SMALL)/(1.0+phid)*exp(-0.06354*0.013*sqr(1.0+phid)/particleRho/stabilise(pow(r[j+2],1.66666),SMALL)/pow(eps,0.66666));

            Birth_b[L] += deltar * betaKernel[0]*bKernel[0]*2.0*fj[j+2] + deltar/2 * (betaKernel[1]*bKernel[1]*2.0*fj[j+3] - betaKernel[0]*bKernel[0]*2.0*fj[j+2]) ;
        }

        bKernel[0] = 0.01031*pow(eps,0.33333)/stabilise(pow(r[L+1],0.66666),SMALL)/(1.0+phid)*exp(-0.06354*0.013*sqr(1.0+phid)/particleRho/stabilise(pow(r[L+1],1.66666),SMALL)/pow(eps,0.66666));
        Death_b[L] = bKernel[0]*fj[L+2];
    }

    scalar gR = growthRate() ;

    forAll(fw, L)
    {
        fj[L+2] = fj[L+2] + dt * (-(1/deltar)*(gR*(fj[L+2] + 0.5*deltar * fr1[L]) - gR*(fj[L+1] + 0.5*deltar * fr0[L])) + Birth_b[L] - Death_b[L]) ; 
        F(L) = fj[L+2] ;
    }

    scalar n3 = 0.0 ;
    scalar n2 = 0.0 ;
    forAll(fw, i)
    {
        n3 += fj[i+2] * pow3(r[i+1]) ;
        n2 += fj[i+2] * sqr(r[i+1]) ;
    }

    d() = n3 / stabilise(n2,SMALL);


    forAll(fw, i)
    {
        fw[i] = deltar * particleRho * pi / 6.0 * pow3(r[i+1]) * fj[i+2] + deltar/2 * (particleRho * pi / 6.0 * pow3(r[i+2]) * fj[i+3] - particleRho * pi / 6.0 * pow3(r[i+1]) * fj[i+2]);
    }

    scalar total_mass = 0.0 ;
    forAll(fw, i)
    {
	total_mass += fw[i] ;
    }

    scalar total_mass_cell(total_mass * cloud.pMesh().cellVolumes()[this->cell()]) ;
    scalar noParticle(total_mass_cell / mass()) ;

    int roundToInt(noParticle + 0.5);
    nParticle() = max(1, roundToInt);

//    Info<< "noparticle = "<< noParticle << endl;
//    Info<< "total_v_cell = "<< total_v_cell << endl;
//    Info<< "volume = "<< volume() << endl;
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::polymerization
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{

    scalar r_const = 1.987 ; // J/mol/K
    scalar mw_m = 100.121 ; // g/mol
    scalar Tp = T() ;

    scalar rho_m = 968.0 - 1.15*(Tp-273.15);
    scalar rho_p = 1212.0 - 0.845*(Tp-273.15);
    scalar volume_f = (rho_p - rho_m)/rho_p;
    M() = M0()*(1-X()) / (1+volume_f*X()) ; 

    scalar f_propa = 1.0 ; // -
    scalar k_d = 1.7e14*exp(-3.0e4 / r_const / Tp) ; // 1/s
    scalar k_trm = 0.0 ; 
    scalar k_td_0 = 9.8e7*exp(-701 / r_const / Tp);
    scalar k_tc_0 = 0.0 ;
    scalar k_p_0 = 4.92e5*exp(-4353/ r_const / Tp);

    scalar g_t = 0.0;
    scalar g_p = 0.0;
    scalar v_f_m = 0.025 + 0.001*(Tp-167);
    scalar v_f_p = 0.025 + 0.00048*(Tp-387);
    scalar v_f_cr = 0.1856 - 2.965e-4*(Tp-273.15);
    scalar M_fraction = M()/(M()+lamda0()+Mom0());//*mw_m/rho_m / (M()/(M()+lamda0()+Mom0())*mw_m/rho_m + (lamda0()+Mom0())/(M()+lamda0()+Mom0())*MWw()/rho_p);
    scalar v_f = M_fraction*v_f_m + (1-M_fraction)*v_f_p;

    if(v_f>v_f_cr)
    {  
        g_t = 0.10575*exp(17.15*v_f - 0.01715*(Tp-273.15));
        g_p = 1.0;
    }
    else
    {
        g_t = 0.23e-5*exp(75*v_f);
        g_p = 0.71e-4*exp(71.53*v_f);
    }

    scalar k_td = k_td_0 * g_t ;
    scalar k_tc = k_tc_0 * g_t ;
    scalar k_p = k_p_0 * g_p ;

    scalarList C(8) ;
    scalarList C0(8) ;
    scalarList preC(8) ;
    scalarList Func(8) ;
    scalarList dFunc(8) ;
    scalar error ;

    C[0] = C0[0] = I2() ;
    C[1] = C0[1] = X() ;
    C[2] = C0[2] = lamda0() ;
    C[3] = C0[3] = lamda1() ;
    C[4] = C0[4] = lamda2() ;
    C[5] = C0[5] = Mom0() ;
    C[6] = C0[6] = Mom1() ;
    C[7] = C0[7] = Mom2() ;

    scalar CM0 = M0() ;
    scalar CM = CM0*(1-C[1]) / (1+volume_f*C[1]) ;
    scalar dxdt = 0.0 ;

    do
    {
        CM = CM0*(1-C[1]) / (1+volume_f*C[1]) ; 
        dxdt = 2*f_propa*k_d*C[0]*(1+volume_f*C[1])/CM0 + (k_p + k_trm)*(1 - C[1])*C[2];       

        Func[0] = C[0] - C0[0] - dt*(-(k_d + volume_f/(1+volume_f*C[1])*dxdt)*C[0]);
        Func[1] = C[1] - C0[1] - dt*(2*f_propa*k_d*C[0]*(1+volume_f*C[1])/CM0 + (k_p + k_trm)*(1-C[1])*C[2]);
        Func[2] = C[2] - C0[2] - dt*(2*f_propa*k_d*C[0] - (k_td + k_tc)*pow(C[2],2.0) - C[2]*volume_f/(1+volume_f*C[1])*dxdt);
        Func[3] = C[3] - C0[3] - dt*(2*f_propa*k_d*C[0] + k_p*CM*C[2] + k_trm*CM*(C[2]-C[3]) - (k_td+k_tc)*C[2]*C[3] - C[3]*volume_f/(1+volume_f*C[1])*dxdt);
        Func[4] = C[4] - C0[4] - dt*(2*f_propa*k_d*C[0] + k_p*CM*(2*C[3]+C[2]) + k_trm*CM*(C[2]-C[4]) - (k_td+k_tc)*C[2]*C[4] - C[4]*volume_f/(1+volume_f*C[1])*dxdt);
        Func[5] = C[5] - C0[5] - dt*(k_trm*CM*C[2] + (k_td+0.5*k_tc)*pow(C[2],2.0) - C[5]*volume_f/(1+volume_f*C[1])*dxdt);
        Func[6] = C[6] - C0[6] - dt*(k_trm*CM*C[3] + (k_td+k_tc)*C[2]*C[3] - C[6]*volume_f/(1+volume_f*C[1])*dxdt);
        Func[7] = C[7] - C0[7] - dt*(k_trm*CM*C[4] + k_td*C[2]*C[4] + k_tc*(C[2]*C[4]+pow(C[3],2.0)) - C[7]*volume_f/(1+volume_f*C[1])*dxdt);

        dFunc[0] = 1 - dt*(-(k_d + volume_f/(1+volume_f*C[1])*dxdt));
        dFunc[1] = 1 - dt*(2*f_propa*k_d*C[0]*volume_f/CM0 - (k_p + k_trm)*C[2]);
        dFunc[2] = 1 - dt*(-2*(k_td + k_tc)*C[2] - volume_f/(1+volume_f*C[1])*dxdt);
        dFunc[3] = 1 - dt*(-k_trm*CM - (k_td+k_tc)*C[2] - volume_f/(1+volume_f*C[1])*dxdt);
        dFunc[4] = 1 - dt*(-k_trm*CM - volume_f/(1+volume_f*C[1])*dxdt);
        dFunc[5] = 1 - dt*(-volume_f/(1+volume_f*C[1])*dxdt);
        dFunc[6] = 1 - dt*(-volume_f/(1+volume_f*C[1])*dxdt);
        dFunc[7] = 1 - dt*(-volume_f/(1+volume_f*C[1])*dxdt);

        forAll(C, i)
        {
            preC[i] = C[i];
            C[i] = C[i] - Func[i]/dFunc[i] ;
        }

        error = 0 ;
        forAll(Func, i)
        {
            error += abs((preC[i] - C[i])/stabilise(preC[i],SMALL)) ;
        }
       
    }while(error > 1e-5);

    I2() = C[0] ;
    X() = C[1] ;
    lamda0() = C[2] ;
    lamda1() = C[3] ;
    lamda2() = C[4] ;
    Mom0() = C[5] ;
    Mom1() = C[6] ;
    Mom2() = C[7] ;
    M() = CM0*(1-C[1]) / (1+volume_f*C[1]) ;
    MWn() = (Mom1()+lamda1())/stabilise((lamda0()+Mom0()),SMALL)*mw_m;
    MWw() = (Mom2()+lamda2())/stabilise((lamda1()+Mom1()),SMALL)*mw_m;
    PDI() = MWw()/stabilise(MWn(),SMALL);
    growthRate() = volume_f*L0()/3*dxdt ;
    rho() = (I2()*242.23 + M()*100.121 + (Mom0()+lamda0())*MWw()) / (I2()*242.23/1330.0 + M()*100.121/rho_m + (Mom0()+lamda0())*MWw()/rho_p);

}

template<class ParcelType>
template<class TrackCloudType>
const Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
) const
{
    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const typename TrackCloudType::forceType& forces = cloud.forces();

    // Momentum source due to particle forces
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, mu);
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, mu);
    const scalar massEff = forces.massEff(p, ttd, mass);

    /*
    // Proper splitting ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Sp()*td.Uc() + Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;
    const scalar bncp = Fncp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaUcp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, acp, bcp
        );
    const vector deltaUncp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, ancp, bncp
        );
    const vector deltaT = deltaUcp + deltaUncp;
    */

    // Shortcut splitting assuming no implicit non-coupled force ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaU = cloud.UIntegrator().delta(U_, dt, acp + ancp, bcp);
    const vector deltaUncp = ancp*dt;
    const vector deltaUcp = deltaU - deltaUncp;

    // Calculate the new velocity and the momentum transfer terms
    vector Unew = U_ + deltaU;

    dUTrans -= massEff*deltaUcp;

    Spu = dt*Fcp.Sp();

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = cloud.pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_)
{}


template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime,
    const label withPBE
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    ttd.switchProcessor = false;
    ttd.keepParticle = true;

    const scalarField& cellLengthScale = cloud.cellLengthScale();
    const scalar maxCo = cloud.solution().maxCo();

    while (ttd.keepParticle && !ttd.switchProcessor && p.stepFraction() < 1)
    {
        // Cache the current position, cell and step-fraction
        const point start = p.position();
        const scalar sfrac = p.stepFraction();

        // Total displacement over the time-step
        const vector s = trackTime*U_;

        // Cell length scale
        const scalar l = cellLengthScale[p.cell()];

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        // Fraction of the displacement to track in this loop. This is limited
        // to ensure that the both the time and distance tracked is less than
        // maxCo times the total value.
        scalar f = 1 - p.stepFraction();
        f = min(f, maxCo);
        f = min(f, maxCo*l/max(small*l, mag(s)));
        if (p.active())
        {
            // Track to the next face
            p.trackToFace(f*s - d, f);
        }
        else
        {
            // At present the only thing that sets active_ to false is a stick
            // wall interaction. We want the position of the particle to remain
            // the same relative to the face that it is on. The local
            // coordinates therefore do not change. We still advance in time and
            // perform the relevant interactions with the fixed particle.
            p.stepFraction() += f;
        }

        const scalar dt = (p.stepFraction() - sfrac)*trackTime;

        // Avoid problems with extremely small timesteps
        if (dt > rootVSmall)
        {
            // Update cell based properties
            p.setCellValues(cloud, ttd);

            if (withPBE == 0)
            {
                p.polymerization(cloud, ttd, dt);
                p.PBE(cloud, ttd, dt);
            }

            p.calcDispersion(cloud, ttd, dt);

            if (cloud.solution().cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(cloud, ttd, dt);
            }

            p.calc(cloud, ttd, dt);
        }

        p.age() += dt;

        if (p.active() && p.onFace())
        {
            cloud.functions().postFace(p, ttd.keepParticle);
        }

        cloud.functions().postMove(p, dt, start, ttd.keepParticle);

        if (p.active() && p.onFace() && ttd.keepParticle)
        {
            p.hitFace(f*s - d, f, cloud, ttd);
        }
    }

    return ttd.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicParcel<ParcelType>::hitPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    // Invoke post-processing model
    cloud.functions().postPatch(p, pp, td.keepParticle);

    // Invoke surface film model
    if (cloud.surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // All interactions done
        return true;
    }
    else if (pp.coupled())
    {
        // Don't apply the patchInteraction models to coupled boundaries
        return false;
    }
    else
    {
        // Invoke patch interaction model
        return cloud.patchInteraction().correct(p, pp, td.keepParticle);
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    TrackCloudType&,
    trackingData&
)
{
    // wall interactions are handled by the generic hitPatch method
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "KinematicParcelIO.C"

// ************************************************************************* //
