/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(dict_, -1),
    rhoMin_(dict_, 0.0),
    rho0_(dict_, 0.0),
    minParcelMass_(dict_, 0.0)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParcelMass_(cp.minParcelMass_)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(dict_, "parcelTypeId", -1),
    rhoMin_(dict_, "rhoMin", 1e-15),
    rho0_(dict_, "rho0"),
    minParcelMass_(dict_, "minParcelMass", 1e-15)
{}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    dr_(0.0),
    noNode_(0.0),
    T_(0.0),
    L0_(0.0),
    growthRate_(0.0),
    I2_(0.0),
    X_(0.0),
    Mom0_(0.0),
    Mom1_(0.0),
    Mom2_(0.0),
    lamda0_(0.0),
    lamda1_(0.0),
    lamda2_(0.0),
    M_(0.0),
    M0_(0.0),
    MWn_(0.0),
    MWw_(0.0),
    PDI_(0.0)
{
    for(int i=0 ; i<40 ; i++)
    {
        F(i) = 0.0 ;
    }
}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli
)
:
    ParcelType(owner, position, celli),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    dr_(0.0),
    noNode_(0.0),
    T_(0.0),
    L0_(0.0),
    growthRate_(0.0),
    I2_(0.0),
    X_(0.0),
    Mom0_(0.0),
    Mom1_(0.0),
    Mom2_(0.0),
    lamda0_(0.0),
    lamda1_(0.0),
    lamda2_(0.0),
    M_(0.0),
    M0_(0.0),
    MWn_(0.0),
    MWw_(0.0),
    PDI_(0.0)
{
    for(int i=0 ; i<40 ; i++)
    {
        F(i) = 0.0 ;
    }
}


template<class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const constantProperties& constProps
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    dTarget_(dTarget0),
    U_(U0),
    rho_(constProps.rho0()),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    dr_(0.0),
    noNode_(0.0),
    T_(0.0),
    L0_(0.0),
    growthRate_(0.0),
    I2_(0.0),
    X_(0.0),
    Mom0_(0.0),
    Mom1_(0.0),
    Mom2_(0.0),
    lamda0_(0.0),
    lamda1_(0.0),
    lamda2_(0.0),
    M_(0.0),
    M0_(0.0),
    MWn_(0.0),
    MWw_(0.0),
    PDI_(0.0)
{
    for(int i=0 ; i<40 ; i++)
    {
        F(i) = 0.0 ;
    }
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::KinematicParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::minParcelMass() const
{
    return minParcelMass_.value();
}


// * * * * * * * KinematicParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::KinematicParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::age() const
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::F(int i) const
{
    return F_[i];
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::dr() const
{
    return dr_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::noNode() const
{
    return noNode_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::T() const
{
    return T_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::L0() const
{
    return L0_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::growthRate() const
{
    return growthRate_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::I2() const
{
    return I2_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::X() const
{
    return X_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Mom0() const
{
    return Mom0_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Mom1() const
{
    return Mom1_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Mom2() const
{
    return Mom2_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::lamda0() const
{
    return lamda0_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::lamda1() const
{
    return lamda1_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::lamda2() const
{
    return lamda2_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::M() const
{
    return M_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::M0() const
{
    return M0_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::MWn() const
{
    return MWn_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::MWw() const
{
    return MWw_;
}
template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::PDI() const
{
    return PDI_;
}

template<class ParcelType>
inline bool& Foam::KinematicParcel<ParcelType>::active()
{
    return active_;
}


template<class ParcelType>
inline Foam::label& Foam::KinematicParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::nParticle()
{
    return nParticle_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::d()
{
    return d_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::dTarget()
{
    return dTarget_;
}

template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::U()
{
    return U_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::rho()
{
    return rho_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::age()
{
    return age_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::tTurb()
{
    return tTurb_;
}

template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb()
{
    return UTurb_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::F(int i)
{
    return F_[i];
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::dr()
{
    return dr_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::noNode()
{
    return noNode_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::T()
{
    return T_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::L0()
{
    return L0_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::growthRate()
{
    return growthRate_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::I2()
{
    return I2_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::X()
{
    return X_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::Mom0()
{
    return Mom0_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::Mom1()
{
    return Mom1_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::Mom2()
{
    return Mom2_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::lamda0()
{
    return lamda0_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::lamda1()
{
    return lamda1_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::lamda2()
{
    return lamda2_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::M()
{
    return M_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::M0()
{
    return M0_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::MWn()
{
    return MWn_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::MWw()
{
    return MWw_;
}
template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::PDI()
{
    return PDI_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::massCell
(
    const trackingData& td
) const
{
    return td.rhoc()*this->mesh().cellVolumes()[this->cell()];
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Re
(
    const trackingData& td
) const
{
    return Re(td.rhoc(), U_, td.Uc(), d_, td.muc());
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Re
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar muc
)
{
    return rhoc*mag(U - Uc)*d/max(muc, rootVSmall);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::We
(
    const trackingData& td,
    const scalar sigma
) const
{
    return We(td.rhoc(), U_, td.Uc(), d_, sigma);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::We
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar sigma
)
{
    return rhoc*magSqr(U - Uc)*d/max(sigma, rootVSmall);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Eo
(
    const trackingData& td,
    const scalar sigma
) const
{
    return Eo(td.g(), rho_, td.rhoc(), U_, d_, sigma);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::Eo
(
    const vector& g,
    const scalar rho,
    const scalar rhoc,
    const vector& U,
    const scalar d,
    const scalar sigma
)
{
    const vector dir = U/max(mag(U), rootVSmall);
    return mag(g & dir)*mag(rho - rhoc)*sqr(d)/max(sigma, rootVSmall);
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::moment(const label no) const
{
    int node = noNode() ;
    scalarList r(node+2) ;
    scalar deltar = dr();
    forAll(r, L)
    {
        r[L] = L * deltar ;
    }
    scalar n1 = 0.0 ;
    scalar n0 = 0.0 ;
    int i = 0;
    
    if(no == 1)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * r[i+1] ;
            n0 += F(i) ;
        }
    }
    else if(no == 2)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * sqr(r[i+1]) ;
            n0 += F(i) ;
        }
    }
    else if(no == 3)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * pow3(r[i+1]) ;
            n0 += F(i) ;
        }
    }
    else if(no == 4)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * pow4(r[i+1]) ;
            n0 += F(i) ;
        }
    }
    else if(no == 5)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * pow5(r[i+1]) ;
            n0 += F(i) ;
        }
    }
    else if(no == 6)
    {
        for(i=0 ; i<node ; i++)
        {
            n1 += F(i) * pow6(r[i+1]) ;
            n0 += F(i) ;
        }
    }
    return n1 / stabilise(n0,SMALL); 
}


// ************************************************************************* //
