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

#include "KinematicParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propertyList_ =
    Foam::KinematicParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::KinematicParcel<ParcelType>::sizeofFields_
(
    sizeof(KinematicParcel<ParcelType>)
  - offsetof(KinematicParcel<ParcelType>, active_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
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
    for(int i=0 ; i<noNode_ ; i++)
    {
        F(i) = 0.0 ;
    }
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
            dr_ = readScalar(is);
            noNode_ = readScalar(is);
            T_ = readScalar(is);
            L0_ = readScalar(is);
            for(int i=0 ; i<noNode_ ; i++)
            {
                F_[i] = readScalar(is);
            }
	    I2_ = readScalar(is);
	    X_ = readScalar(is);
	    Mom0_ = readScalar(is);
	    Mom1_ = readScalar(is);
	    Mom2_ = readScalar(is);
	    lamda0_ = readScalar(is);
	    lamda1_ = readScalar(is);
	    lamda2_ = readScalar(is);
	    M_ = readScalar(is);
	    M0_ = readScalar(is);
	    MWn_ = readScalar(is);
	    MWw_ = readScalar(is);
	    PDI_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicParcel<ParcelType>::KinematicParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::readFields(CloudType& c)
{
    bool write = c.size();

    ParcelType::readFields(c);

    IOField<label> active
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, active);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget
    (
        c.fieldIOobject("dTarget", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U
    (
        c.fieldIOobject("U", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb
    (
        c.fieldIOobject("tTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.fieldIOobject("UTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, UTurb);

    c.checkFieldIOobject(c, UTurb);

    IOField<scalar> dr
    (
        c.fieldIOobject("dr", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, dr);

    IOField<scalar> noNode
    (
        c.fieldIOobject("noNode", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, noNode);

    IOField<scalar> T
    (
        c.fieldIOobject("T", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, T);

    IOField<scalar> L0
    (
        c.fieldIOobject("L0", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, L0);

    IOField<scalarField> F
    (
        c.fieldIOobject("F", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, F);

    IOField<scalar> I2
    (
        c.fieldIOobject("I2", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, I2);
    IOField<scalar> X
    (
        c.fieldIOobject("X", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, X);
    IOField<scalar> Mom0
    (
        c.fieldIOobject("Mom0", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, Mom0);
    IOField<scalar> Mom1
    (
        c.fieldIOobject("Mom1", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, Mom1);
    IOField<scalar> Mom2
    (
        c.fieldIOobject("Mom2", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, Mom2);
    IOField<scalar> lamda0
    (
        c.fieldIOobject("lamda0", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, lamda0);
    IOField<scalar> lamda1
    (
        c.fieldIOobject("lamda1", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, lamda1);
    IOField<scalar> lamda2
    (
        c.fieldIOobject("lamda2", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, lamda2);
    IOField<scalar> M
    (
        c.fieldIOobject("M", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, M);
    IOField<scalar> M0
    (
        c.fieldIOobject("M0", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, M0);
    IOField<scalar> MWn
    (
        c.fieldIOobject("MWn", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, MWn);
    IOField<scalar> MWw
    (
        c.fieldIOobject("MWw", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, MWw);
    IOField<scalar> PDI
    (
        c.fieldIOobject("PDI", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, PDI);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        KinematicParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        p.dr_ = dr[i];
        p.noNode_ = noNode[i];
        p.T_ = T[i];
        p.L0_ = L0[i];
        for(int j=0 ; j<p.noNode() ; j++)
        {
            p.F_[j] = F[i][j];
        } 
        p.I2_ = I2[i];
        p.X_ = X[i];
        p.Mom0_ = Mom0[i];
        p.Mom1_ = Mom1[i];
        p.Mom2_ = Mom2[i];
        p.lamda0_ = lamda0[i];
        p.lamda1_ = lamda1[i];
        p.lamda2_ = lamda2[i];
        p.M_ = M[i];
        p.M0_ = M0[i];
        p.MWn_ = MWn[i];
        p.MWw_ = MWw[i];
        p.PDI_ = PDI[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);
    IOField<scalarField> F(c.fieldIOobject("F", IOobject::NO_READ), np);
    IOField<scalar> dr(c.fieldIOobject("dr", IOobject::NO_READ), np);
    IOField<scalar> noNode(c.fieldIOobject("noNode", IOobject::NO_READ), np);
    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> L0(c.fieldIOobject("L0", IOobject::NO_READ), np);
    IOField<scalar> I2(c.fieldIOobject("I2", IOobject::NO_READ), np);
    IOField<scalar> X(c.fieldIOobject("X", IOobject::NO_READ), np);
    IOField<scalar> Mom0(c.fieldIOobject("Mom0", IOobject::NO_READ), np);
    IOField<scalar> Mom1(c.fieldIOobject("Mom1", IOobject::NO_READ), np);
    IOField<scalar> Mom2(c.fieldIOobject("Mom2", IOobject::NO_READ), np);
    IOField<scalar> lamda0(c.fieldIOobject("lamda0", IOobject::NO_READ), np);
    IOField<scalar> lamda1(c.fieldIOobject("lamda1", IOobject::NO_READ), np);
    IOField<scalar> lamda2(c.fieldIOobject("lamda2", IOobject::NO_READ), np);
    IOField<scalar> M(c.fieldIOobject("M", IOobject::NO_READ), np);
    IOField<scalar> M0(c.fieldIOobject("M0", IOobject::NO_READ), np);
    IOField<scalar> MWn(c.fieldIOobject("MWn", IOobject::NO_READ), np);
    IOField<scalar> MWw(c.fieldIOobject("MWw", IOobject::NO_READ), np);
    IOField<scalar> PDI(c.fieldIOobject("PDI", IOobject::NO_READ), np);
    IOField<scalar> moment1(c.fieldIOobject("Moment1", IOobject::NO_READ), np);
    IOField<scalar> moment2(c.fieldIOobject("Moment2", IOobject::NO_READ), np);
    IOField<scalar> moment3(c.fieldIOobject("Moment3", IOobject::NO_READ), np);
    IOField<scalar> moment4(c.fieldIOobject("Moment4", IOobject::NO_READ), np);
    IOField<scalar> moment5(c.fieldIOobject("Moment5", IOobject::NO_READ), np);
    IOField<scalar> moment6(c.fieldIOobject("Moment6", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();
        scalarField dataF(p.noNode());
        for(int j=0 ; j<p.noNode_ ; j++)
        {
            dataF[j] = p.F(j) ;
        }
        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        F[i] = dataF;
        dr[i] = p.dr();
        noNode[i] = p.noNode();
        T[i] = p.T();
        L0[i] = p.L0();
        I2[i] = p.I2();
        X[i] = p.X();
        Mom0[i] = p.Mom0();
        Mom1[i] = p.Mom1();
        Mom2[i] = p.Mom2();
        lamda0[i] = p.lamda0();
        lamda1[i] = p.lamda1();
        lamda2[i] = p.lamda2();
        M[i] = p.M();
        M0[i] = p.M0();
        MWn[i] = p.MWn();
        MWw[i] = p.MWw();
        PDI[i] = p.PDI();
        moment1[i] = p.moment(1);
        moment2[i] = p.moment(2);
        moment3[i] = p.moment(3);
        moment4[i] = p.moment(4);
        moment5[i] = p.moment(5);
        moment6[i] = p.moment(6);

        i++;
    }

    const bool write = np > 0;

    active.write(write);
    typeId.write(write);
    nParticle.write(write);
    d.write(write);
    dTarget.write(write);
    U.write(write);
    rho.write(write);
    age.write(write);
    tTurb.write(write);
    UTurb.write(write);
    F.write(write);
    dr.write(write);
    noNode.write(write);
    T.write(write);
    L0.write(write);
    I2.write(write);
    X.write(write);
    Mom0.write(write);
    Mom1.write(write);
    Mom2.write(write);
    lamda0.write(write);
    lamda1.write(write);
    lamda2.write(write);
    M.write(write);
    M0.write(write);
    MWn.write(write);
    MWw.write(write);
    PDI.write(write);
    moment1.write(write);
    moment2.write(write);
    moment3.write(write);
    moment4.write(write);
    moment5.write(write);
    moment6.write(write);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb()
            << token::SPACE << p.dr()
            << token::SPACE << p.noNode()
            << token::SPACE << p.T()
            << token::SPACE << p.L0()
            << token::SPACE << p.I2()
            << token::SPACE << p.X()
            << token::SPACE << p.Mom0()
            << token::SPACE << p.Mom1()
            << token::SPACE << p.Mom2()
            << token::SPACE << p.lamda0()
            << token::SPACE << p.lamda1()
            << token::SPACE << p.lamda2()
            << token::SPACE << p.M()
            << token::SPACE << p.M0()
            << token::SPACE << p.MWn()
            << token::SPACE << p.MWw()
            << token::SPACE << p.PDI();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            KinematicParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
