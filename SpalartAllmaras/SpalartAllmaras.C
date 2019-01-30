
#include "SpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(RASModel, SpalartAllmaras, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::chi() const
{
    return rho_*nuTilda_/mu();
}


tmp<volScalarField> SpalartAllmaras::fv1(const volScalarField& chi) const
{
    volScalarField chi3 = pow3(chi);
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmaras::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
    //return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> SpalartAllmaras::fv3
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField chiByCv2 = (1/Cv2_)*chi;

    return
        (scalar(1) + chi*fv1)
       *(1/Cv2_)
       *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
       /pow3(scalar(1) + chiByCv2);
}


tmp<volScalarField> SpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r = min
    (
        nuTilda_
       /(
           max(Stilda, dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
           *sqr(kappa_*d_)
        ),
        scalar(10.0)
    );
    r.boundaryField() == 0.0;

    volScalarField g = r + Cw2_*(pow6(r) - r);

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel
)
:
    RASModel(typeName, rho, U, phi, thermophysicalModel),

    sigmaNut_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Prt_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),

    Cb1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),
    Cv2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),
    Cn1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn1",
            coeffDict_,
            16.0
        )
    ),
    Cn2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn2",
            coeffDict_,
            0.7
        )
    ),
    Cn3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn3",
            coeffDict_,
            0.9
        )
    ),

    StildaModification_
    (
        coeffDict_.lookupOrDefault("StildaModification", false)
    ),

    negativeNuTilda_
    (
        coeffDict_.lookupOrDefault("negativeNuTilda", false)
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    ),

    d_(mesh_)
{
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();

    if (StildaModification_)
    {
        Info<< "    Enabling new Stilda modification" << endl;
    }
    if (negativeNuTilda_)
    {
        Info<< "    Enabling negative nuTilda" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::DnuTildaEff(const volScalarField& chi) const
{
    volScalarField pow3chi = pow(chi, 3);
    volScalarField fn =
        pos(chi) + neg(chi)*(Cn1_ + pow3chi)/(Cn1_ - pow3chi);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "DnuTildaEff", 
            ((rho_*fn*nuTilda_ + mu())/sigmaNut_)
        )
    );
}


tmp<volSymmTensorField> SpalartAllmaras::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k() - (mut_/rho_)*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<volSymmTensorField> SpalartAllmaras::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> SpalartAllmaras::divDevRhoReff(volVectorField& U) const
{
    volScalarField muEff_ = muEff();

    return
    (
      - fvm::laplacian(muEff_, U)
      - fvc::div(muEff_*dev2(fvc::grad(U)().T()))
    );
}


bool SpalartAllmaras::read()
{
    if (RASModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
        Cn1_.readIfPresent(coeffDict());
        Cn2_.readIfPresent(coeffDict());
        Cn3_.readIfPresent(coeffDict());

        StildaModification_.readIfPresent
        (
            "StildaModification", coeffDict()
        );
        negativeNuTilda_.readIfPresent
        (
            "negativeNuTilda", coeffDict()
        );

        return true;
    }
    else
    {
        return false;
    }
}


void SpalartAllmaras::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }

    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*nuTilda_*fv1(chi());
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    if (mesh_.changing())
    {
        d_.correct();
    }

    volScalarField chi = this->chi();
    volScalarField fv1 = this->fv1(chi);

    volTensorField gradU = fvc::grad(U_);
    volSymmTensorField S = symm(gradU);
    volTensorField W = skew(gradU);
    volScalarField Omega = sqrt(2.0)*mag(W);

    volScalarField Stilda =
        fv3(chi, fv1)*Omega
      + fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_);

    // new Stilda modification is used when StildaModification = true
    if (StildaModification_)
    {
        volScalarField Sbar
        (
            fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_)
        );
        Stilda = Omega
               + pos(Cn2_*Omega + Sbar)*Sbar
               + neg(Cn2_*Omega + Sbar)
               *(Omega*(sqr(Cn2_)*Omega + Cn3_*Sbar))
               /max(((Cn3_ - 2.0*Cn2_)*Omega - Sbar),
                    dimensionedScalar("SMALL", Omega.dimensions(), SMALL));
    }

    volScalarField pow3chi = pow(chi, 3);
    volScalarField fn =
        pos(chi) + neg(chi)*(Cn1_ + pow3chi)/(Cn1_ - pow3chi);

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(rho_, nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*rho_*magSqr(fvc::grad(nuTilda_))
      + mu()/rho_*(fvc::grad(rho_)&fvc::grad(nuTilda_))/sigmaNut_
      + fvm::Sp(fn*(fvc::grad(rho_)&fvc::grad(nuTilda_))/sigmaNut_, nuTilda_)
     ==
        pos(nuTilda_)
       *(
        // RHS for nuTilda >= 0
            Cb1_*rho_*Stilda*nuTilda_
          - fvm::Sp(Cw1_*fw(Stilda)*nuTilda_*rho_/sqr(d_), nuTilda_)
        )
      + neg(nuTilda_)
       *(
        // RHS for nuTilda < 0
            Cb1_*rho_*Omega*nuTilda_
          + fvm::Sp(Cw1_*nuTilda_*rho_/sqr(d_), nuTilda_)
        )
    );

    nuTildaEqn().relax();
    solve(nuTildaEqn);
    // bound nuTilda only when negative nuTilda is not active
    if(!negativeNuTilda_)
    {
        // bound nuTilda when negativeNuTilda = false
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }
    nuTilda_.correctBoundaryConditions();
    
    // racalculate fv1 for update of mut
    chi = this->chi();
    fv1 = this->fv1(chi);
    
    // Re-calculate viscosity mag(fv1) is used to ensure that negative nuTilda
    // makes negative mut.
    mut_.internalField() = mag(fv1)*nuTilda_.internalField()*rho_.internalField();
    // bound mut when negative nuTilda is active
    if(negativeNuTilda_)
    {
        bound(mut_, dimensionedScalar("0", mut_.dimensions(), 0.0));
    }
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
