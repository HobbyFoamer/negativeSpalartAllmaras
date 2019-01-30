/*
Description
    Spalart-Allmaras one-eqn mixing-length model for compressible
    external flows.
    Modified SA model which allows negative nuTilda.
    Refered as SA-neg in NASA's documentation
    (https://turbmodels.larc.nasa.gov/spalart.html).
    
    Reference:
    @verbatim
        "Modifications and Clarifications for the Implementation
         of the Spalart-Allmaras Turbulence Model"
        S.R. Allmaras,
        F.T. Johnson,
        P.R. Spalart
        Seventh International Conference on Computational Fluid Dynamics 
        (ICCFD7), 2012, pp. 1-11.
    @endverbatim

    Additional model coefficients correspond to the following:
    @verbatim
        SpalartAllmarasCoeffs
        {
            StildaModification  false;
            negativeNuTilda     false;
            Cn1                 16.0;
            Cn2                 0.7;
            Cn3                 0.9;
        }
    @endverbatim
    
    Additional modifications are activated by switches.
    Therefore this model is the same as the original SpalartAllmaras of
    foam-extend-3.2 unless you activate the modifications in RASProperties.
*/

#ifndef compressibleSpalartAllmaras_H
#define compressibleSpalartAllmaras_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            dimensionedScalar Prt_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;
            
        // Additional model coefficients
            dimensionedScalar Cn1_;
            dimensionedScalar Cn2_;
            dimensionedScalar Cn3_;

        //- Optional flag to activate the new STilda modification
        Switch StildaModification_;

        //- Optional flag to activate the negative nuTilda model
        Switch negativeNuTilda_;

        // Fields

            volScalarField nuTilda_;
            volScalarField mut_;
            volScalarField alphat_;


        //- Wall distance
        wallDist d_;


    // Private member functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Construct from components
        SpalartAllmaras
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel
        );


    //- Destructor
    virtual ~SpalartAllmaras()
    {}


    // Member Functions

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const
/* moved to .C file because DnuTildaEff depends on whether
   nuTilda is pos/neg in SA-neg.
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DnuTildaEff",
                    rho_*nuTilda_/sigmaNut_ + mu()
                )
            );
        }
*/
        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the effective turbulent thermal diffusivity
        virtual tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphat_ + alpha())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epslion",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
