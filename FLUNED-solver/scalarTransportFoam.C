/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvModels.H"
#include "fvConstraints.H"

#include "simpleControl.H"

#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

	#include "CourantNo.H"

  bool includeDecayScalar ( runTime.controlDict().lookupOrDefault<bool>("includeDecayScalar", true) );

    
//    volScalarField Dturbulent = nut/Sct;



    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.name() << nl << endl;

	

        while (simple.correctNonOrthogonal())
        {

            // radio-isotope equation for activation
            fvScalarMatrix TaEqn
            (

               fvm::ddt(Ta)
             + fvm::div(phi, Ta)
	     + fvm::Sp(lambda,Ta)
             - fvm::laplacian(DT, Ta)
	     - fvm::laplacian(Dturbulent, Ta)
             ==
	       Source

            );
            TaEqn.relax();
            TaEqn.solve();

            // radio-isotope equation for inlet flow decay
            fvScalarMatrix TdEqn
            (

               fvm::ddt(Td)
             + fvm::div(phi, Td)
	     + fvm::Sp(lambda,Td)
             - fvm::laplacian(DT, Td)
	     - fvm::laplacian(Dturbulent, Td)

            );
            TdEqn.relax();
            TdEqn.solve();


			if (includeDecayScalar ){
				T = Ta + Td;}
			else {
				T = Ta;
			}

            // residence time equation
            fvScalarMatrix TrEqn
            (

               fvm::ddt(Tr)
             + fvm::div(phi, Tr)
     	     - fvm::laplacian(Dturbulent, Tr)
             ==
	       TrSource
            );
            TrEqn.relax();
            TrEqn.solve();
        }



        runTime.write();
    }



    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
