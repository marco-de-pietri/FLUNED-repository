/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;
	
	#include "CourantNo.H"



    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
	
            // radio-isotope equation
            fvScalarMatrix TEqn
            (

               fvm::ddt(T)
             + fvm::div(phi, T)
	     + fvm::Sp(lambda,T)
             - fvm::laplacian(DT, T)
	     - fvm::laplacian(Dturbulent, T)
             ==
	       Source

            );
            TEqn.relax();
            TEqn.solve();

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
