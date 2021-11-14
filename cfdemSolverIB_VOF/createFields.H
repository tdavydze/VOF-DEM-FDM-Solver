    Info<< "Reading field p\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    /*Info<< "\nCreating outside field\n" << endl;
    volScalarField outside
    (
        IOobject
        (
            "outside",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(0)
    );*/


    //========================
    // drag law modelling
    //========================

    Info<< "Reading field phiIB\n" << endl;
    volScalarField phiIB
    (
        IOobject
        (
            "phiIB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field phiIB\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//========================

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    volScalarField gh
    (
        "gh",
        g & mesh.C()
    );

    surfaceScalarField ghf
    (
        "ghf",
        g & mesh.Cf()
    );


    //========================
    // Two-phase Properties
    //========================

    #include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;

    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    //singlePhaseTransportModel laminarTransport(U, phi);
    mesh.setFluxRequired(pd.name());
    mesh.setFluxRequired(alpha1.name());


    Info<< "\nCreating dummy density field rho = 1\n" << endl;
    volScalarField rho                      // NEW !!
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //mesh//,
        alpha1*(rho1 + (scalar(1) - alpha1) * rho2)
    );
    rho.oldTime();

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd+rho*gh
    );



    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, pd, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    if (pd.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        pd = p - rho*gh;
    }

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
    );

    // Mass flux         // NEW!!!
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi   // Original
    //fvc::interpolate(rho)*phi //phi
    );


//=== dyM ===================

    Info<< "Reading field interFace\n" << endl;
    volScalarField interFace
    (
        IOobject
        (
            "interFace",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        //dimensionedScalar("0", dimensionSet(0, -1, 0, 0, 0), 0.0)
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );

//===========================
