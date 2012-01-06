    forAll (ign.sites(), i)
    {
        const ignitionSite& ignSite = ign.sites()[i];

        if (ignSite.igniting())
        {
            forAll(ignSite.cells(), icelli)
            {
                label ignCell = ignSite.cells()[icelli];
                Info<< "Igniting cell " << ignCell;

                Info<< " state :"
                    << ' ' << b[ignCell]
                    << endl;

                bEqn.diag()[ignSite.cells()[icelli]] +=
                (
                    ignSite.strength()*ignSite.cellVolumes()[icelli]
                   *rho[ignSite.cells()[icelli]]/ignSite.duration()
                )/(b[ignSite.cells()[icelli]] + 0.001);
            }
        }
    }
