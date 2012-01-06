    if (runTime.outputTime())
    {
        OFstream averagingTimeFile
        (
            runTime.path()/runTime.timeName()/"averagingTime.raw"
        );

        if (averagingTimeFile.good())
        {
            averagingTimeFile << averagingTime << endl;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file "
                << runTime.path()/runTime.timeName()/"averagingTime.raw"
                << abort(FatalError);
        }
    }
