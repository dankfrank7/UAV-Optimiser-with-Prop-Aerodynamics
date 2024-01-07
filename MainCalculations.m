    switch Optimisation.Algorithm
        case('MOPSO')
            MOPSOMain
        case('CMOPSO')
            CMOPSOMain
        case('MOEADDE')
            MOEADDEMain
        case('GDE3')
            GDE3Main
        case('NSGA3')
            NSGA3Main
        case('NSGA2')
            NSGA2Main
        case('NMPSO')
            NMPSOMain
        case('MMOPSO')
            MMOPSOMain
        case('LMOCSO')
            LMOCSOMain
        case('C2ODE')
            C2ODEMain
    end