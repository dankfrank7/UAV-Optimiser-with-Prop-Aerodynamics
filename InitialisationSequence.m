switch Optimisation.Algorithm
    case('MOPSO')
        MOPSOInit
    case('CMOPSO')
        MOPSOInit
    case('MOEADDE')
        if nPop > 50
            nPop = 50
        end
        MOEADDEInit
    case('GDE3')
        GDE3Init
    case('NSGA3')
        NSGA3Init
    case('NSGA2')
        NSGA2Init
    case('NMPSO')
        NMPSOInit
    case('MMOPSO')
        MMOPSOInit
    case('LMOCSO')
        LMOCSOInit
    case('C2ODE')
        C2ODEInit
    otherwise
        error('the algorithm you selected is not set up yet')
end

if isfield(Optimisation,'Algorithm2')
    switch Optimisation.Algorithm2
    case('MOPSO')
        MOPSOInit
    case('CMOPSO')
        MOPSOInit
    case('MOEADDE')
        if nPop > 50
            nPop = 50
        end
        MOEADDEInit
    case('GDE3')
        GDE3Init
    case('NSGA3')
        NSGA3Init
    case('NSGA2')
        NSGA2Init
    case('NMPSO')
        NMPSOInit
    case('MMOPSO')
        MMOPSOInit
    case('LMOCSO')
        LMOCSOInit
    case('C2ODE')
        C2ODEInit
    otherwise
        error('the algorithm you selected is not set up yet')
    end
end

if isfield(Optimisation,'Algorithm3')
    switch Optimisation.Algorithm3
    case('MOPSO')
        MOPSOInit
    case('CMOPSO')
        MOPSOInit
    case('MOEADDE')
        if nPop > 50
            nPop = 50
        end
        MOEADDEInit
    case('GDE3')
        GDE3Init
    case('NSGA3')
        NSGA3Init
    case('NSGA2')
        NSGA2Init
    case('NMPSO')
        NMPSOInit
    case('MMOPSO')
        MMOPSOInit
    case('LMOCSO')
        LMOCSOInit
    case('C2ODE')
        C2ODEInit
    otherwise
        error('the algorithm you selected is not set up yet')
    end
end