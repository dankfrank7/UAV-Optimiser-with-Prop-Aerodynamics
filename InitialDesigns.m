switch Optimisation.initialdesign
    case('latin')
        InitialGeneration               = lhsdesign(ceil(nPop/2),nVar,'iterations',100+999*rand);
        InitialGeneration               = [InitialGeneration;1-InitialGeneration];
        if (rem(nPop,2) == 1)
            InitialGeneration           = InitialGeneration(1:end-1,:);
        end
    case('sobol')
        p                               = sobolset(nVar,'Skip',1e3,'Leap',1e2);
        InitialGeneration               = net(p,nPop);
    case('restart')
        load(restartfile)
        InitialGeneration               = reshape(extractfield(Pbest,'Position'),[],size(Pbest,1))';
    otherwise
        error('unknown initialisation method')
        
end
