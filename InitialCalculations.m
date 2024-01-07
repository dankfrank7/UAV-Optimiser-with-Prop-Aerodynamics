% This file contains all the initial calculations for each of the
% implemented algorithms

switch Optimisation.Algorithm
    case('MOPSO')
        Archive     = UpdateArchive(pop,nPop,MOPSO.nGrid);
        Pbest       = pop;
    case('MOEADDE')
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
    case('NSGA3')
        [Z,nPop]    = UniformPoint(nPop,nObj);
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        PopCon      = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        Zmin        = min(PopObj(all(PopCon<=0,2),:),[],1);
    case('NSGA2')
        [~,FrontNo,CrowdDis,~] = EnvironmentalSelectionNSGA2(pop ,nPop);
    case('NMPSO')
        Pbest       = pop;
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Archive     = UpdateArchiveNMPSO(pop(NDSort(PopObj,1)==1),[],nPop);
    case('MMOPSO')
        [W,nPop]    = UniformPoint(nPop,nObj);
        W           = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
        [pop,ind]   = ClassificationMMOPSO(pop,W,Z);
        if ~isempty(ind)
            popn = pop(ind);
            parfor i=1:length(popn)
                switch runcase
                    case('test')
                        popn(i).Cost      = CostFunction(popn(i).Position);
                    case('full')
                        popn(i).Cost       = CostFunction(popn(i).Position,AoA,density,velocity,fields,VarMin,VarMax,i,OptSettings);
                end
            end
            pop(ind) = popn;
        end
        Archive     = UpdateArchiveMMOPSO(pop,nPop);
    case('LMOCSO')
        [V,nPop]    = UniformPoint(nPop,nObj);
        pop         = EnvironmentalSelectionLMOCSO(pop,V,(it/MaxIt)^2);
    case('C2ODE')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);
        Gbest       = pop;
    case('GDE3')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);

end

if isfield(Optimisation,'Algorithm2')
    switch Optimisation.Algorithm2
    case('MOPSO')
        Archive     = UpdateArchive(pop,nPop,MOPSO.nGrid);
        Pbest       = pop;
    case('MOEADDE')
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
    case('NSGA3')
        [Z,nPop]    = UniformPoint(nPop,nObj);
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        PopCon      = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        Zmin        = min(PopObj(all(PopCon<=0,2),:),[],1);
    case('NSGA2')
        [~,FrontNo,CrowdDis,~] = EnvironmentalSelectionNSGA2(pop ,nPop);
    case('NMPSO')
        Pbest       = pop;
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Archive     = UpdateArchiveNMPSO(pop(NDSort(PopObj,1)==1),[],nPop);
    case('MMOPSO')
        [W,nPop]    = UniformPoint(nPop,nObj);
        W           = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
        [pop,ind]   = ClassificationMMOPSO(pop,W,Z);
        if ~isempty(ind)
            popn = pop(ind);
            parfor i=1:length(popn)
                switch runcase
                    case('test')
                        popn(i).Cost      = CostFunction(popn(i).Position);
                    case('full')
                        popn(i).Cost       = CostFunction(popn(i).Position,AoA,density,velocity,fields,VarMin,VarMax,i,OptSettings);
                end
            end
            pop(ind) = popn;
        end
        Archive     = UpdateArchiveMMOPSO(pop,nPop);
    case('LMOCSO')
        [V,nPop]    = UniformPoint(nPop,nObj);
        pop         = EnvironmentalSelectionLMOCSO(pop,V,(it/MaxIt)^2);
    case('C2ODE')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);
        Gbest       = pop;
    case('GDE3')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);

    end
end

if isfield(Optimisation,'Algorithm3')
    switch Optimisation.Algorithm3
    case('MOPSO')
        Archive     = UpdateArchive(pop,nPop,MOPSO.nGrid);
        Pbest       = pop;
    case('MOEADDE')
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
    case('NSGA3')
        [Z,nPop]    = UniformPoint(nPop,nObj);
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        PopCon      = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        Zmin        = min(PopObj(all(PopCon<=0,2),:),[],1);
    case('NSGA2')
        [~,FrontNo,CrowdDis,~] = EnvironmentalSelectionNSGA2(pop ,nPop);
    case('NMPSO')
        Pbest       = pop;
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Archive     = UpdateArchiveNMPSO(pop(NDSort(PopObj,1)==1),[],nPop);
    case('MMOPSO')
        [W,nPop]    = UniformPoint(nPop,nObj);
        W           = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
        PopObj      = reshape(extractfield(pop,'Cost'),[],size(pop,1))';
        Z           = min(PopObj,[],1);
        [pop,ind]   = ClassificationMMOPSO(pop,W,Z);
        if ~isempty(ind)
            popn = pop(ind);
            parfor i=1:length(popn)
                switch runcase
                    case('test')
                        popn(i).Cost      = CostFunction(popn(i).Position);
                    case('full')
                        popn(i).Cost       = CostFunction(popn(i).Position,AoA,density,velocity,fields,VarMin,VarMax,i,OptSettings);
                end
            end
            pop(ind) = popn;
        end
        Archive     = UpdateArchiveMMOPSO(pop,nPop);
    case('LMOCSO')
        [V,nPop]    = UniformPoint(nPop,nObj);
        pop         = EnvironmentalSelectionLMOCSO(pop,V,(it/MaxIt)^2);
    case('C2ODE')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);
        Gbest       = pop;
    case('GDE3')
        PopCon    = reshape(extractfield(pop,'Constraint'),[],size(pop,1))';
        ConV        = sum(max(0, PopCon), 2);
        Var0        = max(ConV);

    end
end