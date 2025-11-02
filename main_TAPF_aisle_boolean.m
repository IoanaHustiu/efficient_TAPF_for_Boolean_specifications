close all;
clear;
clc;
fclose('all');

addpath(['.' filesep 'functions']);
addpath(['.' filesep 'maps']);

n_exp = 20; %number of experiments to be performed
N_robots = 100; %number of robots for experiments

B = load_map('warehouse-10-20-10-2-1.map');
[robotPts, scenTbl] = loadAllScens('warehouse-10-20-10-2-1-random-.scen', 25);

% Convert coordinates from (0,0 top-left) → (0,0 bottom-left)
[H, W] = size(B);

for i = 1:numel(robotPts)
    p = robotPts{i};
    % y_cartesian = (height - 1) - y_original
    p(:,2) = (H - 1) - p(:,2);
    robotPts{i} = p;
end

%% --- Extrae todos los puntos de inicio y final
starts = cellfun(@(p) p(1,:), robotPts, 'UniformOutput', false);
goals  = cellfun(@(p) p(2,:), robotPts, 'UniformOutput', false);

starts = vertcat(starts{:});
goals  = vertcat(goals{:});

% --- Selección según coordenada x
mask_init = (starts(:,1) < 25);
mask_goal = (goals(:,1)  > 28) & (goals(:,1)  < 130);

% --- Puntos iniciales y finales candidatos
initPts_all = starts(mask_init, :);
goalPts_all = goals(mask_goal, :);

% --- Asegurar que no hay puntos repetidos
[~, ia] = unique(initPts_all, 'rows');
initPts = num2cell(initPts_all(ia, :), 2);

[~, ig] = unique(goalPts_all, 'rows');
goalPts = num2cell(goalPts_all(ig, :), 2);

% --- Eliminar solapamiento entre conjuntos (por si acaso)
initMat = vertcat(initPts{:});
goalMat = vertcat(goalPts{:});

[~, ia] = setdiff(initMat, goalMat, 'rows', 'stable');
[~, ig] = setdiff(goalMat, initMat, 'rows', 'stable');

initPts = initPts(ia);
goalPts = goalPts(ig);

fprintf('Generated %d unique initial points and %d unique final points.\n', ...
        numel(initPts), numel(goalPts));

% Flip B vertically to match Cartesian coordinates
B = flipud(B);

T = build_topology_from_B(B);                     % T.adj (NxN, sparse), T.free, T.size, ...
T.map2D = B;
% Extract indices of free cells
fullAdj = T.adj;
N       = size(fullAdj, 1);
if isfield(T, 'rem_cells') && ~isempty(T.rem_cells)
    freeIdx = T.rem_cells(:);
else
    error('Invalid T structure: T.rem_cells is missing or empty.');
end
% Compute (x,y) coordinates of the center of each free cell
for k = 1:numel(freeIdx)
    [r, c] = ind2sub([H, W], freeIdx(k));
    T.centr{k} = [c - 0.5, r - 0.5];   % <-- center, not corner
end

% Subgraph on free cells (sparse, binary, symmetric, zero-diagonal)
adj = spones(fullAdj(freeIdx, freeIdx));         % ensure 0/1
adj = adj | adj.';                               % force symmetry (undirected)
adj = adj - diag(diag(adj));                     % clear diagonal

%% --- Petri net from adjacency on free cells ---
[Pre, Post] = construct_PN(adj);
[nplaces, ntrans] = size(Pre);

fprintf(1,"\nThe Petri net has %i places and %i transitions.\n",nplaces,ntrans);

% Index mappings
invMap          = zeros(N, 1);                   % original -> reduced (0 if obstacle)
invMap(freeIdx) = 1:numel(freeIdx);
fwdMap          = freeIdx;                       % reduced -> original (linear index in B)


plot_animation = 0;%input("Do you want to plot the environment and the trajectories? (1 - yes, 0 - no)\n");

flag_ILP = 1;%input("Do you want to solve also the ILP formulation? This might take a while... (1 - yes, 0 - no)\n");

%%
for i = 1 : numel(N_robots)
    N_r = N_robots(i); % Set the number of robots for the current iteration
    for max_reg = 1:10 %1 : min(10,floor(numel(goalPts)/N_r))
        fprintf(1,'\nMaximum number of destinations per robot %i',max_reg);
        success = 0; % Initialize success counter for the current number of robots
        for exp=1:n_exp
            fprintf(1,"\n\n=======================================\n");
            fprintf(1,"Experiment number %i (%i robots)\n",exp,N_r);
            fprintf(1,"=======================================\n");

            %Generate Boolean goal
            el = randi([1 max_reg],1,N_r);
            N_p = sum(el);
            At = zeros(N_r,sum(el));
            bt = -ones(N_r,1);

            At(1,1:el(1)) = -ones(1,el(1));
            for j = 2:N_r
                At(j,sum(el(1:j-1))+1:sum(el(1:j))) = -ones(1,el(j));
            end

            rng('shuffle');                                  % different sample each run
            selIdxStart      = randperm(numel(initPts), N_r);
            selectedStart = initPts(selIdxStart);
            selIdxFin      = randperm(numel(goalPts), N_p);
            selectedFin = goalPts(selIdxFin);

            fprintf('Selected %i start and %i goal points.\n', N_r, N_p);

            [m0, idxStart] = initial_marking_multi_boolean(selectedStart, B, invMap, nplaces);
            matt = cell2mat(selectedFin);
            rows = matt(:,2) + 1;
            cols = matt(:,1) + 1;
            lin  = sub2ind([H, W], rows, cols);
            idxFin = double(invMap(lin));

            T.props = idxFin;

            if plot_animation
                plot_environment_new_SG_boolean(selectedStart, selectedFin, T.map2D);
            end

            [optVal, flag] = solve_LPs_collision_avoidance_boolean_CM(Post,Pre,At,bt,m0,T,flag_ILP,plot_animation);
            if flag, success = success + 1; end

            sim(exp).optim = optVal;
            sim(exp).flag  = flag;
            sim(exp).m0    = m0;
            sim(exp).selectedFin    = selectedFin;
            sim(exp).At    = At;
            sim(exp).T     = T;
            sim(exp).success = flag;
            sim(exp).max_reg = max_reg;
            sim(exp).formula = el;
        end
        fprintf(1,'\nSolved %i from %i experiments.',success,n_exp);
        save(sprintf('simulations_aisle_boolean_%drobots_%dN_p.mat', N_r, max_reg), 'sim', '-v7.3');
        fprintf(1,'\n');
        clear sim;
    end
end

dataDir = pwd;
files = dir(fullfile(dataDir, 'simulations_aisle_boolean_100robots_*N_p.mat'));
if isempty(files)
    warning('No se encontraron ficheros simulations_aisle_boolean_100robots_*N_p en %s', dataDir);
    T = table(); return;
end

%    expr = 'simulations_TAPF_(\d+)robots\.mat';
expr = 'simulations_aisle_boolean_100robots_(\d*)N_p';
results = [];

for k = 1:numel(files)
    fname = files(k).name;
    fpath = fullfile(files(k).folder, fname);

    tok = regexp(fname, expr, 'tokens', 'once');
    if isempty(tok), fprintf('Saltando %s (patrón)\n', fname); continue; end
    nRobots = str2double(tok{1});

    S = load(fpath);
    if ~isfield(S,'sim') || ~isstruct(S.sim)
        fprintf('Saltando %s (sin struct "sim")\n', fname); continue;
    end
    sim = S.sim;
    nExp = numel(sim);

    % Arrays por experimento
    rLP    = nan(nExp,1);   % runtimeLP por experimento (suma de partes válidas)
    rILP   = nan(nExp,1);   % runtimeILP por experimento (suma de partes válidas)
    svals  = nan(nExp,1);   % sbar por experimento
    costs  = nan(nExp,1);   % coste por experimento
    succ   = nan(nExp,1);   % éxito global

    for i = 1:nExp
        % Éxito global sim(i).success
        if isfield(sim(i),'success')
            succ(i) = toScalar(sim(i).success);
        else
            succ(i) = NaN;
        end

        hasOpt = isfield(sim(i),'optim') && isstruct(sim(i).optim);

        % ---------------- runtimeLP ----------------
        part_sum = 0; part_cnt = 0;
        if hasOpt && isfield(sim(i).optim,'LP1') && isstruct(sim(i).optim.LP1)
            lp1 = sim(i).optim.LP1;
            if isfield(lp1,'exitflag') && toScalar(lp1.exitflag)==1 ...
                    && isfield(lp1,'runtime') && ~isnan(succ(i)) && succ(i) == 1
                part_sum = part_sum + toScalar(lp1.runtime);
                part_cnt = part_cnt + 1;
            end
        end
        if isfield(sim(i).optim,'MILP') && isstruct(sim(i).optim.MILP) ...
                && isfield(sim(i).optim.MILP,'success') && toScalar(sim(i).optim.MILP.success)==1 ...
                && isfield(sim(i).optim.MILP,'runtime_all')  && isfield(sim(i).optim.MILP,'runtime') ...
                && isfield(lp1,'cellCapacity') && toScalar(lp1.cellCapacity) > 1 && ~isnan(succ(i)) ...
                && succ(i) == 1
            if (toScalar(sim(i).optim.MILP.runtime_all) == 0)
                part_sum = part_sum + toScalar(sim(i).optim.MILP.runtime);
            else
                part_sum = part_sum + toScalar(sim(i).optim.MILP.runtime_all);
            end
            part_cnt = part_cnt + 1;
        end
        if part_cnt > 0, rLP(i) = part_sum; end

        % ---------------- runtimeILP ----------------
        part_sum = 0; part_cnt = 0;
        if hasOpt && isfield(sim(i).optim,'ILP1') && isstruct(sim(i).optim.ILP1) ...
                && isfield(sim(i).optim,'LP1')  && isstruct(sim(i).optim.LP1) ...
                && ~isnan(succ(i)) && succ(i) == 1
            ilp1 = sim(i).optim.ILP1; lp1 = sim(i).optim.LP1;
            if isfield(lp1,'exitflag') && toScalar(lp1.exitflag)==1 ...
                    && isfield(ilp1,'runtime') && ~isnan(succ(i)) && succ(i) == 1
                part_sum = part_sum + toScalar(ilp1.runtime);
                part_cnt = part_cnt + 1;
            end
        end
        if isfield(sim(i).optim,'MILP') && isstruct(sim(i).optim.MILP) ...
                && isfield(sim(i).optim.MILP,'success') && toScalar(sim(i).optim.MILP.success)==1 ...
                && isfield(sim(i).optim,'ILP2') && isstruct(sim(i).optim.ILP2) ...
                && isfield(sim(i).optim.ILP2,'runtime_all') && isfield(ilp1,'cellCapacity') ...
                && toScalar(ilp1.cellCapacity) > 1 && ~isnan(succ(i)) && succ(i) == 1
            if toScalar(sim(i).optim.ILP2.runtime_all) > 0
                part_sum = part_sum + toScalar(sim(i).optim.ILP2.runtime_all);
            else
                part_sum = part_sum + toScalar(sim(i).optim.ILP2.runtime);
            end
            part_cnt = part_cnt + 1;
        end
        if part_cnt > 0, rILP(i) = part_sum; end

        % ---------------- sbar (por experimento) ----------------
        sval = NaN;
        if isfield(sim(i).optim,'MILP') && isstruct(sim(i).optim.MILP) ...
                && isfield(sim(i).optim.MILP,'success') && toScalar(sim(i).optim.MILP.success)==1 ...
                && isfield(sim(i).optim.MILP,'num_intermediate') && ~isnan(succ(i)) && succ(i) == 1
            sval = toScalar(sim(i).optim.MILP.num_intermediate);
        elseif hasOpt && isfield(sim(i).optim,'LP1') && isstruct(sim(i).optim.LP1) ...
                && isfield(sim(i).optim.LP1,'exitflag') && toScalar(sim(i).optim.LP1.exitflag)==1 ...
                && isfield(sim(i).optim.LP1,'cellCapacity') && ~isnan(succ(i)) && succ(i) == 1
            sval = toScalar(sim(i).optim.ILP1.cellCapacity);
        end
        svals(i) = sval;

        % ---------------- cost (por experimento) ----------------
        c = NaN;
        if hasOpt && isfield(sim(i).optim,'MILP') && isstruct(sim(i).optim.MILP) ...
                && isfield(sim(i).optim.MILP,'success') && toScalar(sim(i).optim.MILP.success)==1 ...
                && isfield(sim(i).optim.MILP,'cost') && ~isnan(succ(i)) && succ(i) == 1
            c = toScalar(sim(i).optim.MILP.cost);
        elseif hasOpt && isfield(sim(i).optim,'LP1') && isstruct(sim(i).optim.LP1) ...
                && isfield(sim(i).optim.LP1,'exitflag') && toScalar(sim(i).optim.LP1.exitflag)==1 ...
                && isfield(sim(i).optim.LP1,'cost') && isfield(sim(i).optim.LP1,'cellCapacity') ...
                && toScalar(sim(i).optim.LP1.cellCapacity)==1 && ~isnan(succ(i)) && succ(i) == 1
            c = toScalar(sim(i).optim.LP1.cost);
        end
        costs(i) = c;
    end

    % ---- Agregar fila por nº de robots ----
    row = struct();
    row.Np      = nRobots;
    row.runtimeLP   = mean(rLP,  'omitnan');
    row.runtimeILP  = mean(rILP, 'omitnan');
    row.cost        = mean(costs,'omitnan');

    valid_s = svals(~isnan(svals));
    if isempty(valid_s)
        row.sbar_mean = NaN; row.sbar_min = NaN; row.sbar_max = NaN;
    else
        row.sbar_mean = double(mean(valid_s));
    end

    if all(isnan(succ))
        row.SR_percent = NaN;
    else
        row.SR_percent = 100 * mean(succ == 1, 'omitnan');
    end

    results = [results; row]; %#ok<AGROW>
end

if isempty(results)
    T = table();
    warning('No se pudo construir la tabla (¿faltan datos válidos?)');
    return;
end

T = struct2table(results);
T = sortrows(T, 'Np');

disp(T);
outCSV = fullfile(dataDir, 'summary_TAPF.csv');
try
    writetable(T, outCSV);
    fprintf('Tabla guardada en: %s\n', outCSV);
catch ME
    warning('No se pudo guardar el CSV: %s', ME.message);
end


function x = toScalar(v)
if isempty(v), x = NaN; return; end
if isnumeric(v)
    if isscalar(v), x = double(v); else, x = double(v(1)); end
elseif islogical(v)
    x = double(v(1));
else
    try x = double(v(1)); catch, x = NaN; end
end
end
