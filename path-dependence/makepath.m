%% =============================================================================
% Description:
%{
    TITLE: MAKING FULL PATHS (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function creates a Paths structure containing σ-T points for
         all paths based on the specified "corner points" in
         "ARUN_predict.m". Full σ-T paths are constructed by making linear
         segments between the corner points].

    INPUTS:
        - all_paths: Cell array. Each cell is [N x 2] endpoints (T,σ).
        - n_leg: Number of points per leg.

    OUTPUTS:
        - Paths structure containing Leg_# tables and FullTable table

    NOTES:
        - Straight line legs
        - Value of n_leg is preserved so that the total number of points in
          a path is equal to n_leg times the number of legs in the path.
%}
%  =============================================================================
function Paths = makepath(all_paths, n_leg)

    % === Convert all_paths to a structure ===
    legpoints = struct();
    for p = 1:length(all_paths)
        legpoints.(sprintf('Path_%d', p)) = all_paths{p};
    end

    % === Initialize output structure ===
    Paths = struct();

    % === Loop over each path ===
    path_names = fieldnames(legpoints);
    for i = 1:length(path_names)
        path_name = path_names{i};
        endpoints = legpoints.(path_name);

        % === Number of legs ===
        n_l = size(endpoints, 1) - 1;

        % === Initialize full path storage ===
        allpoints = NaN(n_l * n_leg, 2); % Full table size
        k = 1;

        %% BUILD FULL PATHS
        for j = 1:n_l
            T_s = endpoints(j, 1); % Starting temperature
            St_s = endpoints(j, 2); % Starting stress
            T_e = endpoints(j+1, 1); % Final temperature
            St_e = endpoints(j+1, 2); % Final stress

            % === Ensure total n points = n_leg*(# legs) ===
            n_points = n_leg;
            if j > 1
                n_points = n_leg + 1;
            end
            T = linspace(T_s,  T_e,  n_points)';
            St = linspace(St_s, St_e, n_points)';
            if j > 1
                T(1) = [];
                St(1) = [];
            end

            % === Store legs ===
            Paths.(path_name).(sprintf('Leg_%d', j)) = table(T, St, ...
                'VariableNames', {'Temperature', 'Stress'} );

            % === Store Full Path ===
            m = numel(T);
            allpoints(k:k+m-1, :) = [T, St];
            k = k + m; % Build on previous legs
        end

        % === Store FullTable ===
        Paths.(path_name).FullTable = array2table( ...
            allpoints, 'VariableNames', {'Temperature', 'Stress'} );
    end
end
