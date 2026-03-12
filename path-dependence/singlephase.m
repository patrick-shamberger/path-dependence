%%  ============================================================================
% Description:
%{
    TITLE: DEFINING SINGLE PHASE PROPERTIES (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function takes the raw ISHC data and specified boundary/α points
         and finds the single phase material properties of interest 
         (compliance(σ) & α(σ)). It also stores the full reference point of
         each single phase region (T,σ,ε).]

    INPUTS:
        - InitialInputs structure for ISHC data.
        - MatProps structure for phase boundaries and reference points.
        - Plot (boolean logic) to declare if figures associated with this
          fuction should be shown (1 - YES, 0 - NO).

    OUTPUTS:
        - InitialInputs structure.
        - MatProps structure after updating with strain values for each of
          the single phase reference points and the single phase material
          properties.

    NOTES:
%}
%  =============================================================================
function [InitialInputs, MatProps] = singlephase(InitialInputs, MatProps, Plot)

    % === Extract Inputs ===
    ISHC = InitialInputs.ISHCdata;
    St = InitialInputs.Stresses(:);
    CTE_i = InitialInputs.CTEIndices(:);

    % === Reference points ===
    RefA = MatProps.Reference.Austenite;
    RefM = MatProps.Reference.Martensite; 

    % === Storage for CTE fits ===
    n = numel(ISHC);
    CTEA = cell(n,1); % Austenite
    CTEM = cell(n,1); % Martensite

    %% CTE fits
    for i = 1:n
        eps = ISHC{i}(:,3); % Strain [%]
        T = ISHC{i}(:,2) + 273; % Temperature [K]

        % === Split indices ===
        [~, ic_end] = min(T);
        ic = 1:ic_end; % Cooling indices
        ih = (ic_end+1):numel(T); % Heating indices

        % === CTE bounds ===
        St_i = St(i);
        [A_bounds, M_bounds] = ctebounds(CTE_i, St_i);

        % === Austenite CTE from cooling ===
        A_k = (T(ic) >= A_bounds(1)) & (T(ic) <= A_bounds(2));
        A_i = ic(A_k);
        CTEA{i} = polyfit(T(A_i), eps(A_i), 1);

        % === Martensite CTE from heating ===
        M_k = (T(ih) >= M_bounds(1)) & (T(ih) <= M_bounds(2));
        M_i = ih(M_k);
        CTEM{i} = polyfit(T(M_i), eps(M_i), 1);
    end

    % === Collect CTE(σ) values === 
    CTEA_all = NaN(n,1);
    CTEM_all = NaN(n,1);
    for i = 1:n
        CTEA_all(i) = CTEA{i}(1);
        CTEM_all(i) = CTEM{i}(1);
    end

    %% Fit α(σ) relationships (A & M)
    CTEA_p = polyfit(St, CTEA_all, 1); % Linear Fit
    CTEM_p = polyfit(St, CTEM_all, 2); % Quadratic Fit

    % === Store α(σ) ===
    MatProps.CTE.Austenite = struct('Type',...
        'Linear', 'Coefficients', CTEA_p);
    MatProps.CTE.Martensite = struct('Type',...
        'Quadratic', 'Coefficients', CTEM_p);

    %% Plot ε–T & α lines
    if Plot == 1
        figure('Name','Strain vs Temperature with CTE Lines');
        hold on;
        xlabel('Temperature (K)'); 
        ylabel('Strain (%)');
        cmap = lines(n);
        for i = 1:n
            eps = ISHC{i}(:,3); % [%]
            T = ISHC{i}(:,2) + 273; % [K]
            plot(T, eps, '-', 'Color', cmap(i,:));
            x_range = linspace(min(T), max(T), 200);
            plot(x_range, polyval(CTEA{i}, x_range), '-r');
            plot(x_range, polyval(CTEM{i}, x_range), '-r');
        end
        hold off;

        stress_fit = linspace(min(St), max(St), 200);

        % === Austenite ===
        figure('Name','CTE vs Stress - Austenite');
        plot(St, CTEA_all, 'o', 'LineWidth', 2); 
        hold on;
        plot(stress_fit, polyval(CTEA_p, stress_fit), '-', 'LineWidth', 1);
        xlabel('Stress (MPa)'); 
        ylabel('CTE (Austenite) (%/K)'); 
        hold off;

        % === Martensite ===
        figure('Name','CTE vs Stress - Martensite');
        plot(St, CTEM_all, 's', 'LineWidth', 2); 
        hold on;
        plot(stress_fit, polyval(CTEM_p, stress_fit), '-', 'LineWidth', 1);
        xlabel('Stress (MPa)'); 
        ylabel('CTE (Martensite) (%/K)'); 
        hold off;
    end

    %% Austenite Compliance
    RA_T = RefA.Temperature;
    AStrains = NaN(n,1);

    for i = 1:n
        eps = ISHC{i}(:,3); % [%]
        T = ISHC{i}(:,2) + 273; % [K]

        % === Cooling section ===
        [~, ic_end] = min(T);
        ic = 1:ic_end;
        T_c = T(ic);
        eps_c = eps(ic);

        % === Austenite strains at reference temperature ===
        dist = abs(T_c - RA_T);
        mindist = min(dist);
        j = find(dist == mindist, 1, 'last');
        AStrains(i) = eps_c(j);
    end

    % === Store 5 MPa reference strain ===
    MatProps.Reference.Austenite.Strain = AStrains(1);

    % === Elastic compliance (Austenite) ===
    Astrain_p = polyfit(St, AStrains, 1);
    Acomp_p = Astrain_p(1);
    MatProps.Compliance.Austenite = struct('Type','Linear', ...
                                           'StressStrain',Astrain_p, ...
                                           'Compliance',Acomp_p);

    stress_fit = linspace(min(St), max(St), 200);
    if Plot == 1
        figure('Name','Austenite Strain vs Stress');
        hold on;
        plot(St, AStrains, 'o-', 'LineWidth', 2);
        plot(stress_fit, polyval(Astrain_p, stress_fit), 'k--');
        xlabel('Stress (MPa)');
        ylabel(sprintf('Strain at %.2f K (Austenite)', RA_T));
        hold off;
    end

    %% === Martensite "Compliance" ===
    RM_T = RefM.Temperature;
    MStrains = NaN(n,1);

    for i = 1:n
        eps = ISHC{i}(:,3); % [%]
        T = ISHC{i}(:,2) + 273; % [K]

        % === Heating section ===
        [~, ic_end] = min(T);
        ih = (ic_end+1):numel(T);
        T_h = T(ih);
        eps_h = eps(ih);

        % === Martensite strains at reference temperature ===
        dist = abs(T_h - RM_T);
        mindist = min(dist);
        j = find(dist == mindist, 1, 'last');
        MStrains(i) = eps_h(j);
    end

    % === Store 5 MPa reference strain ===
    MatProps.Reference.Martensite.Strain = MStrains(1);

    % === Cubic compliance (Martensite) ===
    Mstrain_p = polyfit(St, MStrains, 3);
    p = Mstrain_p;
    dStrain_dSt = [3*p(1), 2*p(2), p(3)]; % Derivative coeffs

    MatProps.Compliance.Martensite = struct('Type','Cubic', ...
                                            'StressStrain',Mstrain_p, ...
                                            'Compliance',dStrain_dSt);

    if Plot == 1
        figure('Name','Martensite Strain vs Stress');
        hold on;
        plot(St, MStrains, 's-', 'LineWidth', 2);
        plot(stress_fit, polyval(Mstrain_p, stress_fit), 'k--');
        xlabel('Stress (MPa)');
        ylabel(sprintf('Strain at %.2f K (Martensite)', RM_T));
        hold off;
    end
end

%% Helper functions
function [A_bounds, M_bounds] = ctebounds(CTE_i, St_i)
    k = find([CTE_i.Stress] == St_i, 1, 'first');
    A = CTE_i(k).Austenite; 
    M = CTE_i(k).Martensite;  
    A_bounds = [A.LowerT, A.UpperT];
    M_bounds = [M.LowerT, M.UpperT];
end
