%% =============================================================================
% Description:
%{
    TITLE: ACCOUNTING FOR PATH-DEPENDENCE (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)    

    SUMMARY:
        [Applies chosen method for accounting for path-depenence and calculates
         final strain value for each path.]

    INPUTS:
        - Paths structure containing all paths.
        - MatProps structure containing relevant material properties.
        - Plotting structure that contains grid size.
        - Path-dependence method (1 - Approach 1, 2 - Approach 2).
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - Adjusted final strain value located within the
          Paths structure for each path.

    NOTES:
        - Approach 1 = Experimental ε_final parametrization on σ_avg
        - Approach 2 = Self-similarity
%}
%  =============================================================================
function Paths = pathdep(Paths, MatProps, Plotting, Type, Plot)
    % === Extract Required Inputs ===
    Meps = MatProps.Compliance.Martensite.StressStrain; % ISHC Strain
    StG = Plotting.Grid.StressRange;
    St_min = min(StG);
    St_max = max(StG);

    % === Import Reorientation Experiment ===
    Sys = Plotting.Materialsys;
    folder_cur = fileparts(mfilename('fullpath'));
    folder_data = fullfile(folder_cur, 'Experimental Data', 'Detwinning', Sys);
    file = [Sys '_Reversible_Stress.xlsx'];
    fpath = fullfile(folder_data, file);

    if isfile(fpath)

        % === Loading Leg of Experiment ===
        loadd = readmatrix(fpath, 'Sheet', 3);
        eps_exp = loadd(:,3); % [%]
        St_exp = loadd(:,4); % [MPa]

        % === Fit General Logistic Function ===
        % ε(σ) = d + a / [1 + exp(-b*(x-c))]^n
        Sig = @(p, x) p(4) + p(1) ./ (1 + exp(-p(2) * (x - p(3)))).^p(5);
        p0 = [5, 0.03, mean(St_exp), min(eps_exp), 1]; % Initial Guess
        p = nlinfit(St_exp, eps_exp, Sig, p0);

        % === Plot Fit ===
        px = linspace(min(St_exp), max(St_exp), 100);
        py = Sig(p, px);
        if Plot == 1
            figure('Name', 'Sigmoid Fit');
            hold on;
            plot(St_exp, eps_exp,'bo');
            plot(px, py, 'r-');
            xlabel('Stress (MPa)');
            ylabel('Strain (%)');
            hold off;
        end

        % === Align Strain Values ===
        Sig_e0 = Sig(p, St_min);
        ISHC_e0 = polyval(Meps, St_min);
        d_e0 = ISHC_e0 - Sig_e0;
        if abs(d_e0) ~= 0
            fprintf('ISHC and Reorientation strains differ by %.4f %%\n', d_e0);
        end
        ASig = @(x) Sig(p, x) + d_e0;

        %% Approach 1: Experimental parameterization at σ_max
        if Type == 1
            % === Max and Min Strains ===
            x1 = 5; % [MPa]
            x2 = 200; % [MPa] Experimental σ_max for type III experiment
            y1 = ASig(x2); % [%]
            y2 = polyval(Meps, x2); % [%]

            % === Parameterization ===
            %{ 
               This uses a logarithmic function to describe the type III 
               experimental behavior at 200 MPa, but it should be changed
               to use a function that best describes the experimental
               behavior at the given σ_max pf interest.
            %}
            
            % y = a*log(x) + b
            a = (y2 - y1) / (log(x2) - log(x1));
            b = y1 - a*log(x1);

            % === Assign ε_final(σ_avg) ===
            pns = fieldnames(Paths);
            pns = pns(contains(pns, 'Path_'));

            for i = 1:length(pns)
                pn = pns{i};
                St_avg = Paths.(pn).Avg_Stress;
                eps_f = a*log(St_avg) + b;
                Paths.(pn).FinalStrain = eps_f;
            end
        elseif Type == 2
            %% Approach 2: Phenomenological self-similarity
            e_min = ASig(St_max); % Type II strain at reversible stress
            e_max = polyval(Meps, St_max); % Type I strain at reversible stress
    
            % === Sigmoid Strain Bounds ===
            eS_min = ASig(St_min);
            eS_max = ASig(St_max);
    
            %{
                This uses a logarithmic function to describe the type-III 
                experimental behavior at σ_rev, but it should be changed to use a 
                function that best describes the experimental behavior at σ_rev for
                a given material of interest.
            %}
    
            % y = a*log(x) + b
            a = (e_max - e_min) / (log(St_max) - log(St_min));
            b = e_min - a*log(St_min);

            pns = fieldnames(Paths);
            pns = pns(contains(pns, 'Path_'));
            for i = 1:length(pns)
                pn = pns{i};

                % === Max Stress of Path ===
                PSt = Paths.(pn).FullTable.Stress;
                PSt_max = max(PSt);

                % === Assign ε_final(σ_avg) ===
                St_avg = Paths.(pn).Avg_Stress;

                if St_avg == St_min
                    eps_f = ASig(PSt_max);
                else
                    ys = polyval(Meps, St_avg); % Start strain mapped from ISHC
                    yf = a*log(St_avg) + b;

                    % === Reparameterize Stress Interval ===
                    if St_avg == St_max
                        sigma = St_max;
                    else
                        sigma = St_min + (PSt_max - St_avg) / ... 
                            (St_max - St_avg) * (St_max - St_min);
                    end

                    % === Evaluate on Mapped Stress ===
                    sig_map = ASig(sigma);

                    % === Normalize and Scale Strain ===
                    sig_n = (sig_map - eS_min) / (eS_max - eS_min);
                    eps_f = ys + (yf - ys) * sig_n;
                end
                Paths.(pn).FinalStrain = eps_f;
            end
        else
            error('Type must be 1 or 2.');
        end
    else
        %% Not enough information to establish a path-dependence
        fprintf(['Path dependence can not be established for %s because ' ...
                 'the detwinning experiment file was not found\n'], Sys);
        
        pns = fieldnames(Paths);
        pns = pns(contains(pns, 'Path_'));

        for i = 1:length(pns)
            pn = pns{i};
            Paths.(pn).FinalStrain = NaN; % NaN final strain
        end
    end
end