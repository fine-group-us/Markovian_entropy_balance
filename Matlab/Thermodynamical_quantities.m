
clearvars; clc;
%% Description-------------------------------------------------------
%This script computes various thermodynamic quantities (Wm,Ix,Iu,I_after,I_before, W_FP,S_hk) from the raw data
%produced by the Fortran simulations of the two-state flashing ratchet. It can also compute phase diagrams by 
%finding intersections of Wm with its bounds.
%It is designed to perform sweeps over one of the main parameters: (Fext),(deltax), (deltat, Nstep, nummed).
%It produces an Excel file with the results of the sweep.
%It requires the functions S_hk_phase.m and cycle_average.m to be in the same folder.
%% Helper: safe log (log(0) -> 0)
safe_log = @(x) (x == 0).*0 + (x ~= 0).*log(max(x, realmin));

%% Fixed parameters
L      = 1;
beta   = 1;
F      = 0.1;     % height of the potential
offset = F/2;     % offset for the potential

% From simulation 
Nsim   = 1e5; % number of simulated particles
Nbins  = 600; % number of spatial bins
deltax = 0;   % measurement error  
a      = 1;   % position of potential minimum
M      = 10;  % size of the control memory
h      = 2*L/Nbins;  % spatial step size
deltat = 5e-4;  % time step size
Nstep  = 1000;  % steps per cycle
nummed = 100;   % number of cycles

nummed_transit = nummed/2; % number of transient cycles to discard
nummed_steady  = nummed - nummed_transit;  % number of steady-state cycles

deltatm = Nstep*deltat; % duration of one cycle
nt      = Nstep*nummed; % total number of time steps
nr      = Nstep*nummed_transit; % transient steps
ht      = deltatm/Nstep; % time step size for cycle-averaged quantities

%% Parameter sweep (choose ONE of these — keep others commented)
% --- Sweep in Fext ---
paramv1 = (0:0.05:0.9); % external force values     

% --- Sweep in deltax (examples) ---
% paramv1 = [0,0.025,0.05,0.075,0.10,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4]; 
% paramv1 = 0:0.1:0.9;                % for Fext = 0
% paramv1 = [0,0.05,0.1:0.05:0.9,0.99]; % for Fext > 0

% --- Sweep in deltat / Nstep / nummed (example template) ---
% paramv1 = [1e-5,2.5e-5,5e-5,7.5e-5,1e-4,2.5e-4,5e-4];         % deltat
% paramv2 = [100,100,100,100,100,1000,1000];                    % Nstep
% paramv3 = [1000,1000,1000,1000,1000,100,100];                 % nummed

%% Preallocate outputs
nv = numel(paramv1);
WFPv   = zeros(nv,1);
Wmv    = zeros(nv,1);
deltaIv= zeros(nv,1);
deltaFv= zeros(nv,1);
Icv    = zeros(nv,1);
Shkv   = zeros(nv,1);
Ixv    = zeros(nv,1);
Iuv    = zeros(nv,1);

%% String helpers (consistent with Fortran folder layout)
Nsim_str    = sprintf('%.1E', Nsim);
step_str    = sprintf('%.1E', Nstep);
F_str       = sprintf('%.1E', F);
deltax_str0 = sprintf('%.1E', deltax);
deltatm_str = sprintf('%.1E', deltatm);
nummed_str  = sprintf('%.1E', nummed);
Nbins_str   = num2str(Nbins);
a_str       = sprintf('%.1E', a);
M_str       = num2str(M);

baseDir = fullfile('FORTRAN','DATA_phases','Wm', ...
                   ['V0_', F_str], ...
                   ['Deltax_', deltax_str0]);

%% Verify that all required files exist
fprintf('Checking required files...\n');
for iv = 1:nv
    % Select parameter:----------------------------
    Fext = paramv1(iv);
    % deltax = paramv1(iv); % uncomment if sweeping deltax
    % deltat = paramv1(iv); % uncomment if sweeping deltat
    % Nstep  = paramv2(iv); % uncomment if sweeping Nstep
    % nummed = paramv3(iv); % uncomment if sweeping nummed
    % deltatm = Nstep*deltat; % update cycle duration if needed
    % nr      = Nstep*nummed_transit; % update transient steps if needed
    % nummed_steady  = nummed - nummed_transit;  % update steady cycles if needed
    % ht      = deltatm/Nstep; % update ht if needed

    Fext_str = sprintf('%.1E', Fext);
    %deltax_str = sprintf('%.1E', deltax);
    %deltatm_str = sprintf('%.1E', deltatm);
    %nummed_str  = sprintf('%.1E', nummed);
    %step_str    = sprintf('%.1E', Nstep);
    %----------------------------------------------
    file_tag = [Nsim_str, deltatm_str, nummed_str, step_str, M_str, Nbins_str, Fext_str, a_str, filesep];

    f_pc1 = fullfile(baseDir, file_tag, 'PC1.dat');
    if ~isfile(f_pc1)
        fprintf('MISSING: case %g — %s\n', paramv1(iv), f_pc1);
    end
end

%% Main loop over parameter values
for iv = 1:nv
    % Select parameter:----------------------------
    Fext = paramv1(iv);
    % deltax = paramv1(iv); % uncomment if sweeping deltax
    % deltat = paramv1(iv); % uncomment if sweeping deltat
    % Nstep  = paramv2(iv); % uncomment if sweeping Nstep
    % nummed = paramv3(iv); % uncomment if sweeping nummed
    % deltatm = Nstep*deltat; % update cycle duration if needed
    % nr      = Nstep*nummed_transit; % update transient steps if needed
    % nummed_steady  = nummed - nummed_transit;  % update steady cycles if needed
    % ht      = deltatm/Nstep; % update ht if needed
    %---------------------------------------------
    fprintf('Case %d / %d (param = %g)\n', iv, nv, paramv1(iv));

    % Build file names
    Fext_str = sprintf('%.1E', Fext);
    file_tag = [Nsim_str, deltatm_str, nummed_str, step_str, M_str, Nbins_str, Fext_str, a_str, filesep];

    f_pc1 = fullfile(baseDir, file_tag, 'PC1.dat');
    f_pc0 = fullfile(baseDir, file_tag, 'PC0.dat');
    f_pc  = fullfile(baseDir, file_tag, 'P_C.dat');
    f_hp  = fullfile(baseDir, file_tag, 'H_p.dat');

    % Load data from simulation
    % (Using load since files are plain whitespace-separated)
    Pc1      = load(f_pc1);
    Pc0      = load(f_pc0);
    Pc       = load(f_pc);
    H_prima  = load(f_hp);

    % Filter rows to avoid NaNs (legacy Fortran output behaviour)
    keepRows = Nstep*(nummed_steady - 2);
    Pc1 = Pc1(1:keepRows, :);
    Pc0 = Pc0(1:keepRows, :);

    % Cycle-averaged probabilities
    F1      = Nsim * Pc1;                        % frequency matrices
    F0      = Nsim * Pc0;
    F1ciclo = cycle_average(F1, Nstep);
    F0ciclo = cycle_average(F0, Nstep);

    P1ciclo = F1ciclo / (Nsim * (nummed_steady - 2)); % NOTE: change if Fortran is fixed
    P0ciclo = F0ciclo / (Nsim * (nummed_steady - 2));
    Pcciclo = mean(Pc);

    % Indices and forces per piecewise potential
    ia      = floor(a/h);          % index corresponding to x = a
    Fleft0  = -F / a;              % V(x,0) left force
    Fright0 =  F / (2*L - a);      % V(x,0) right force
    Fleft1  =  F / a;              % V(x,1) left force
    Fright1 = -F / (2*L - a);      % V(x,1) right force

    % --- Work contributed by Fext over one cycle (W_FP) ---
    WFPciclo = 0;
    for m = 1:Nstep
        % left region 1..ia
        for i = 1:ia
            WFPciclo = WFPciclo + Fext*ht*(Fleft1*P1ciclo(m,i) + Fleft0*P0ciclo(m,i));
        end
        % right region (ia+1)..Nbins
        for i = (ia+1):Nbins
            WFPciclo = WFPciclo + Fext*ht*(Fright1*P1ciclo(m,i) + Fright0*P0ciclo(m,i));
        end
    end
    WFPciclo = WFPciclo + Fext^2 * deltatm;

    % --- Protocol work per cycle (W_m), free energy-like Fciclo,
    %     and markovian information (I_after, I_before) ---
    Wmciclo  = 0;
    Fciclo   = 0;
    I_after  = 0;
    I_before = 0;

    % Left region (x < a)
    for i = 1:ia
        P1a = P1ciclo(Nstep,   i);  P1b = P1ciclo(Nstep-1, i);
        P0a = P0ciclo(Nstep,   i);  P0b = P0ciclo(Nstep-1, i);

        x_i    = i*h;
        V1left = F/a*(a - x_i) - offset;   % V1(x,0)
        V0left = F/a*(x_i)     - offset;   % V0(x,0)

        % Wm(x < a)
        Wmciclo = Wmciclo + V1left*(P1a - P1b) + V0left*(P0a - P0b);

        % F(x < a) “free energy” term
        Fciclo = Fciclo ...
               + P1a*(V1left + safe_log(P1a/h)) + P0a*(V0left + safe_log(P0a/h)) ...
               - P1b*(V1left + safe_log(P1b/h)) - P0b*(V0left + safe_log(P0b/h));

        % ΔI = I+ - I-
        Px_plus = P1a + P0a;

        if P1a == 0
            I_after = I_after + P0a * log( P0a / (Px_plus * (1 - Pcciclo)) );
        elseif P0a == 0
            I_after = I_after + P1a * log( P1a / (Px_plus * Pcciclo) );
        else
            I_after = I_after ...
                    + P1a * log( P1a / (Px_plus * Pcciclo) ) ...
                    + P0a * log( P0a / (Px_plus * (1 - Pcciclo)) );
        end

        if P1b == 0
            I_before = I_before + P0b * log( P0b / (Px_plus * (1 - Pcciclo)) );
        elseif P0b == 0
            I_before = I_before + P1b * log( P1b / (Px_plus * Pcciclo) );
        else
            I_before = I_before ...
                     + P1b * log( P1b / (Px_plus * Pcciclo) ) ...
                     + P0b * log( P0b / (Px_plus * (1 - Pcciclo)) );
        end
    end

    % Right region (x > a)
    for i = (ia+1):Nbins
        P1a = P1ciclo(Nstep,   i);  P1b = P1ciclo(Nstep-1, i);
        P0a = P0ciclo(Nstep,   i);  P0b = P0ciclo(Nstep-1, i);

        x_i     = i*h;
        V1right = F/(2*L - a)*(x_i - a)   - offset; % V1(x,0)
        V0right = F/(2*L - a)*(2*L - x_i) - offset; % V0(x,0)

        % Wm(x > a)
        Wmciclo = Wmciclo + V1right*(P1a - P1b) + V0right*(P0a - P0b);

        % F(x > a)
        Fciclo = Fciclo ...
               + P1a*(V1right + safe_log(P1a/h)) + P0a*(V0right + safe_log(P0a/h)) ...
               - P1b*(V1right + safe_log(P1b/h)) - P0b*(V0right + safe_log(P0b/h));

        % ΔI
        Px_plus = P1a + P0a;

        if P1a == 0
            I_after = I_after + P0a * log( P0a / (Px_plus * (1 - Pcciclo)) );
        elseif P0a == 0
            I_after = I_after + P1a * log( P1a / (Px_plus * Pcciclo) );
        else
            I_after = I_after ...
                    + P1a * log( P1a / (Px_plus * Pcciclo) ) ...
                    + P0a * log( P0a / (Px_plus * (1 - Pcciclo)) );
        end

        if P1b == 0
            I_before = I_before + P0b * log( P0b / (Px_plus * (1 - Pcciclo)) );
        elseif P0b == 0
            I_before = I_before + P1b * log( P1b / (Px_plus * Pcciclo) );
        else
            I_before = I_before ...
                     + P1b * log( P1b / (Px_plus * Pcciclo) ) ...
                     + P0b * log( P0b / (Px_plus * (1 - Pcciclo)) );
        end
    end

    % Channel-information term Sc and Ic
    Sc  = -Pcciclo*log(Pcciclo) - (1 - Pcciclo)*log(1 - Pcciclo);
    Ic  = I_after - Sc + H_prima;

    % Housekeeping entropy over a cycle
    xvec = linspace(-L, L, Nbins);
    Shk  = 0;
    for m = 1:Nstep
        rho1 = @(x) interp1(xvec, P1ciclo(m,:)/h, x, "linear");
        rho0 = @(x) interp1(xvec, P0ciclo(m,:)/h, x, "linear");
        Shk  = Shk + ht * S_hk_phase(Fext, a, offset, F, Nbins, rho1, rho0);
    end



    %% Compute Ix - Iu (steady-state analytic densities)
    V1 = @(x)  F*(abs(x) - L/2)/L - Fext*x;
    V0 = @(x) -F*(abs(x) - L/2)/L - Fext*x;

    % Unnormalized steady densities fs1, fs0 -- built via two-step integration+interp
    fs1_raw = @(x) exp(-beta*V1(x)) + ...
        (exp(-beta*V1(-L)) - exp(-beta*V1(L))) ./ ...
        integral(@(u) exp(-beta*(V1(L) - V1(u))), -L, L, 'RelTol', 1e-6) .* ...
        integral(@(y) exp(-beta*(V1(x) - V1(y))), -L, x, 'RelTol', 1e-10);

    fs0_raw = @(x) exp(-beta*V0(x)) + ...
        (exp(-beta*V0(-L)) - exp(-beta*V0(L))) ./ ...
        integral(@(u) exp(-beta*(V0(L) - V0(u))), -L, L, 'RelTol', 1e-6) .* ...
        integral(@(y) exp(-beta*(V0(x) - V0(y))), -L, x, 'RelTol', 1e-10);

    % Tabulate then interpolate for robust numeric integration on [-L, L]
    vals1 = zeros(1, Nbins);
    vals0 = zeros(1, Nbins);
    for ii = 1:Nbins
        z = xvec(ii);
        vals1(ii) = fs1_raw(z);
        vals0(ii) = fs0_raw(z);
    end
    fs1 = @(x) interp1(xvec, vals1, x, "linear");
    fs0 = @(x) interp1(xvec, vals0, x, "linear");

    % Normalize
    norm1 = integral(@(x) fs1(x), -L, L);
    norm0 = integral(@(x) fs0(x), -L, L);
    fs1 = @(x) fs1(x) ./ norm1;   % P(x | c=1)
    fs0 = @(x) fs0(x) ./ norm0;   % P(x | c=0)

    % I_x (before) and I_u (after) using last two time slices of the cycle
    Ixciclo = 0;
    Iuciclo = 0;
    for ii = 1:Nbins
        Ixciclo = Ixciclo - P1ciclo(Nstep-1, ii)*log(fs1(xvec(ii))) ...
                           - P0ciclo(Nstep-1, ii)*log(fs0(xvec(ii)));
        Iuciclo = Iuciclo - P1ciclo(Nstep,   ii)*log(fs1(xvec(ii))) ...
                           - P0ciclo(Nstep,   ii)*log(fs0(xvec(ii)));
    end

    %% Store results
    WFPv(iv)   = WFPciclo;
    Wmv(iv)    = Wmciclo;
    deltaIv(iv)= I_after - I_before;
    deltaFv(iv)= Fciclo;
    Icv(iv)    = Ic;
    Shkv(iv)   = Shk;
    Ixv(iv)    = Ixciclo;
    Iuv(iv)    = Iuciclo;
end

%% Save table to Excel (single file for the sweep)
T = table(Wmv(:), WFPv(:), deltaIv(:), Icv(:), Shkv(:), Ixv(:), Iuv(:), ...
    'VariableNames', {'Wmv','WFPv','deltaIv','Icv','Shkv','Ixv','Iuv'});

% One filename for the swept parameter (uses the *last* Fext if you swept Fext).
% If you want one file per case, move writetable() inside the loop and
% add an index or the parameter value to the filename.
filename = sprintf('Fext_%.2f_deltax_%.4f_deltam_%.4f_a_%d_V_%d.xlsx', Fext, deltax, deltatm, a, F);
writetable(T, filename);

%% Phase diagram--------------------------------------
    % Uncomment to find intersections of Wm with its bounds
    % paramv = paramv1; % choose the swept parameter  

    % W_f= @(x)interp1(paramv,Wmv+WFPv,x,'linear');
    % deltax0_Wm=fzero(Wmv_f,0.4)

    % Cota_basica_f= @(x)interp1(paramv,-deltaIv,x,'linear');
    % deltax0_Cota_basica=fzero(Cota_basica_f,0.4)

    % Cota_Shk_f= @(x)interp1(paramv,+Shkv-deltaIv,x,'linear');
    % deltax0_Cota_Shk=fzero(Cota_Shk_f,0.5)

    % %------------------
    % Cota_basica_c_f= @(x)interp1(paramv,-Icv,x,'linear');
    % deltax0_Cota_basica_c=fzero(Cota_basica_c_f,0.5)

    % Cota_Shk_c_f= @(x)interp1(paramv,Shkv-Icv,x,'linear');
    % deltax0_Cota_Shk_c=fzero(Cota_Shk_c_f,0.8)

%%------------------------------------------------------------------
%% Plots
paramv = paramv1;

% Figure 1
figure(1); clf; hold on; box on;
plot(paramv, Wmv + WFPv, 'k-',  'LineWidth', 2);
plot(paramv, Shkv - deltaIv, 'k--', 'LineWidth', 2);
plot(paramv, Shkv - Icv,     'b--', 'LineWidth', 2);
plot(paramv, -deltaIv, 'k:', 'LineWidth', 2);
plot(paramv, -Icv,     'b:', 'LineWidth', 2);
yline(0,'--');
ylabel('$\langle W \rangle$','fontsize',20,'interpreter','latex');
xlabel('$\Delta x$','fontsize',20,'interpreter','latex');
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
set(0,'DefaultAxesFontName','Times New Roman');
ax = gca;
ax.FontSize = 35;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 2;
ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';
set(0,'DefaultAxesFontName','Times New Roman');
shg;


