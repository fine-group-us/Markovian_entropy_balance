%% COMPUTE I AND Iu AS FUNCTIONS OF Δt_m (f = 0, Δx = 0)
clearvars; clc;

% Parameters
L     = 1;
F     = 10;
beta  = 1;
Nbins = 600;

xvec = linspace(-L, L, Nbins);
xL   = -L/2;
xR   =  L/2;

deltatmv = linspace(5e-3, 1, 5);
Ixv      = zeros(1, numel(deltatmv));
Iuv      = zeros(1, numel(deltatmv));
Wv       = zeros(1, numel(deltatmv));

% Safe log to avoid log(0)
eps_log  = 1e-30;
safe_log = @(z) log(max(z, eps_log));

% Numerical tolerances
RelTol = 1e-8;  AbsTol = 1e-12;

% Loop over Δt_m values
for k = 1:numel(deltatmv)
    deltatm = deltatmv(k);

    % Initial condition for P_x (will be overwritten by fixed-point iteration)
    P_x = @(x) 1/(2*L);
    tol = inf;

    % Fixed-point iteration to obtain stationary P_x
    maxIter = 200;                    % safety cap
    iter    = 0;
    P_x_v   = zeros(1, Nbins);        % preallocate

    while tol > 1e-2 && iter < maxIter
        iter = iter + 1;

        for in = 1:Nbins   % integral over x_p on the grid
            x    = xvec(in);
            k1   = @(xp) K_C_1(x, xp, deltatm, F);
            k0   = @(xp) K_C_0(x, xp, deltatm, F);

            % Region split: c=1 for x in [xL, xR], else c=0
            P_x_v(in) = ...
                integral(@(xp) k1(xp) .* P_x(xp),  xL, xR, 'RelTol',RelTol,'AbsTol',AbsTol) + ...
                integral(@(xp) k0(xp) .* P_x(xp), -L,  xL, 'RelTol',RelTol,'AbsTol',AbsTol) + ...
                integral(@(xp) k0(xp) .* P_x(xp),  xR,  L, 'RelTol',RelTol,'AbsTol',AbsTol);
        end

        % Renormalize to ensure ∫ P_x = 1
        Z = trapz(xvec, P_x_v);
        if Z > 0, P_x_v = P_x_v / Z; end

        % Convergence check against current P_x evaluated on xvec
        tol  = max(abs(P_x_v - P_x(xvec)));
        P_x  = @(x) interp1(xvec, P_x_v, x, "spline", "extrap");
    end

    %% Compute I_x
    aux_v = zeros(1, Nbins);
    for in = 1:Nbins
        x  = xvec(in);
        k1 = @(xp) K_C_1(x, xp, deltatm, F);
        k0 = @(xp) K_C_0(x, xp, deltatm, F);

        t1 = integral(@(xp) P_x(xp).*k1(xp).*safe_log(k1(xp)), xL, xR, 'RelTol',RelTol,'AbsTol',AbsTol);
        t0 = integral(@(xp) P_x(xp).*k0(xp).*safe_log(k0(xp)), -L, xL, 'RelTol',RelTol,'AbsTol',AbsTol) ...
           + integral(@(xp) P_x(xp).*k0(xp).*safe_log(k0(xp)),  xR,  L, 'RelTol',RelTol,'AbsTol',AbsTol);

        aux_v(in) = t1 + t0;
    end
    aux = @(x) interp1(xvec, aux_v, x, "spline", "extrap");
    Ix  = - integral(@(x) aux(x), -L, L, 'RelTol',RelTol,'AbsTol',AbsTol);

    %% Compute I_u
    % Use the “reversed-argument” kernels inside the log
    auxu_v = zeros(1, Nbins);
    for in = 1:Nbins
        x    = xvec(in);
        k1   = @(xp) K_C_1(x,  xp, deltatm, F);
        k0   = @(xp) K_C_0(x,  xp, deltatm, F);
        k1_r = @(xp) K_C_1(xp, x,  deltatm, F);
        k0_r = @(xp) K_C_0(xp, x,  deltatm, F);

        u1 = integral(@(xp) P_x(xp).*k1(xp).*safe_log(k1_r(xp)), xL, xR, 'RelTol',RelTol,'AbsTol',AbsTol);
        u0 = integral(@(xp) P_x(xp).*k0(xp).*safe_log(k0_r(xp)), -L, xL, 'RelTol',RelTol,'AbsTol',AbsTol) ...
           + integral(@(xp) P_x(xp).*k0(xp).*safe_log(k0_r(xp)),  xR,  L, 'RelTol',RelTol,'AbsTol',AbsTol);

        auxu_v(in) = u1 + u0;
    end
    auxu = @(x) interp1(xvec, auxu_v, x, "spline", "extrap");
    Iu   = - integral(@(x) auxu(x), -L, L, 'RelTol',RelTol,'AbsTol',AbsTol);

    % Store
    Ixv(k) = Ix;
    Iuv(k) = Iu;

    %% Joint distributions at t_k^+ (after measurement)
    Indicatriz_1 = @(x) (x >= xL & x <= xR);
    P_x_1_a = @(x) P_x(x) .* Indicatriz_1(x);
    P_x_0_a = @(x) P_x(x) .* ~Indicatriz_1(x);

    %% Joint distributions at t_k^- (before measurement)
    P_x_1_b_v = zeros(1, Nbins);
    P_x_0_b_v = zeros(1, Nbins);
    for iv = 1:Nbins
        xv = xvec(iv);
        P_x_1_b_v(iv) = integral(@(xp) K_C_1(xv, xp, deltatm, F).*P_x_1_a(xp), -L, L, 'RelTol',RelTol,'AbsTol',AbsTol);
        P_x_0_b_v(iv) = integral(@(xp) K_C_0(xv, xp, deltatm, F).*P_x_0_a(xp), -L, L, 'RelTol',RelTol,'AbsTol',AbsTol);
    end
    P_x_1_b = @(x) interp1(xvec, P_x_1_b_v, x, "spline", "extrap");
    P_x_0_b = @(x) interp1(xvec, P_x_0_b_v, x, "spline", "extrap");

    %% Measurement work W
    Vx1 = @(x)  F*(abs(x) - L/2)/L;
    Vx0 = @(x) -F*(abs(x) - L/2)/L;

    W = integral(@(x) Vx1(x).*(P_x_1_a(x) - P_x_1_b(x)) + ...
                       Vx0(x).*(P_x_0_a(x) - P_x_0_b(x)), -L, L, 'RelTol',RelTol,'AbsTol',AbsTol);
    Wv(k) = W;
end

%% Plot
figure(1); clf; hold on; box on;
plot(deltatmv, Wv,          'k-', 'LineWidth', 2.5);
plot(deltatmv, -(Ixv - Iuv),'ro', 'MarkerSize', 10, 'LineWidth', 2);

ylabel('$\beta \langle W \rangle,\ \overline{\mathcal{I}_u-\mathcal{I}_{\vec{x}}}$', ...
       'fontsize', 30, 'interpreter','latex');
xlabel('$\Delta t_m$', 'fontsize', 35, 'interpreter','latex');
set(gca,'FontSize',35,'TickLabelInterpreter','latex');
set(0,'DefaultAxesFontName','Times New Roman');
ax = gca; ax.LineWidth = 2; ax.Box = 'on'; ax.XColor = 'k'; ax.YColor = 'k';
shg;
