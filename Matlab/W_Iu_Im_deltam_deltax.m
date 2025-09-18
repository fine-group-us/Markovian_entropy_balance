%% COMPUTE I AND Iu AS A FUNCTION OF DELTATM (f=0, deltax=0)

clear all;
% Safe log: use a tiny epsilon to avoid log(0). You can also use realmin.
eps_log  = 1e-30;
safe_log = @(z) log(max(z, eps_log));

% Parameters
L     = 1;
F     = 5;
beta  = 1;
Nbins = 600;

xvec = linspace(-L, L, Nbins);
xL   = -L/2;
xR   =  L/2;

% Sweep of measurement interval Δt_m
deltatmv = linspace(5e-3, 0.5, 20);

Ixv     = zeros(1, numel(deltatmv));
Iuv     = zeros(1, numel(deltatmv));
Wv      = zeros(1, numel(deltatmv));
Imasv   = zeros(1, numel(deltatmv));
Imenosv = zeros(1, numel(deltatmv));

% Measurement kernel (protocol): Θ(x|c)
deltax = 0.8;
Theta1 = @(x) interp1(linspace(-L, L, 1000), Protocol_Wm(deltax, 1, F, F/2), x, "linear");
Theta0 = @(x) 1 - Theta1(x);

% Iterate over Δt_m values
for k = 1:numel(deltatmv)
    fprintf('Δt_m case %d / %d\n', k, numel(deltatmv));
    deltatm = deltatmv(k);

    % Initial condition for P_x (will be overwritten by fixed-point iteration)
    P_x  = @(x) 1/(2*L);
    tol  = inf;

    % Fixed-point iteration to obtain stationary P_x
    maxIter = 200;                % safety cap
    iter    = 0;
    P_x_v   = zeros(1, Nbins);    % preallocate

    while tol > 1e-2 && iter < maxIter
        iter = iter + 1;

        for in = 1:Nbins   % loop over x on the evaluation grid
            x      = xvec(in);
            k_c_1  = @(xp) K_C_1(x, xp, deltatm, F);
            k_c_0  = @(xp) K_C_0(x, xp, deltatm, F);

            % P_x^{new}(x) = ∫ [k1 * P_x * Θ1 + k0 * P_x * Θ0] dxp
            P_x_v(in) = integral(@(xp) k_c_1(xp).*P_x(xp).*Theta1(xp), -L, L, ...
                                  'RelTol',1e-8,'AbsTol',1e-12) ...
                       + integral(@(xp) k_c_0(xp).*P_x(xp).*Theta0(xp), -L, L, ...
                                  'RelTol',1e-8,'AbsTol',1e-12);
        end

        % Optional: renormalize to ensure ∫ P_x = 1
        Z = trapz(xvec, P_x_v);
        if Z > 0, P_x_v = P_x_v / Z; end

        % Convergence check against current P_x evaluated on xvec
        Px_old = P_x(xvec);
        tol    = max(abs(P_x_v - Px_old));
        fprintf('  iter %3d  tol = %.3e\n', iter, tol);

        % Update P_x from tabulated P_x_v
        P_x = @(x) interp1(xvec, P_x_v, x, "spline", "extrap");
    end

    %% Compute Ix
    aux_v = zeros(1, Nbins);
    for in = 1:Nbins
        x     = xvec(in);
        k_c_1 = @(xp) K_C_1(x, xp, deltatm, F);
        k_c_0 = @(xp) K_C_0(x, xp, deltatm, F);

        mix = @(xp) k_c_1(xp).*Theta1(xp) + k_c_0(xp).*Theta0(xp);

        term1 = integral(@(xp) P_x(xp).*k_c_1(xp).*Theta1(xp).*safe_log(mix(xp)), -L, L, ...
                         'RelTol',1e-8,'AbsTol',1e-12);
        term0 = integral(@(xp) P_x(xp).*k_c_0(xp).*Theta0(xp).*safe_log(mix(xp)), -L, L, ...
                         'RelTol',1e-8,'AbsTol',1e-12);

        aux_v(in) = term1 + term0;
    end
    aux = @(x) interp1(xvec, aux_v, x, "spline", "extrap");
    Ix  = - integral(@(x) aux(x), -L, L, 'RelTol',1e-8,'AbsTol',1e-12);

    %% Compute Iu
    % (be careful: we need log of the kernel in the "reverse" argument order)
    auxu_v = zeros(1, Nbins);
    for in = 1:Nbins
        x        = xvec(in);
        k_c_1    = @(xp) K_C_1(x, xp, deltatm, F);
        k_c_0    = @(xp) K_C_0(x, xp, deltatm, F);
        k_c_1_r  = @(xp) K_C_1(xp, x, deltatm, F);  % reversed args for log k
        k_c_0_r  = @(xp) K_C_0(xp, x, deltatm, F);

        term1 = integral(@(xp) P_x(xp).*k_c_1(xp).*Theta1(xp).*safe_log(k_c_1_r(xp)), -L, L, ...
                         'RelTol',1e-8,'AbsTol',1e-12);
        term0 = integral(@(xp) P_x(xp).*k_c_0(xp).*Theta0(xp).*safe_log(k_c_0_r(xp)), -L, L, ...
                         'RelTol',1e-8,'AbsTol',1e-12);

        auxu_v(in) = term1 + term0;
    end
    auxu = @(x) interp1(xvec, auxu_v, x, "spline", "extrap");
    Iu   = - integral(@(x) auxu(x), -L, L, 'RelTol',1e-8,'AbsTol',1e-12);

    % Store results
    Ixv(k) = Ix;
    Iuv(k) = Iu;

    %% Joint distributions at t_k^+ (after measurement)
    P_x_1_a = @(x) P_x(x).*Theta1(x);
    P_x_0_a = @(x) P_x(x).*Theta0(x);

    %% Joint distributions at t_k^- (before measurement)
    P_x_1_b_v = zeros(1, Nbins);
    P_x_0_b_v = zeros(1, Nbins);
    for iv = 1:Nbins
        xv = xvec(iv);
        P_x_1_b_v(iv) = integral(@(xp) K_C_1(xv, xp, deltatm, F).*P_x_1_a(xp), -L, L, ...
                                 'RelTol',1e-8,'AbsTol',1e-12);
        P_x_0_b_v(iv) = integral(@(xp) K_C_0(xv, xp, deltatm, F).*P_x_0_a(xp), -L, L, ...
                                 'RelTol',1e-8,'AbsTol',1e-12);
    end
    P_x_1_b = @(x) interp1(xvec, P_x_1_b_v, x, "spline", "extrap");
    P_x_0_b = @(x) interp1(xvec, P_x_0_b_v, x, "spline", "extrap");

    %% Channel probabilities p(c)
    pc1 = integral(@(x) P_x_1_a(x), -L, L, 'RelTol',1e-10,'AbsTol',1e-14);
    pc0 = integral(@(x) P_x_0_a(x), -L, L, 'RelTol',1e-10,'AbsTol',1e-14);

    %% Measurement work W
    Vx1 = @(x)  F*(abs(x) - L/2)/L;
    Vx0 = @(x) -F*(abs(x) - L/2)/L;

    W = integral(@(x) Vx1(x).*(P_x_1_a(x) - P_x_1_b(x)) + ...
                       Vx0(x).*(P_x_0_a(x) - P_x_0_b(x)), -L, L, ...
                 'RelTol',1e-8,'AbsTol',1e-12);
    Wv(k) = W;

    %% ΔI_m = I^+ - I^-
    Imasv(k) = integral(@(x) P_x_1_a(x).*safe_log(Theta1(x)/pc1) + ...
                              P_x_0_a(x).*safe_log(Theta0(x)/pc0), -L, L, ...
                        'RelTol',1e-8,'AbsTol',1e-12);

    Imenosv(k) = integral(@(x) P_x_1_b(x).*safe_log(P_x_1_b(x)./(P_x(x)*pc1)) + ...
                               P_x_0_b(x).*safe_log(P_x_0_b(x)./(P_x(x)*pc0)), -L, L, ...
                         'RelTol',1e-8,'AbsTol',1e-12);
end

%% PLOTS
figure(10); clf; hold on; box on;
plot(deltatmv, Wv, 'k-', 'LineWidth', 2.5);
plot(deltatmv, -(Ixv - Iuv), 'ro', 'MarkerSize', 10, 'LineWidth', 2.5);
plot(deltatmv, -(Imasv - Imenosv), 'k--', 'LineWidth', 2.5);

ylabel('$\beta \langle W \rangle$ and bounds', 'fontsize', 30, 'interpreter', 'latex');
xlabel('$\Delta t_m$', 'fontsize', 35, 'interpreter', 'latex');
set(gca, 'FontSize', 35, 'TickLabelInterpreter','latex');
set(0, 'DefaultAxesFontName', 'Times New Roman');
ax = gca; ax.LineWidth = 2; ax.Box = 'on'; ax.XColor = 'k'; ax.YColor = 'k';
shg;

   