%% Wm_I_Iu_vs_Deltax  (Delta t_m -> infinity)
% Computes W (measurement work), I_x, I_u, and Delta I_m as functions of Δx
hold on

%% Parameters
L     = 1;
F     = 10;
beta  = 1;
Nbins = 600;

xvec  = linspace(-L, L, Nbins);
xL    = -L/2;
xR    =  L/2;

deltaxv   = linspace(0, 0.5, 10);
Ixv       = zeros(1, numel(deltaxv));
Iuv       = zeros(1, numel(deltaxv));
Wv        = zeros(1, numel(deltaxv));
DeltaImv  = zeros(1, numel(deltaxv));

%% Safe log to avoid log(0)
eps_log  = 1e-30;
safe_log = @(z) log(max(z, eps_log));

%% Stationary conditional densities P_s(x | c)
% Potentials defined on [-L, L]
V1 = @(x)  F*(abs(x) - L/2)/L;
V0 = @(x) -F*(abs(x) - L/2)/L;

% Unnormalized steady densities (two-step formula, then tabulate)
fs1_raw = @(x) exp(-beta*V1(x)) + ...
    (exp(-beta*V1(-L)) - exp(-beta*V1(L))) ./ ...
     integral(@(u) exp(-beta*(V1(L) - V1(u))), -L, L, 'RelTol',1e-6) .* ...
     integral(@(y) exp(-beta*(V1(x) - V1(y))), -L, x, 'RelTol',1e-10);

fs0_raw = @(x) exp(-beta*V0(x)) + ...
    (exp(-beta*V0(-L)) - exp(-beta*V0(L))) ./ ...
     integral(@(u) exp(-beta*(V0(L) - V0(u))), -L, L, 'RelTol',1e-6) .* ...
     integral(@(y) exp(-beta*(V0(x) - V0(y))), -L, x, 'RelTol',1e-10);

% Tabulate for robust numeric integration on [-L, L]
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
Z1  = integral(@(x) fs1(x), -L, L);
Z0  = integral(@(x) fs0(x), -L, L);
fs1 = @(x) fs1(x) ./ Z1;     % P(x|c=1)
fs0 = @(x) fs0(x) ./ Z0;     % P(x|c=0)

%% Control probabilities p(c)
pc1 = integral(@(x) fs0(x), xL, xR) / (1 - integral(@(x) fs1(x) - fs0(x), xL, xR));
pc0 = 1 - pc1;

%% Stationary marginal P_s(x)
P_s = @(x) pc1*fs1(x) + pc0*fs0(x);

%% Numerical integration tolerances (adjust if needed)
RelTol_outer = 1e-2;  AbsTol_outer = 1e-2;
RelTol_inner = 1e-2;  AbsTol_inner = 1e-2;

%% Sweep over Δx
for i = 1:numel(deltaxv)
    fprintf('Δx case %d / %d : %.4f\n', i, numel(deltaxv), deltaxv(i));
    deltax = deltaxv(i);
    % Measurement kernel Θ(x|c)
    % Protocol_Wm is external and returns Θ1 sampled on linspace(-L, L, 1000)
    Theta1 = @(x) interp1(linspace(-L, L, 1000), Protocol_Wm(deltax, 1, F, F/2), x, "linear");
    Theta0 = @(x) 1 - Theta1(x);

    %% I_x (information with mixed kernel)
    % Integrand wrappers (scalar-by-scalar) to keep integrals stable
    G = @(x,xp) fs1(xp).*Theta1(x) + fs0(xp).*Theta0(x);  % common log argument

    Integrand1 = @(x,xp) fs1(xp).*Theta1(x).*P_s(x) .* safe_log(G(x,xp));
    Integrand0 = @(x,xp) fs0(xp).*Theta0(x).*P_s(x) .* safe_log(G(x,xp));

    % Iterated integration: first over xp | fixed x, then over x
    inner1 = @(x) arrayfun(@(xx) integral(@(xp) Integrand1(xx,xp), -L, L, ...
                                          'RelTol',RelTol_inner,'AbsTol',AbsTol_inner, ...
                                          'ArrayValued', true), x);
    I1 = integral(inner1, -L, L, 'RelTol',RelTol_outer, 'AbsTol',AbsTol_outer, 'ArrayValued', true);

    inner0 = @(x) arrayfun(@(xx) integral(@(xp) Integrand0(xx,xp), -L, L, ...
                                          'RelTol',RelTol_inner,'AbsTol',AbsTol_inner, ...
                                          'ArrayValued', true), x);
    I0 = integral(inner0, -L, L, 'RelTol',RelTol_outer, 'AbsTol',AbsTol_outer, 'ArrayValued', true);

    Ixv(i) = - I1 - I0;

    %% I_u (using steady conditional densities directly)
    Iuv(i) = - integral(@(x) P_s(x).*Theta1(x).*log(fs1(x)), -L, L, 'RelTol',1e-8,'AbsTol',1e-12) ...
             - integral(@(x) P_s(x).*Theta0(x).*log(fs0(x)), -L, L, 'RelTol',1e-8,'AbsTol',1e-12);

    %% Measurement work W
    Wv(i) = integral(@(x) V1(x).*(P_s(x).*Theta1(x) - pc1*fs1(x)) + ...
                          V0(x).*(P_s(x).*Theta0(x) - pc0*fs0(x)), -L, L, ...
                     'RelTol',1e-8,'AbsTol',1e-12);

    %% ΔI_m = I^+ - I^- (channel information change)
    Imas  = integral(@(x) P_s(x).*Theta1(x).*safe_log(Theta1(x)/pc1), -L, L, 'RelTol',1e-8,'AbsTol',1e-12) ...
          + integral(@(x) P_s(x).*Theta0(x).*safe_log(Theta0(x)/pc0), -L, L, 'RelTol',1e-8,'AbsTol',1e-12);

    Imenos = integral(@(x) (pc1*fs1(x)).*safe_log( fs1(x) ./ P_s(x) ), -L, L, 'RelTol',1e-8,'AbsTol',1e-12) ...
           + integral(@(x) (pc0*fs0(x)).*safe_log( fs0(x) ./ P_s(x) ), -L, L, 'RelTol',1e-8,'AbsTol',1e-12);

    DeltaImv(i) = Imas - Imenos;
end

%% Plot
figure(11); clf; hold on; box on;
plot(deltaxv, -(Ixv - Iuv), 'ro', 'MarkerSize', 10, 'LineWidth', 2.5);
plot(deltaxv, Wv,          'k-', 'LineWidth', 2.5);
plot(deltaxv, -DeltaImv,   'k:', 'LineWidth', 2.5);

% Shade a vertical band (example values)
xinf = 0.06;     % lower x bound for V0=5 (adjust to your case)
xsup = 0.43;     % upper x bound for V0=5
yl   = [-0.5, 1.5];

fill([xinf xsup xsup xinf], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], ...
     'EdgeColor','none','FaceAlpha',0.2);
xline(xinf, 'k--', 'LineWidth', 1);
xline(xsup, 'k--', 'LineWidth', 1);

ylabel('$\beta \langle W \rangle$ and bounds','fontsize',35,'interpreter','latex');
xlabel('$\Delta x$','fontsize',35,'interpreter','latex');
set(gca,'FontSize',35,'TickLabelInterpreter','latex');
set(0,'DefaultAxesFontName','Times New Roman');
ax = gca; ax.LineWidth = 2; ax.Box = 'on'; ax.XColor = 'k'; ax.YColor = 'k';
shg;
