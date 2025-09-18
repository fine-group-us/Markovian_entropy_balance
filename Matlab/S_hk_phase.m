function Shk2=S_hk_phase(Fext,a,offset,F,Nbins,g1,g2) 
%% Computing S_hk
%This script compute: S_hk=sum_c\intdx J_s^2/P_s^2P
%Parameters
beta=1;
L=1;
xvec=linspace(-L,L,Nbins);
% Defining the potential in the interval [-L,L]
%V_1
V1= @(x) ((x+L)<a).*(F/a*((a-(x+L))))+((x+L)>a).*(F/(2*L-a).*((x+L)-a))-Fext*x-offset;
%V_0
V0= @(x) ((x+L)<a).*(F/a*(x+L))+((x+L)>a).*(F/(2*L-a).*(2*L-(x+L)))-Fext*x-offset;
% Computing rho_s(x|c=1) and rho_s(x|c=0) (without normalization constant)
fs1=@(x) exp(- beta*V1(x))+(exp(-beta*V1(-L))-exp(-beta*V1(L)))./(integral(@(x) exp( -beta*(V1(L)-V1(x))),-L,L, 'RelTol', 1e-10))*integral(@(y) exp(- beta*(V1(x)-V1(y))),-L,x, 'RelTol', 1e-10);
fs0=@(x) exp(- beta*V0(x))+(exp(-beta*V0(-L))-exp(-beta*V0(L)))./(integral(@(x) exp( -beta*(V0(L)-V0(x))),-L,L, 'RelTol', 1e-10))*integral(@(y) exp(- beta*(V0(x)-V0(y))),-L,x, 'RelTol', 1e-10);
valor1=zeros(1,Nbins);
valor0=zeros(1,Nbins);
for i=1:Nbins
    z=xvec(i);
    valor1(i)=fs1(z);
    valor0(i)=fs0(z);
end
fs1=@(x) interp1(xvec,valor1,x,"linear");
fs0=@(x) interp1(xvec,valor0,x,"linear");
% Computing J_1 and J_0 (without normalization constant, since it cancels out in Shk)
J_1= (exp(-V1(L))-exp(-V1(-L)))./(integral(@(x) exp(-(V1(L)-V1(x))),-L,L));
J_0= (exp(-V0(L))-exp(-V0(-L)))./(integral(@(x) exp(-(V0(L)-V0(x))),-L,L));
% Computing Shk
Shk2=integral(@(x) (J_1./fs1(x)).^2.*g1(x),-L,L)+integral(@(x) (J_0./fs0(x)).^2.*g2(x),-L,L);
end
