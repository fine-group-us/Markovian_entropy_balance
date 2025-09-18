program langevine_evol_x
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! ======= External routines you already provide elsewhere =======
  real ran3
  external :: protocol_Wm       ! call protocol_Wm(x, deltax, L, p, xl, xr)
  external :: probability      ! call probability(x, L, Nsim, nbins, PCx)
  external :: binary_fun        ! call binary_fun(chain(:), M, PDF)

  ! ======= Simulation sizes & parameters =======
  integer*8, parameter :: Nsim=100000, Nstep=int(10**3), nummed=100
  integer*8, parameter :: nummed_transit=int(dble(nummed)/2.d0)
  integer*8, parameter :: nummed_steady=nummed-nummed_transit
  integer*8, parameter :: nt=Nstep*nummed, nbins=600, M=10
  real*8,    parameter :: deltat=5*10.d0**(-4)
  real*8,    parameter :: deltatm=deltat*dble(Nstep)
  real*8,    parameter :: Fext=0.d0*10.d0**(0), V0=0.1d0, a=1.0d0
  real*8,    parameter :: gamma=1.d0, beta=1.d0, D=1.d0, L=2.d0
  real*8,    parameter :: deltax=0.d0, tol0=10.d0**(-100)
  real*8,    parameter :: ht=deltatm/dble(Nstep)

  ! ======= State, accumulators, and helpers =======
  integer*8 :: i, j, k, n, iseed
  real*8    :: t1, t2, t3, xl, xr, pi, p, Fin, Feff, u1, u2, norm
  real*8    :: Sseqb, Sseqs, PC
  logical   :: control(Nsim), chainb(Nsim,M), chains(Nsim,M-1)
  real*8    :: x(Nsim), PC1x(nbins), PC0x(nbins)
  integer*8 :: auxiliar(Nsim)
  double precision :: PDFb(2**M), PDFs(2**(M-1))
  double precision :: IM(nummed_steady), I_after(nummed_steady), Sc(nummed_steady)
  double precision :: normalizacionb, normalizacions, Wext, Wint

  ! ======= Output names =======
  character(len=100) :: folder_name,  file_name_parameters,  file_name_pc1,  file_name_pc0
  character(len=100) :: file_name_pc, file_name_H_p, file_name_SM
  character(len=10)  :: deltatm_str, M_str, Nsim_str, nummed_str, Fext_str, a_str, Nbins_str
  character(len=10)  :: Nstep_str, V0_str, deltax_str

  ! ======= Initialization =======
  iseed = -2345678
  pi    = 4.d0*atan(1.d0)

  ! Characteristic times 
  t1 = L**2/(2*D)
  t2 = gamma*L/(V0/a+Fext)
  t3 = gamma*L/((L-a)/V0+Fext)

  ! Strings for folder naming 
  write(deltatm_str, '(ES10.1)') deltatm
  write(M_str,       '(I10)')    M
  write(Nsim_str,    '(ES10.1)') real(Nsim)
  write(Nbins_str,   '(I10)')    nbins
  write(Nstep_str,   '(ES10.1)') real(Nstep)
  write(nummed_str,  '(ES10.1)') real(nummed)
  write(a_str,       '(ES10.1)') real(a)
  write(Fext_str,    '(ES10.1)') Fext
  write(V0_str,      '(ES10.1)') V0
  write(deltax_str,  '(ES10.1)') deltax

  Nsim_str    = adjustl(Nsim_str)
  Nstep_str   = adjustl(Nstep_str)
  Nbins_str   = adjustl(Nbins_str)
  M_str       = adjustl(M_str)
  deltatm_str = adjustl(deltatm_str)
  nummed_str  = adjustl(nummed_str)
  a_str       = adjustl(a_str)
  Fext_str    = adjustl(Fext_str)
  V0_str      = adjustl(V0_str)
  deltax_str  = adjustl(deltax_str)

  write(folder_name, '(A,A,A,A,A,A,A,A,A,A,A,A,A)') 'DATA_phases/Wm/V0_',trim(V0_str), &
       '/Deltax_',trim(deltax_str),'/',trim(Nsim_str),trim(deltatm_str),trim(nummed_str),  &
       trim(Nstep_str), trim(M_str), trim(Nbins_str),trim(Fext_str), trim(a_str)

  call system('mkdir -p ' // trim(folder_name))

  file_name_parameters = trim(folder_name) // '/Param.txt'
  file_name_pc1        = trim(folder_name) // '/PC1.dat'
  file_name_pc0        = trim(folder_name) // '/PC0.dat'
  file_name_pc         = trim(folder_name) // '/P_C.dat'
  file_name_H_p        = trim(folder_name) // '/H_p.dat'
  file_name_SM         = trim(folder_name) // '/S_M.dat'

  open(23,file=file_name_parameters,status="unknown")
  open(24,file=file_name_pc1,status="unknown")
  open(25,file=file_name_pc0,status="unknown")
  open(26,file=file_name_pc,status="unknown")
  open(27,file=file_name_H_p,status="unknown")
  open(28,file=file_name_SM,status="unknown")

  write(23,*) 'Nsim= ', Nsim, ' Nstep= ',Nstep, ' delta t= ', deltat, ' delta t_m= ', deltatm, &
              ' F_ext= ', Fext, ' V0= ', V0, ' a= ', a, ' Nummed= ', nummed, ' Nbins= ', nbins, ' M= ', M

  ! Reset PDFs / normalizations
  do k=1,2**M
     PDFb(k)=0.d0
  enddo
  normalizacionb=0.d0

  do k=1,2**(M-1)
     PDFs(k)=0.d0
  enddo
  normalizacions=0.d0

  ! Initialize positions uniformly in [0,L) and control flags
  do k=1,Nsim
     x(k)       = L*ran3(iseed)
     control(k) = .FALSE.
     auxiliar(k)= 0_8
  enddo

  ! Zero accumulators
  do n=1,nummed_steady
     IM(n)=0.d0; I_after(n)=0.d0; Sc(n)=0.d0
  enddo
  Wext=0.d0; Wint=0.d0

  ! ======= Main temporal loop =======
  do i=1,nt

     ! Reset histograms at each time step
     do k=1,nbins
        PC1x(k)=0.d0
        PC0x(k)=0.d0
     enddo
     PC=0.d0

     ! Potential region limits
     xl=a/2.d0
     xr=(L+a)/2.d0

     do n=1,Nsim

        if (mod(i,Nstep)==0) then
           ! ---- Measurement step ----
           auxiliar(n)=auxiliar(n)+1_8

           call protocol_Wm(x(n), deltax, L, p, xl, xr)

           u1 = ran3(iseed)
           if (u1 .GE. p) then
              control(n) = .FALSE.
              if (x(n)<a) then
                 Fin=-V0/a
              else
                 Fin= V0/(L-a)
              endif
           else
              control(n) = .TRUE.
              if (x(n)<a) then
                 Fin= V0/a
              else
                 Fin=-V0/(L-a)
              endif
           endif

           ! Build joint P(x,C) and marginal P(C) using x BEFORE evolution
           if (i>nummed_transit*Nstep) then
              if (control(n)) then
                 call probability(x(n), L, Nsim, nbins, PC1x)
                 PC=PC+1.d0/dble(Nsim)
              else
                 call probability(x(n), L, Nsim, nbins, PC0x)
              endif
           endif

           ! One Euler–Maruyama step with Gaussian noise (Box–Muller via ran3)
           u1 = ran3(iseed)
           u2 = ran3(iseed)
           do while (u1==0.d0)
              u1 = ran3(iseed)
           enddo
           norm = sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)

           Feff=Fin+Fext
           x(n)=x(n)+Feff*deltat+norm*sqrt(2.d0)*sqrt(deltat)
           x(n)=modulo(x(n),L)

           ! ---- Build chain PDFs at steady state ----
           if ((auxiliar(n) .GE. nummed_transit) .AND. ((auxiliar(n)-M) < nummed_transit)) then
              chainb(n, auxiliar(n)-nummed_transit) = control(n)
           endif

           if ((auxiliar(n)-M) .GE. nummed_transit) then
              do j=1,(M-1)
                 chains(n,j)=chainb(n,j+1)
                 chainb(n,j)=chainb(n,j+1)
              enddo
              chainb(n,M)=control(n)

              call binary_fun(chainb(n,:), M,   PDFb)
              normalizacionb=normalizacionb+1.d0
              call binary_fun(chains(n,:), M-1, PDFs)
              normalizacions=normalizacions+1.d0
           endif

        else
           ! ---- Free evolution (no measurement) ----
           if (.NOT. control(n)) then
              if (x(n)<a) then
                 Fin=-V0/a
              else
                 Fin= V0/(L-a)
              endif
           else
              if (x(n)<a) then
                 Fin= V0/a
              else
                 Fin=-V0/(L-a)
              endif
           endif

           ! Build P(x,C) and P(C) during free evolution (steady-state only)
           if (i>nummed_transit*Nstep) then
              if (control(n)) then
                 call probability(x(n), L, Nsim, nbins, PC1x)
                 PC=PC+1.d0/dble(Nsim)
              else
                 call probability(x(n), L, Nsim, nbins, PC0x)
              endif
           endif

           u1 = ran3(iseed)
           u2 = ran3(iseed)
           do while (u1==0.d0)
              u1 = ran3(iseed)
           enddo
           norm = sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)

           Feff=Fin+Fext
           x(n)=x(n)+Feff*deltat+norm*sqrt(2.d0)*sqrt(deltat)
           x(n)=modulo(x(n),L)

        endif

     enddo  ! trajectories loop

     ! ---- Write histograms once in steady regime ----
     if (i>(nummed_transit*Nstep)) then
        write(24,199) PC1x
        write(25,199) PC0x
        write(26,199) PC
     endif
     !---- Progress report ----
     if ( floor(dble(i)/dble(nt)*100.d0) > floor(dble(i-1)/dble(nt)*100.d0) ) then
        print*, floor(dble(i)/dble(nt)*100.d0), '%'
     endif

  enddo  ! temporal loop

  ! ======= Sequence entropies =======
  if (normalizacionb>0.d0) PDFb=PDFb/normalizacionb
  if (normalizacions>0.d0) PDFs=PDFs/normalizacions

  Sseqb=0.d0
  do k=1,2**M
     if (PDFb(k)>0.d0) Sseqb = Sseqb - PDFb(k)*log(PDFb(k))
  enddo

  Sseqs=0.d0
  do k=1,2**(M-1)
     if (PDFs(k)>0.d0) Sseqs = Sseqs - PDFs(k)*log(PDFs(k))
  enddo

  write(27,199) Sseqb - Sseqs
  write(28,199) Sseqb

199 format(1000000000(e14.7,2x))
end program langevine_evol_x
