!========Subrutine for building the probabilty==============
! This subroutine builds the probability distribution
!-----------------------------------------------------------
      subroutine probability(x,L,Nsim, nbins,Px)

      integer :: i,k,n
      integer, intent(in) :: nbins, Nsim
      real*8, intent(in) :: x,L
      double precision, intent(out):: Px(nbins)
      real*8 ::  delta_L
      
  
      delta_L=L/dble(nbins)
      n=int(min(floor(dble(x)/dble(delta_L)), nbins-1))
      Px(n+1)=Px(n+1)+1.d0/dble(Nsim)

      
     
      end subroutine probability

      
