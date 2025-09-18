!!========Subrutine for division and binary-to-decimal conversion==============

      subroutine binary_fun(chain,M, PDF)

      integer :: i,k
      integer, intent(in) :: M
      logical, intent(in) :: chain(M)
      double precision, intent(out):: PDF(2**M)
      real*8 :: two(M), biyec
      
      do k=1,M
         two(k)=2.d0**(dble(k-1))
      enddo
     
    
         biyec=1.d0+sum(-dble(chain)*two)
         PDF(biyec)=PDF(int(biyec))+1.d0
     
         
     
      end subroutine binary_fun

      
      
∫∫