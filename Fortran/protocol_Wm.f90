!========Subrutine for the Wm protocol==============

! This subroutine calculates the probability p of
! switching the potential state, Theta(c=1|x), in the
! Wm protocol, given the current position x of the
! particle, the potential parameters (L, a, V0), the
! the potential state (on/off), and the measurement
! error (deltax).

!---------------------------------------------------

      subroutine protocol_Wm(x, deltax,L,p,xl,xr)

      integer :: i,k,n
      real*8, intent(in) :: x,L, deltax,xl,xr
      real*8, intent(out)::p
      real*8::  ll, lu
      
      lu=modulo(x+deltax,L)
      ll=modulo(x-deltax,L)
      if (ll.LE. lu) then
         if (lu .LE.xl) then
            p=0
         elseif ((xl.LE.lu).AND.(lu.LE.xr)) then
            if (ll .LE. xl) then
               p=(lu-xl)/(2.d0*deltax)
            else
               p=(lu-ll)/(2.d0*deltax)
            endif
            else !lu > xr
               if (ll.LE.xl) then
                    p=(xr-xl)/(2.d0*deltax);
                elseif ((xl.LE.ll).AND.(ll.LE.xr)) then
                    p=(xr-ll)/(2.d0*deltax);
                else !ll>xr
                    p=0.d0;
                endif
           endif
      else !ll>lu
         if (ll.LE.xl) then
                p=(xr-xl)/(2.d0*deltax);
            elseif ((xl.LE.ll).AND.(ll.LE.xr))then
                if (lu.LE.xl) then
                    p=(xr-ll)/(2.d0*deltax);
                else !lu>xl
                    p= (xr-xl+lu-ll)/(2*deltax);
                endif
            else !ll>xr
                if (lu.LE.xl)then
                    p=0.d0;
                elseif ((xl.LE.lu).AND.(lu.LE.xr)) then
                    p=(lu-xl)/(2.d0*deltax);
                else !lu>xr
                    p=(xr-xl)/(2.d0*deltax);
                endif
            endif
        endif
      
     
      end subroutine protocol_Wm
