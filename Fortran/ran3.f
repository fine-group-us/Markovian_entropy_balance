C====================================================================
!     Subrutine Ran3 random number
!====================================================================
!       PROGRAM: ran3.f
!       TYPE   : function
!       PURPOSE: generate random numbers
!       VERSION: 17 June 94
!       COMMENT: Initialize idum with negative integer
!======================================================================
        real function ran3(idum)
        integer mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp,iff,idum
        double precision fac
        Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
        Dimension MA(55)
        save
        if (idum.lt.0.or.iff.eq.0) then
	   iff=1
	   mj=mseed-iabs(idum)
	   mj=mod(mj,mbig)
	   ma(55)=mj
	   mk=1
	   do 11 i=1,54
	      ii=mod(21*i,55)
	      ma(ii)=mk
	      mk=mj-mk
	      if (mk.lt.mz) mk=mk+mbig
	      mj=ma(ii)
 11	   continue
	   do 13 k=1,4
	      do 12 i=1,55
		 ma(i)=ma(i)-ma(1+mod(i+30,55))
		 if (ma(i).lt.mz) ma(i)=ma(i)+mbig
 12	      continue
 13	   continue
	   inext=0
	   inextp=31
	   idum=1
	end if
	inext=inext+1
	if (inext.eq.56) inext=1
	inextp=inextp+1
	if (inextp.eq.56) inextp=1
	mj=ma(inext)-ma(inextp)
	if (mj.lt.mz) mj=mj+mbig
	ma(inext)=mj
	ran3=mj*fac
	return
	end
