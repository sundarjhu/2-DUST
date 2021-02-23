!!$---------------------------------------------------------------------
FUNCTION MIPCUB(n, x, xarr, yarr, ipos)
!!$---------------------------------------------------------------------
  USE PRCSN
  IMPLICIT NONE
!!$---------------------------------------------------------------------
!!$  (M)onotonic (I)nterpolation by (P)iecewise (CUB)ic functions in one
!!$  dimension - see M. Steffen (1990), A&A 349, p.447
!!$
!!$  For given fields, XARR(1:n) and YARR(1:n), interpolation in YARR
!!$  is done at the abscisa X, which is located in the open interval
!!$  [XARR(ipos), XARR(ipos+1) ).  The position, IPOS, of X has to be 
!!$  found by some means and has to be given as an input value.  The 
!!$  result of the interpolation is returned.  The function quits when
!!$  X is out of range specified by XARR or in case of n < 3.
!!$---------------------------------------------------------------------
  INTEGER,   INTENT(IN)               :: n
  REAL(PRC), INTENT(IN), DIMENSION(n) :: xarr, yarr
  REAL(PRC), INTENT(IN)               :: x
  INTEGER,   INTENT(IN)               :: ipos
!!$---------------------------------------------------------------------
  REAL(PRC) :: mipcub,hl,hm,hr,sl,sm,sr,pm,pr,ydm,ydr,a,b,c,d,dx
!!$---------------------------------------------------------------------
  if (n >= 3) then
     if (ipos == 1) then
        hm  = xarr(ipos+1) - xarr(ipos)
        hr  = xarr(ipos+2) - xarr(ipos+1)
        sm  = (yarr(ipos+1) - yarr(ipos)) / hm
        sr  = (yarr(ipos+2) - yarr(ipos+1)) / hr
        pm  = sm * (1.0_PRC + hm / (hm+hr)) - sr * hm / (hm+hr)
        pr  = (sm * hr + sr * hm ) / (hm+hr)
        ydm = (sign(1.0_PRC,pm) + sign(1.0_PRC,sm)) &
             * min( abs(sm), 0.5_PRC * abs(pm) )
        ydr = (sign(1.0_PRC,sm) + sign(1.0_PRC,sr)) &
             * min( abs(sm), abs(sr), 0.5_PRC * abs(pr) )
     elseif (ipos == (n-1)) then
        hl  = xarr(ipos) - xarr(ipos-1)
        hm  = xarr(ipos+1) - xarr(ipos)
        sl  = (yarr(ipos) - yarr(ipos-1)) / hl
        sm  = (yarr(ipos+1) - yarr(ipos)) / hm
        pm  = (sl * hm + sm * hl) / (hl+hm)
        pr  = sm * (1.0_PRC + hm / (hl + hm)) - sl * hm / (hl+hm)
        ydm = (sign(1.0_PRC,sl) + sign(1.0_PRC,sm)) &
             * min( abs(sl), abs(sm), 0.5_PRC * abs(pm) )
        ydr = (sign(1.0_PRC,pr) + sign(1.0_PRC,sm)) &
             * min( abs(sm), 0.5_PRC * abs(pr) )
     else
        hl  = xarr(ipos) - xarr(ipos-1)
        hm  = xarr(ipos+1) - xarr(ipos)
        hr  = xarr(ipos+2) - xarr(ipos+1)
        sl  = (yarr(ipos) - yarr(ipos-1)) / hl
        sm  = (yarr(ipos+1) - yarr(ipos)) / hm
        sr  = (yarr(ipos+2) - yarr(ipos+1)) / hr
        pm  = (sl * hm + sm * hl) / (hl+hm)
        pr  = (sm * hr + sr * hm ) / (hm+hr)
        ydm = (sign(1.0_PRC,sl) + sign(1.0_PRC,sm)) &
             * min( abs(sl), abs(sm), 0.5d0 * abs(pm) )
        ydr = (sign(1.0_PRC,sm) + sign(1.0_PRC,sr)) &
             * min( abs(sm), abs(sr), 0.5_PRC * abs(pr) )
     endif
     a      = (ydm + ydr - 2.0_PRC*sm) / hm**2.0_PRC
     b      = (3.0_PRC*sm - 2.0_PRC*ydm - ydr) / hm
     c      = ydm
     d      = yarr(ipos)
     dx     = x - xarr(ipos)
     mipcub = a * dx**3.0_PRC + b * dx**2.0_PRC + c * dx + d
  else 
     if (ipos < 1) then
        mipcub = yarr(1)
     elseif ((ipos == n) .and. (x == xarr(n))) then
        mipcub = yarr(n)
     elseif ((ipos >= n) .and. (x /= xarr(n))) then
        mipcub = yarr(n)
     else
        mipcub = 0.0_PRC
     end if
  end if
end FUNCTION MIPCUB
