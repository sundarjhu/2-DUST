!!$---------------------------------------------------------------------
program mapgrid
!!$---------------------------------------------------------------------
!!$  Reads output (user specified) from mapspec and transform polar
!!$  grid format of mapspec output into cartesian grid format of 
!!$  user specified size and pixel scale.
!!$---------------------------------------------------------------------
  USE CONST
  implicit none
!!$---------------------------------------------------------------------
  character     :: ckflg
  character(40) :: infile,outfile1,outfile2,dummy,dummy1,dummy2
  integer :: i,j,k,l,m,n,iflag,nl,ntmp
  integer :: ierr,nwav,nlam,nphi,nrad,imxsize,imysize,xpix,ypix,nlen,nout
  integer :: xpixmax,xpixmin,ypixmax,ypixmin,count1,count2
  integer :: num,kstart,lstart,xsearchbox,ysearchbox,maxnum,minnum
  real(prc) :: dum,xmax,xmin,ymax,ymin,xpixscale,ypixscale,p,q,factor
  real(prc) :: tmpds2,tmpval1,tmpval2,tmpds22,sum1,sum2,wgt,thwidth
  real(prc) :: pixarea,rval,pval,dum1,dum2,oldds2
  integer, allocatable, dimension(:)   :: numwav
  integer, allocatable, dimension(:,:) :: npts1,flag1
  real(prc), allocatable, dimension(:)   :: rad,phi,val1,val2,ds2,lam
  real(prc), allocatable, dimension(:,:) :: flux,tau,xx,yy,fluxpix, &
       taupix,imagebox1,imagebox2,area,totarea,imagebox3,imagebox4,pix
!!$ fitsio variables
  integer oim,status,blocksize,naxes(2)
  integer :: nr,nq,mxstep,dflag,sflag,bflag,n0,nzone,ierror
  real(prc) :: vspace,rstar,tstar,distance,velocity,a,b,c,d,e,f,tau0, &
       f1,f2,f3
  real(prc), allocatable, dimension(:) :: rlayer      
!!$---------------------------------------------------------------------
!!$  Read input file
!!$---------------------------------------------------------------------
  open(unit=13,file='getmap.dat',status='old',iostat=iflag)
  if (iflag == 0) then
     read(13,'(A30)') infile
     read(13,*) imxsize, imysize
     read(13,*) xpixscale, ypixscale
     read(13,*) factor
     pixarea    = xpixscale*ypixscale ! pixel area in arcsec^2
     xsearchbox = ANINT(imxsize*factor)
     ysearchbox = ANINT(imysize*factor)
     read(13,*) nout
     allocate(numwav(nout))
     do i=1,nout
        read(13,*) numwav(i)
     end do
     close(13)
     open(unit=10,file=infile,status='old',iostat=ierr)
     if (ierr /= 0) then
        write(*,'("  Error opening the file!  Try again! ")')
        stop
     end if
  else     
     write(*,'("  ")')
     write(*,'("      *** MAPGRID ***    ")')
     write(*,'("  ")')
     do
        write(*,'(" Enter the input map file name >> ")',advance='no')
        read(*,'(A30)') infile
        open(unit=10,file=infile,status='old',iostat=ierr)
        if (ierr == 0) exit
        write(*,'("  Error opening the file!  Try again! ")')
     end do
  end if
  if (iflag /= 0) write(*,'(" READING INPUT FILE...")')   
  read(10,*) dum
  read(10,*) nwav
  do i=1,nwav
     read(10,*) dum
     if (iflag /= 0) then
        write(*,'("  ",i2,": ",ES12.6," microns")') i,dum
     end if
  end do
  if (iflag /= 0) then
     do
!!$        write(*,'(" Output to .dat (y) or just .fits (N) => ")',advance='no')
!!$        read(*,'(A)') datflag
        write(*,'(" Select multiple wavelengths? ")')
        write(*,'("  (Y / n) >> ")',advance='no') 
        read(*,'(A)') ckflg
        if (CKFLG == ' ' .OR. CKFLG == "y" .OR. CKFLG == "Y") then
           do
              write(*,'(" How many wavelengths? ")')
              write(*,'("  (number) >> ")',advance='no') 
              read(*,*) ntmp
              if (0 < ntmp .AND. ntmp <= nwav) then
                 nout = ntmp
                 exit
              else
                 write(*,'(" Invalid choice!  Select the number &
                      &again. ")')
                 write(*,'("  ")')
              end if
           end do
           exit
        else if (CKFLG == "n" .OR. CKFLG == "N") then    
           nout = 1
           exit
        end if
     end do
     allocate(numwav(nout))
     do i=1,nout
        do 
           write(*,'(" Select a wavelength by the number ")')
           write(*,'("  (number) >> ")',advance='no')
           read(*,*) ntmp
           if (0 < ntmp .AND. ntmp <= nwav) then
              numwav(i) = ntmp
              exit
           else
              write(*,'(" Invalid choice!  Select the number again. ")')
              write(*,'("   ")')
           end if
        end do
     end do
  end if
  read(10,*) nrad,nphi
  close(10)
  allocate(lam(nout))
  allocate(rad(0:nrad))
  allocate(phi(nphi))
  allocate(flux(nrad,nphi))
  allocate(tau(nrad,nphi))
  allocate(xx(nrad,nphi))
  allocate(yy(nrad,nphi))
  allocate(area(nrad,nphi))
  thwidth = pi/nphi
!!$---------------------------------------------------------------------
!!$  Partially done reading input file
!!$---------------------------------------------------------------------

  OPEN(UNIT=30,FILE='datafiles.dat',STATUS='OLD',IOSTAT=IERROR)
  READ(30,*) IERROR
  READ(30,'(A)') dummy
  CLOSE(30)
  nlen=len_trim(dummy)
  WRITE(dummy(NLEN+1:NLEN+4),'(''.dat'')')
  nlen=len_trim(dummy)
  OPEN(UNIT=4,file=dummy,STATUS='old',IOSTAT=IERROR) 
  IF (IERROR /= 0) THEN
     WRITE(*,14) dummy(1:nlen)
     WRITE(*,'(" skipping header information")')
     GO TO 111
  END IF
 14    FORMAT(' Error: File ',A,' does not exist! Try again. ')
  READ(4,*) NR,NQ          
  READ(4,*) MXSTEP,VSPACE           
  READ(4,*) DFLAG,SFLAG,BFLAG
  READ(4,*) RSTAR,TSTAR      
  READ(4,*) DISTANCE         
  READ(4,*) VELOCITY         
  READ(4,*) A,B,C,D,E,F
  READ(4,*) TAU0,N0          
  READ(4,*) F1,F2,F3         
  READ(4,*) NZONE            
  ALLOCATE(RLAYER(0:NZONE), STAT=IERROR)
  IF (NZONE /= 1) THEN
     READ(4,*) (RLAYER(I),I=1,NZONE-1)
  ENDIF
  RLAYER(0)     = 1.0_PRC ! = Rmin
  RLAYER(NZONE) = F2      ! = Rmax
  CLOSE(UNIT=4,STATUS='keep')

 111  continue

!!$---------------------------------------------------------------------
!!$  Set output image size and pixel scales
!!$---------------------------------------------------------------------
     if (iflag /= 0) then
        write(*,'(" Enter output image SIZE in even integer pixels ")')
        write(*,'("  (xsize ysize) >> ")',advance='no') 
        read(*,*) imxsize, imysize
        write(*,'("  ")')
        write(*,'("   Output images will have ",i5," x ",i5," pixels.")') &
             imxsize,imysize
        write(*,'("  ")')
        write(*,'(" Enter output image PIXEL SCALE in arcsec/pix ")')
        write(*,'("  (xpixscale ypixscale) >> ")',advance='no')
        read(*,*) xpixscale, ypixscale
        write(*,'("  ")')
        write(*,'("   Output images will cover ",es10.4," x ",es10.4, &
             &" arcsec^2.")') imxsize*xpixscale,imysize*ypixscale
        write(*,'("  ")')
        pixarea = xpixscale*ypixscale ! pixel area in arcsec^2
     end if
!!$---------------------------------------------------------------------
!!$  Done setting output image size and pixel scales
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Define a serach box
!!$---------------------------------------------------------------------
!!$      An "imagebox" has (xpixmin:xpixmax,ypixmin:ypixmax) pixels. 
!!$      Thus, the total number of pixels in the imagebox is AT MOST 
!!$      twice as many as the final image for interpolating purposes.
!!$---------------------------------------------------------------------
     if (iflag /= 0) then
        factor = 0.3_prc
        write(*,'(" Enter the search box size relative to final image ")')
        write(*,'("  (factor = 0.3 -> box is 30% of the image, max = 1) [",&
             &f6.4,"] >> ")',advance='no') factor
        read(*,'(A)') dummy
        if (dummy == ' ') then
           factor = factor
        else
           read(dummy,*) factor
           if (factor > 1.0_prc) factor = 1.0_prc
        end if
        xsearchbox = ANINT(imxsize*factor)
        ysearchbox = ANINT(imysize*factor)
        write(*,'("  ")')
        write(*,'("   Searchbox of [",i4,", ",i4,"] will be used.")') &
             xsearchbox,ysearchbox
        write(*,'("  ")')
     end if
!!$---------------------------------------------------------------------
!!$  Done defining a serach box
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Repeat the following steps to generate the desired number of maps
!!$---------------------------------------------------------------------
  do nl=1,nout
     nlam = numwav(nl)
     rad  = 0.0_prc
     phi  = 0.0_prc
     flux = 0.0_prc
     tau  = 0.0_prc
     xx   = 0.0_prc
     yy   = 0.0_prc
     area = 0.0_prc
     open(unit=10,file=infile,status='old',iostat=ierr)
     do i=1,nwav+3
        read(10,'(A30)') dummy
     end do
     do i=1,nphi
        do j=1,nrad
           if (j /= 1) then
              read(10,*) rval,dum
              rad(j) = rval 
           else
              read(10,*) rval,pval
              rad(j) = rval
              phi(i) = pval
           end if
           dum = rad(j-1)
           area(j,i) = (rval*rval-dum*dum)*thwidth*0.5_prc ! polargrid area
!!$---------------------------------------------------------------------
!!$     Read values only at the specified lambda and convert the polar 
!!$     coordinates of the grid, (rad(j), phi(i)), into the cartesian 
!!$     coordinates, (xx(j,i), yy(j,i)).  Also, convert intensity from
!!$     erg s^-1 cm^-2 Hz^-1 sr^-1 to mJy arcsec^-2.
!!$---------------------------------------------------------------------
!!$         xx(j,i) runs from -Rmax to Rmax (in arcsec)
!!$         yy(j,i) runs from     0 to Rmax (in arcsec)
!!$---------------------------------------------------------------------
           do k=1,nwav
              if (k == nlam) then
                 read(10,*) lam(nl),tau(j,i),flux(j,i)
                 flux(j,i) = flux(j,i)*2.35044305391E+15_prc 
                 xx(j,i) = rad(j) * sin(phi(i))
                 yy(j,i) = rad(j) * cos(phi(i))
              else 
                 read(10,*) dum,dum,dum
              end if
           end do
        end do
     end do
     close(10)
!!$---------------------------------------------------------------------
!!$  Find the maximum and minimum values along the cartesian grids
!!$---------------------------------------------------------------------
     xmax=MAXVAL(xx)
     xmin=MINVAL(xx)
     ymax=MAXVAL(yy)
     ymin=MINVAL(yy)
!!$  write(55,'(4es22.14)') xmax,xmin,ymax,ymin
!!$  do i=1,nphi
!!$     do j=1,nrad
!!$        write(55,'(2es22.14)') xx(j,i),yy(j,i)
!!$     end do
!!$  end do
!!$---------------------------------------------------------------------
!!$  Done reading input file
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Transfer flux/optical depth values in the polar grids into the 
!!$  corresponding cartesian grids.
!!$---------------------------------------------------------------------
!!$    Determine the corresponding grids
!!$---------------------------------------------------------------------
!!$      For a given polar grid, figure out where in the Cartesian grid 
!!$      it would fall on.  
!!$
!!$      Then, for each Cartesian grid, we do bookkeeping of;
!!$        imagebox1 : sum of intensities    (mJy arcsec^-2)
!!$        imagebox2 : sum of optical depths
!!$        imagebox3 : sum of fluxes         (mJy)
!!$        imagebox4 : sum of optical depths * area
!!$        npts1     : number of polar grid values in the cartesian grid
!!$        flag1     : = 1 if pixel has a value and = 0 if empty 
!!$        totarea   : total polar grid area tjat has fallen onto the
!!$                    cartesian grid
!!$---------------------------------------------------------------------
     xpixmax    = imxsize + ANINT(xsearchbox*0.5_prc)
     xpixmin    = -ANINT(xsearchbox*0.5_prc)
     ypixmax    = imysize + ANINT(ysearchbox*0.5_prc)
     ypixmin    = -ANINT(ysearchbox*0.5_prc)
     allocate(imagebox1(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(imagebox2(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(imagebox3(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(imagebox4(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(npts1(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(flag1(xpixmin:xpixmax,ypixmin:ypixmax))
     allocate(totarea(xpixmin:xpixmax,ypixmin:ypixmax))
     imagebox1 = 0.0_prc
     imagebox2 = 0.0_prc
     imagebox3 = 0.0_prc
     imagebox4 = 0.0_prc
     npts1     = 0
     flag1     = 0
     totarea   = 0.0_prc
     do i=1,nphi
        do j=1,nrad
           do k=1,2 ! for +/- regions of y space
              p    = xx(j,i)/xpixscale + imxsize*0.5_prc
              q    = (-1.0)**k * yy(j,i)/ypixscale + imysize*0.5_prc
              xpix = ANINT(p-0.5)
              ypix = ANINT(q-0.5)
              if ((xpix>=xpixmin .AND. xpix<=xpixmax) .AND. &
                   (ypix>=ypixmin .AND. ypix<=ypixmax)) then
                 imagebox1(xpix,ypix) = imagebox1(xpix,ypix) + flux(j,i)
                 imagebox3(xpix,ypix) = imagebox3(xpix,ypix) + &
                      flux(j,i)*area(j,i)
                 imagebox2(xpix,ypix) = imagebox2(xpix,ypix) + tau(j,i)
                 imagebox4(xpix,ypix) = imagebox4(xpix,ypix) + &
                      tau(j,i)*area(j,i)
                 npts1(xpix,ypix)     = npts1(xpix,ypix)     + 1
                 totarea(xpix,ypix)   = totarea(xpix,ypix)   + area(j,i)
                 flag1(xpix,ypix)     = 1
              end if
           end do
        end do
     end do
!!$---------------------------------------------------------------------
!!$  Done transferring polar grid values to the Cartesian grid
!!$---------------------------------------------------------------------

!!$---------------------------------------------------------------------
!!$  Determine the pixel values in the Cartesian grid
!!$---------------------------------------------------------------------
!!$    Loop through all pixels in the output image
!!$      (1) if values exist then take average
!!$      (2) if value does not exist then interpolates from nearby
!!$          (a) find valued pixels in the search box
!!$          (b) sort them in order of the distance (from the closest)
!!$          (c) use minimum number of vicinity pixels to interpolate
!!$---------------------------------------------------------------------
     count1 = 0
     count2 = 0
     allocate(pix(imxsize,imysize))
     allocate(fluxpix(imxsize,imysize))
     allocate(taupix(imxsize,imysize))
     maxnum = xsearchbox*ysearchbox
     allocate(ds2(0:maxnum))
     allocate(val1(maxnum))
     allocate(val2(maxnum))
     ds2    = float(xsearchbox*ysearchbox/4)
     ds2(0) = 0.0_prc 
     do i=1,imxsize
        do j=1,imysize
           if (flag1(i,j) == 1) then
!!$---------------------------------------------------------------------
!!$     When the current cartesian grid has polar grid values, calculate
!!$     area-averaged intensity and optical depth for that grid.
!!$---------------------------------------------------------------------
              fluxpix(i,j) = imagebox3(i,j)/totarea(i,j)
              taupix(i,j)  = imagebox4(i,j)/totarea(i,j)
!!$           write(51,*) i,j,fluxpix(i,j),taupix(i,j)
!!$---------------------------------------------------------------------
!!$     When the current cartesian grid has no polar grid value, try
!!$     interpolating the pixel value using nearby non empty pixels.
!!$---------------------------------------------------------------------
           else
              val1   = 0.0_prc
              val2   = 0.0_prc
              count1 = count1  + 1
              kstart = xpixmin + i
              lstart = ypixmin + j
              num    = 0
!!$---------------------------------------------------------------------
!!$        loop through within the "searchbox"
!!$---------------------------------------------------------------------
!!$          tmpds2: distance (squared) between the current pixel and 
!!$                  the present searchbox pixel
!!$          ds2   : distances (squared) for the searchbox pixels to be
!!$                  used in interporation (max is half the searchbox
!!$                  size)
!!$          val1  : area-averaged flux for the searchbox pixel to be
!!$                  used in interporation
!!$          val2  : area-averaged opt. depth for the searchbox pixel
!!$                  to be used in interporation
!!$---------------------------------------------------------------------
              do k=kstart+1,kstart+xsearchbox
                 do l=lstart+1,lstart+ysearchbox
                    if (flag1(k,l) == 1) then
!!$---------------------------------------------------------------------
!!$              if a value exists in the searchbox get the distance
!!$              from the current pix to the current searchbox
!!$---------------------------------------------------------------------
                       tmpds2 = (i-k)*(i-k)+(j-l)*(j-l)
                       if (tmpds2 < ds2(maxnum)) then
!!$---------------------------------------------------------------------
!!$                 if the pixel is within the allowed search radius
!!$                 get flux/opt. depth values from the searchbox pixel
!!$---------------------------------------------------------------------
                          m       = num + 1 ! temporaliry update counter
                          ds2(m)  = tmpds2
                          val1(m) = imagebox3(k,l) / totarea(k,l)
                          val2(m) = imagebox4(k,l) / totarea(k,l)
                          num     = m       ! update counter
                          if (m > 1) then
!!$---------------------------------------------------------------------
!!$                    if there are more than one values to be used for
!!$                    interpolation, sort the ds2, val2, and val2
!!$                    arrays in order of distance
!!$---------------------------------------------------------------------
                             do ! sort the ds2, val1, and val2 arrays
                                oldds2  = ds2(m-1) 
                                if (tmpds2 >= oldds2) then
                                   exit
                                else
                                   tmpds22   = ds2(m)
                                   tmpval1   = val1(m)
                                   tmpval2   = val2(m)
                                   ds2(m)    = ds2(m-1)
                                   val1(m)   = val1(m-1)
                                   val2(m)   = val2(m-1)
                                   ds2(m-1)  = tmpds22
                                   val1(m-1) = tmpval1
                                   val2(m-1) = tmpval2
                                   m         = m - 1
                                end if
                             end do
                          end if
                       end if
                    end if
                 end do
              end do
!!$---------------------------------------------------------------------
!!$        The nearby, non-empty pixels are found.  Do interpolation.
!!$---------------------------------------------------------------------
!!$          num  = number of non-empty pixels in the vicinity
!!$          val1 = flux values in the nearby pixels (area-averaged)
!!$          val2 = tau values in the nearby pixels (area-averaged)
!!$          ds2  = distance squared
!!$---------------------------------------------------------------------
              if (num /= 0) then ! when non-empty pixels are found
                 sum1 = 0.0_prc
                 sum2 = 0.0_prc
                 wgt  = 0.0_prc
!!$---------------------------------------------------------------------
!!$           calculate the interpolated pixel value giving weights
!!$           based on the distance from the current pixel
!!$
!!$           however, the maximum number of values used in the 
!!$           interpolation is limited by the number of vicinity cells
!!$           found from (1,1) within search radius "ds2".
!!$---------------------------------------------------------------------
                 if (i == 1 .AND. j == 1) minnum = num
                 do k=1,minnum
                    dum  = 1.0_prc / ds2(k)
                    sum1 = sum1 + val1(k) * dum
                    sum2 = sum2 + val2(k) * dum
                    wgt  = wgt + dum
!!$                 write(51,*) val1(k),val2(k),dum
                 end do
                 fluxpix(i,j) = sum1 / wgt
                 taupix(i,j)  = sum2 / wgt
!!$              write(51,*) i,j,fluxpix(i,j),taupix(i,j)
              else ! if there is still no pixel value
                 count2 = count2 + 1
!!$              write(51,*) i,j,fluxpix(i,j),taupix(i,j)
              end if
           end if
        end do
     end do
!!$---------------------------------------------------------------------
!!$  Display summary and output
!!$---------------------------------------------------------------------
     if (iflag /= 0) then
        write(*,'(" Interpolation summary ")')
        write(*,'("  ")')
        write(*,'("   Number of pixels:                   ",i8)') &
             imxsize*imysize
        write(*,'("   Number of zero value pixels before: ",i8)') count1
        write(*,'("   Number of zero value pixels after:  ",i8)') count2
     end if
     num = LEN_TRIM(infile)
     outfile1 = infile
     outfile2 = infile
     write(dummy1,'(g9.4)')lam(nl)
     write(dummy2,'(f5.3)')xpixscale
     dummy=dummy1(1:len_trim(dummy1))//".dat"
!!$     dummy=dummy1(1:len_trim(dummy1))//"_"//dummy2(1:len_trim(dummy2))//".dat"
     write(outfile1(num-3:num+15),'("f_",A)') dummy(1:len_trim(dummy))
     write(outfile2(num-3:num+15),'("t_",A)') dummy(1:len_trim(dummy))
!!$---------------------------------------------------------------------
!!$  Output flux is in mJy arcsec^-2
!!$---------------------------------------------------------------------
!!$    previously outputs are written in ascii
!!$---------------------------------------------------------------------
!!$     if (datflag.eq.'y' .or. datflag.eq.'Y') then
!!$       open(unit=11,file=outfile1)
!!$       open(unit=12,file=outfile2)
!!$       do i=1,imxsize
!!$          do j=1,imysize
!!$             write(11,*) fluxpix(i,j)
!!$             write(12,*) taupix(i,j)
!!$          end do
!!$       end do
!!$       close(11)
!!$       close(12)
!!$     endif
!!$---------------------------------------------------------------------
!!$    now output to fits
!!$---------------------------------------------------------------------
     outfile1=outfile1(1:len_trim(outfile1)-3)//"fits"
     outfile2=outfile2(1:len_trim(outfile2)-3)//"fits"
     status=0
     naxes(1)=imxsize
     naxes(2)=imysize
     num=naxes(1)*naxes(2)
     !
     call ftgiou(oim,status)
     !
     call ftopen(oim,outfile1,1,blocksize,status)
     if (status .eq. 0)then !file was opened;  so now delete it 
        call ftdelt(oim,status)
     else if (status .eq. 103)then ! file doesn't exist
        status=0
        call ftcmsg
     else ! there was some other error; delete anyway
        status=0
        call ftcmsg
        call ftdelt(oim,status)
     end if
!!$
     call ftinit(oim,outfile1,1,status)
     call ftphpr(oim,.true.,-32,2,naxes,0,1,.true.,status)
     call ftpkyd(oim,'SECPPIX',xpixscale,3,'arcsec per pixel',status)
     call ftpkyd(oim,'LAMBDA',lam(nl),3,'wavelength of image',status)
     call ftpkys(oim,'TYPE','Flux','mJy/arcsec^2',status)
!!$
     if (ierror==0) then
        call ftpkyj(oim,'NRAD',nr,' ',status)
        call ftpkyj(oim,'NQ',nq,' ',status)
        call ftpkyj(oim,'MXSTEP',mxstep,' ',status)
        call ftpkyd(oim,'VSPACE',vspace,3,' ',status)
        call ftpkyj(oim,'DFLAG',dflag,' ',status)
        call ftpkyj(oim,'SFLAG',sflag,' ',status)
        call ftpkyj(oim,'BFLAG',bflag,' ',status)
        call ftpkyd(oim,'RSTAR',rstar,3,' ',status)
        call ftpkyd(oim,'TSTAR',tstar,3,' ',status)
        call ftpkyd(oim,'DISTANCE',distance,3,' ',status)
        call ftpkyd(oim,'VELOCITY',velocity,3,' ',status)
        call ftpkyd(oim,'A',a,3,' ',status)
        call ftpkyd(oim,'B',b,3,' ',status)
        call ftpkyd(oim,'C',c,3,' ',status)
        call ftpkyd(oim,'D',d,3,' ',status)
        call ftpkyd(oim,'E',e,3,' ',status)
        call ftpkyd(oim,'F',f,3,' ',status)
        call ftpkyd(oim,'TAU0',tau0,3,' ',status)
        call ftpkyd(oim,'RMIN',f1,3,' ',status)
        call ftpkyd(oim,'RMAX',f2,3,' ',status)
        call ftpkyd(oim,'RSW',f3,3,' ',status)
        call ftpkyj(oim,'N0',n0,' ',status)
        call ftpkyj(oim,'NZONE',nzone,' ',status)
        do i=0,nzone,1
           write(dummy(1:7),'("RLAYER",1i1)') i
           call ftpkyd(oim,dummy(1:7),rlayer(i),3,' ',status)
        enddo
     endif
!!$
     do i=1,imxsize
        do j=1,imysize
           pix(j,i)=fluxpix(i,j)
        end do
     end do
     call ftpprd(oim,1,1,num,pix,status)
     call ftclos(oim,status)
     call ftfiou(oim,status)
!!$
     call ftgiou(oim,status)
!!$
     call ftopen(oim,outfile2,1,blocksize,status)
     if (status .eq. 0)then !file was opened;  so now delete it 
        call ftdelt(oim,status)
     else if (status .eq. 103)then ! file doesn't exist
        status=0
        call ftcmsg
     else ! there was some other error; delete anyway
        status=0
        call ftcmsg
        call ftdelt(oim,status)
     end if
!!$
     call ftinit(oim,outfile2,1,status)
     call ftphpr(oim,.true.,-32,2,naxes,0,1,.true.,status)
     call ftpkyd(oim,'SECPPIX',xpixscale,3,'arcsec per pixel',status)
     call ftpkyd(oim,'LAMBDA',lam(nl),3,'wavelength of image',status)
     call ftpkys(oim,'TYPE','Tau','optical depth',status)
     do i=1,imxsize
        do j=1,imysize
           pix(j,i)=taupix(i,j)
        end do
     end do
     call ftpprd(oim,1,1,num,pix,status)
     call ftclos(oim,status)
     call ftfiou(oim,status)
!!$
     num = LEN_TRIM(outfile1)
     if (iflag /= 0) then
        write(*,'(" ")')
        write(*,'("  Output flux map is ",A)') outfile1(1:num) 
        write(*,'("  Output tau map is  ",A)') outfile2(1:num) 
        write(*,'(" ")')
     end if
!!$
     deallocate(imagebox1,imagebox2,imagebox3,imagebox4,npts1,flag1, &
          totarea,pix,fluxpix,taupix,ds2,val1,val2)
  end do
  if (iflag /= 0) write(*,'(" DONE! ")')
!!$---------------------------------------------------------------------
!!$  Done making maps!
!!$---------------------------------------------------------------------
  deallocate(rad,phi,flux,tau,xx,yy,area,numwav,lam)
end program mapgrid
