c	subroutine MRNBC(indxx)
c  Programme to read atmospheric data from different GCM files, correct for bias and write them into a single file
	integer ndaymax,nsmax,monmax,ndmax
	integer nyrmax,nvarmax,nnmax,mvmax,nsplmax
	parameter (ndaymax=366,nsmax=6,monmax=13,ndmax=31)
	parameter (nyrmax=150,nvarmax=50,mvmax=31,nygmax=50)
	parameter (nnmax=nyrmax*ndaymax,mmmax=nyrmax*mvmax)
	parameter (nsplmax=4)
	character*200 frawc,frawf,fbcc,fbcf
	character*200 fyle,frec,fref
	real*4 gcmc(nvarmax,nyrmax,monmax,ndmax)
	real*4 gcmf(nvarmax,nyrmax,monmax,ndmax)
	real*4 rec(nvarmax,nyrmax,monmax,ndmax)
	real*4 a(300),miss
	real*4 phlwr(nvarmax),phupr(nvarmax),thres(nvarmax)
	real*4 ref(nvarmax,nyrmax,monmax,ndmax)
	integer nday(monmax),isn,idays(monmax,4)
        integer isj(nsmax,monmax),ij(monmax),leap(4),
     1        igg(nvarmax),ilimit(nvarmax),indxx
c JRP
c       integer isj(nsmax,monmax),ij(monmax),leap(4),indxx
c       real*4  igg(nvarmax),ilimit(nvarmax)

c	leap(1), LEAP years options for GCM current climate raw data
c	leap(2), LEAP years options for GCM future climate raw data
c	leap(3), LEAP years options for Observed current climate data
c	leap(4), LEAP years options for Observed future climate data
c	leap()=0, normal 365/366 days in a year format
c	leap()=1, fixed days in a month, as specified in the data  
c	ilimit()=0 no threshold below/above adjustment
c	ilimit()>0 appy threshold below/above adjustment using the threshold given
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	nout=12
	leap=0
	iseed=-1341
c       Removed write staement for parallel implementation
c       write(*,*)'    Running bias correction'
c	read necessary information from the file
	open(1,file='basic.dat',status='old')
	read(1,*)
	Read(1,*)
	read(1,*)nycur, nst
	nsc=nst-1
	read(1,*)
	read(1,'(a)')frec
	Read(1,*)
	Read(1,*)
	read(1,*)nyfut, nst
	nsf=nst-1
	read(1,*)
	read(1,'(a)')fref
 	read(1,*)
	read(1,*)
	read(1,*)ngcur,nst
	nsgc=nst-1
	read(1,*)
	read(1,'(a)')frawc
	read(1,*)
	read(1,'(a)')fbcc
	read(1,*)
	read(1,*)
	read(1,*)ngfut,nst
	read(1,*)
	read(1,'(a)')frawf
	read(1,*)
	read(1,'(a)')fbcf
	nsgf=nst-1
	call check_blank(frec)
	call check_blank(fref)
	call check_blank(frawc)
	call check_blank(fbcc)
	call check_blank(frawf)
	call check_blank(fbcf)

	read(1,*)
	read(1,*)nvar
	read(1,*)
	read(1,*)itime
	read(1,*)
	read(1,*)miss
	if(itime.eq.0)then
	read(1,*)
	read(1,*)iband
	if(iband.gt.15)stop 'Maximum iband value allowed is 15'
c	enter information whether data follows a usual leap year format(0), or a fixed days in a month format (1)
	read(1,*)
	read(1,*)(leap(i),i=1,4)
	endif
c	read physical lower and upper limits on the variables, time aggregation criteria and a threshold indicator
c	if threshold is greater than zero, the program corrects number of occurrences below and above that threshold for that variable
	read(1,*)
	read(1,*)
	do k=1,nvar
	read(1,*)i1,phlwr(k),phupr(k),igg(k),ilimit(k)
	if(ilimit(k).gt.0)then
	backspace(unit=1)
	read(1,*)i1,phlwr(k),phupr(k),igg(k),ilimit(k),thres(k)
	else
	thres(k)=0.0
	endif
	if(thres(k).eq.0.0)thres(k)=0.00000001
	enddo
c	read additional information about number of days in a year for each dataset
	if(itime.eq.0)then
	isum=0.0
	do i=1,4
	isum=isum+leap(i)
	enddo
	if(isum.gt.0)then
c	enter information as data follows a 'fixed days in a month'format (for datasets with leap years opions day entries are ignored)
	read(1,*)
	do j=1,nout
	read(1,*)(idays(j,l),l=1,4)
	enddo
	endif
	endif
c	
	close(unit=1)

	open(22,file=fbcc)
	open(23,file=fbcf)
	natm=nvar

	call day(nday,monmax)


c	read atmospheric data 
	call readdata(ngcur,ngfut,nycur,nyfut,nday,nout,nvar,nsc,
	1              nsf,nsgc,nsgf,leap,idays,frec,fref,
     2              frawc,frawf,rec,ref,gcmc,gcmf,itime,
	3              nvarmax,nyrmax,monmax,ndmax)
c	write(*,*)' Finished data reading and start threshold adjustment'

c	adjust the lower values to have same percentage of values below or above daily/monthly threshold
	call adj_thres(rec,gcmc,gcmf,nvar,nsc,nycur,nsgc,ngcur,ngfut,
     1           nsgf,miss,itime,nout,leap,idays,thres,ilimit,nday,
     2           nvarmax,nyrmax,monmax,ndmax,nnmax)
c	write(*,*)' Finished threshold adjustment & start bias corerction'

c	standardise the GCM data by using the reanalysis mean and sd
	call standardise(nsc,nsgc,nycur,ngcur,nvar,iband,itime,
     1                 ngfut,nsgf,leap,idays,nout,
     2                 miss,phlwr,phupr,thres,ilimit,rec,gcmc,gcmf,
     4                 nvarmax,nyrmax,monmax,ndmax,nsmax,nnmax)

c	write(*,*)' Finished Bias correction'


	do i=1,ngcur
	ii=0
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,1)
	if(leap(1).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgc,i,nd)
	endif
	endif
	do l=1,nd
	if(j.eq.2.and.l.eq.29)then
	do k=1,nvar
	gcmc(k,i,j,l)=gcmc(k,i,j,l-1)
	enddo
	endif
	ii=ii+1
	do k=1,nvar
	if(ilimit(k).gt.0)then
	if(gcmc(k,i,j,l).le.0.1*thres(k).and.gcmc(k,i,j,l).gt.0.0)
     1                 gcmc(k,i,j,l)=0.0
	endif

c       Removed write staements for parallel implementation
c       if(gcmc(k,i,j,l).lt.phlwr(k))
c    1  write(*,*)'Cur Climate variable takes low value',k,i,j,l
cc	if(gcmc(k,i,j,l).lt.phlwr(k))write(*,*)gcmc(k,i,j,l)
c       Removed write statement for parallel implementation
c       if(gcmc(k,i,j,l).gt.phupr(k))
c    1  write(*,*)'Cur climate variable takes high value',k,i,j,l
cc	if(gcmc(k,i,j,l).gt.phlwr(k))write(*,*)gcmc(k,i,j,l)

	if(gcmc(k,i,j,l).lt.phlwr(k))gcmc(k,i,j,l)=phlwr(k)
	if(gcmc(k,i,j,l).gt.phupr(k))gcmc(k,i,j,l)=phupr(k)
	a(k)=gcmc(k,i,j,l)
	enddo
	if(itime.eq.0)write(22,134)i+nsgc,j,l,(a(k),k=1,natm)
	if(itime.ne.0)write(22,136)i+nsgc,j,(a(k),k=1,natm)
	enddo
	enddo
	enddo
	close(unit=22)
c	for future climate
	do i=1,ngfut
	ii=0
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,2)
	if(leap(2).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgf,i,nd)
	endif
	endif
	do l=1,nd
	ii=ii+1
	if(j.eq.2.and.l.eq.29)then
	do k=1,nvar
	gcmf(k,i,j,l)=gcmf(k,i,j,l-1)
	enddo
	endif

	do k=1,nvar
	if(ilimit(k).gt.0)then
	if(gcmf(k,i,j,l).le.0.1*thres(k).and.gcmf(k,i,j,l).gt.0.0)
     1                 gcmf(k,i,j,l)=0.0
	endif

c       Removed writre statements for parallel implementation
c       if(gcmf(k,i,j,l).lt.phlwr(k))
c     1  write(*,*)'Fut Climate variable takes low value',k,i,j,l
cc	if(gcmf(k,i,j,l).lt.phlwr(k))write(*,*)gcmf(k,i,j,l)
c       if(gcmf(k,i,j,l).gt.phupr(k))
c     1  write(*,*)'Fut climate variable takes high value',k,i,j,l
cc	if(gcmf(k,i,j,l).gt.phlwr(k))write(*,*)gcmf(k,i,j,l)

	if(gcmf(k,i,j,l).lt.phlwr(k))gcmf(k,i,j,l)=phlwr(k)
	if(gcmf(k,i,j,l).gt.phupr(k))gcmf(k,i,j,l)=phupr(k)
	a(k)=gcmf(k,i,j,l)
	enddo
	if(itime.eq.0)write(23,134)i+nsgf,j,l,(a(k),k=1,natm)
	if(itime.ne.0)write(23,136)i+nsgf,j,(a(k),k=1,natm)
	enddo
	enddo
	enddo
	close(unit=23)

c*************************************************
c	write(*,*)' Calculating statistics of BIAS CORRECTED data'
c	fyle='statbcc.dat'
c	ilpg=1
c	ilph=3
c	call stats(rec,gcmc,nvar,nsc,nycur,nsgc,ngcur,miss,fyle,isn,
c	1          itime,ij,isj,nout,leap,idays,igg,ilph,ilpg,thres,ilimit,
c     2          nnmax,nsmax,nvarmax,nyrmax,monmax,ndmax,nsplmax)
c	fyle='statbcf.dat'
c	ilpg=2
c	ilph=4
c	call stats(ref,gcmf,nvar,nsf,nyfut,nsgf,ngfut,miss,fyle,isn,
c	1          itime,ij,isj,nout,leap,idays,igg,ilph,ilpg,thres,ilimit,
c     2          nnmax,nsmax,nvarmax,nyrmax,monmax,ndmax,nsplmax)
c*********************************************
c      Removed write statement for parallel implementation
c       write(*,*)'    Bias correction is over'
 136	format(i5,i3,300g18.8)
 134	format(i5,2i3,300g18.8)
 132	format(' yr mon day',70i12)
c	indxx=1
c	return
	end
c************************************************
	subroutine check_blank(fre)
	character*200 fre,fra
        fra=' '
	do i=1,200
        if(fre(i:i).eq.' '.or.fre(i:i).eq.'	')goto 10
c       if(fre(i:i).eq.' '.or.fre(i:i).eq.'     ')goto 10
	ii=i
	goto 100
  10    continue
	enddo

100     j=0
	do i=ii,200,1
	j=j+1
	fra(j:j)=fre(i:i)
	enddo
        fre=' '
	do i=1,200
	fre(i:i)=fra(i:i)
	enddo
	return
	end

c**************************************************
	subroutine readdata(ngcur,ngfut,nycur,nyfut,nday,nout,nvar,
	1                    nsc,nsf,nsgc,nsgf,leap,idays,frec,
     2                    fref,frawc,frawf,rec,ref,gcmc,gcmf,itime,
	3                    nvarmax,nyrmax,monmax,ndmax)

	character*280 head,frec*200,fref*200,frawc*200,frawf*200
	real*4 gcmc(nvarmax,nyrmax,monmax,ndmax)
	real*4 gcmf(nvarmax,nyrmax,monmax,ndmax)
	real*4 rec(nvarmax,nyrmax,monmax,ndmax)
	real*4 ref(nvarmax,nyrmax,monmax,ndmax)
	real*4 avr(nvarmax),ano(nvarmax)
	real*4 b(300)

	integer nday(monmax),leap(4),idays(monmax,4)

	avr=0.0
	ano=0.0
c	read RAW calibration data
	open(1,file=frawc,status='old')
c	read(1,'(A)')HEAD
c	write(22,'(A)')HEAD
 41	read(1,*,err=99,end=19)i1,i2
	if(i1.lt.nsgc+1)goto 41
	backspace(unit=1)
	do i=1,ngcur
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,1)
	if(leap(1).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgc,i,nd)
	endif
	endif
	do l=1,nd
	if(itime.eq.0)read(1,*,err=99,end=19)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(1,*,err=99,end=19)i1,i2,(b(k),k=1,nvar)
	if(itime.eq.0)then
	if(i+nsgc.ne.i1.or.j.ne.i2.or.l.ne.i3)
     1 write(*,*)'Check date of raw calibration data',i1,i2,i3
	if(i+nsgc.ne.i1.or.j.ne.i2.or.l.ne.i3)stop
	else
	if(i+nsgc.ne.i1.or.j.ne.i2)write(*,*)'Check date rawc',i1,i2
	if(i+nsgc.ne.i1.or.j.ne.i2)stop
	endif
	do k=1,nvar
	if(b(k).gt.-990)avr(k)=avr(k)+b(k)
	if(b(k).gt.-990)ano(k)=ano(k)+1
	gcmc(k,i,j,l)=b(k)
	enddo
	enddo
	enddo
	enddo
19 	if(i.lt.ngcur)
     1 write(*,*)'Only a part of current raw data is read'
	close(unit=1)
	do k=1,nvar
	if(ano(k).gt.0)avr(k)=avr(k)/ano(k)
	enddo
	do i=1,ngcur
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,1)
	if(leap(1).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgc,i,nd)
	endif
	endif
	do l=1,nd
	do k=1,nvar
	if(gcmc(k,i,j,l).lt.-990)gcmc(k,i,j,l)=avr(k)
	enddo
	enddo
	enddo
	enddo

c	read RAW - to be bias corrected data
	open(1,file=frawf,status='old')
c	read(1,'(A)')HEAD
c	write(23,'(A)')HEAD
 42	read(1,*,err=99,end=119)i1,i2
	if(i1.lt.nsgf+1)goto 42
	backspace(unit=1)
	avr=0.0
	ano=0.0
	do i=1,ngfut
	ii=0
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,2)
	if(leap(2).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgf,i,nd)
	endif
	endif
	do l=1,nd
	ii=ii+1
	if(itime.eq.0)read(1,*,err=99,end=119)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(1,*,err=99,end=119)i1,i2,(b(k),k=1,nvar)
	if(itime.eq.0)then
	if(i+nsgf.ne.i1.or.j.ne.i2.or.l.ne.i3)
     1        write(*,*)'Check date Rawf',i1,i2,i3
	if(i+nsgf.ne.i1.or.j.ne.i2.or.l.ne.i3)stop
	else
	if(i+nsgf.ne.i1.or.j.ne.i2)write(*,*)'Check date Rawf',i1,i2
	if(i+nsgf.ne.i1.or.j.ne.i2)stop
	endif
	do k=1,nvar
	if(b(k).gt.-990)avr(k)=avr(k)+b(k)
	if(b(k).gt.-990)ano(k)=ano(k)+1
	gcmf(k,i,j,l)=b(k)
	enddo
	enddo
	enddo
	enddo
119 	if(i.lt.ngfut)
	1 write(*,*)'Only a part of future raw data is read'
 	close(unit=1)
	do k=1,nvar
	if(ano(k).gt.0)avr(k)=avr(k)/ano(k)
	enddo
	do i=1,ngfut
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,2)
	if(leap(2).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsgf,i,nd)
	endif
	endif
	do l=1,nd
	do k=1,nvar
	if(gcmf(k,i,j,l).lt.-990)gcmf(k,i,j,l)=avr(k)
	enddo
	enddo
	enddo
	enddo
c	read reanalysis calibration period data
	AVR=0.0
	ano=0.0
	open(1,file=frec,status='old')
43	read(1,*,err=99,end=121)i1,i2
	if(i1.lt.nsc+1)goto 43
	backspace(unit=1)
	do i=1,nycur
	ii=0
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,3)
	if(leap(3).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsc,i,nd)
	endif
	endif
	do l=1,nd
	ii=ii+1
	if(itime.eq.0)read(1,*,err=98,end=121)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(1,*,err=98,end=121)i1,i2,(b(k),k=1,nvar)
	if(itime.eq.0)then
	if(i+nsc.ne.i1.or.j.ne.i2.or.l.ne.i3)
	1  write(*,*)'Check date Obsc',i1,i2,i3,i+nsc,j,l
	if(i+nsc.ne.i1.or.j.ne.i2.or.l.ne.i3)stop
	else
	if(i+nsc.ne.i1.or.j.ne.i2)write(*,*)'Check date Obsc',i1,i2
	if(i+nsc.ne.i1.or.j.ne.i2)stop
	endif
221	do k=1,nvar
	if(b(k).gt.-990)avr(k)=avr(k)+b(k)
	if(b(k).gt.-990)ano(k)=ano(k)+1
	rec(k,i,j,l)=b(k)
	enddo
	enddo
	enddo
	enddo
121 	if(i.lt.nycur)
	1 write(*,*)'Only a part of current observed data is read'
	close(unit=1)
	do k=1,nvar
	if(ano(k).gt.0)avr(k)=avr(k)/ano(k)
	enddo
	do i=1,nycur
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,3)
	if(leap(3).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsc,i,nd)
	endif
	endif
	do l=1,nd
	do k=1,nvar
	if(rec(k,i,j,l).lt.-990)rec(k,i,j,l)=avr(k)
	enddo
	enddo
	enddo
	enddo

c	read reanalysis validation period data
	avr=0.0
	ano=0.0
	open(1,file=fref,status='old')
44	read(1,*,err=99,end=122)i1,i2
	if(i1.lt.nsf+1)goto 44
	backspace(unit=1)
	do i=1,nyfut
	ii=0
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,4)
	if(leap(4).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsf,i,nd)
	endif
	endif
	do l=1,nd
	ii=ii+1
	if(itime.eq.0)read(1,*,err=97,end=122)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(1,*,err=97,end=122)i1,i2,(b(k),k=1,nvar)
	if(itime.eq.0)then
	if(i+nsf.ne.i1.or.j.ne.i2.or.l.ne.i3)write(*,*)
	1 'Check date Obsf',i1,i2,i3,i+nsf,j,l
	if(i+nsf.ne.i1.or.j.ne.i2.or.l.ne.i3)stop
	else
	if(i+nsf.ne.i1.or.j.ne.i2)write(*,*)'Check date Obsf',i1,i2
	if(i+nsf.ne.i1.or.j.ne.i2)stop
	endif
222	do k=1,nvar
	if(b(k).gt.-990)avr(k)=avr(k)+b(k)
	if(b(k).gt.-990)ano(k)=ano(k)+1
	ref(k,i,j,l)=b(k)
	enddo
	enddo
	enddo
	enddo
122 	if(i.lt.nyfut)
	1 write(*,*)'Only a part of future observed data is read'
	close(unit=1)
	do k=1,nvar
	if(ano(k).gt.0)avr(k)=avr(k)/ano(k)
	enddo
	do i=1,nyfut
	do j=1,nout
	if(itime.ne.0)then
	nd=1
	else
	nd=idays(j,4)
	if(leap(4).eq.0)then
	nd=nday(j)
	if(j.eq.2)call daycount(nsf,i,nd)
	endif
	endif
	do l=1,nd
	do k=1,nvar
	if(ref(k,i,j,l).lt.-990)ref(k,i,j,l)=avr(k)
	enddo
	enddo
	enddo
	enddo

	return
99	backspace(unit=1)
	read(1,'(a)')head
	write(*,*)nvar
	write(*,'(a)')head
	stop
98	backspace(unit=1) 
	read(1,'(a)')head
	do i1=1,280
	if(head(i1:i1).eq.'N'.and.head(i1+2:i1+2).eq.'N')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'n'.and.head(i1+2:i1+2).eq.'N')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'N'.and.head(i1+2:i1+2).eq.'n')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'n'.and.head(i1+2:i1+2).eq.'n')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo

	if(itime.eq.0)read(head,*)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(head,*)i1,i2,(b(k),k=1,nvar)
	do k=1,nvar
	if(b(k).eq.-99.0)b(k)=-999.0
	enddo
	goto 221
97	backspace(unit=1) 
	read(1,'(a)')head
	do i1=1,280
	if(head(i1:i1).eq.'N'.and.head(i1+2:i1+2).eq.'N')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'n'.and.head(i1+2:i1+2).eq.'N')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'N'.and.head(i1+2:i1+2).eq.'n')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo
	do i1=1,280
	if(head(i1:i1).eq.'n'.and.head(i1+2:i1+2).eq.'n')then
	head(i1:i1)='-'
	head(i1+1:i1+1)='9'
	head(i1+2:i1+2)='9'
	endif
	enddo

	if(itime.eq.0)read(head,*)i1,i2,i3,(b(k),k=1,nvar)
	if(itime.ne.0)read(head,*)i1,i2,(b(k),k=1,nvar)
	do k=1,nvar
	if(b(k).eq.-99.0)b(k)=-999.0
	enddo
	goto 222
	return
	end
c**************************************************************
