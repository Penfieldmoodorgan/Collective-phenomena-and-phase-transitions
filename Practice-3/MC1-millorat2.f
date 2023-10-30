c pràctica 2, problema 2 Fenòmens. Algorisme Metropolis per a 
c simular model dIsing 2D. 

      program practica
      implicit none
      integer*4 seed,l,i,j,k,imc,ipas
      parameter(l=32)
      integer*4 pbc(0:l+1)
      integer*2 s(1:l,1:l), de1
      integer*8 suma,de,mctot,n
      real*8 mag,magne,genrand_real2,ene,energ,temp,enebis
      real*8 w(-8:8)
      character*29 nom

      seed=234567
      mctot=10000
      temp=1.5d0 !Temperatura
c      temp=1.3d0 
c      temp=3.6d0 
      

c Dimensió de la matriu s

      n=l**2

c Genera vector pbc que implementa less condicions periòdiques de contorn

      pbc(0)=l
      pbc(l+1)=1

      do k=1,l
            pbc(k) = k
      end do

c genera vector w 

      do de=-8,8
            w(de)=dexp(-dfloat(de)/temp)
      end do

c genera cadena nombres aleatoris

      call init_genrand(seed) 

c omple la matriu inicial amb valors de spin iguals a +1, -1

      do i=1,l
            do j=1,l
                  if(genrand_real2().ge.0.5d0) then
                        s(i,j)=-1
                  else
                        s(i,j)=1
                  end if
            end do
      end do

c Calcula i la magnetització inicial i la energia inicial

      mag=magne(s,l)
      ene=energ(s,l,pbc)
      enebis=energ(s,l,pbc)

      nom="SIM-L-032-TEMP-1500-MCTOT-10K"
      open(UNIT=12,FILE=nom//".dat")

      imc=0

      write(12,*) imc,ene,enebis,mag

c Bucle que repeteix fins a mctot passes 
      do imc=1,mctot
c Bucle per a una passa de mc
            do ipas=1,n
            call mc(s,l,w,temp,ene,pbc)
            end do
            enebis=energ(s,l,pbc)
            mag=magne(s,l)
            write(12,*) imc,ene,enebis,mag
      end do

      close(12)

      stop
      end program

c-----------------------------------------------------------------------
c                        SUBRUTINES I FUNCIONS     
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Fa un pas Monte Carlo
c-----------------------------------------------------------------------

      subroutine mc(s,l,w,temp,ene,pbc)
      implicit none
      real*8 temp,valor,prob,delta,genrand_real2,genrand_real3,ene
      integer*2 s(1:l,1:l),canvi
      integer*8 de,suma
      real*8 w(-8:8)
      integer*4 pbc(0:l+1),i,j,l

c Escull un spin s(i,j) a l'atzar escollint valors i,j a l'atzar. 
c Emprem genrand_real2 que genera un # aleatori en l'interval real [0,1) 

      i=int(1.d0+l*genrand_real2())
      j=int(1.d0+l*genrand_real2())

c Proposem canvi -s(i,j)
      
      canvi=-s(i,j)

c Variació energia associada al canvi

      suma=s(i,pbc(j+1))+s(i,pbc(j-1))+
     +s(pbc(i+1),j)+s(pbc(i-1),j)

      de=2*s(i,j)*suma


c  si canvi energia menor que zero acceptem canvi,si canvi energia major
c o igual que zero acceptem canvi amb una probabilitat e^-de/T
c probabilitat acceptar canvi p=e^-de/T, prob no acceptar 1-p

      if (de.le.0) then
            s(i,j)=canvi
            ene=ene+de
      else 
            valor=genrand_real3()
            prob=w(de)
            if(valor.lt.prob) then
                  s(i,j)=canvi
                  ene=ene+de
            end if
      end if

      return
      end subroutine

c-----------------------------------------------------------------------
c     Obre un arxiu "P1.configuration.dat" i escriu sequencialment en 
c     una columna les parelles i,j corresponents als spins s(i,j)= +1
c-----------------------------------------------------------------------
      subroutine writeconfig(s,l)
      implicit none 
      integer*4 i,j,l
      integer*2 S(1:l,1:l)
      
      open(1,file='configuration.dat',status='unknown')

      do i=1,l
            do j=1,l
                  if(s(i,j).eq.1) then
                  write(1,*) i,j
                  end if
            end do
      end do

      close(1)

      return
      end

c-----------------------------------------------------------------------
c     funcio magnet
c-----------------------------------------------------------------------     
      
      function magne(s,l)
      implicit none
      integer*4 l,i,j
      integer*2 s(1:l,1:l)
      real*8 sum,magne

      sum=0.0d0

      do i=1,l
            do j=1,l
                  sum=sum+s(i,j)
            end do
      end do      

      magne=sum

      return
      end

c-----------------------------------------------------------------------
c     Calcula l'energia h = - SUM(SiSj) per i,j primers veïns
c-----------------------------------------------------------------------

      real*8 function energ(s,l,pbc)
      integer*4 l,i,j
      integer*4 pbc(0:l+1)
      integer*2 s(1:l,1:l)
      real*8 e

      e=0.0d0

      do i=1,l
            do j=1,l
                  e=e-s(i,j)*s(i,pbc(j+1))-s(i,j)*s(pbc(i+1),j)
c                  write(*,*) e
            end do
      end do  
      
      energ=e

      return
      end

c-----------------------------------------------------------------------
c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     initialize by an array with array-length
c     init_key is the array for initializing keys
c     key_length is its length
c-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
c
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)
     &           +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0x7fffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1]-real-interval
c-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1) with 53-bit resolution
c-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
