                 program main

                IMPLICIT REAL*8 (A-H,O-T,V-Z)

                write(6,*) 'eneutrino gev'
                read(5,*)eingev
                write(6,*)' coslepton, xmagev/1.049, RPA (0/1)'
                read(5,*)  cost, xx,ipol
                 xmagev = xx*1.049d0
                write(6,*) ' position as a fraction of RMax'
                read(5,*) rmaxfrac


       delta = (eingev-0.025)/100.d0
       do it = 0,100
          tmugev = it*delta
          
          call amunu_calc(eingev,tmugev,cost,xmagev,ipol,
     f      rmaxfrac)
          
c     Print Q2, form factors, tensor elts
c     No header is printed because this file is not meant to be read. A
c     script will label and plot the results

       end do
          
 13             format ('#', 1x, a,1x,a)
 11               format(3(e14.7,1x),i1,1x,e14.7)
 12               format('#', 1x, a,1x,a)
 2               format('#', 1x, a,1x,e9.3)
 1               format('#', 1x, a,1x,'&',1x,a)

               stop
               end
******************************************************************************
      subroutine amunu_calc(eingev,tmugev,cost,xmagev,ipol,
     f        rmaxfrac)
        IMPLICIT REAL*8 (A-H,O-T,V-Z)
        IMPLICIT COMPLEX*16 (U)
        DIMENSION WD(5) ! funciones de respuesta
        dimension dr(2000),drp(2000),df0(2000),vcd(2000)
        
        COMMON/NXRO/NXRO
        COMMON/VFPI/DROP,DRON
        COMMON/CTEROPN/DXP,DA0P,DRO0P,DXN,DA0N,DRO0N

        COMMON/NCXRO/NCXRO
        COMMON/CTEROPNC/DNCXP,DNCA0P,DRONC0P



        COMMON/ROFEROS/KLAVE
        COMMON/NUC/DZZ,DAA

        common/datos/dpi,hbarc,GF0,DMNU,DMA
        common/datos2/dmneutrino,dmlepton,dmi,dmf,
     f                coscabibbo

        common/calculo/ic
        common/vcfto/vcd
        common/importantisimo/xma


        if (ipol.ne.0 .and. ipol.ne.1) stop'wrong ipol'
        ic = 0
        DPI=3.141592653589793D0
        hbarc = 197.3269602d0
        xuma = 931.494013 / hbarc
        dmneutrino = 0.d0 / hbarc
        dmuon = 105.658357d0 / hbarc ! masa muon
        dmelectron = 0.510998902d0 / hbarc ! masa electron         
        dmnu = 940.d0 / hbarc
        coscabibbo=0.974d0
        GF0= 1.1664d-11*hbarc**2 ! units fm^2
        xma= xmagev*1.d3/ hbarc
     
c


         READ(4,*)DZZ,DAA,KLAVE,DNCXP,DNCA0P,DNCXN,DNCA0N
         DMA = DAA*xuma ! masa del nucleo.
         DNN = DAA-DZZ
         read(4,*) qvalue,lepton
        
        if (lepton.eq.1) then 
        dmlepton = dmuon 

        elseif (lepton.eq.0) then
        dmlepton = dmelectron 
          else
          stop'lepton incorrecto' 
        endif
         
        dmi = dmneutrino
        dmf = dmlepton
c       ieta = 1 


ccccccc         read (4,*) mn,mnr 
           mn=2
           mnr=1
c         READ(4,*)ILIN
         read(4,*)ieta

         rewind(4)          
 
 


                
c********************************************************************
c********************************************************************
         ilin=2
c********************************************************************
c*******************************************************************



        nxro=0 ! normaliza densidad
        ncxro=0 ! normaliza densidad

 
c         lc = 0 ! 0 desconvoluciona cualquier otro numero no
c     CONVOLUTION TURNED OFF
         lc = 1
         call convolucion (dncxp,dnca0p,dncxn,dnca0n,klave,lc,
     f   dxp,da0p,dxn,da0n)

        CALL XRO(0.D0)
 
             IF (KLAVE.EQ.0) THEN
              RMAX= MAX(DXP,DXN)+9.25*MAX(DA0P,DA0N)
             ELSEIF (KLAVE.EQ.1) THEN
              RMAX=DSQRT(20.D0)*MAX(DXN,DXP)
             ENDIF




             if (ipol.eq.1) then
c calculo potencial coulomb tamanyo finito

               call dsg20r(0.d0,rmax,mnr,dr,nr)

                alpha = 1.d0/137.036d0
 
                do ir=1,nr
                r = dr(ir)
             
              CALL DSG20R (0.D0,R,5,DRP,NRP)

                  do irp = 1,nrP
                   rp=drp(irp)
                   
                   df0 (irp) = (rp**2)*densq(rp)/r
c                  write(10,*)rp,densq(rp)
                  enddo                   

              CALL DRG20R (0.D0,R,5,df0,f1)

              CALL DSG20R (R,RMAX,5,DRP,NRP)

                  do irp = 1,nrP
                   rp=drp(irp)
                   
                   df0 (irp) = rp*densq(rp)
c                  write(10,*)rp,densq(rp)
                  enddo                  

              CALL DRG20R (R,RMAX,5,df0,f2)

               vcd(ir)=-alpha*ieta*4.d0*dpi*(f1+f2)
c               write(40,*)r,vcd(ir)*hbarc
               enddo

                endif

               cosaux= cost
           theta=acos(cosaux)!60.d0*dpi/180.d0


              eout = dmf+(tmugev*1.d3)/hbarc 

             POUT=DSQRT(EOUT**2-DMF**2)



                 ein = eingev*1d3/hbarc
                 PIN=DSQRT(EIN**2-DMI**2)
                 Q0=EIN-eout-qvalue/hbarc


                 if (q0.ge.0) then                   
             DQ=DSQRT(PIN**2+POUT**2-2.D0*DCOS(THETA)*PIN*POUT)


             call secciondif(q0,dq,theta,ein,pin,eout,pout,ilin,ieta,
     f        rmax,mn,mnr,wd,sig,ipol,rmaxfrac)

             else
                axx=0.
                azz=0.
                aoz=0.
                a00=0.
                axy=0.
             end if

100          FORMAT(A)


        return
        END




           subroutine secciondif(q0,dq,theta,ein,pin,eout,pout,ilin,
     f     ieta,rmax,mn,mnr,wd,sig,ipol,rmaxfrac)
        IMPLICIT REAL*8 (A-H,O-T,V-Z)
        IMPLICIT COMPLEX*16 (U)
 
        dimension wd(5)
        common/datos/dpi,hbarc,GF0,DMNU,DMA
        common/datos2/dmneutrino,dmlepton,dmi,dmf,
     f                coscabibbo
 
             Q2=Q0**2-DQ**2

            factor = pout*eout*dma*(coscabibbo*GF0)**2/(dpi**2)
             cte = dmlepton**2/(eout*(eout+pout))
             sin12 = sin(theta/2)  
             cos12 = dsqrt(1.d0 - sin12**2)
             coseno = cos12**2-sin12**2             

       CALL WSELF(ILIN,Q0,dq,coseno,Q2,RMAX,MN,mnr,WD,ipol,eout,pin,
     f  ieta,rmaxfrac)

             return
             end






         subroutine WSELF(ILIN,Q0,dq,coseno,Q2,RMAX,MN,mnr,WD,ipol,
     f eout,pin,ieta,rmaxfrac)


******* TODAS LAS UNIDADES EN FMS; Q0 Y DQ DEFINEN EL4-MOMENTO DEL FOTON
******* RMAX DEFINE EL RANGO DE INTEGRACION Y MN EL NUMERO DE PTOS.
******* GAUSS. FINALMENTE WL Y WT SON LAS FUNCIONES DE ESTRUCTURA.
*******         IF (ILIN.EQ.1) THEN
*******          'LINHARD NORMAL'
*******         ELSEIF (ILIN.EQ.2) THEN
*******          'LINHARD RELATIVISTA'
*******         ELSEIF (ILIN.EQ.3) THEN
*******          'LINHARD F.SPECTRALES'
*******         ENDIF
*******         Unidades FMs y qvector va en la direccion z
 
              IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
              IMPLICIT COMPLEX*16 (U)
              complex*16 cunuc, delta_lind

              DIMENSION vcd(2000),DR(2000),DFXX(2000),DFZZ(2000)
              DIMENSION DF0Z(2000),DFXY(2000),DF00(2000)
              DIMENSION tulin(0:3),rulin(0:3,0:3),WD(5)

 
              COMMON/VFPI/DROP,DRON
              common/datos/dpi,hbarc,GF0,DMNU,DMA
              common/datos3/xmpi,gp
              common/densidad/dro,dro0
              common/cunuc/cunuc
              common/vcfto/vcd
              common/fsi/xlindfs

        common/datos2/dmneutrino,dmlepton,dmi,dmf,
     f                coscabibbo
        common/importantisimo/xma

       if (ilin.ne.2) stop'no preparado'

              cte_r= -4.d0*dpi/(8.d0*dpi*dmnu**2)



              CALL DSG20R (0.D0,RMAX,MNR,DR,NR)
                
 
c              DO IR=1,NR

              iflag=0
              q0old = q0
              dqold = dq
              q2old = q2
              vc= 0.d0
              qth = 0.d0



c               R=DR(IR)
              R=rmaxfrac*RMAX

              CALL xro(R)


c            correccion del potencial coulombiano

              DKFP=(3.D0*DPI*DPI*DROP)**(1.D0/3.D0)
              DKFN=(3.D0*DPI*DPI*DRON)**(1.D0/3.D0)

              if (ipol.eq.1) then 
                 vc = vcd(ir)
 
              epro = dsqrt(dmnu**2+dkfp**2)              
              eneu = dsqrt(dmnu**2+dkfn**2)              
              qth = ieta*(epro-eneu)
              endif
                 

              q0 = q0+qth

               eoutlocal = eout-vc
               if ((eout-vc-dmf).lt.0) then
                 iflag=1 ! asi vale cero la respuesta
                 poutlocal = 0.d0
               else
               poutlocal = dsqrt(eoutlocal**2-dmf**2)
               endif 
 
               pout = dsqrt(eout**2-dmf**2)
               dq =DSQRT(PIN**2+POUTlocal**2-2.D0*coseno*PIN*POUTlocal) 
               fema=eoutlocal*poutlocal/(eout*pout)

              q2=(q0**2-dq**2)



* CONSTANTS
               
              
c              xmassproton = 938.3d0/hbarc
c              xmassneutron = 939.6d0/hbarc
              UYI=(0.D0,1.D0)
              DMAGP=2.792847D0
              DMAGN=-1.913043D0
              FFAC=0.843D0*1.D3/hbarc
              xln=5.6d0
              if (ipol.eq.2) xln = 0.d0            
              GQ=1.D0/(1.D0-Q2/(FFAC*FFAC))**2
              tau = -q2/4.d0/dmnu**2
cccccccccccccccccc              xma = 1.049*1.D3/hbarc
              ga = 1.257d0
              xmpi = 139.57d0 / hbarc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc             
              f1v = 0.5d0*gq*(1.d0+tau*dmagp
     f         -dmagn*xln*tau**2/(1.d0+xln*tau))/(1.d0+tau)
              xmuf2v =0.5d0*gq*(-1.d0+dmagp
     f         -dmagn*(1.d0+xln*tau+tau)/(1.d0+xln*tau))/(1.d0+tau)


              GAQ = ga /(1.d0-Q2/xma**2)**2 
              GPQ = 2.d0*dmnu*gaq/(xmpi**2-q2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             if (ipol.eq.1) then
                C0POL=380.D0/hbarc
                fp0in= 0.33d0
                fp0ex=0.45d0 
                gp=0.63d0
                VT=VTP(Q0,DQ)
                VL=VLP(Q0,DQ)
             endif
****






             if (ieta.eq.1) then
              xlind=DUR(q0,DQ,DKFn,dkfp,DMNU,t0,r00)
               elseif (ieta.eq.-1) then
              xlind=DUR(q0,DQ,DKFp,dkfn,DMNU,t0,r00)
              endif

 
               if (ipol.eq.1) then
                    if (iflag.eq.1) then
                 xlind=0.d0
                 t0 = 0.d0
                 r00 = 0.d0
c                 write(6,*)'vc juega un papel', q0*hbarc
                 endif
               endif
            tulin(0) = t0
            tulin(1) = 0.d0
            tulin(2) = 0.d0
            tulin(3) = (0.5d0*q2 + q0*tulin(0))/dq 

            rulin(0,0) = r00
            rulin(0,1) = 0.d0
            rulin(0,2) = 0.d0
            rulin(0,3) = (0.5d0*q2*tulin(0) + q0*rulin(0,0))/dq 
 
            rulin(1,0) = rulin(0,1)
            rulin(2,0) = rulin(0,2)
            rulin(3,0) = rulin(0,3)
            trazaerre =r00-dmnu**2
            errezz = (q2**2+4.d0*r00*q0**2+4.d0*q2*q0*t0)/(4.d0*dq**2)
            gamma = (trazaerre-errezz)/2.d0
            delta = (-trazaerre+3.d0*errezz)/2.d0/dq**2
            rulin(1,1) = gamma
            rulin(1,2) = 0.d0
            rulin(1,3) = 0.d0
            rulin(2,2) = gamma
            rulin(2,3) = 0.d0
            rulin(3,3) = gamma+delta*dq**2

              do i = 2,3
               do j=1,i-1
                rulin(i,j) = rulin(j,i)
               enddo
              enddo  



               if (ipol.eq.1) then
                if (iflag.eq.1) then
                 xlind=0.d0
                 t0 = 0.d0
                 r00 = 0.d0
c                 write(6,*)'vc juega un papel', q0*hbarc
                 endif
               endif  


              fact=1.d0
              facl=1.d0
              f00=1.d0

            if (ipol.eq.1) then

              DKF=(1.5D0*DPI*DPI*DRO)**(1.D0/3.D0)


              irel= 1 ! relativista  


              if (abs(xlind).gt.1.d-10) then  

              unuc_j = ulinrel(q0,dq,dkf)

              u_zero=q0+uyi*0.0d0
              udel=delta_lind(u_zero,dq,dro,dkf)

              else ! no es que valga cero, pero como xlind es cero da igual
               unuc_j = 0.d0+uyi*0.d0
               udel = 0.d0 +uyi*0.d0
              endif

              cunuc= unuc_j
              ufuli = unuc_j+udel          
             DELTA=dreal(ufuli-cunuc)





          FACT=abs(1.D0/(1.d0-VT*UFULI))**2
          FACL=abs(1.D0/(1.d0-VL*UFULI))**2
          delfact=(1.d0-vt*delta)**2
          delfacl=(1.d0-vl*delta)**2
c          fact=fact*delfact  
c          facl=facl*delfacl  


          FPRIMA0=(DRO/dro0*fp0in+(1.D0-DRO/dro0)*fp0ex)*C0POL
          F00=ABS(1.D0/(1.D0-FPRIMA0*CUNUC))**2 

 
           endif


          call amunu_pol(axx,azz,a0z,a00,axy,tulin,rulin,q2,q0,dq,f1v,
     f                  xmuf2v,gaq,gpq,fact,facl,f00,ipol)

          if (abs(xlind).gt.1.d-10) then  
             hbarc = 197.3269602d0
             write(61,15) -q2*hbarc**2*1d-6,q0*hbarc*1d-3,dq*hbarc*1d-3,
     f            axx*hbarc**2*1d-6,azz*hbarc**2*1d-6,a0z*hbarc**2*1d-6,
     f            a00*hbarc**2*1d-6,axy*hbarc**2*1d-6,
     f            fact,facl,f00,(eout-dmf)*hbarc*1d-3,
     f            xlind*hbarc**2*1d-6,
     f            f1v,xmuf2v,gaq,gpq*1d3/hbarc
          end if
 15       format(e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,
     f         e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,e15.7,e15.7)


            RETURN
            END
 

          subroutine amunu_pol(axx,azz,a0z,a00,axy,tulin,rulin,q2,q0
     f    ,dq,f1v,xmuf2v,gaq,gpq,fact,facl,f00,ipol)


              IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
              DIMENSION tulin(0:3),rulin(0:3,0:3)
              common/datos/dpi,hbarc,GF0,DMNU,DMA


              axx = 16.d0*f1v**2*(2.d0*rulin(1,1)-q2/2.d0)+ 
     f    2.d0*q2*xmuf2v**2*(-4.d0*fact-4.d0*rulin(1,1)/dmnu**2) +
     f    4.d0*gaq**2*(2.d0*rulin(1,1)-(q2/2.d0-2.d0*fact*dmnu**2))-
     f    16.d0*f1v*xmuf2v*fact*q2

              azz = 16.d0*f1v**2*(2.d0*rulin(3,3)+2.d0*dq*tulin(3)
     f                            -q2/2.d0)+ 
     f    2.d0*q2*xmuf2v**2*(-4.d0-4.d0*rulin(3,3)/dmnu**2-
     f    4.d0*dq*tulin(3)/dmnu**2-dq**2*(4.d0/q2+1.d0/dmnu**2)) +
     f    4.d0*gaq**2*(2.d0*rulin(3,3)+2.d0*dq*tulin(3)
     f                 -(q2/2.d0-2.d0*facl*dmnu**2))-
     f    (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*dq**2-
     f    16.d0*f1v*xmuf2v*(q2+dq**2)


              a0z = 16.d0*f1v**2*((2.d0*rulin(0,3)+tulin(0)*dq)*f00
     f                          +tulin(3)*q0)+
     f    2.d0*q2*xmuf2v**2*(-4.d0*rulin(0,3)/dmnu**2-
     f    2.d0*(dq*tulin(0)+q0*tulin(3))/dmnu**2
     f            -dq*q0*(4.d0/q2+1.d0/dmnu**2)) +
     f    4.d0*gaq**2*((2.d0*rulin(0,3)+dq*tulin(0))*facl+q0*tulin(3))-
     f    (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*dq*q0-
cc     f    16.d0*f1v*xmuf2v*f00*dq*q0
     f    16.d0*f1v*xmuf2v*dq*q0

cc              axy = 16.d0*gaq*(xmuf2v+f1v)*fact*(-dq*tulin(0)
cc     f         +q0*tulin(3))

              axy = 16.d0*gaq*(xmuf2v+f1v)*(-dq*tulin(0)*fact
     f         +q0*tulin(3))

              a00 = 16.d0*f1v**2*(2.d0*rulin(0,0)*f00+2.d0*q0*tulin(0)
     f                            +q2/2.d0)+ 
     f    2.d0*q2*xmuf2v**2*(4.d0-4.d0*rulin(0,0)/dmnu**2-
     f    4.d0*q0*tulin(0)/dmnu**2-q0**2*(4.d0/q2+1.d0/dmnu**2)) +
     f    4.d0*gaq**2*(2.d0*rulin(0,0)+2.d0*q0*tulin(0)
     f                 +(q2/2.d0-2.d0*dmnu**2))-
     f    (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*q0**2-
     f    16.d0*f1v*xmuf2v*(-q2+q0**2)*f00



             return
             end
                       




        SUBROUTINE XRO(DR)
        IMPLICIT REAL*8 (A-H,O-T,V-Z)
        IMPLICIT COMPLEX*16 (U)
        DIMENSION DID(2000),DFD(2000)
        COMMON/VFPI/DROP,DRON
        COMMON/ROFEROS/KLAVE
        COMMON/NXRO/NXRO
        COMMON/NUC/DZZ,DAA
        common/densidad/dro,dro0
C
C DATOS PARA LA FUNCION DENSIDAD DE FERMI DE PROTONES Y NEUTRONES
        COMMON/CTEROPN/DXP,DA0P,DRO0P,DXN,DA0N,DRO0N
C DROPR... densidad protones /dro0p
C DRONR....densidad neutrones /dro0n
C KLAVE=0 densidad de fermi . KLAVE=1 densidad de oscilador
        DATA DPI/3.141592653589793D0/
 
        IF(NXRO.EQ.0) THEN
         IF (KLAVE.EQ.0) THEN
C Calculamos las ctes. de normalizacion de las densidades de neutrones y
C  protones normalizadas a DZZ y DAA-DZZ, respectivamente
        	DRMAX1P=5.D0*DXP
        	CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
        	 DO I=1,NPP
        	DI=DID(I)
        	DFD(I)=DI*DI/(DEXP((DI-DXP)/DA0P)+1.D0)
        	 END DO
        	CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
        	DRO0P=DZZ/(DRESP*4.D0*DPI)
        	DRMAX1N=5.D0*DXN
        	CALL DSG20R(0.D0,DRMAX1N,10,DID,NPP)
        	 DO I=1,NPP
        	DI=DID(I)
        	DFD(I)=DI*DI/(DEXP((DI-DXN)/DA0N)+1.D0)
        	 END DO
        	CALL DRG20R(0.D0,DRMAX1N,10,DFD,DRESN)
        	DRO0N=(DAA-DZZ)/(DRESN*4.D0*DPI)
        ELSE ! oscilador
        	 DRMAX1P=8.D0*DXP
        	CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
        	 DO I=1,NPP
        	DI=DID(I)
        	DFD(I)=DI*DI*(1.D0+DA0P*(DI/DXP)**2)*DEXP(-(DI/DXP)**2)
        	 END DO
        	CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
        	DRO0P=DZZ/(DRESP*4.D0*DPI)
        	DRMAX1N=8.D0*DXN
        	CALL DSG20R(0.D0,DRMAX1N,10,DID,NPP)
        	 DO I=1,NPP
        	DI=DID(I)
        	DFD(I)=DI*DI*(1.D0+DA0N*(DI/DXN)**2)*DEXP(-(DI/DXN)**2)
        	 END DO
        	CALL DRG20R(0.D0,DRMAX1N,10,DFD,DRESN)
        	DRO0N=(DAA-DZZ)/(DRESN*4.D0*DPI)
        ENDIF
 
        NXRO=1
 
        END IF
 
        IF(DR.GT.5.D0*DXN.AND.KLAVE.EQ.0)THEN
        DROPR=0.D0
        DRONR=0.D0
        DROP=0.D0
        DRON=0.D0
        dro  = drop+dron
        dro0 = dro0p+dro0n

        RETURN
        END IF
 
        IF(DR.GT.8.D0*DXN.AND.KLAVE.EQ.1)THEN
        DROPR=0.D0
        DRONR=0.D0
        DROP=0.D0
        DRON=0.D0
        dro  = drop+dron
        dro0 = dro0p+dro0n

        RETURN
        END IF
 
 
                IF(KLAVE.EQ.0)THEN
C densidad reducida de fermi
        DROPR=1.D0/(DEXP((DR-DXP)/DA0P)+1.D0)
        DRONR=1.D0/(DEXP((DR-DXN)/DA0N)+1.D0)
        DROP=DRO0P*DROPR
        DRON=DRO0N*DRONR
 
                ELSE
C densidad reducida oscilador
        DROPR=(1.D0+DA0P*(DR/DXP)**2)*DEXP(-(DR/DXP)**2)
        DRONR=(1.D0+DA0N*(DR/DXN)**2)*DEXP(-(DR/DXN)**2)
        DROP=DRO0P*DROPR
        DRON=DRO0N*DRONR
                END IF
        dro  = drop+dron
        dro0 = dro0p+dro0n

        RETURN
        END
 



        double precision function densq(DR)
        IMPLICIT REAL*8 (A-H,O-T,V-Z)
        IMPLICIT COMPLEX*16 (U)
        DIMENSION DID(2000),DFD(2000)

        COMMON/ROFEROS/KLAVE
        COMMON/NCXRO/NCXRO
        COMMON/CTEROPNC/DXP,DA0P,DRO0P
        COMMON/NUC/DZZ,DAA

C DROPR... densidad protones /dro0p
C KLAVE=0 densidad de fermi . KLAVE=1 densidad de oscilador
        DATA DPI/3.141592653589793D0/
 
        IF(NCXRO.EQ.0) THEN
         IF (KLAVE.EQ.0) THEN
C Calculamos las ctes. de normalizacion de las densidades de 
c protones normalizada a DZZ.
        DRMAX1P=5.D0*DXP
        CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
         DO I=1,NPP
        DI=DID(I)
        DFD(I)=DI*DI/(DEXP((DI-DXP)/DA0P)+1.D0)
         END DO
        CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
        DRO0P=DZZ/(DRESP*4.D0*DPI)

        ELSE ! oscilador
         DRMAX1P=8.D0*DXP
        CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
         DO I=1,NPP
        DI=DID(I)
        DFD(I)=DI*DI*(1.D0+DA0P*(DI/DXP)**2)*DEXP(-(DI/DXP)**2)
         END DO
        CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
        DRO0P=DZZ/(DRESP*4.D0*DPI)

        ENDIF
 
        NCXRO=1
 
        END IF
 
        IF(DR.GT.5.D0*DXP.AND.KLAVE.EQ.0)THEN
        densq=0.D0
        RETURN
        END IF
 
        IF(DR.GT.8.D0*DXP.AND.KLAVE.EQ.1)THEN
        densq=0.D0
        RETURN
        END IF
 
 
                IF(KLAVE.EQ.0)THEN
C densidad reducida de fermi
        DROPR=1.D0/(DEXP((DR-DXP)/DA0P)+1.D0)
        densq=DRO0P*DROPR

 
                ELSE
C densidad reducida oscilador
        DROPR=(1.D0+DA0P*(DR/DXP)**2)*DEXP(-(DR/DXP)**2)
        densq=DRO0P*DROPR

                END IF

        RETURN
        END
 



         subroutine convolucion (dncxp,dnca0p,dncxn,dnca0n,klave,lc,
     f   dxp,da0p,dxn,da0n)
         IMPLICIT REAL*8 (A-H,O-T,V-Z)
         DATA DPI/3.141592653589793D0/

         PI2 = DPI*DPI
* DESCONVOLUCION DE LOS PARAMETROS
*
      if (lc.eq.0) then
 
        IF (KLAVE.EQ.0) THEN
         DXP=DNCXP+3.45D0*DNCXP/(15.*DNCXP**2+7*PI2*DNCA0P**2)
         DA0P=DSQRT((DNCXP**3+PI2*DXP*DNCA0P**2-DXP**3)/(DXP*PI2))
         DXN=DNCXN+3.45D0*DNCXN/(15.*DNCXN**2+7*PI2*DNCA0N**2)
         DA0N=DSQRT((DNCXN**3+PI2*DXN*DNCA0N**2-DXN**3)/(DXN*PI2))
        ELSE
         DXPP=0.46D0
* DXPP ES IGUAL A 2/3 DEL RADIO CUADRATICO MEDIO DEL PROTON EN FM2
         DXP=DSQRT(DNCXP**2-DXPP)
 
         DXXP=DNCA0P*(DNCXP**2)/((1.D0+1.5D0*DNCA0P)*DXP**2)
 
         DA0P=2.D0*DXXP/(2.D0-3.D0*DXXP)
 
         DXN=DSQRT(DNCXN**2-DXPP)
 
         DXXN=DNCA0N*(DNCXN**2)/((1.D0+1.5D0*DNCA0N)*DXN**2)
 
         DA0N=2.D0*DXXN/(2.D0-3.D0*DXXN)
 
 
        ENDIF


        else !lc.ne.0

        dxp=dncxp
        da0p=dnca0p
        dxn=dncxn
        da0n=dnca0n

        endif
      
        return
        end




      FUNCTION DUR (Q0,DQ,PF1,PF2,DMNU,t0,r00)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)


      IF (Q0.LT.0) THEN
c       WRITE(6,*)'EN LA LLAMADA A DUR,Q0 ES NEGATIVA'
       DUR=0.
       t0=0.d0
       r00 =0.d0 
       return
      ENDIF
 
      DPI=3.141592653589793D0
      CTE=-0.5D0*DMNU*DMNU/DPI/DQ
      Q2=Q0**2-DQ**2
      IF (Q2.gE.0) THEN
       DUR=0.
       t0=0.d0
       r00 =0.d0 
       RETURN
      endif
      EF1=DSQRT(DMNU**2+PF1**2)
      EF2=DSQRT(DMNU**2+PF2**2)
      IF ((EF2-Q0).GT.EF1) THEN
       DUR=0.
       t0=0.d0
       r00 =0.d0 
       RETURN
      ENDIF
      a=0.5d0*(-q0+dq*sqrt(1.d0-4.d0*dmnu**2/q2))
   
      IF (a.GT.EF1) THEN
       DUR=0.
       t0=0.d0
       r00 =0.d0
       RETURN
      ENDIF
         
      EINF=MAX(DMNU,EF2-Q0,a)

      t0  = 0.5d0*(ef1+einf)
      r00 = (ef1**2+einf**2+ef1*einf)/3.d0 
      DUR=(ef1-einf)*CTE
      RETURN
      END


 
      FUNCTION VTP (A0,A1)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      common/datos/dpi,hbarc,GF0,DMNU,DMA
      common/datos3/xmpi,gp
       
      DFM2=0.08*4.D0*DPI/(XMPI*XMPI)
      A2=A0**2-A1**2
      XL2=(2500.d0/hbarc)**2
      CRO=2.0
      DMRO2=(770.d0/hbarc)**2
      AFORMA=CRO*((XL2-DMRO2)/(XL2-A2))**2
      VTP=DFM2*(AFORMA*A1**2/(A2-DMRO2)+GP)
      RETURN
      END


      FUNCTION VLP (A0,A1)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      common/datos/dpi,hbarc,GF0,DMNU,DMA
      common/datos3/xmpi,gp
      
      DMPI2=xmpi*xmpi 
      DFM2=0.08*4.D0*DPI/DMPI2
      A2=A0**2-A1**2
      XL2=(1200.d0/hbarc)**2
      
      AFORMA=((XL2-DMPI2)/(XL2-A2))**2
      VLP=DFM2*(AFORMA*A1**2/(A2-DMPI2)+GP)
      RETURN
      END


*********************************************************************
*      FUNCIONES DE LINHARD
**********************************************************************
       SUBROUTINE ULIND(QZR,Q,XKF,CUFUN)
c	subroutine contained obsolete if statements that would not c	compile
       RETURN
       END





        function delta_lind (q_zero,q_mod,rho,k_fermi)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   complex Lindhard function for symmetric nuclear matter:
c                    from Appendix of
c                    E.Oset et al Phys. Rept. 188:79, 1990
c                    formula A.4 
c
c   input variables: 
c     q_zero [fm^-1] : Energy
c     q_mod  [fm^-1] : Momentum
c     rho    [fm^3]  : Nuclear density
c     k_fermi[fm^-1] : Fermi momentum
c
c   All variables are real*8
c
c   output variable: 
c     delta_lind [fm^-2]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ATTENTION!!!
c Only works properly for real q_zero,
c if q_zero has an imaginary part calculates the L. function
c assuming Gamma= 0.
c Therefore this subroutine provides two different functions
c depending on whether q_zero is real or not!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        implicit none
        real*8 q_mod,rho,k_fermi,m,rq_zero,gamma,s,srot,wr
        real*8 fdel_f,mpi,gammap,sp,srotp,qcm,qcmp,md,pi
        complex*16 delta_lind,q_zero,z,zp,i,pzeta,pzetap
c
c m = 939/197.3,md = 1232./179.3, mpi = 139./197.3 
c
        data md,m,i/6.2433d0,4.7592d0,(0.d0,1.d0)/
        data mpi,pi/0.7045d0,3.14159265358979323846264338d0/
c
c  f*/f = 2.13 --> f*^2/4pi = .36
c  NOTE: fdel = fdel_f because f =1.0
c

        fdel_f = 2.13d0
        
        wr = md-m
        
        
        gamma = 0.d0 
        gammap = 0.d0
        if(imag(q_zero)**2.lt.1.d-36)then
          rq_zero = dreal(q_zero)
          s = m**2+rq_zero**2-q_mod**2+
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(s.gt.(m+mpi)**2)then
           srot = sqrt(s)
           qcm = sqrt(s**2+mpi**4+m**4-2.d0*(s*mpi**2+s*m**2+
     &    (mpi*m)**2)) /(2.d0*srot)
           gamma = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcm**3/srot*(m+sqrt(m**2+qcm**2))/mpi**2  
         endif          
          sp = m**2+rq_zero**2-q_mod**2-
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(sp.gt.(m+mpi)**2)then
           srotp = sqrt(sp)
           qcmp=sqrt(sp**2+mpi**4+m**4-2.d0*(sp*mpi**2+sp*m**2+
     &                    (mpi*m)**2))/(2.d0*srotp)
           gammap = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcmp**3/srotp*(m+sqrt(m**2+qcmp**2))/mpi**2  
          endif          
        endif
              
        z=md/(q_mod*k_fermi)*(q_zero-q_mod**2/(2.d0*md)
     &                         -wr +i*gamma/2.d0)

        zp=md/(q_mod*k_fermi)*(-q_zero-q_mod**2/(2.d0*md)
     &                          -wr +i*gammap/2.d0)
c
c care with limit cases
c
        if(abs(z).gt.50.d0)then
            pzeta =  2.d0/(3.d0*z) +2.d0/(15.d0*z**3)
        else if(abs(z).lt.1.d-2)then
            pzeta =  2.d0*z -2.d0/3.d0*z**3 -i*pi/2.d0*(1.d0 - z**2) 
        else
            pzeta =  z +  (1.d0-z**2) * log((z+1.d0)/(z-1.d0))/2.d0
        endif
        
        if(abs(zp).gt.50.d0)then
            pzetap =  2.d0/(3.d0*zp) +2.d0/(15.d0*zp**3)
        else if(abs(zp).lt.1.d-2)then
            pzetap =  2.d0*zp -2.d0/3.d0*zp**3 -i*pi/2.d0*(1.d0 - zp**2) 
        else
            pzetap =  zp + (1.d0-zp**2) * log((zp+1.d0)/(zp-1.d0))/2.d0
        endif


       delta_lind = 2.d0/3.d0 * rho * md/(q_mod*k_fermi) * ( 
     &        pzeta +pzetap) * fdel_f **2    
     
        return
        end
	
****************************************************************************
C    *******************************************************************GAU00010
C               SUBROUTINAS DE INTEGRACION NUMERICA POR GAUSS           GAU00020
C    *******************************************************************GAU00030
      SUBROUTINE SG20R(A,B,N,X,NP)                                      GAU00040
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               GAU00050
      DIMENSION Y(10),X(2000)                                           GAU00060
      DATA Y/.9931285991,.9639719272,.9122344282,.8391169718,           GAU00070
     F .7463319064,.6360536807,.5108670019,.3737060887,                 GAU00080
     F .2277858511,.0765265211/                                         GAU00090
      NP=20*N                                                           GAU00100
      DINT=(B-A)/FLOAT(N)                                               GAU00110
      DELT=DINT*0.5                                                     GAU00120
      ORIG=A-DELT                                                       GAU00130
      I1=-20                                                            GAU00140
      DO 1 I=1,N                                                        GAU00150
      ORIG=ORIG+DINT                                                    GAU00160
      DORIG=ORIG+ORIG                                                   GAU00170
      I1=I1+20                                                          GAU00180
      I2=I1+21                                                          GAU00190
      DO 2 J=1,10                                                       GAU00200
      J1=I1+J                                                           GAU00210
      J2=I2-J                                                           GAU00220
      X(J1)=ORIG-DELT*Y(J)                                              GAU00230
 2    X(J2)=DORIG-X(J1)                                                 GAU00240
 1    CONTINUE                                                          GAU00250
      RETURN                                                            GAU00260
      END                                                               GAU00270
      SUBROUTINE RG20R(A,B,N,CF,CRES)                                   GAU00280
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               GAU00290
      DIMENSION W(10),CF(2000)                                          GAU00300
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00310
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00320
     F .1527533871/                                                     GAU00330
      CR=(0.,0.)                                                        GAU00340
      I1=-20                                                            GAU00350
      DO 1 I=1,N                                                        GAU00360
      I1=I1+20                                                          GAU00370
      I2=I1+21                                                          GAU00380
      DO 2 J=1,10                                                       GAU00390
      J1=I1+J                                                           GAU00400
      J2=I2-J                                                           GAU00410
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00420
 1    CONTINUE                                                          GAU00430
      CRES=CR*0.5*(B-A)/FLOAT(N)                                        GAU00440
      RETURN                                                            GAU00450
      END                                                               GAU00460
                                                                        GAU00470
                                                                        GAU0048 0
C    *******************************************************************GAU00490
C               SUBROUTINAS DE INTEGRACION NUMERICA POR GAUSS           GAU00500
C    *******************************************************************GAU00510
      SUBROUTINE DSG20R(A,B,N,X,NP)                                     GAU00520
      IMPLICIT REAL*8 (A-H,O-Z)                                         GAU00530
      DIMENSION Y(10),X(2000)                                           GAU00540
      DATA Y/.9931285991,.9639719272,.9122344282,.8391169718,           GAU00550
     F .7463319064,.6360536807,.5108670019,.3737060887,                 GAU00560
     F .2277858511,.0765265211/                                         GAU00570
      NP=20*N                                                           GAU00580
      DINT=(B-A)/DBLE(N)                                                GAU00590
      DELT=DINT*0.5D0                                                   GAU00600
      ORIG=A-DELT                                                       GAU00610
      I1=-20                                                            GAU00620
      DO 1 I=1,N                                                        GAU00630
      ORIG=ORIG+DINT                                                    GAU00640
      DORIG=ORIG+ORIG                                                   GAU00650
      I1=I1+20                                                          GAU00660
      I2=I1+21                                                          GAU00670
      DO 2 J=1,10                                                       GAU00680
      J1=I1+J                                                           GAU00690
      J2=I2-J                                                           GAU00700
      X(J1)=ORIG-DELT*Y(J)                                              GAU00710
 2    X(J2)=DORIG-X(J1)                                                 GAU00720
 1    CONTINUE                                                          GAU00730
      RETURN                                                            GAU00740
      END                                                               GAU00750
      SUBROUTINE DRG20R(A,B,N,CF,CRES)                                  GAU00760
      IMPLICIT REAL*8 (A-H,O-Z)                                         GAU00770
      DIMENSION W(10),CF(2000)                                          GAU00780
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00790
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00800
     F .1527533871/                                                     GAU00810
      CR=(0.D0,0.D0)                                                    GAU00820
      I1=-20                                                            GAU00830
      DO 1 I=1,N                                                        GAU00840
      I1=I1+20                                                          GAU00850
      I2=I1+21                                                          GAU00860
      DO 2 J=1,10                                                       GAU00870
      J1=I1+J                                                           GAU00880
      J2=I2-J                                                           GAU00890
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00900
 1    CONTINUE                                                          GAU00910
      CRES=CR*0.5D0*(B-A)/DBLE(N)                                       GAU00920
      RETURN                                                            GAU00930
      END                                                               GAU00940


      SUBROUTINE CRG20R(A,B,N,CF,CRES)                                  GAU00760
      IMPLICIT REAL*8 (A-H,O-Z)
      complex*16 CF,cres,CR 
      DIMENSION W(10),CF(2000)                                          GAU00780
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00790
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00800
     F .1527533871/                                                     GAU00810
      CR=(0.D0,0.D0)                                                    GAU00820
      I1=-20                                                            GAU00830
      DO 1 I=1,N                                                        GAU00840
      I1=I1+20                                                          GAU00850
      I2=I1+21                                                          GAU00860
      DO 2 J=1,10                                                       GAU00870
      J1=I1+J                                                           GAU00880
      J2=I2-J                                                           GAU00890
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00900
 1    CONTINUE                                                          GAU00910
      CRES=CR*0.5D0*(B-A)/DBLE(N)                                       GAU00920
      RETURN                                                            GAU00930
      END                                                               GAU00940
****************************************************************************
c          PROGRAM TEST
c             IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
c             IMPLICIT COMPLEX*16 (U)
c   1    write(6,*)'kf,q0'
c        read(5,*)dkf,q0

c        do dq=q0+0.01D0,10.d0,.1d0
c        irel=0
c        call unucjuan(irel,Q0,DQ,DKF,unuc0)
c
c        irel=1
c        call unucjuan(irel,Q0,DQ,DKF,unuc1)
c
c        irel=1
c        call unucjuanb(irel,Q0,DQ,DKF,unucx)
c        
c        unuc2=ulinrel(q0,dq,dkf)
c  
c        write(27,111)dq,unuc0,unuc1,unucx,unuc2

c        enddo
         
c  111    format(11g15.4)
c        end







c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NUCLEON relativistic Lindhard Function
c Same normalization as ULIN
c Real part
c taken from Eur.Phys.J.A25:299-318,2005 (Barbaro et al)
c Eq. 61
c
c Im. part: Juan. 
c
c INPUT: Real*8 
c  q0:Energy   [fm]
c  qm: modulus 3mom [fm]
c  kf: Fermi mom [fm]
c
c OUTPUT: Complex*16 [fm]
c
c USES: RULINRELX, DUR_J
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          COMPLEX*16 FUNCTION ULINREL(q0,qm,kf)
          IMPLICIT REAL*8 (A-T,V-Z)
          IMPLICIT COMPLEX*16 (U)
c range validity
          if(q0.gt.qm)then
           write(6,*)' Rlinrel is not valid for q0 > q'
           write(6,*)' q0=',q0,' q=',qm
           stop
          endif
          m = 940.d0/197.3269602d0
          uyi=(0.d0,1.d0)
c
          RLINREL=RulinrelX(q0,qm,kf)+RulinrelX(-q0,qm,kf)
          ximag= 2.d0*DUR_j(q0,qm,kf,kf,m)
          ULINREL=RLINREL+uyi*ximag 

          RETURN
          END
C
          REAL*8 function RulinrelX(q0,qm,kf)
c fm, only real part, only |q0|<|q|
          IMPLICIT REAL*8 (A-T,V-Z)
          IMPLICIT COMPLEX*16 (U)
c
          m = 940.d0/197.3269602d0
          pi=3.1415926535d0
          pi2=pi**2
          uy=(0.d0,1.d0)
c
       ef=sqrt(m**2+kf**2)
       q2=q0**2-qm**2
       ds=Sqrt(1.d0-(4.d0*m**2)/q2)
c
       L1=Log((kf + ef)/m)
       uL2=Log(Abs((ef + q0 - Sqrt(m**2+(kf-qm)**2))/
     &             (ef + q0 - Sqrt(m**2 + (kf + qm)**2)))) + 
     &  Log(Abs((ef + q0 + Sqrt(m**2 + (kf - qm)**2))/
     &          (ef + q0 + Sqrt(m**2 + (kf + qm)**2))))

       uL3=Log(Abs(((2*kf + q0*ds)**2-qm**2)/
     &             ((2*kf - q0*ds)**2-qm**2))) + 
     &  Log(Abs(((kf-ef*ds)**2 - (4*m**4*qm**2)/q2**2)/
     &          ((kf+ef*ds)**2 - (4*m**4*qm**2)/q2**2)))
c
       RlinrelX = -L1/(16.d0*pi2)+uL2*(2.d0*ef+q0)/(32.d0*pi2*qm)
     &    -uL3*ds/(64.d0*Pi2)
c factor definition
       RulinrelX=RlinrelX*16.*m**2
       return
       end
c
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
c  JUAN REL.+norel LINDHARD, SLOW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             subroutine unucjuanb(irel,Q0,DQ,DKF,unuc_jb)               
              IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
              IMPLICIT COMPLEX*16 (U)
              dimension pd(4000), xmud(4000), xmud2(4000)
              dimension fmu(4000), fp(4000)
              dimension fmu2(4000)
c irel=0 no relativista
c irel=1 relativista 
         ener(irel,dmnu,pmon2)=dfloat(1-irel)*(dmnu+0.5d0*pmon2/dmnu) +
     f       dfloat(irel)*dsqrt(dmnu**2+pmon2)
       

             if (irel.ne.1 .and. irel.ne.0) stop'wrong irel'

              call dsg20r(0.d0, dkf, 20, pd,np)
              call dsg20r(-1.d0, 1.d0, 20, xmud2,nmu2)

              dkf2= dkf**2
              dq2=dq**2 
              uyi= (0.d0,1.d0)
              dmnu = 940.d0 / 197.3269602d0
              epsilon = (0.5d0/ 197.3269602d0)*dq*197.3269602d0/650.
              pi2= 3.14159265d0**2


              xnuc_j= 0.d0

c contrib. q y -q

              do ij=1,2


              do ip=1,np
               p= pd(ip)
               p2= p**2
               ep=  ener(irel,dmnu,p2)              
               xmu0= (dkf2-p2-dq2)/(2.d0*p*dq)


               if (xmu0.le.1.d0) then
  
                xmu_min= max(-1.d0, xmu0)

               call dsg20r(xmu_min, 1.d0, 20, xmud,nmu)

                 do imu=1,nmu
                  xmu = xmud(imu)

                  pq2=p2+dq2+2.d0*p*dq*xmu  


                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep-epq+uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep-epq+uyi*epsilon
                    
                   xreal_uaux= dreal(1.d0/uaux)

                  if (irel.eq.1) fmu(imu) =xreal_uaux*dmnu/epq
                  if (irel.eq.0) fmu(imu) =xreal_uaux
c enddo xmu
                 enddo


                call DRG20R(xmu_min, 1.d0, 20,fmu,resmu)

                else
                  resmu= 0.d0
               endif


                   if (irel.eq.1) then !  nueva contrib
       if (q0.gt.2.d0*dmnu)stop'contribuc imag, de antipart'

                 do imu=1,nmu2
                  xmubis = xmud2(imu)

                  pq2=p2+dq2+2.d0*p*dq*xmubis  
                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep+epq-uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep+epq-uyi*epsilon
                    
                   xreal_uaux= dreal(1.d0/uaux)

                   fmu2(imu) =xreal_uaux*dmnu/epq

c enddo xmubis
                 enddo


                call DRG20R(-1.d0, 1.d0, 20,fmu2,resmu2)

c nueva contri
                  endif                

               if (irel.eq.1) fp(ip) = (resmu-resmu2)*p2*dmnu/ep 
               if (irel.eq.0) fp(ip) = resmu*p2
c enddo p
              enddo  

              call DRG20R(0.d0, dkf, 20,fp,resp)


              xnuc_j= xnuc_j+resp/pi2

c enddo ij
               enddo
 

              if (irel.eq.1) ximag= 2.d0*DUR_j(q0,DQ,DKF,dkf,DMNU)
              if (irel.eq.0) ximag= 2.d0*DUIM_j(q0,DQ,DKF,dkf,DMNU)

               unuc_jb=xnuc_j+uyi*ximag 

               return
               end
c
c
      FUNCTION DUIM_j(W,Q,K1,K2,DMNU)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z)
      DOUBLE PRECISION KF,L2,K1,K2,NUM
 
      PI=3.141592653589793D0
      L2=(DMNU*W/Q-Q/2.d0)**2
      NUM=K2**2-2*DMNU*W
        ef1 = dmnu + k1**2 / (2.d0*dmnu)         
      IF((L2.LE.NUM).AND.(NUM.LE.K1**2))THEN
        DUIM_j=-DMNU/(4.d0*PI*Q)*(K1**2-K2**2+2.d0*DMNU*W)
      ELSEIF((NUM.LE.L2).AND.(L2.LT.K1**2))THEN
        DUIM_j=-DMNU/(4.d0*PI*Q)*(K1**2-L2)
      ELSE
        DUIM_j=0.d0
       ENDIF
      RETURN
      END
c
c
      FUNCTION DUR_j (Q0,DQ,PF1,PF2,DMNU)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION DE(2000),DFE(2000)

c      IF (Q0.LT.0) THEN
c       WRITE(6,*)'EN LA LLAMADA A DUR,Q0 ES NEGATIVA'
c       DUR=0.
c       return
c      ENDIF
 
      DPI=3.141592653589793D0
      CTE=-0.5D0*DMNU*DMNU/DPI/DQ
      Q2=Q0**2-DQ**2
      IF (Q2.gE.0) THEN
       DUR_j=0.
       RETURN
      endif
      EF1=DSQRT(DMNU**2+PF1**2)
      EF2=DSQRT(DMNU**2+PF2**2)
      IF ((EF2-Q0).GT.EF1) THEN
       DUR_j=0.
       RETURN
      ENDIF
      a=0.5d0*(-q0+dq*sqrt(1.d0-4.d0*dmnu**2/q2))
      IF (a.GT.EF1) THEN
       DUR_j=0.
       RETURN
      ENDIF
         
      EINF=MAX(DMNU,EF2-Q0,a)
      DUR_j=(ef1-einf)*CTE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccc
c Juan Rel. anterior / additional approximation
cccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine unucjuan(irel,Q0,DQ,DKF,unuc_j)
       
              IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
              IMPLICIT COMPLEX*16 (U)
              dimension pd(4000), xmud(4000)
              dimension fmu(4000), fp(4000)
c irel=0 no relativista
c irel=1 relativista 
         ener(irel,dmnu,pmon2)=dfloat(1-irel)*(dmnu+0.5d0*pmon2/dmnu) +
     f       dfloat(irel)*dsqrt(dmnu**2+pmon2)
       

             if (irel.ne.1 .and. irel.ne.0) stop'wrong irel'
              call dsg20r(0.d0, dkf, 20, pd,np)
              dkf2= dkf**2
              dq2=dq**2 
              uyi= (0.d0,1.d0)
              dmnu = 940.d0 / 197.3269602d0
              epsilon = (0.5d0/ 197.3269602d0)*dq*197.3269602d0/650.
              pi2= 3.14159265d0**2

             xnuc_j= 0.d0
c contrib. q y -q
              do ij=1,2
              do ip=1,np
               p= pd(ip)
               p2= p**2
               xmu0= (dkf2-p2-dq2)/(2.d0*p*dq)

               if (xmu0.le.1.d0) then  
                xmu_min= max(-1.d0, xmu0)

               call dsg20r(xmu_min, 1.d0, 20, xmud,nmu)
                 do imu=1,nmu
                  xmu = xmud(imu)

                  pq2=p2+dq2+2.d0*p*dq*xmu  

                  ep=  ener(irel,dmnu,p2)              
                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep-epq+uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep-epq+uyi*epsilon
                    
                   xreal_uaux= dreal(1.d0/uaux)

                  if (irel.eq.1) fmu(imu) =xreal_uaux*dmnu/epq
                  if (irel.eq.0) fmu(imu) =xreal_uaux
c enddo xmu
                 enddo

                call DRG20R(xmu_min, 1.d0, 20,fmu,resmu)

                else
                  resmu= 0.d0
               endif

               if (irel.eq.1) fp(ip) = resmu*p2*dmnu/ep 
               if (irel.eq.0) fp(ip) = resmu*p2
c enddo p
              enddo  
              call DRG20R(0.d0, dkf, 20,fp,resp)
              xnuc_j= xnuc_j+resp/pi2
c enddo ij
               enddo
              if (irel.eq.1) ximag= 2.d0*DUR_j(q0,DQ,DKF,dkf,DMNU)
              if (irel.eq.0) ximag= 2.d0*DUIM_j(q0,DQ,DKF,dkf,DMNU)

               unuc_j=xnuc_j+uyi*ximag 

               return
               end


