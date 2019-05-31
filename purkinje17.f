! 21 March 2007, version of purkinje.f to run on a node.
c 13 March 2007, begin construction of Purkinje cell model, starting with bask.f
c Sources include Llinas & Sugimori; Roth & Hausser for cell structure and passive
c parameters; de Schutter & Bower & Miyasho et al. for many of the active properties;
c Roth & Hausser for anomalous rectifier; Martina et al for data on delayed rectifier.
c 1 soma cylinder (level 1), 24 smooth dendritic compartments, #2 - 25.
c compartments 2 & 3 are dendritic shaft (level 2), 4 - 25 rest of smooth dendrites (level 3.
c There are 44 "canonical treelets" for the spiny dendrites; 2 treelets attach to
c each of the level 3 smooth dendritic compartments.
c A spiny treelet has 12 compartments: 4 along the main root, then 4 in each of 2 symmetrical
c branches.  The spiny dendrites are compartments 26 - 553.
c Axon has 6 5-micron compartments, #554 - 559, level 0.

c Active conductances: gnaf (m cubed h) - axon will use shifted version. 
c                      gnap (m cubed)
c                      gcaP (m) P channels
c                      gcaT (m h) T channels
c                      gcaR (m h) R channels
c                      gar  (m)   anomolous rectifier
c                      gKDR (m fourth) delayed rectifier, leaving out inactivation
c                      gKM  (m) M current, also called persistent K
c                      gKA  (m fourth x h)  A current
c                      gKD  (m h)  D current
c                      gKC  - try my original structure
c                      gKAHP - try my original structure.
c  12 Active conductances in all

        PROGRAM PURKINJE

       INTEGER J1, I, J, K, L, O, ISEED, K1, thisno
       REAL*8  TIMTOT, Z, Z1, Z2, curr(559), c(559), DT, time
       REAL*8 TIMER, gettime, time1, time2

c CINV is 1/C, i.e. inverse capacitance
       real*8 v(559), chi(559), cinv(559), gL(559), 
     x gNaF(559), gNaP(559), gCaP(559), gCaT(559), gCaR(559),
     x gar(559), gKDR(559), gKM(559), gKA(559), gKD(559),
     x gKC(559), gKAHP(559)

        real*8 jacob(559,559), betchi(559), gam(0:559,0:559)

        real*8
     x  mnaf(559), hnaf(559),
     x  mnap(559),
     x  mcap(559),
     x  mcat(559), hcat(559),
     x  mcar(559), hcar(559),
     x  mar (559),
     x  mkdr(559),
     x  mkm(559),
     x  mka(559), hka(559),
     x  mkd(559), hkd(559),
     x  mkc(559),
     x  mkahp(559)


       real*8  gampa(559),gnmda(559),ggaba_a(559),cafor(559)

       real*8 kdr_shift

       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_nap(0:640),betam_nap(0:640),dalpham_nap(0:640),
     X   dbetam_nap(0:640),

     X alpham_cap(0:640),betam_cap(0:640),dalpham_cap(0:640),
     X   dbetam_cap(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_car(0:640),betam_car(0:640),dalpham_car(0:640),
     X   dbetam_car(0:640),
     X alphah_car(0:640),betah_car(0:640),dalphah_car(0:640),
     X   dbetah_car(0:640),

     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640),

     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_kd(0:640), betam_kd(0:640),dalpham_kd(0:640) ,
     X   dbetam_kd(0:640),
     X alphah_kd(0:640), betah_kd(0:640), dalphah_kd(0:640),
     X   dbetah_kd(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640) 
       real*8 vL,vk,vna,var,vca,vgaba_a

        INTEGER NEIGH(559,6), NNUM(559)

c START EXECUTION PHASE
          include 'mpif.h'
          call mpi_init (info)
          call mpi_comm_rank(mpi_comm_world, thisno, info)
          call mpi_comm_size(mpi_comm_world, nodes , info)
          time1 = gettime()

        time1 = gettime()

c       kdr_shift =  7.d0
        kdr_shift = 10.d0

       CALL   PURKINJE_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_nap, betam_nap, dalpham_nap, dbetam_nap,

     X    alpham_cap, betam_cap, dalpham_cap, dbetam_cap,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_car, betam_car, dalpham_car, dbetam_car,
     X    alphah_car, betah_car, dalphah_car, dbetah_car,

     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,
     X    kdr_shift,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_kd , betam_kd , dalpham_kd , dbetam_kd ,
     X    alphah_kd , betah_kd , dalphah_kd , dbetah_kd ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc) 

          

        CALL PURKINJE_MAJ (GL,GAM,
     X    gNaF, gNaP, gCaP, gCaT, gCaR, gar, gKDR,
     X    gKM, gKA, gKD, gKC, gKAHP,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)

          do i = 1, 559
             cinv(i) = 1.d0 / c(i)
          end do

c                 gnaf = 0.d0
c                 gnap = 0.d0
c                 gkdr = 0.d0
c                 gka = 0.d0
c                 gkd = 0.d0
c                 gkm = 0.d0
c                 gkc = 0.d0
c                 gkahp = 0.d0
c                 gcat = 0.d0
c                 gcaP = 0.d0
c                 gcaR = 0.d0
c                 gar = 0.d0

c Below introduces somatic shunt
c            gL(1) = 100.d0 * gL(1)
             gL(1) =   5.d0 * gL(1)

c below 4 statements disconnect axon from soma.
c            gam(1,554) = 0.d0
c            gam(554,1) = 0.d0
c            jacob(1,554) = 0.d0
c            jacob(554,1) = 0.d0

        MG = 1.0d0
C  IN MILLIMOLAR

        VL = -80.d0
        VK =  - 85.d0
        VNA = 45.d0
        VCA = 135.d0
        VAR = -30.d0
        VGABA_A = -75.d0

        DT = .00060d0
c       DT = .00200d0

        TIMTOT = 1000.00d0
c       TIMTOT = 500.d0 * dt

c ? initialize membrane state variables?
        v = VL

        k1 = idnint (4.d0 * (v(1) + 120.d0))

      hnaf = alphah_naf(k1)/(alphah_naf(k1)+betah_naf(k1))
      hka = alphah_ka(k1)/(alphah_ka(k1)+betah_ka(k1))
      hkd = alphah_kd(k1)/(alphah_kd(k1)+betah_kd(k1))
      hcat=alphah_cat(k1)/(alphah_cat(k1)+betah_cat(k1))
      hcar=alphah_car(k1)/(alphah_car(k1)+betah_car(k1))

       O = 0
       TIME = 0.d0
       TIMER = 0.d0 ! for periodic spikes

2000      TIME = TIME + DT
          O = O + 1
          IF (TIME .GT. TIMTOT) GO TO 2001

          if (mod(O,8000).eq.0) timer = timer + 4.8d0

c         curr(1) = 1.50d0
c         curr(1) = 0.75d0

c         IF  (MOD(O,100).EQ. 0 )                      THEN
c         IF  ( i.eq.i          )                      THEN
c           WRITE (6,904) TIME, v(1), v(14), chi(14),
c    X V(25), v(554), v(557), v(558), v(559), curr(559)
c    X ,curr(14)
c         ENDIF
904       FORMAT(2X,F8.3,2X,8f8.2,f8.4,f9.4)

           if (time.le. 50.d0) then
          curr(1) = -1.00d0
           else if (time.le. 950.d0) then
          curr(1) = 1.50d0
           else
          curr(1) = -1.00d0
           endif

c           if (time.le. 50.d0) then
c         curr(14) = -0.1d0
c           else if (time.le.450.d0) then
c         curr(14) = 0.15d0
c           else
c         curr(14) = -0.1d0
c           endif

c Single antidromic spike
c          if (time.le.20.d0) then
c         curr(559) =  0.00d0
c          else if (time.le.20.50d0) then
c         curr(559) = 0.40d0 * (time-20.d0) / 0.5d0
c          else if (time.le.21.d0) then
c      curr(559) = 0.4d0 - 0.40d0 * (time-20.5d0) / 0.5d0
c          else
c         curr(559) = 0.d0
c          endif

c Periodic antidromic spikes
c            if (mod(O, 8000).le.900) then
c         curr(559) = 0.4d0 * (time-timer) / 0.54d0
c            else if (mod(O, 8000).lt.1800) then
c      curr(559) = 0.4d0 -0.4d0*(time-timer-0.54d0)/0.54d0
c            else
c          curr(559) = 0.0d0
c            endif

c            if (mod(O, 8000).le.1300) then
c         curr(559) = 0.5d0
c            else
c          curr(559) = 0.15d0
c            endif

909                  CONTINUE

         CALL   PURKINJE_INT (V,CHI,CINV,
     X mnaf,hnaf,mnap,mcap,mcat,hcat,mcar,hcar,mar,
     x mkdr,mkm,mka,hka,mkd,hkd,mkc,mkahp,
     x dt,neigh,nnum,jacob,mg,
     x vL,vk,vna,var,vca,vgaba_a,betchi,gam,gL,
     x gnaf,gnap,gcaP,gcaT,gCaR,gar,
     x gKDR,gKM,gKA,gKD,gKC,gKAHP,     
     x gampa,gnmda,ggaba_a,
     x O,time,
     X    alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_nap, betam_nap, dalpham_nap, dbetam_nap,

     X    alpham_cap, betam_cap, dalpham_cap, dbetam_cap,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_car, betam_car, dalpham_car, dbetam_car,
     X    alphah_car, betah_car, dalphah_car, dbetah_car,

     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,

     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_kd , betam_kd , dalpham_kd , dbetam_kd ,
     X    alphah_kd , betah_kd , dalphah_kd , dbetah_kd ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    cafor,curr)


              GO TO 2000
2001          CONTINUE

         time2 = gettime()
         write(6,309) time2 - time1
309      format(' Elapsed time = ',f8.0,' secs')


1000    CONTINUE
        call mpi_finalize(info)
        END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE PURKINJE_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_nap, betam_nap, dalpham_nap, dbetam_nap,

     X    alpham_cap, betam_cap, dalpham_cap, dbetam_cap,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_car, betam_car, dalpham_car, dbetam_car,
     X    alphah_car, betah_car, dalphah_car, dbetah_car,

     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,
     X    kdr_shift,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_kd , betam_kd , dalpham_kd , dbetam_kd ,
     X    alphah_kd , betah_kd , dalphah_kd , dbetah_kd ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc) 

      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_nap(0:640),betam_nap(0:640),dalpham_nap(0:640),
     X   dbetam_nap(0:640),
  

     X alpham_caP(0:640),betam_caP(0:640),dalpham_caP(0:640),
     X   dbetam_caP(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_car(0:640),betam_car(0:640),dalpham_car(0:640),
     X   dbetam_car(0:640),
     X alphah_car(0:640),betah_car(0:640),dalphah_car(0:640),
     X   dbetah_car(0:640),

     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640),

     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_kd(0:640), betam_kd(0:640),dalpham_kd(0:640) ,
     X   dbetam_kd(0:640),
     X alphah_kd(0:640), betah_kd(0:640), dalphah_kd(0:640),
     X   dbetah_kd(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640) 

       real*8  kdr_shift

C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION

       DO 1, I = 0, 640
          V = dble (I)
          V = (V / 4.d0) - 120.d0

c gNaF
c T. Miyasho et al. 2001 
           alpham_naf(i) = 35.d0 / (0.d0 + dexp((v + 5.d0)/(-10.d0))) 
           betam_naf(i) =  7.d0 / (0.d0 + dexp((v+65.d0)/20.d0))

c T. Miyasho et al. 2001
            alphah_naf(i) = 0.225d0 / (1.d0 + dexp((v+80.d0)/10.d0)) 
            betah_naf(i) =  7.5d0 / (0.d0 + dexp((v-3.d0)/(-18.d0)))

c T. Miyasho et al. 2001
            alpham_nap(i) = 200.d0  / (1.d0 + dexp((v-18.d0)/(-16.d0))) 
            betam_nap(i) =  25.d0 / (1.d0 + dexp((v+58.d0)/8.d0))

c P channels, T. Miyasho et al. 2001
            alpham_cap(i) = 8.5d0 / (1.d0 + dexp((v-8.d0)/(-12.5d0)))  
            betam_cap(i) =  35.d0 / (1.d0 + dexp((v+74.d0)/14.5d0))

c T-current, from T. Miyasho et al., 2001
            alpham_cat(i) = 2.6d0 / (1.d0 + dexp((v+21.d0)/(-8.d0))) 
            betam_cat(i) =  0.18d0 / (1.d0 + dexp((v+40.d0)/4.d0))
            alphah_cat(i) = 0.0025d0 / (1.d0 + dexp((v+40.d0)/8.d0))
            betah_cat(i) =  0.19d0 / (1.d0 + dexp((v+50.d0)/(-10.d0)))

c R-current, from T. Miyasho et al., 2001
            alpham_car(i) = 2.6d0 / (1.d0 + dexp((v+7.d0)/(-8.d0))) 
            betam_car(i) =  0.18d0/ (1.d0 + dexp((v+26.d0)/4.d0))
            alphah_car(i) = 0.0025/ (1.d0 + dexp((v+32.d0)/8.d0))
            betah_car(i) =  0.19d0/ (1.d0 + dexp((v+42.d0)/(-10.d0)))

c h-current (anomalous rectifier), Roth and Hausser 
           alpham_ar(i) = 0.00063d0 * dexp (-0.063d0 * (v + 73.2d0))
           betam_ar(i) =  0.00063d0 * dexp(0.079d0 * (v +73.2d0))

c delayed rectifier, non-inactivating. Start with mine (e.g. bask.f)
c Perhaps modify after comparison with Martina et al., 2003
           v = v + kdr_shift
c      minf = 1.d0/(1.d0+dexp((-v-27.d0)/11.5d0))
       minf = 1.d0/(1.d0+dexp((-v-20.d0)/11.5d0)) ! elevate threshold
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998 - for bask.f
           v = v - kdr_shift

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c A current: T. Miyasho et al., 2001
           alpham_ka(i) = 1.4d0 / (1.d0 + dexp((v+27.d0)/(-12.d0))) 
           betam_ka(i) =  0.49d0/ (1.d0 + dexp((v+30.d0)/4.d0))
           alphah_ka(i) = 0.0175d0/(1.d0 + dexp((v+50.d0)/8.d0))
           betah_ka(i) = 1.3 / (1.d0 + dexp((v+13.d0)/(-10.d0)))

c D current: T. Miyasho et al., 2001
           alpham_kd(i) = 8.5d0 / (1.d0 + dexp((v+17.d0)/(-12.5d0))) 
           betam_kd(i) =  35.0d0/ (1.d0 + dexp((v+99.d0)/14.5d0))
           alphah_kd(i) = 0.0015d0/(1.d0 + dexp((v+89.d0)/8.d0))
           betah_kd(i) = 0.0055d0 / (1.d0 + dexp((v+83.d0)/(-8.d0)))

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif
c Speed-up of C kinetics here.
          alpham_kc(i) = 2.d0 * alpham_kc(i)
           betam_kc(i) = 2.d0 *  betam_kc(i)

1        CONTINUE

         do 2, i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_nap(i) = (alpham_nap(i+1)-alpham_nap(i))/.25d0
      dbetam_nap(i) = (betam_nap(i+1)-betam_nap(i))/.25d0

      dalpham_cap(i) = (alpham_cap(i+1)-alpham_cap(i))/.25d0
      dbetam_cap(i) = (betam_cap(i+1)-betam_cap(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_car(i) = (alpham_car(i+1)-alpham_car(i))/.25d0
      dbetam_car(i) = (betam_car(i+1)-betam_car(i))/.25d0
      dalphah_car(i) = (alphah_car(i+1)-alphah_car(i))/.25d0
      dbetah_car(i) = (betah_car(i+1)-betah_car(i))/.25d0

      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0

      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_kd(i) = (alpham_kd(i+1)-alpham_kd(i))/.25d0
      dbetam_kd(i) = (betam_kd(i+1)-betam_kd(i))/.25d0
      dalphah_kd(i) = (alphah_kd(i+1)-alphah_kd(i))/.25d0
      dbetah_kd(i) = (betah_kd(i+1)-betah_kd(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
2      CONTINUE
       END

       SUBROUTINE PURKINJE_MAJ (GL,GAM,
     X    gNaF, gNaP, gCaP, gCaT, gCaR, gar, gKDR,
     X    gKM, gKA, gKD, gKC, gKAHP,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)

c Conductances: see main program.

c Note VAR = equil. potential for anomalous rectifier.

c 1 soma cylinder (level 1), 24 smooth dendritic compartments, #2 - 25.
c compartments 2 & 3 are dendritic shaft (level 2), 4 - 25 rest of smooth dendrites (level 3.
c There are 44 "canonical treelets" for the spiny dendrites; 2 treelets attach to
c each of the level 3 smooth dendritic compartments (22 of them).
c A spiny treelet has 12 compartments: 4 along the main root, then 4 in each of 2 symmetrical
c branches.  The spiny dendrites are compartments 26 - 553, level 4.
c Axon has 6 5-micron compartments, #554 - 559, level 0.

c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

        REAL*8 C(559),GL(559),GAM(0:559,0:559)
        REAL*8 gnaf(559), gnap(559), gcap(559), gcat(559), gcar(559)
        REAL*8 gar(559), gKDR(559), gKM(559), gKA(559), gKD(559),
     X   gKC(559), gKAHP(559) 
        REAL*8 JACOB(559,559),RI_SD,RI_AXON,RM_SD,RM_AXON,CDENS
        INTEGER LEVEL(559)
        REAL*8 GNAF_DENS(0:4), GNAP_DENS(0:4) 
        REAL*8 GCAP_DENS(0:4), GCAT_DENS(0:4), GCAR_DENS(0:4) 
        REAL*8 GAR_DENS(0:4)
        REAL*8 GKDR_DENS(0:4), GKM_DENS(0:4), GKA_DENS(0:4)
        REAL*8 GKD_DENS(0:4), GKC_DENS(0:4), GKAHP_DENS(0:4)

        REAL*8 RES, RINPUT, z1, z2, z, somaar
        REAL*8 RSOMA, PI, BETCHI(559), CAFOR(559)
        REAL*8 RAD(559), LEN(559), GAM1, GAM2, ELEN(559)
        REAL*8 RIN, D(559), AREA(559), RI
        INTEGER NEIGH(559,6), NNUM(559)
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS

        RI_SD = 115.d0 ! Roth & Hausser
c       RM_SD = 125000.d0 ! Roth & Hausser
c       RM_SD = 100000.d0 ! close to Roth & Hausser
        RM_SD =  50000.d0 ! close to Miyasho et al.
        RI_AXON = 100.d0
c       RM_AXON = 125000.d0 ! NOTE same as for SD
        RM_AXON =   2000.d0 
        CDENS = 0.8d0 ! Roth & Hausser

        PI = 3.14159d0

c See T. Miyasho et al., Brain Research 2001, for initial densities, but
c there will be changes.  e.g. Ca spikes have voltage-dependent fast-repol.
c so I will make dendritic gKDR larger and gKC smaller.
c ALSO, since there is now axon with voltage-shifted gNaF activation, I will
c try much smaller gNaF somatic density (Miyasho et al. have 10,000 ohm-cm-sqd)
c gNaF
           gnaf_dens(0) = 3500.d0
           gnaf_dens(1) = 5000.d0
c          gnaf_dens(1) = 7000.d0
           gnaf_dens(2) =  10.d0
           gnaf_dens(3) =   0.d0
           gnaf_dens(4) =   0.d0
c gNaP
           gnap_dens(0) = 0.1d0
           gnap_dens(1) = 5.0d0
           gnap_dens(2) = 1.0d0
           gnap_dens(3) = 0.0d0
           gnap_dens(4) = 0.d0
c gCaP
           gcaP_dens(0) = 0.0d0
           gcaP_dens(1) = 0.0d0
           gcaP_dens(2) = 0.0d0
           gcaP_dens(3) = 5.0d0
           gcaP_dens(4) = 5.0d0
c gCaT
           gcaT_dens(0) = 0.0d0
           gcaT_dens(1) = 0.0d0
           gcaT_dens(2) = 0.5d0
           gcaT_dens(3) = 1.5d0
           gcaT_dens(4) = 1.5d0
c gCaR
           gcaR_dens(0) = 0.0d0
           gcaR_dens(1) = 0.0d0
           gcaR_dens(2) = 0.0d0
           gcaR_dens(3) = 8.0d0
           gcaR_dens(4) = 8.0d0

c gar - Miyasho only have Ih at soma, but many princ cells have dendritic Ih
           gar_dens(0) = 0.0d0
           gar_dens(1) = 0.005d0
           gar_dens(2) = 0.005d0
           gar_dens(3) = 0.005d0
           gar_dens(4) = 0.005d0

c gKDR
c          gKDR_dens(0) = 3500.0d0
           gKDR_dens(0) = 1000.0d0
c          gKDR_dens(1) = 7000.0d0
           gKDR_dens(1) = 1000.0d0
           gKDR_dens(2) = 0.5d0
           gKDR_dens(3) = 0.5d0
           gKDR_dens(4) = 0.5d0
c gKM
           gKM_dens(0) = 1.00d0
           gKM_dens(1) = 1.00d0
           gKM_dens(2) = 1.00d0
           gKM_dens(3) = 0.04d0
           gKM_dens(4) = 0.04d0
c gKA
           gKA_dens(0) =  1.0d0
           gKA_dens(1) = 15.0d0
           gKA_dens(2) = 80.0d0
           gKA_dens(3) = 80.0d0
           gKA_dens(4) = 80.0d0
c gKD
           gKD_dens(0) =  0.0d0
           gKD_dens(1) =  0.0d0
           gKD_dens(2) = 80.0d0
           gKD_dens(3) = 80.0d0
c          gKD_dens(3) = 200.0d0
           gKD_dens(4) = 80.0d0
c          gKD_dens(4) = 200.0d0
c gKC
           gKC_dens(0) =  0.0d0
           gKC_dens(1) =  0.0d0
           gKC_dens(2) = 25.0d0
           gKC_dens(3) = 25.0d0
           gKC_dens(4) = 25.0d0
c gKAHP
         gkahp_dens(0) =   0.d0
         gkahp_dens(1) =   0.d0
         gkahp_dens(2) =   0.d0
         gkahp_dens(3) = 1.60d0
         gkahp_dens(4) = 1.60d0

        WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(P)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 4
          WRITE (6,9990) I, gnaf_dens(i), gcaR_dens(i), gkdr_dens(i),
     X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1

        level(2) = 2
        level(3) = 2

        do i = 4, 25
         level(i) = 3
        end do

        do i = 26, 553
         level(i) = 4
        end do

        do i = 554, 559
         level(i) = 0
        end do

c connectivity of axon
        nnum(554) = 2
        nnum(555) = 2
        nnum(556) = 2
        nnum(557) = 2
        nnum(558) = 2
        nnum(559) = 1
         neigh(554,1) =  1
         neigh(554,2) = 555
         neigh(555,1) = 554
         neigh(555,2) = 556
         neigh(556,1) = 555
         neigh(556,2) = 557
         neigh(557,1) = 556
         neigh(557,2) = 558
         neigh(558,1) = 557
         neigh(558,2) = 559
         neigh(559,1) = 558

c connectivity of SD part
          nnum(1) = 2  ! SOMA
          neigh(1,1) = 554
          neigh(1,2) =  2

c Now proximal smooth dendrite: no spiny treelets attached
          nnum(2) = 2
          neigh(2,1) = 1
          neigh(2,2) = 3

          nnum(3) = 3
          neigh(3,1) = 2
          neigh(3,2) = 4
          neigh(3,3) = 5

c Now distal smooth dendrite (level 3): 2 spiny treelets to each
          nnum(4) = 6
          neigh(4,1) = 3
          neigh(4,2) = 5
          neigh(4,3) = 12
          neigh(4,4) = 18
          neigh(4,5) = 26 ! base of spiny treelet
          neigh(4,6) = 38 ! base of spiny treelet

          nnum(5) = 5
          neigh(5,1) = 3
          neigh(5,2) = 4
          neigh(5,3) = 6 
          neigh(5,4) = 50 ! base of spiny treelet
          neigh(5,5) = 62 ! base of spiny treelet

          nnum(6) = 4
          neigh(6,1) = 5
          neigh(6,2) = 7
          neigh(6,3) = 74
          neigh(6,4) = 86

          nnum(7) = 4
          neigh(7,1) = 6
          neigh(7,2) = 8
          neigh(7,3) = 98
          neigh(7,4) = 110

          nnum(8) = 4
          neigh(8,1) = 7
          neigh(8,2) = 9
          neigh(8,3) = 122
          neigh(8,4) = 134

          nnum(9) = 4
          neigh(9,1) = 8
          neigh(9,2) = 10
          neigh(9,3) = 146
          neigh(9,4) = 158

          nnum(10) = 4
          neigh(10,1) = 9
          neigh(10,2) = 11
          neigh(10,3) = 170
          neigh(10,4) = 182

          nnum(11) = 3
          neigh(11,1) = 10
          neigh(11,2) = 194
          neigh(11,3) = 206

          nnum(12) = 5
          neigh(12,1) = 4 
          neigh(12,2) = 18 
          neigh(12,3) = 13  
          neigh(12,4) = 218 
          neigh(12,5) = 230 

          nnum(13) = 4
          neigh(13,1) = 12
          neigh(13,2) = 14 
          neigh(13,3) = 242 
          neigh(13,4) = 254 

          nnum(14) = 4
          neigh(14,1) = 13
          neigh(14,2) = 15 
          neigh(14,3) = 266 
          neigh(14,4) = 278 

          nnum(15) = 4
          neigh(15,1) = 14
          neigh(15,2) = 16 
          neigh(15,3) = 290 
          neigh(15,4) = 302 

          nnum(16) = 4
          neigh(16,1) = 15
          neigh(16,2) = 17 
          neigh(16,3) = 314 
          neigh(16,4) = 326 

          nnum(17) = 3
          neigh(17,1) = 16
          neigh(17,2) = 338
          neigh(17,3) = 350 

          nnum(18) = 5
          neigh(18,1) = 4 
          neigh(18,2) = 12 
          neigh(18,3) = 19  
          neigh(18,4) = 362 
          neigh(18,5) = 374 

          nnum(19) = 4
          neigh(19,1) = 18
          neigh(19,2) = 20 
          neigh(19,3) = 386 
          neigh(19,4) = 398 

          nnum(20) = 4
          neigh(20,1) = 19
          neigh(20,2) = 21 
          neigh(20,3) = 410 
          neigh(20,4) = 422 

          nnum(21) = 4
          neigh(21,1) = 20
          neigh(21,2) = 22 
          neigh(21,3) = 434 
          neigh(21,4) = 446 

          nnum(22) = 4
          neigh(22,1) = 21
          neigh(22,2) = 23 
          neigh(22,3) = 458 
          neigh(22,4) = 470 

          nnum(23) = 4
          neigh(23,1) = 22
          neigh(23,2) = 24 
          neigh(23,3) = 482 
          neigh(23,4) = 494 

          nnum(24) = 4
          neigh(24,1) = 23
          neigh(24,2) = 25 
          neigh(24,3) = 506 
          neigh(24,4) = 518 

          nnum(25) = 3
          neigh(25,1) = 24
          neigh(25,2) = 530
          neigh(25,3) = 542 

c Attachments of bases of treelets: one side to a smooth compartment, the other to next more distal part of treelet
          do i = 26, 530, 24
           nnum(i) = 2
           neigh(i,1) = (i - 26)/24 + 4
           neigh(i,2) = i + 1
          end do

          do i = 38, 542, 24
           nnum(i) = 2
           neigh(i,1) = (i - 38)/24 + 4
           neigh(i,2) = i + 1
          end do

c Now connect up the rest of the treelets
           do i = 27, 543, 12
            nnum(i) = 2
            neigh(i,1) = i-1
            neigh(i,2) = i+1
           end do

           do i = 28, 544, 12
            nnum(i) = 2
            neigh(i,1) = i-1
            neigh(i,2) = i+1
           end do

           do i = 29, 545, 12
            nnum(i) = 3
            neigh(i,1) = i-1
            neigh(i,2) = i+1
            neigh(i,3) = i+5
           end do

           do i = 30, 546, 12
            nnum(i) = 3
            neigh(i,1) = i-1
            neigh(i,2) = i+1
            neigh(i,3) = i+4
           end do

           do i = 31, 547, 12
            nnum(i) = 2
            neigh(i,1) = i-1
            neigh(i,2) = i+1
           end do

           do i = 32, 548, 12
            nnum(i) = 2
            neigh(i,1) = i-1
            neigh(i,2) = i+1
           end do

           do i = 33, 549, 12
              nnum(i) = 1
              neigh(i,1) = i-1
           end do

           do i = 34, 550, 12
              nnum(i) = 3
              neigh(i,1) = i-5
              neigh(i,2) = i-4
              neigh(i,3) = i+1
           end do

           do i = 35, 551, 12
              nnum(i) = 2
              neigh(i,1) = i-1
              neigh(i,2) = i+1
           end do

           do i = 36, 552, 12
              nnum(i) = 2
              neigh(i,1) = i-1
              neigh(i,2) = i+1
           end do

           do i = 37, 553, 12
              nnum(i) = 1
              neigh(i,1) = i-1
           end do

         DO 332, I = 1, 559
           WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
     X NEIGH(I,5),NEIGH(I,6)
3330     FORMAT(2X,I5,I5,I5,I5,I5,I5,I5)
332      CONTINUE
          DO 858, I = 1,559
           DO 858, J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
            DO 859, L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
859         CONTINUE
             IF (IT.EQ.0) THEN
              WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
             ENDIF
858       CONTINUE

c length and radius of axonal compartments
          do i = 554, 559
c           len(i) = 5.d0
            len(i) = 10.d0
          end do
          rad(554) = 0.750d0
          rad(555) = 0.700d0
          rad(556) = 0.650d0
          rad(557) = 0.600d0
          rad(558) = 0.550d0
          rad(559) = 0.5d0

c  length and radius of SD compartments
          len(1) = 29.d0
          rad(1) = 9.d0 ! SOMA

c Smooth dendrites
          do i = 2, 25
           len(i) = 15.d0
          end do
c Possibly lengthen dendritic shaft, to isolate (relatively) rest of dendrites.
          len(2) = 30.d0
          len(3) = 30.d0

          rad(2) =   2.25d0  ! proximal shaft
          rad(3) =   2.25d0  ! proximal shaft
          rad(2) =   1.80d0  ! proximal shaft ! What happens if this is narrower?
          rad(3) =   1.80d0  ! proximal shaft

           do i = 4, 11
             rad(i) = 1.8d0 ! a bit smaller
           end do
           do i = 12, 25
             rad(i) = 1.42d0
           end do

c Spiny dendrites
          do i = 26, 553
            len(i) = 25.d0
          end do

          do j = 26, 29
           do i = 0, 516, 12
             rad(i + j) = 0.75d0
           end do
          end do

          do j = 30, 37
           do i = 0, 516, 12
             rad(i + j) = 0.60d0
           end do
          end do


        WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
        DO 920, I = 1, 559
920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 559
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
          K = LEVEL(I)
         if (k .eq. 4) area(i) = 3.d0 * area(i)   ! spine correction
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAP(I) = GCAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GCAR(I) = GCAR_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKD(I) = GKD_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
         DO I = 2, 25
           Z = Z + AREA(I)
1019     END DO    

         z1 = 0.d0
         do i = 26, 553
           z1 = z1 + area(i)
         end do
         WRITE(6,1020) area(1),Z, z1
1020     FORMAT(2X,' SOMA AREA ',F7.0,
     X 'smooth dend ',f7.0,' spiny dends c spines ',f7.0)

        DO 140, I = 1, 559
        DO 140, K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )

         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
140     CONTINUE
c gam computed in microS

        DO 299, I = 1, 559
c299       BETCHI(I) = .20d0 ! NOTE HOW FAST
299       BETCHI(I) = .80d0 ! NOTE HOW FAST
        BETCHI( 1) =  .05d0

        DO 300, I = 1, 559
c300     D(I) = 1.D-4
c300     D(I) = 3.D-2
300     D(I) = 6.D-2
        DO 301, I = 1, 559
         IF (LEVEL(I).EQ.1) D(I) = 3.D-2
301     CONTINUE

       DO 160, I = 1, 559
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
C     NOTE CORRECTION

        do 200, i = 1, 559
200     C(I) = 1000.d0 * C(I)
C     TO GO FROM MICROF TO NF.

      DO 909, I = 1, 559
       JACOB(I,I) = - GL(I)
      DO 909, J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
             WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, 559
          do j = 1, 559
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

       DO 500, I = 1, 559
        WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F8.5)
500     CONTINUE
        END


       SUBROUTINE PURKINJE_INT (V,CHI,CINV,
     X mnaf,hnaf,mnap,mcap,mcat,hcat,mcar,hcar,mar,
     x mkdr,mkm,mka,hka,mkd,hkd,mkc,mkahp,
     x dt,neigh,nnum,jacob,mg,
     x vL,vk,vna,var,vca,vgaba_a,betchi,gam,gL,
     x gnaf,gnap,gcaP,gcaT,gCaR,gar,
     x gKDR,gKM,gKA,gKD,gKC,gKAHP,     
     x gampa,gnmda,ggaba_a,
     x O,time,
     X    alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_nap, betam_nap, dalpham_nap, dbetam_nap,

     X    alpham_cap, betam_cap, dalpham_cap, dbetam_cap,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_car, betam_car, dalpham_car, dbetam_car,
     X    alphah_car, betah_car, dalphah_car, dbetah_car,

     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar,

     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_kd , betam_kd , dalpham_kd , dbetam_kd ,
     X    alphah_kd , betah_kd , dalphah_kd , dbetah_kd ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    cafor,curr)

c CINV is 1/C, i.e. inverse capacitance
       real*8 v(559), chi(559), cinv(559)
       real*8 mnaf(559), hnaf(559), mnap(559), mkdr(559),
     x mka(559),hka(559),mkm(559),mkc(559),mkahp(559),
     x mkd(559),hkd(559),mcap(559),mcar(559),hcar(559),
     x mcat(559),hcat(559),mar(559),jacob(559,559),betchi(559),
     x gam(0:559,0:559),gL(559),
     x gnaf(559),gnap(559),gkdr(559),gka(559),gkd(559),
     x gkm(559),gkc(559),gkahp(559),gcat(559),gar(559),
     x gcap(559),gcar(559),
     x gampa(559),gnmda(559),ggaba_a(559),cafor(559) 

       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_nap(0:640),betam_nap(0:640),dalpham_nap(0:640),
     X   dbetam_nap(0:640),

     X alpham_cap(0:640),betam_cap(0:640),dalpham_cap(0:640),
     X   dbetam_cap(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_car(0:640),betam_car(0:640),dalpham_car(0:640),
     X   dbetam_car(0:640),
     X alphah_car(0:640),betah_car(0:640),dalphah_car(0:640),
     X   dbetah_car(0:640),

     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640),

     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_kd(0:640), betam_kd(0:640),dalpham_kd(0:640) ,
     X   dbetam_kd(0:640),
     X alphah_kd(0:640), betah_kd(0:640), dalphah_kd(0:640),
     X   dbetah_kd(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640) 

        real*8 fastna_shift

c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(559), fchi(559),fmnaf(559),fhnaf(559),fmkdr(559),
     x fmka(559),fhka(559),fmkd(559),fhkd(559),fmnap(559),
     x fmkm(559),fmkc(559),fmkahp(559),
     x fmcat(559),fhcat(559),fmar(559),
     x fmcap(559),fmcar(559),fhcar(559)

c below are for calculating the partial derivatives
       real*8 dfv_dv(559,559), dfv_dchi(559), dfv_dmnaf(559),
     x  dfv_dhnaf(559),dfv_dmnap(559),
     x  dfv_dmkdr(559),dfv_dmka(559),dfv_dhka(559),
     x  dfv_dmkd(559),dfv_dhkd(559),dfv_dmkm(559),dfv_dmkc(559),
     xdfv_dmkahp(559),dfv_dmcat(559),dfv_dhcat(559),dfv_dmcap(559),
     x  dfv_dmar(559),dfv_dmcar(559),dfv_dhcar(559)

        real*8 dfchi_dv(559), dfchi_dchi(559),
     x dfmnaf_dmnaf(559), dfmnaf_dv(559),dfhnaf_dhnaf(559),
     x dfhnaf_dv(559),
     x dfmnap_dmnap(559), dfmnap_dv(559),
     x dfmkdr_dmkdr(559),dfmkdr_dv(559),
     x dfmka_dmka(559),dfmka_dv(559),dfhka_dhka(559),dfhka_dv(559),
     x dfmkd_dmkd(559),dfmkd_dv(559),dfhkd_dhkd(559),dfhkd_dv(559),
     x dfmkm_dmkm(559),dfmkm_dv(559),dfmkc_dmkc(559),dfmkc_dv(559),
     x dfmcap_dmcap(559),dfmcap_dv(559),
     x dfmcat_dmcat(559),dfmcat_dv(559),dfhcat_dhcat(559),
     x dfhcat_dv(559),
     x dfmcar_dmcar(559),dfmcar_dv(559),dfhcar_dhcar(559),
     x dfhcar_dv(559),
     x dfmar_dmar(559),dfmar_dv(559),dfmkahp_dchi(559),
     x dfmkahp_dmkahp(559), dt2, outrcd(20), time

         REAL*8 dt,mg,vL,vk,vna,var,vca,vgaba_a,curr(559),Z,Z1,Z2
         INTEGER O, K0, K1, K2, NEIGH(559,6), NNUM(559)
       REAL*8 OPEN(559),gamma(559),gamma_prime(559)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(559), alpham_ahp_prime(559)
       REAL*8 gna_tot(559),gk_tot(559),gca_tot(559),gar_tot(559)

       DO 301, I = 1, 559
          FV(I) = -GL(I) * (V(I) - VL) * cinv(i)
          DO 302, J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K) - V(I)) * cinv(i)
301    CONTINUE


c       CALL FNMDA (V, OPEN, MG)

      DO 421, I = 1, 559
421    FV(I) = FV(I) + ( CURR(I)
     X   - gampa(I) * V(I) 
c    X   - (gampa(I) + open(i) * gnmda(I))*V(I)
     X   - ggaba_a(I)*(V(I)-Vgaba_a) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV

       do i = 1, 559
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i))
        if (chi(i).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
       end do

      DO 88, I = 1, 559
       gna_tot(i) = gnaf(i) * (mnaf(i)**3) * hnaf(i) +
     x     gnap(i) * (mnap(i)**3)
       gk_tot(i) = gkdr(i) * (mkdr(i)**4) +
     x             gka(i)  * (mka(i)**4) * hka(i) +
     x             gkd(i)  * (mkd(i)**4) * hkd(i) +
     x             gkm(i)  * mkm(i) +
     x             gkc(i)  * mkc(i) * gamma(i) +
     x             gkahp(i)* mkahp(i)
       gca_tot(i) = gcat(i) *  mcat(i) * hcat(i) + ! NOTE CHANGE TO 1st power m
     x              gcaP(i) *  mcaP(i) +
     x              gcaR(i) *  mcaR(i) * hcaR(i) 

       gar_tot(i) = gar(i) * mar(i)


88     FV(I) = FV(I) - ( gna_tot(i) * (v(i) - vna)
     X  + gk_tot(i) * (v(i) - vK)
     X  + gca_tot(i) * (v(i) - vCa)
     X  + gar_tot(i) * (v(i) - var) ) * cinv(i)

         do i = 1, 559
         do j = 1, 559
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X    + ggaba_a(i) + gampa(i)
     X   + open(i) * gnmda(I) )
          endif
         end do
         end do

          do i = 1, 559
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i) *
     x                     gamma_prime(i) * (v(i)-vK)
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i)**2) *
     X    (gnaf(i) * hnaf(i) + gnap(i)) * (v(i) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i)**3) *
     X                    (v(i) - vna)
        dfv_dmnap(i) = -3.d0 * cinv(i) * (mnap(i)**2) *
     X    (gnap(i) + gnap(i)) * (v(i) - vna)

        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i)**3)
     X                   * (v(i) - vK)
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i)**3) *
     X                   hka(i) * (v(i) - vK)
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i)**4) *
     X                    (v(i) - vK)
        dfv_dmkd(i)  = -4.d0 * cinv(i) * gkd(i) * (mkd(i)**3) *
     X                   hkd(i) * (v(i) - vK)
        dfv_dhkd(i)  = - cinv(i) * gkd(i) * (mkd(i)**4) *
     X                    (v(i) - vK)

        dfv_dmkm(i)  = - cinv(i) * gkm(i) * (v(i) - vK)
        dfv_dmkc(i)  = - cinv(i) * gkc(i) * gamma(i) * (v(i)-vK)
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i) - vK)

        dfv_dmcat(i)  = -1.d0 * cinv(i) * gcat(i) * 1.d0    *
     X                    hcat(i) * (v(i) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * mcat(i) *
     X                  (v(i) - vCa)
        dfv_dmcar(i)  = -1.d0 * cinv(i) * gcar(i) * 1.d0    *
     X                    hcar(i) * (v(i) - vCa)
        dfv_dhcar(i) = - cinv(i) * gcar(i) * mcar(i) *
     X                  (v(i) - vCa)
        dfv_dmcap(i) = -1.d0 * cinv(i) * gcap(i) * 1.d0    *
     X                      (v(i) - vCa)

        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i) - var)
          end do

         do i = 1, 559
          fchi(i) = - cafor(i) * gca_tot(i) * (v(i) - vca)
     x       - betchi(i) * chi(i)
          dfchi_dv(i) = - cafor(i) * gca_tot(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, 559
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i),0.01d0) ! original from baskint
        alpham_ahp(i) = dmin1(6.0d-4 * chi(i),0.30d0) ! speed this up
        if (chi(i).le.500.d0) then
          alpham_ahp_prime(i) = 6.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, 559
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i))
c    x                  -.001d0 * mkahp(i)
     x                  -.060d0 * mkahp(i) ! speed this up significantly
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i))
       end do

          do i = 1, 559

       K1 = IDNINT ( 4.d0 * (V(I) + 120.d0) ) ! For SD
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c            fastNa_shift =  2.0d0 ! For axon
             fastNa_shift =  6.0d0 ! For axon
c            fastNa_shift =  0.0d0 ! For axon
       K0 = IDNINT ( 4.d0 * (V(I)+      fastNa_shift+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0


           if (i.lt.554) then
            k2 = k1
           else
            k2 = k0
           endif
        fmnaf(i) = alpham_naf(k2) * (1.d0 - mnaf(i)) -
     X              betam_naf(k2) * mnaf(i)
        fhnaf(i) = alphah_naf(k1) * (1.d0 - hnaf(i)) -
     X              betah_naf(k1) * hnaf(i)
        fmnap(i) = alpham_nap(k1) * (1.d0 - mnap(i)) -
     X              betam_nap(k1) * mnap(i)

        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i)) -
     X              betam_kdr(k1) * mkdr(i)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i)) -
     X              betam_ka (k1) * mka(i)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i)) -
     X              betah_ka (k1) * hka(i)
        fmkd(i)  = alpham_kd (k1) * (1.d0 - mkd(i)) -
     X              betam_kd (k1) * mkd(i)
        fhkd(i)  = alphah_kd (k1) * (1.d0 - hkd(i)) -
     X              betah_kd (k1) * hkd(i)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i)) -
     X              betam_km (k1) * mkm(i)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i)) -
     X              betam_kc (k1) * mkc(i)

        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i)) -
     X              betam_cat(k1) * mcat(i)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i)) -
     X              betah_cat(k1) * hcat(i)
        fmcar(i) = alpham_car(k1) * (1.d0 - mcar(i)) -
     X              betam_car(k1) * mcar(i)
        fhcar(i) = alphah_car(k1) * (1.d0 - hcar(i)) -
     X              betah_car(k1) * hcar(i)
        fmcap(i) = alpham_cap(k1) * (1.d0 - mcap(i)) -
     X              betam_cap(k1) * mcap(i)

        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i)) -
     X              betam_ar (k1) * mar(i)

       dfmnaf_dv(i) = dalpham_naf(k2) * (1.d0 - mnaf(i)) -
     X                  dbetam_naf(k2) * mnaf(i)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i)) -
     X                  dbetah_naf(k1) * hnaf(i)
       dfmnap_dv(i) = dalpham_nap(k1) * (1.d0 - mnap(i)) -
     X                  dbetam_nap(k1) * mnap(i)

       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i)) -
     X                  dbetam_kdr(k1) * mkdr(i)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i)) -
     X                  dbetam_ka(k1) * mka(i)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i)) -
     X                  dbetah_ka(k1) * hka(i)
       dfmkd_dv(i)  = dalpham_kd(k1) * (1.d0 - mkd(i)) -
     X                  dbetam_kd(k1) * mkd(i)
       dfhkd_dv(i)  = dalphah_kd(k1) * (1.d0 - hkd(i)) -
     X                  dbetah_kd(k1) * hkd(i)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i)) -
     X                  dbetam_km(k1) * mkm(i)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i)) -
     X                  dbetam_kc(k1) * mkc(i)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i)) -
     X                  dbetam_cat(k1) * mcat(i)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i)) -
     X                  dbetah_cat(k1) * hcat(i)
       dfmcar_dv(i) = dalpham_car(k1) * (1.d0 - mcar(i)) -
     X                  dbetam_car(k1) * mcar(i)
       dfhcar_dv(i) = dalphah_car(k1) * (1.d0 - hcar(i)) -
     X                  dbetah_car(k1) * hcar(i)
       dfmcap_dv(i) = dalpham_cap(k1) * (1.d0 - mcap(i)) -
     X                  dbetam_cap(k1) * mcap(i)

       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i)) -
     X                  dbetam_ar(k1) * mar(i)

       dfmnaf_dmnaf(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmnap_dmnap(i) =  - alpham_nap(k1) - betam_nap(k1)

       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmkd_dmkd(i)  =   - alpham_kd (k1) - betam_kd (k1)
       dfhkd_dhkd(i)  =   - alphah_kd (k1) - betah_kd (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)

       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcar_dmcar(i) =  - alpham_car(k1) - betam_car(k1)
       dfhcar_dhcar(i) =  - alphah_car(k1) - betah_car(k1)
       dfmcap_dmcap(i) =  - alpham_cap(k1) - betam_cap(k1)

       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, 559
          v(i) = v(i) + dt * fv(i)
           do j = 1, 559
        v(i) = v(i) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i) = v(i) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmkd(i)  * fmkd(i)
     X          + dfv_dhkd(i)  * fhkd(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcar(i)  * fmcar(i)
     X          + dfv_dhcar(i) * fhcar(i)
     X          + dfv_dmcap(i) * fmcap(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i) = chi(i) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))

        mnaf(i) = mnaf(i) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        hnaf(i) = hnaf(i) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mnap(i) = mnap(i) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))

        mkdr(i) = mkdr(i) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i) =  mka(i) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i) =  hka(i) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mkd(i) =  mkd(i) + dt * fmkd(i) + dt2 *
     X   (dfmkd_dmkd(i) * fmkd(i) + dfmkd_dv(i) * fv(i))
        hkd(i) =  hkd(i) + dt * fhkd(i) + dt2 *
     X   (dfhkd_dhkd(i) * fhkd(i) + dfhkd_dv(i) * fv(i))
        mkm(i) =  mkm(i) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i) =  mkc(i) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i) = mkahp(i) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i) =  mcat(i) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i) =  hcat(i) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcar(i) =  mcar(i) + dt * fmcar(i) + dt2 *
     X   (dfmcar_dmcar(i) * fmcar(i) + dfmcar_dv(i) * fv(i))
        hcar(i) =  hcar(i) + dt * fhcar(i) + dt2 *
     X   (dfhcar_dhcar(i) * fhcar(i) + dfhcar_dv(i) * fv(i))
        mcap(i) =  mcap(i) + dt * fmcap(i) + dt2 *
     X   (dfmcap_dmcap(i) * fmcap(i) + dfmcap_dv(i) * fv(i))

        mar(i) =   mar(i) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
         end do


c      IF (MOD(O,75).EQ.0) THEN
       IF (MOD(O,150).EQ.0) THEN
           OUTRCD(1) = TIME
           OUTRCD(2) = v(1)
           outrcd(3) = v(14)
           outrcd(4) = V(17)    
           outrcd(5) = v( 22)
           outrcd(6) = v(555)
           outrcd(7) = v(559)
           outrcd(8) = chi(1)
       outrcd(9) =  chi(14)
       outrcd(10)=  chi(22)
       outrcd(11)=  chi(17)
           outrcd(12) = curr(1)
           outrcd(13) = curr(559)
           outrcd(14) = curr( 14)
           outrcd(15) = gna_tot(1) * (v(1)-vna) * cinv(1)
           outrcd(16) = gk_tot(1) * (v(1)-vk) * cinv(1)
           outrcd(17) = gca_tot( 14) * (v( 14)-vca)*cinv( 14)
           outrcd(18) = gk_tot( 14) * (v( 14)-vk)*cinv( 14)
         OPEN(11,FILE='PURK17.OU')
         WRITE (11,FMT='(18F10.3)') (OUTRCD(I),I=1,18)
C900      FORMAT (100A4)
       ENDIF


              END

               SUBROUTINE FNMDA (VSTOR,OPEN,MG)
               REAL*8 VSTOR(559), OPEN(559)
               REAL*8 A, BB1, BB2, VM, A1, A2, B1, B2, MG
c modify so that potential is absolute and not relative to
c  "rest"
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)
C  TO DETERMINE VOLTAGE-DEPENDENCE OF NMDA CHANNELS
           DO 1, I = 1, 559
c          VM = VSTOR(I) - 60.
           VM = VSTOR(I)
           A1 = DEXP(-.016d0*VM - 2.91d0)
           A2 = 1000.d0 * MG * DEXP (-.045d0 * VM - 6.97d0)
           B1 = DEXP(.009d0*VM + 1.22d0)
           B2 = DEXP(.017d0*VM + 0.96d0)
        OPEN(I)     = 1.d0/(1.d0 + (A1+A2)*(A1*BB1 + A2*BB2) /
     X   (A*A1*(B1+BB1) + A*A2*(B2+BB2))  )
C  FROM JAHR & STEVENS, EQ. 4A
C               DO 124, J = 1, 19
C          OPEN(J) = 1./(1.+.667* EXP(-0.07*(VSTOR(J)-60.)) )
C  FROM CHUCK STEVENS
1               CONTINUE
                        END

