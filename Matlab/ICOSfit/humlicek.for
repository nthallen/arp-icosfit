      SUBROUTINE HUMDEV ( N, X, Y, K, L, DKDX, DKDY )

*     To calculate the Faddeeva function and partial derivatives of the Voigt function for y>=0

* Arguments
      INTEGER N                                                         ! IN   Number of points
      REAL    X(0:N-1)                                                  ! IN   Input x array
      REAL    Y                                                         ! IN   Input y value y>= 0
      REAL    K(0:N-1)                                                  ! OUT  Voigt array
      REAL    L(0:N-1)                                                  ! OUT  Imaginary part
      REAL    DKDX(0:N-1)                                               ! OUT  dVoigt/dx array
      REAL    DKDY(0:N-1)                                               ! OUT  dVoigt/dy array

* Constants
      DOUBLE PRECISION                DRTPI                             ! 1/SQRT(pi)
      REAL        RRTPI
      PARAMETER ( RRTPI = 0.56418958, DRTPI = 0.5641895835477563D0 )
      REAL        Y0,       Y0PY0,         Y0Q                          ! for CPF12 algorithm
      PARAMETER ( Y0 = 1.5, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0  )
      REAL C(0:5), S(0:5), T(0:5)
      SAVE C,      S,      T 
*     SAVE preserves values of C, S and T (static) arrays between procedure calls
      DATA C / 1.0117281,     -0.75197147,       0.012557727,
     &         0.010022008,   -0.00024206814,    0.00000050084806 /
      DATA S / 1.393237,       0.23115241,      -0.15535147,
     &         0.0062183662,   0.000091908299,  -0.00000062752596 /
      DATA T / 0.31424038,     0.94778839,       1.5976826,
     &         2.2795071,      3.0206370,        3.8897249 /


* Local variables
      INTEGER I, J                                                      ! Loop variables
      INTEGER RGB, RGC, RGD                                             ! y polynomial flags
      REAL ABX, XQ, YQ                                                  ! |x|, x^2 and y^2
      REAL XLIMA, XLIMB, XLIMC, XLIM4                                   ! x on region boundaries
      REAL MT(0:5), MQ(0:5), PT(0:5), PQ(0:5)                           ! Temporary variables
      REAL XP(0:5), XM(0:5), YP(0:5), YM(0:5), MF(0:5), PF(0:5)
      REAL YP2Y0, YPY0, YPY0Q, YF1, YF2, MFQ, PFQ, D,  U, DUDY, DVDY
      REAL A0, B1, C0, C2, D0, D1, D2, E0, E2, E4, F1, F3, F5
      REAL G0, G2, G4, G6, H0, H2, H4, H6, P0, P2, P4, P6, P8
      REAL Q1, Q3, Q5, Q7, R0, R2, W0, W2, W4, Z0, Z2, Z4, Z6, Z8
      DOUBLE PRECISION DB

***** Start of executable code ****************************************

      RGB = 1
      RGC = 1     
      RGD = 1
      YQ  = Y*Y                                                         ! y^2
      XLIMA = 146.7 - Y                                                 ! Region A boundary
      XLIMB = 24.0 - Y                                                  ! Region B boundary
      XLIMC = 7.4 - Y                                                   ! Region C boundary
      XLIM4 = 18.1*Y + 1.65                                             ! CPF12 I-II boundary
*.....
      DO I = 0, N-1                                                     ! Loop over all points
       ABX = ABS ( X(I) )                                               ! |x|
       XQ  = ABX*ABX                                                    ! x^2
       IF ( ABX .GT. XLIMA ) THEN                                       ! Region A
        D       = 1.0 / (XQ + YQ)
        D1      = D*RRTPI
        K(I)    = D1*Y
        L(I)    = D1*X(I)
        D1      = D1*D
        DKDX(I) = -D1*(Y+Y)*X(I)
        DKDY(I) = D1*(XQ-YQ)
       ELSEIF ( ABX .GT. XLIMB ) THEN                                   ! Region B
        IF ( RGB .NE. 0 ) THEN                                          ! First point in Region B
         RGB = 0
         A0 = YQ + 0.5                                                  ! Region A y-dependents
         B1 = YQ - 0.5
         D0 = A0*A0                                                     ! y^4 + y^2 + 0.25
         D2 = B1 + B1                                                   ! 2y^2 - 1
         C0 = 1.5   + YQ*(1.0 - D2)                                     ! 1.5 + 2y^2 - 2y^4
         C2 = A0 + A0                                                   ! 2y^2 + 1
         R0 = 0.125 + YQ*(0.25 - YQ*(0.5 + YQ))
         R2 = 0.25  + YQ*(5.0  + YQ)
        ENDIF
        D       = 1.0 / (D0 + XQ*(D2 + XQ))
        D1      = RRTPI*D
        K(I)    = D1*(A0 + XQ)*Y
        L(I)    = D1*(B1 + XQ)*X(I)
        D1      = D1*D
        DKDX(I) = D1*X(I)*Y*(C0 - (C2 + XQ)*(XQ+XQ))
        DKDY(I) = D1*(R0 - XQ*(R2 - XQ*(B1 + XQ)))
       ELSE                                                             ! Not Region A
        IF ( ABX .GT. XLIMC ) THEN                                      ! Region C
         IF ( RGC .NE. 0 ) THEN                                         ! First point in Region C
          RGC = 0
          H0 =  0.5625 + YQ*( 4.5  + YQ*(10.5 + YQ*(6.0 + YQ)) )        ! Region B y-dependents
          H2 = -4.5    + YQ*( 9.0  + YQ*( 6.0 + YQ* 4.0))
          H4 = 10.5    - YQ*( 6.0  - YQ*  6.0)
          H6 = -6.0    + YQ*  4.0
          W0 =  1.875  + YQ*(24.25 + YQ*(27.5 + YQ* 7.0))
          W2 =  5.25   + YQ*( 3.0  + YQ* 15.0)
          W4 = -4.5    + YQ*  9.0
          F1 = -1.875  + YQ*( 5.25 + YQ*( 4.5 + YQ))
          F3 =  8.25   - YQ*( 1.0  - YQ*  3.0)
          F5 = -5.5    + YQ*  3.0
          E0 = Y*(1.875 + YQ*( 8.25 + YQ*( 5.5 + YQ)))
          E2 = Y*(5.25  + YQ*( 1.0  + YQ*  3.0))
          E4 = Y*0.75*H6
          G0 = Y*(  9.0 + YQ*(42.0  + YQ*(36.0 + YQ* 8.0)))
          G2 = Y*( 18.0 + YQ*(24.0  + YQ* 24.0))
          G4 = Y*(-12.0 + YQ* 24.0)
          G6 = Y*   8.0
         ENDIF
         U = E0 + XQ*(E2 + XQ*(E4 + XQ*Y))
         D = 1.0 / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
         K(I) = RRTPI*D*U
         L(I) = RRTPI*D*X(I)*(F1 + XQ*(F3 + XQ*(F5 + XQ)))
         DUDY = W0 + XQ*(W2 + XQ*(W4 + XQ))
         DVDY = G0 + XQ*(G2 + XQ*(G4 + XQ*G6))
         DKDY(I) = RRTPI*D*(DUDY - D*U*DVDY)
        ELSEIF ( ABX .LT. 0.85 ) THEN                                   ! Region C
         IF ( RGD .NE. 0 ) THEN                                         ! First point in Region D
          RGD = 0
          Z0 =     272.1014 + Y*(1280.829   + Y*(2802.870  + Y*(3764.966! Region C y-dependents
     &        + Y*(3447.629 + Y*(2256.981   + Y*(1074.409
     &        + Y*(369.1989 + Y*(  88.26741 + Y*(  13.3988 + Y)))))))))
          Z2 =      211.678 + Y*( 902.3066  + Y*(1758.336  + Y*(2037.310
     &        + Y*(1549.675 + Y*( 793.4273  + Y*( 266.2987
     &                      + Y*(  53.59518 + Y*    5.0)))))))
          Z4 =     78.86585 + Y*( 308.1852  + Y*( 497.3014 + Y*(479.2576
     &        + Y*(269.2916 + Y*(  80.39278 + Y*   10.0)))))
          Z6 =     22.03523 + Y*(  55.02933 + Y*(  92.75679
     &                      + Y*(  53.59518 + Y*   10.0)))
          Z8 =     1.496460 + Y*(  13.39880 + Y*    5.0)
          P0 =     153.5168 + Y*( 549.3954  + Y*( 919.4955 + Y*(946.897
     &        + Y*(662.8097 + Y*( 328.2151  + Y*( 115.3772
     &        + Y*(27.93941 + Y*(   4.264678+ Y*    0.3183291))))))))
          P2 =    -34.16955 + Y*(  -1.322256+ Y*( 124.5975 + Y*(189.773
     &                      + Y*( 139.4665  + Y*(  56.81652
     &                      + Y*(  12.79458 + Y*    1.2733163))))))
          P4 =     2.584042 + Y*(  10.46332 + Y*(  24.01655
     &        + Y*(29.81482 + Y*(  12.79568 + Y*    1.9099744))))
          P6 =  -0.07272979 + Y*(   0.9377051
     &        + Y*(4.266322 + Y*    1.273316))
          P8 = 0.0005480304 + Y*    0.3183291
          Q1 =     173.2355 + Y*( 508.2585  + Y*( 685.8378
     &        + Y*(557.5178 + Y*( 301.3208  + Y*( 111.0528
     &        + Y*( 27.6294 + Y*(   4.26413 + Y*    0.3183291)))))))        
          Q3 =     18.97431 + Y*( 100.7375  + Y*( 160.4013 + Y*(130.8905
     &        + Y*(55.88650 + Y*(  12.79239 + Y*    1.273316)))))
          Q5 =     7.985877 + Y*(  19.83766 + Y*(  28.88480
     &                      + Y*(  12.79239 + Y*    1.909974)))
          Q7 =    0.6276985 + Y*(   4.26413 + Y*    1.273316)
         ENDIF
         U    = 1.7724538*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
         D    = 1.0 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8 + XQ)))))
         K(I) = D*U
         L(I) = 1.7724538*D*X(I)*
     &                 (Q1 + XQ*(Q3 + XQ*(Q5 + XQ*(Q7 + XQ*0.3183291))))
         DB = DBLE (X(I))*DBLE(L(I)) + DBLE(Y)*DBLE(K(I)) - DRTPI       ! Double precision
         DKDY(I) = SNGL ( DB + DB )                                     ! Single precision
        ELSE                                                            ! Use CPF12
         YPY0  = Y + Y0
         YPY0Q = YPY0*YPY0
         K(I) = 0.0
         L(I) = 0.0
         DO J = 0, 5
          MT(J) = X(I) - T(J)
          MQ(J) = MT(J)*MT(J)
          MF(J) = 1.0 / (MQ(J) + YPY0Q)
          XM(J) = MF(J)*MT(J)
          YM(J) = MF(J)*YPY0
          PT(J) = X(I) + T(J)
          PQ(J) = PT(J)*PT(J)
          PF(J) = 1.0 / (PQ(J) + YPY0Q)
          XP(J) = PF(J)*PT(J)
          YP(J) = PF(J)*YPY0
          L(I) = L(I) + C(J)*(XM(J)+XP(J)) + S(J)*(YM(J)-YP(J))
         ENDDO
         IF ( ABX .LE. XLIM4 ) THEN                                     ! Humlicek CPF12 Region I
          YF1 = YPY0 + YPY0
          YF2 = YPY0Q + YPY0Q
          DKDY(I) = 0.0
          DO J = 0, 5
           MFQ = MF(J)*MF(J)
           PFQ = PF(J)*PF(J)
           K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
           DKDY(I) = DKDY(I)
     &              + C(J)*( MF(J) + PF(J) - YF2*(MFQ + PFQ) )
     &              + S(J)*YF1*( MT(J)*MFQ - PT(J)*PFQ )
          ENDDO

         ELSE                                                           ! Humlicek CPF12 Region II
          YP2Y0 = Y + Y0PY0
          DO J = 0, 5
           K(I) = K(I)
     &           + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YP2Y0*XM(J))
     &             / (MQ(J)+Y0Q)
     &           + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YP2Y0*XP(J))
     &             / (PQ(J)+Y0Q)
          ENDDO
          K(I) = Y*K(I) + EXP ( -XQ )
          DB = DBLE (X(I))*DBLE(L(I)) + DBLE(Y)*DBLE(K(I)) - DRTPI      ! Double precision
          DKDY(I) = SNGL ( DB + DB )
         ENDIF
        ENDIF
        DB = DBLE(Y)*DBLE(L(I)) - DBLE(X(I))*DBLE(K(I))                 ! Double precision
        DKDX(I) = SNGL ( DB + DB )
       ENDIF                                                            ! Not region A
      ENDDO
*.....
      END
