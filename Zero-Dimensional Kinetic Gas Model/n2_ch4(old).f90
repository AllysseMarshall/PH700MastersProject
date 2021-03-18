program n2_ch4
 use params
 use reaction_scheme
!
!    20.5. 2016, D. Trunec
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! --- solution of kinetic equation for nitrogen and methane afterglow
!     units cm,s
!     O2, CO, CO2 added
!     13.11. 2017, D. Trunec 
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! link n2_ch4.f90 reaction_scheme.f90 params.f90 radau5.f decsol.f dc_decsol.f
!
        implicit none
        real(DP) ::  rpar(3)
        integer  ::  ipar(1)

! --- NUMBER OF EQUATIONS
        integer,parameter :: neqn=NUM_SPECIES
! --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        integer, parameter :: ND=NEQN,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20
        integer :: IJAC,MLJAC,MUJAC,MLMAS,MUMAS,IOUT,IMAS,IDID,ITOL
        
        real(DP) :: Y(ND),WORK(LWORK)
        integer :: IWORK(LIWORK)
! --- DECLARATIONS
        integer :: i,j,n
        real(DP) :: atol,rtol,t,tend,h
        real(DP) :: p_CH4,p_CO2,p_H2,Temp

  open(12,file='konc.dat')

  call compute_rates(IPAR,RPAR)
  call initiate_reaction_scheme(Y,IPAR,RPAR)
!  open(14,file='kr.txt')
!  do i=1,num_reactions
!     write(14,110) i,rate(i)
!    end do
!110 format (I5,1PE12.2)
!  close(14)

! --- DIMENSION OF THE SYSTEM
        N=NEQN
! --- COMPUTE THE JACOBIAN
        IJAC=0
! --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
! --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
! --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
! --- ENDPOINT OF INTEGRATION
        TEND=1.0e1_dp
! --- REQUIRED TOLERANCE
        RTOL=1.0e-12_dp
        ATOL=1.0_dp*RTOL
        ITOL=0
! --- INITIAL STEP SIZE
        H=1.0e-12_dp	
! --- SET DEFAULT VALUES
        DO I=1,LIWORK
           IWORK(I)=0
           WORK(I)=0.0_dp
        END DO

! --- INITIAL VALUES   indexy pro reakce_new_19.txt
        T=0.0_dp
!
        DO I=1,N
           Y(I)=0.0_dp
        END DO  
!
!      percentage of CH4 and CO2 in the mixture 
!
        p_CH4=0.02_dp
!
        p_CO2=0.01_dp
		
!       
!  N2
        y(86)=(1.0_dp-p_CH4-p_CO2)*n_N2
!
!  2% of nitrogen is dissociated
!
!  N(4S)
        y(85)=2.0_dp*0.02_dp*y(86)
!  N2(A) 
        y(87)=1.0e10_dp
!
!  all methane is dissociated
!
!  CH3
        y(14)=0.89_dp*p_CH4*n_N2	!should this not be 0.89, instead of 0.089?
!  CH2
        y(13)=0.1_dp*p_CH4*n_N2
!  CH
        y(6)=0.01_dp*p_CH4*n_N2
!  C
!        y(7)= 0.02_dp*0.1*n_N2	
!  H
        y(2)=y(14)+2.0_dp*y(13)+3.0_dp*y(6)*4.0_dp*y(7)
! NH3
!        y(96)=1.0e-3_dp*n_N2		
!
!  10% of CO2 is dissociated
!
! CO2
        y(185)=0.9*p_CO2*n_N2
! O(3P)
        y(172)=0.05*p_CO2*n_N2 !changed from 0.1
! CO    
        y(184)=0.05*p_CO2*n_N2	!changed from 0.1
! Wall
	y(236)=1e10_dp	!SSD value changed from 1e20_dp
!
! --- CALL OF THE SUBROUTINE RADAU5
        CALL RADAU5(N,DERIVS,T,Y,TEND,H, &
                        RTOL,ATOL,ITOL,  &
                        F_JAC,IJAC,MLJAC,MUJAC, &
                        DERIVS,IMAS,MLMAS,MUMAS, &
                        SOLOUT,IOUT, &
                        WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)


! --- PRINT STATISTICS

 open(1, file = 'data.dat', POSITION ='APPEND')
        write(1,92) T_g, y(3),y(15),y(17),y(21),y(86),y(96),y(103),y(105),y(184),y(188),y(196),y(206),y(220)
 92     format('',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4,',',1PE12.4)

close(1)

contains

SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, &
                                      RPAR,IPAR,IRTRN)
 
 INTEGER :: NR, LRC,N
 DOUBLE PRECISION :: X,Y(N),CONT(LRC),xold
 DOUBLE PRECISION :: RPAR(:)
 INTEGER :: IPAR(:)
 INTEGER :: IRTRN
 integer :: I
 write(12,*) x,(Y(I), I=1,N)
! write(12,110) x,y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8), &
!              y(9),y(10),y(11),y(12),y(13),y(14),y(15),  &
!              y(16),y(17),y(18),y(19),y(20),y(21),y(22),y(23), &
!              y(24),y(25),y(26),y(27),y(28),y(29),y(30),y(31), &
!              y(32),y(33),y(34),y(35),y(36),y(37)
!110     format(1PE20.10,37(1PE14.4E3))

END SUBROUTINE

end program n2_ch4

