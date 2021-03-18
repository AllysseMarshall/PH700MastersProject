module params
! for he_n2_o2.f90
 implicit none
 integer,parameter  :: DP=kind(1.0d0)

 real(DP),parameter :: pi=3.141592653589793238462643383279502884197_dp		
 real(DP),parameter :: e=1.60217646e-19_dp    ! elementary charge [C] 
 real(DP),parameter :: kb=1.3806503e-23_dp      ! m2 kg s-2 K-1
 real(DP),parameter :: amu=1.660538e-27_dp      ! m2 kg s-2 K-1
 real(DP),parameter :: mHe=4.002602_dp*amu      ! Mass of helium atom [kg]
 real(DP),parameter :: mAr=39.948_dp*amu        ! Mass of argon atom [kg]
 real(DP),parameter :: mO=15.9994_dp            ! Mass of oxygen atom [amu]
 real(DP),parameter :: mN=14.0067_dp            ! Mass of nitrogen atom [amu]
 real(DP),parameter :: me=9.10938188e-31_dp    ! Mass of electron [kg] 
 real(DP),parameter :: eps0=8.8541878176e-12_dp ! F/m (or C2 N-1 m-2)  
 real(DP),parameter :: p=90_dp             ! pressure of neutral gas [Pa]
 real(DP),parameter :: T_g=170_dp             ! temperature of neutral gas [K]
 real(DP),parameter :: T_e=170_dp             ! temperature of electrons [K] 
 real(DP),parameter :: n_N2=p/(kb*T_g)/1.0e6 _dp    ! nitrogen concentration [cm^-3]
 real(DP),parameter :: A = exp(4.0)				! for freezing eqn, A determines smoothness of transition from gas to frozen (higher is steeper, lower is more gradual)
 real(DP),parameter :: Gam = 0_dp				
 real(DP),parameter :: Uptake = 1.7e-4_dp 		!uptake coefficient
 real(DP),parameter :: SSD = n_N2*3e-17_dp	!Surface Space Density for Aerosols
 real(DP),parameter :: rad = 9e-5_dp					!radius of aerosol in cm
 real(DP),parameter :: Kc = 11427.35_dp*(rad**2)*Uptake*(T_g**(0.5))*SSD		
 real(DP),parameter :: Fc = 0_dp		!freezing constant (1=freezing, 0=no freezing)
end module params
