!---------------------------------------------------------------------------------------------------------------------------------
! T.O.C.
!  1. program     ck_calc
!  2. subroutine  read_gas_data
!  3. subroutine  abs_line_data
!  4. subroutine  abs_coeff_calc
!  5. subroutine  gk_calc
!  6. subroutine  get_gaussian_abscissae
!  7. subroutine  quicksort
!  8. subroutine  swap
!  9. subroutine  swap_mask
! 10. function    chi
! 11. subroutine  gruszka
! 12. function    xk1
! 13. subroutine  read_baranov_data
! 14. subroutine  read_h2o_continuum_data
! 15. subroutine  read_co2_continuum_data
! 16. subroutine  bilinear
! 17. function    norm_fac
! 18. function    f_vvw
! 19. subroutine  continuum_calc
! 20. function    chi_h2o
! 21. subroutine  humlik
!---------------------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------------------------------
! 1.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Program to create the reference files for the c-k method
!c-mm
!c-mm This is the main program which generates very accurate g(k) 
!c-mm    distributions for a specified gas and range of wavenumbers for
!c-mm    a specified range of pressures and temperatures.  These g(k) 
!c-mm    distributions are output to files, to be read by the program 
!c-mm    which merges the g(k) distribution into just a few k terms.
!c-mm *NOTE*  All units here are mks*, including wavenumber (m-1).  
!c-mm    Items read in from HITRAN are converted to mks, generally by
!c-mm    multiplying by 100.  Be aware of this peculiarity, especially
!c-mm    for measures of wavenumber, which appear 100x larger than they are.  
!c-mm This version of the code uses dynamically allocated arrays where 
!c-mm    possible.  Read comments for further information.
!-----------------------------------------------------------------------------------------------------------------------------------
PROGRAM ck_calc

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE abs_coeff_calc(line_dat, range_min, range_max, gas,    &
                kvals, vstep, numvals, numlines, temp, maxwid,         &
                too_long, failed_wvn)

       IMPLICIT NONE
       INTEGER, INTENT(IN   ) ::             &
          numlines,                          &
          numvals,                           &
          gas                             
       REAL*8, INTENT(IN   ) ::              &
          line_dat(:,:),                     & 
          range_min,                         &
          range_max,                         &
          vstep,                             &
          temp,                              &
          maxwid
       REAL*8, INTENT(  OUT) ::              &
          kvals(:)
       LOGICAL, INTENT(  OUT) ::             &
          too_long
       REAL*8, INTENT(  OUT) ::              &
          failed_wvn
     END SUBROUTINE abs_coeff_calc

     SUBROUTINE humlik(line_dat, range_min, range_max, gas, kvals,     & 
                vstep, numvals, numlines, temp, maxwid, too_long,      &
                failed_wvn,kval_dim)

       IMPLICIT NONE
       INTEGER, INTENT(IN   ) ::             &
          numlines,                          &
          numvals,                           &
          gas,                               &
          kval_dim
       REAL*8, INTENT(IN   ) ::              &
          line_dat(:,:),                     & 
          range_min,                         &
          range_max,                         &
          vstep,                             &
          temp,                              &
          maxwid
       REAL*8, INTENT(  OUT) ::              &
          kvals(:)
       LOGICAL, INTENT(  OUT) ::             &
          too_long
       REAL*8, INTENT(  OUT) ::              &
          failed_wvn
     END SUBROUTINE humlik

     SUBROUTINE continuum_calc(range_min, range_max, gas, k_self_cont, &
          k_frn_cont,vstep, numvals, temp,pratio,wn_arr,temp_arr,      &
          dim_self_arr,dim_frn_arr,wn_CIA_arr,temp_CIA_arr,            &
          dim_CIA_arr,kval_dim)

       IMPLICIT NONE
       INTEGER, PARAMETER ::         &
          nS = 3000,                 &
          nT = 9                    
       INTEGER, INTENT(IN   ) ::     &
          gas,                       & ! HITRAN number of gas of interest
          kval_dim
       REAL*8, INTENT(IN   ) ::      &
          range_min,                 & ! The min wavenumber of spectral range of g(k) distribution
          range_max,                 & ! The max wavenumber of spectral range of g(k) distribution
          vstep,                     & ! Amount (in m-1) by which the spectrum is stepped over
          temp                         ! Temperature used to calculate chi and scale density
       REAL*8, INTENT(IN   ), OPTIONAL ::                   &
          pratio                       ! Pressure used to scale density.  In bars
       REAL*8, INTENT(IN   ), DIMENSION(nS), OPTIONAL ::    &
          wn_CIA_arr
       REAL*8, INTENT(IN   ), DIMENSION(nT), OPTIONAL ::    &
          temp_CIA_arr
       REAL*8, INTENT(IN   ), DIMENSION(nS,nT), OPTIONAL :: &
          dim_CIA_arr
       REAL*8, INTENT(IN   ), DIMENSION(:,:), OPTIONAL ::   &
          dim_self_arr,                                     &
          dim_frn_arr                
       REAL*8, INTENT(IN   ), DIMENSION(:), OPTIONAL ::     & 
          wn_arr                      
       REAL*8, INTENT(IN   ), DIMENSION(:), OPTIONAL ::     &
          temp_arr                      
       INTEGER, INTENT(  OUT) ::     &
          numvals                      ! Number of values of the absorption coefficient calculated
       REAL*8, INTENT(  OUT) ::      &
          k_self_cont(:),          & ! Array to store self continuum values across the spectral range
          k_frn_cont(:)              ! Array to store foreign continuum values across the spectral range
     END SUBROUTINE continuum_calc

     SUBROUTINE gk_calc(kvals,numvals,gkdat,kk,num_terms)
       IMPLICIT NONE
       INTEGER, INTENT(IN   ) ::             &
          numvals,                           &
          kk,                                &
          num_terms
       REAL*8, INTENT(IN   ) ::              &
          kvals(:)
       REAL*8, INTENT(  OUT) ::              &
          gkdat(num_terms)	       
     END SUBROUTINE gk_calc
  END INTERFACE


  INTEGER, PARAMETER ::                      &
     t_num = 18,                             & ! Number of temperatures being calculated over in the spectral band
     p_num = 26,                             & ! Number of pressures being calculated.
     spectral_bands = 14,                    & ! Number of spectral bands (both IR and solar combined)
     ! The maximum number of k-terms in any one band  12-9-08.  Trying 32 for Gauss-double quad.
     num_terms = 32,                         & 
     max_num_gases = 3,                      & ! Maximum number of gases we'll ever use.
     num_fixed_gases = 2,                    & ! The number of gases in the experiment that have a fixed mixing ratio
     num_variable_gases = 1,                 & ! The number of gases in the experiment that have a variable mixing ratio
     num_mix_ratio = 8,                      & ! Number of mixing ratios being calculated for a given spectral band 
     num_hitran_gases = 36,                  & ! Number of HITRAN gases that are available (based upon hitran_gas_data)
     nS = 3000,                              & ! Number of wavenumbers in Baranov data set
     nT = 9,                                 & ! Number of temperatures in Baranov data set
     nS_h2o = 2001,                          & ! Number of wavenumbers in the water vapor continuum data set
     nT_h2o = 72,                            & ! Number of temperatures in the water vapor continuum data set
     nS_co2 = 5001,                          & ! Number of wavenumbers in the carbon dioxide continuum data set
     nT_co2 = 71,                            & ! Number of temperatures in the carbon dioxide continuum data set
     ! Index (of margin_list) from which to start the water calculations.  H2O is always 25 cm-1
     starting_h2o_margin_index = 5,          & 
     starting_margin_index = 1                 ! Index (of margin list) from which to start all other gas calculations.

  REAL*8, PARAMETER ::                       &
     t_init = 50.0D0,                        & ! Initial T for constructing k-values
     delta_t = 20.0D0                          ! Stride in T for constructing k-values

  ! This derived data type contains all the raw data information from the HITRAN database
  TYPE :: hitrandata                            
     REAL*8, DIMENSION(80000,6) ::           &
       rawdat1                                  ! The entire HITRAN database only has 63197 H2O lines and 62913 CO2 lines.
     INTEGER, DIMENSION(80000)  ::           &
       rawdat2                                  ! This bound (80000) will ensure we are larger than anything we might get.
     INTEGER                    ::           &
       num_rawdat
  END TYPE hitrandata

  TYPE(hitrandata) ::                        &
       raw_data(max_num_gases,8)                ! Each gas/isotope stores its own information.
  
  INTEGER                                    &
     i, j, k, kk, jj, x, y,                  & ! Loop variables
     length_count,                           & ! Number of letters in output_dir, to ease concatenation of filename
     mol_no,                                 & ! Molecular number : the HITRAN label
     iso_no,                                 & ! The isotope number : the HITRAN label
     numvals,                                & ! Number of values of the absorption coefficient calculated
     gas,                                    & ! HITRAN number of gas of interest (H2O=1, CO2=2)
     instat,                                 & ! Variable used for determining end of HITRAN file
     kval_dim,                               & ! Calculated value for the number of kvals in a band.  Used to allocate array sizes
     gas_index(max_num_gases),               & ! The HITRAN index of the gases involved in the calculation
     p_start_index,			                 & ! The index of the press array which we are starting from
     p_end_index,			                 & ! The index of the press array which we are ending at
     t_start_index,                          & ! The index of the temperature array which we are starting from
     t_end_index,                            & ! The index of the temperature array which we are ending at
     spec_index,			                 & ! The index of the spectral band array which we are calculating for
     mix_ratio_start_index,		             & ! The index of the mixratio array which we are starting from
     mix_ratio_end_index,		             & ! The index of the mixratio array which we are ending at
     number_of_isotopes(num_hitran_gases),   & ! The number of isotopes of each gas
     number_of_active_isotopes(num_hitran_gases),   & ! The number of isotopes of each gas we are looking at (often =1)
     margin_index(max_num_gases),             &   ! Counting index for choosing the appropriate band margin
     num_gases ! number of gases we actually use.

  REAL*8                                     & 
     range_min,                              & ! The minimum value of the spectral range for the corresponding g(k) file.
     range_max,                              & ! The maximum value of the spectral range for the corresponding g(k) file.
     press(p_num),                           & ! Array to store all pressures (in Pa) at which g(k) calculation done.  
     wn(spectral_bands+1),                   & ! Array to store the wavenumber boundaries
     layer_temp(t_num),                      & ! Array to store all temperatures (kelvin) for which g(k) calculation is done
     pratio,temp,                            & ! Pressure and temperature at which calculation is done
     vstep,                                  & ! Amount (in m-1) by which the spectrum is stepped over
     gkdat(num_terms),                       & ! Array for all g(k) distribution data
     line_cen,                               & ! Wavenumber of line center (cm-1 from HITRAN, outputted as m-1)
     line_stren,                             & ! Line strength (cm-1/(molec.cm-2) from HITRAN, outputted as m/kg)
     temp1,                                  & ! Temporary number, of no real meaning
     abhw,                                   & ! Air broadened halfwidth (cm-1/atm from HITRAN,outputted as m-1/atm)
     sbhw,				                     & ! Self-broadened halfwidth (cm-1/atm from HITRAN, outputted as m-1/atm)
     low_ener,                               & ! Lower state energy (cm-1 from HITRAN, outputted as m-1)
     abhw_coeff,                             & ! Coefficient of temperature dependence of air broadened halfwidth
     kvals_at_one_p(t_num,p_num,num_terms),  & ! Array to temporarily hold kvals for all temperatures at one pressure for output
     q_coef(num_hitran_gases,4,8),	         & ! The polynomial coefficients to calculate TIPS for all gases/isotopes
     q_ref(num_hitran_gases,8),	             & ! Reference values of Q(T=296)
     iso_mass(num_hitran_gases,8),	         & ! Mass of individual isotopes in amu
     iso_fraction(num_hitran_gases,8),       & ! Fraction of total gas abundance that is of specific isotopologue
     wn_CIA_arr(nS),                         & ! Array of wavenumbers for Baranov CIA calcuation
     temp_CIA_arr(nT),                       & ! Array of temperatures for Baranov CIA calculation
     dim_CIA_arr(nS,nT),                     & ! Array of values for Baranov CIA calculation
     wn_h2o_arr(nS_h2o),                     & ! Array of wavenumbers for h2o continuum calculation
     temp_h2o_arr(nT_h2o),                   & ! Array of temperatures for h2o continuum calculation
     dim_h2o_self_arr(nS_h2o,nT_h2o),        & ! Array of values for h2o self continuum calculation
     dim_h2o_frn_arr(nS_h2o,nT_h2o),         & ! Array of values for h2o foreign continuum calculation
     wn_co2_arr(nS_co2),                     & ! Array of wavenumbers for co2 continuum calculation
     temp_co2_arr(nT_co2),                   & ! Array of temperatures for co2 continuum calculation
     dim_co2_arr(nS_h2o,nT_h2o),             & ! Array of values for co2 continuum calculation
     mixratio_list(max_num_gases),           & ! The list of mixing ratios being used for a single loop over jj    
     variable_gas_mixratio(num_mix_ratio, num_variable_gases), & ! Variable gases have the selection of mixing ratios stored here.
     fixed_gas_fixed_mixratio(num_fixed_gases)                   ! Fixed gases have their fixed values here.  This is the analog to mixratio.
  
  REAL*8 ::                                  &
     start,                                  & ! Timing variables
     finish,                                 &
     failed_wvn

  !G added for specifying channel widths
  REAL*8 ::                                  &
     spec_full_width,                        &
     spec_half_width
  !G new logical switch :: if spec_centers == T, then we have specified spectral band CENTERS; if F, then we have specified ENTRIES!
  LOGICAL ::                                 &
     spec_centers

  !G replacing hard-coded margin within which to search for nearby lines whose wings may affect a given band's absorption
  REAL*8 ::                                  &
     margin(max_num_gases),                  &
     margin_list(5),                         &
     vstep_min

  REAL*8, ALLOCATABLE ::                     &
     wn_arr(:),                              & ! Array to store the continuum input wavelengths
     temp_arr(:),                            & ! Array to store the continuum input temperatures
     dim_self_arr(:,:),                      & ! Array to store the self continuum input data
     dim_frn_arr(:,:),                       & ! Array to store the foreign continuum input data
     kvals_single_gas(:),                    & ! Array to store the k values for a single gas, to get summed in variable kvals
     k_self_cont(:),                         & ! Array to store the self continuum values for a single gas.
     k_frn_cont(:),                          & ! Array to store the foreign continuum values for a single gas
     combined_continuum(:),                  & ! Array to store the combined continuum terms, scaled for the abundance of the gas.
     kvals(:),                               & ! Array to store net absorption coefficient values across the spectral range
     line_dat(:,:,:,:)                         ! Array to store basic line data read from file.

  CHARACTER ::                               & 
     hit_data_dir*80,                        & ! The name of the HITRAN database
     hit_data_file*80,                       & ! The name of the HITRAN file
     gk_filename*80,                         & ! Array of filenames which store precise g(k) data (60 files maximum)
     output_dir*80,                          & ! The absolute directory where the output files are going
     op_file*80,                             & ! The full filename, including suffix, of the output file
     op_margin_file*80,                      & ! The full filename, including suffix, of the margin output file
     gas_indices_names*80,                      & ! The gas indices names used in the output file name, from namelist
     mixratiolabel(num_mix_ratio)*8,         & ! A means of labeling the output files with mixing ratio amount
     wnlabel(spectral_bands+1)*7,            & ! A means of labeling the output files with wavelength value
     temp_char*3,                            &
     press_char*2,                           &
     form1*32,                               &
     temp_label*7,                           &
     gas_label*2,                            &
     zero_pad*1

  LOGICAL                                     &
     anylines,	      		              & ! If there are no lines of any gas of interest in the spectral bin, we need to
						!    skip over that spectral bin.  No need to sort through null data.
     margin_for_new_pressure,                 & ! Flag to reset the margin index for a new pressure loop
     too_long                                   ! Flag to indicate that for a gas, time for calculations is too long--reduce margin

   !-----------------------------------------------------------------------------------------
   ! Read in a namelist from the kdm.nml file.
   namelist /input/ gas_index, gas_indices_names, spec_index, fixed_gas_fixed_mixratio
   open(10,file='kdm.nml')
   read(10,nml=input)
   close(10)
   !-----------------------------------------------------------------------------------------

  DATA number_of_isotopes / 4, 8, 5, 5, 6, 3, 3, 3, 2, 1, 2, 1, 3, 1, 2, 2, 1, 2,  &
                            4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 1, 1, 3, 1, 1, 1, 2, 1 /

!c-mm We may choose to ignore most of the more trivial isotopologues.  To avoid looping over these isotopes (in the y loop)
!c-mm    unnecessarily (and to avoid calculating continuum for these non-existant lines), we set the number of active isotopes
!c-mm    to a value reflecting the number of isotopes we want to work on (for CO2 and H2O, this is 1.  For SO2, this is 2.)
  DATA number_of_active_isotopes / 1, 1, 5, 5, 6, 3, 3, 3, 2, 1, 2, 1, 3, 1, 2, 2, 1, 2,  &
                            4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 1, 1, 3, 1, 1, 1, 2, 1 /

  !c-mm Specify pressures here and make sure there are as many as "p_num" specified above!
!c-mm  DATA press/ 0.000100, 0.000158, 0.000251, 0.000398, 0.000631, 0.001000, 0.001580, 0.002510, 0.003980, 0.006310, 0.010000, &
!c-mm              0.015800, 0.025100, 0.039800, 0.063100, 0.100000, 0.158000, 0.251000, 0.398000, 0.631000, 1.000000, 1.580000, &
!c-mm              2.510000, 3.980000, 6.310000, 10.00000, 15.80000, 25.10000, 39.80000, 63.10000, 100.0000, 158.0000, 251.0000, &
!c-mm              398.0000, 631.0000, 1000.000, 1580.000, 2510.000, 3980.000, 6310.000, 10000.00, 15800.00, 25100.00, 39800.00, &
!c-mm	      63100.00, 100000.0, 158000.0, 251000.0, 398000.0, 631000.0, 1000000. /     

  DATA press/ 0.000100, 0.000251, 0.000631, 0.001580, 0.003980, 0.010000, &
              0.025100, 0.063100, 0.158000, 0.398000, 1.000000,  &
              2.510000, 6.310000, 15.80000, 39.80000, 100.0000, 251.0000, &
              631.0000, 1580.000, 3980.000, 10000.00, 25100.00, &
	      63100.00, 158000.0, 398000.0, 1000000. /     

  !G Can either specify spectral band edges or centers.  
  !G   + If edges, then there are spectral_bands+1 entries in wn and the bands WILL NOT overlap
  !G   + If centers, then one must also specify a band width.  There will be only spectral_bands entries in wn and the bands CAN overlap. 
!c-mm  DATA wn/ 10.0D0, 166.667D0, 416.667D0, 625.370D0, 710.020D0, 833.333D0, 1250.0D0, &
!c-mm           2222.22D0, 3087.37D0, 4030.63D0, 5370.57D0, 7651.11D0, 12500.0D0, 25000.0D0, 41666.67D0 /
  DATA wn/ 10.0D0, 166.667D0, 416.667D0, 625.370D0, 710.020D0, 833.333D0, 1250.0D0, &
           2222.22D0, 3014.00D0, 3020.00D0, 5370.57D0, 7651.11D0, 12500.0D0, 25000.0D0, 41666.67D0 /

!c-mm  DATA variable_gas_mixratio/ 0.0D0, 1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2, 1.0D-1 /
  DATA variable_gas_mixratio/ 1.0D-25, 1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2, 1.0D-1 /
  !G here we are only using zero water
  !c-mm  DATA variable_gas_mixratio/ 0.0D-7, 0.0D-6, 0.0D-5, 0.0D-4, 0.0D-3, 0.0D-2 /

!c-mm  DATA mixratiolabel/ '0.00D-07', '1.00D-07', '1.00D-06', '1.00D-05', '1.00D-04', '1.00D-03', '1.00D-02', '1.00D-01' /
  DATA mixratiolabel/ '1.00D-25', '1.00D-07', '1.00D-06', '1.00D-05', '1.00D-04', '1.00D-03', '1.00D-02', '1.00D-01' /
  !G here we are only using zero water
  !c-mm  DATA mixratiolabel/ '0.00D-07', '0.00D-06', '0.00D-05', '0.00D-04', '0.00D-03', '0.00D-02' /

  !G These indices are set for CO2 then CH4, then H2O
  !DATA gas_index/ 2, 6, 1 /      ! For CO2 + CH4 + H2O atmosphere
  !Soto 2/11/2020: gas_index is now read in from a namelist file.

  !G CO2 mixing ratio is set
  !DATA fixed_gas_fixed_mixratio/ 0.9532D0, 1.00D-9 /   !  Just the fixed gases
  !Soto 2/11/2020: fixed_gas_fixed_mixratio is now read in from a namelist file.

  !c-mm Ideally, we would like to go out at least 500 cm-1 from each line, but this is computationally impractical.
  !c-mm    We set up a list of potential margin widths (500 cm-1, 250 cm-1, 100 cm-1, 50 cm-1, 25 cm-1) and, starting
  !c-mm    from 500 cm-1, perform the calculations.  If the amount of time it takes to calculate the k-coefficient
  !c-mm    for a single p,T exceeds some chosen threshold, we bail, and try again with the next value on the list.
  !c-mm    In practice, for Mars, I've determined the viable margins to use for each band and each pressure, so I've
  !c-mm    hardwired these values (You'll see places with a 'margin_list(5)' hardwired.
  DATA margin_list/ 5.0D4, 2.5D4, 1.0D4, 5.0D3, 2.5D3 /   ! Margin steps we will potentially apply

  !G specify TES spectral channel widths
  spec_centers = .FALSE.
  spec_full_width = 10.58D0    ! width in cm^-1
  spec_half_width = spec_full_width / 2.0D0

!c-mm  If you want to hardwire these values without a runscript, uncomment these lines and comment out the NAMELIST line
  p_start_index = 1
  p_end_index = 26
  t_start_index = 1
  t_end_index = 18
!c-Soto  spec_index = 9
  mix_ratio_start_index = 1
  mix_ratio_end_index = 8

!c-mm  NAMELIST/indata/p_start_index,p_end_index,t_start_index,t_end_index,spec_index,mix_ratio_start_index,mix_ratio_end_index
!c-mm  READ(5,indata)

  output_dir = './out/'
  length_count=LEN_TRIM(output_dir)

  !G prepare for formatting below -- fill out strings with numbers (Fortran's equivalent of num2str)
  WRITE(temp_char,'(I3)') t_num
  form1='('//TRIM(temp_char)//'e14.7)'  ! This generates a format statement of 71e14.7 (for t_num=71)
  WRITE(zero_pad,'(I1)') 0

  WRITE(press_char,'(I2)') p_start_index
  IF (p_start_index < 10) press_char=zero_pad//TRIM(ADJUSTL(press_char))
  

  !G construct wavenumber label array based on whether *CENTERS* or *EDGES*
  IF ( spec_centers ) THEN
     DO i = 1, spectral_bands
        IF (wn(i) < 1.0D2) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//zero_pad//zero_pad//temp_label
        ELSE IF (wn(i) >= 1.0D2 .AND. wn(i) < 1.0D3) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//zero_pad//temp_label
        ELSE IF (wn(i) >= 1.0D3 .AND. wn(i) < 1.0D4) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//temp_label
        ELSE IF (wn(i) >= 1.0D4) THEN
           WRITE(temp_label,'(F7.1)') wn(i)
           wnlabel(i) = temp_label
        ENDIF
        WRITE(0,*) 'wnlabel= ',wnlabel(i)
     ENDDO
  ELSE
     DO i = 1, spectral_bands+1
        IF ( wn(i) < 1.0D2 ) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//zero_pad//zero_pad//TRIM(ADJUSTL(temp_label))
        ELSE IF (wn(i) >= 1.0D2 .AND. wn(i) < 1.0D3) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//zero_pad//TRIM(ADJUSTL(temp_label))
        ELSE IF (wn(i) >= 1.0D3 .AND. wn(i) < 1.0D4) THEN
           WRITE(temp_label,'(F7.2)') wn(i)
           wnlabel(i) = zero_pad//TRIM(ADJUSTL(temp_label))
        ELSE IF (wn(i) >= 1.0D4) THEN
           WRITE(temp_label,'(F7.1)') wn(i)
           wnlabel(i) = temp_label
        ENDIF
        WRITE(0,*) 'wnlabel= ',wnlabel(i)
     ENDDO
  ENDIF
  
  !G This is pre-processing step to accessing the full HITRAN database.  read_gas_data requires the file "hitran_gas_data"
  CALL read_gas_data(q_coef,q_ref,iso_mass,iso_fraction)

  !c-mm  Only need to read in CIA and CO2 continuum data if CO2 is one of the gases (which it should always be for Mars)
  IF (ANY(gas_index == 2)) THEN
     !c-mm This is a pre-processing step to calculating CIA
     CALL read_baranov_data(wn_CIA_arr,temp_CIA_arr,dim_CIA_arr)
     !c-mm This is a pre-processing step to read in CO2 continuum using MT_CKD data rather than CIA data. 
     !c-mm    Should only choose one or the other.  Based on Wordsworth et al., (2010), I choose the CIA data
     !c-mm     CALL read_co2_continuum_data(wn_co2_arr,temp_co2_arr,dim_co2_arr)
  ENDIF

  !c-mm This is a pre-processing step to read in water vapor self and foreign continuum.  Only do so if the particular
  !c-mm    execution of the code includes water vapor as a variable gas AND it has a non-zero mixing ratio.
  IF (ANY(gas_index == 1).AND.ANY(variable_gas_mixratio(mix_ratio_start_index:mix_ratio_end_index,1) > 0.0D0)) THEN
     CALL read_h2o_continuum_data(wn_h2o_arr,temp_h2o_arr,dim_h2o_self_arr,dim_h2o_frn_arr)
  ENDIF

  !c-mm Specify the HITRAN data file directory root
  hit_data_dir='./ck/'

  press = press/1.0D+05                 ! Pressure converted into bars (from Pa)

  !-----------------------------------------------------------------------------------------
  !c-mm BIG LOOP generating k-tables

  !c-mm "jj" is index over mixing ratios
  DO jj=mix_ratio_start_index,mix_ratio_end_index			 ! Here, I am writing to mixratio_list only those values

     DO i=1,num_fixed_gases						 !    for the particular jj loop that we are in.  This
        mixratio_list(i)=fixed_gas_fixed_mixratio(i)		         !    includes the appropriate index for the variable
     ENDDO								 !    gas from variable_gas_mixratio, as well as the

     DO i=1,num_variable_gases					         !    fixed values from fixed_gas_fixed_mixratio.
        mixratio_list(num_fixed_gases+i)=variable_gas_mixratio(jj,i)     ! This is used later on wherever the gas mixing ratios
     ENDDO								 !    are required (In two spots below).

     IF (SUM(mixratio_list) > 1.0) mixratio_list(1)=mixratio_list(1)-(SUM(mixratio_list)-1.0)  ! Reduces CO2 mixing ratio if there is
                                                                                               !   an overabundance of other species.

     WRITE(0,*), 'jj = ', jj, ' out of ', (mix_ratio_end_index - mix_ratio_start_index)+1, ' mixing ratios'

     DO kk = spec_index, spec_index
        raw_data(:,:)%num_rawdat=0
        !c-mm  Constructing temperature indices
        DO i=t_start_index,t_end_index
           layer_temp(i) = t_init + REAL(i-1)*delta_t
        ENDDO
        
        !G If-structure to handle filenames differently whether specifying spectral bin CENTERS or EDGES
        IF ( spec_centers ) THEN
           !gk_filename=output_dir(1:length_count)//'h2o+ch4+co2.'//mixratiolabel(jj)//'.'//TRIM(wnlabel(kk))
           gk_filename=output_dir(1:length_count)//TRIM(gas_indices_names)//'.'//mixratiolabel(jj)//'.'//TRIM(wnlabel(kk))
        ELSE
           !gk_filename=output_dir(1:length_count)//'h2o+ch4+co2.'//mixratiolabel(jj)//'.'//TRIM(wnlabel(kk))//'-'//TRIM(wnlabel(kk+1))
           gk_filename=output_dir(1:length_count)//TRIM(gas_indices_names)//'.'//mixratiolabel(jj)//'.'//TRIM(wnlabel(kk))//'-'//TRIM(wnlabel(kk+1))
        ENDIF

        !G if CENTERS, then +/- spec_half_width
        IF ( spec_centers ) THEN
           range_min = wn(kk) - spec_half_width
           range_max = wn(kk) + spec_half_width
           !c-mm           print*, 'kk = ', kk, ' /', spectral_bands, ' :: Center = ',wn(kk)
           WRITE(0,*), 'kk = ', kk, ' /', spec_index, ' :: Center = ',wn(kk)
        !G if EDGES, then assign range to edges
        ELSE
           range_min = wn(kk)
           range_max = wn(kk+1)
           !c-mm           print*, 'kk = ', kk, ' /', spectral_bands, ' ::  Edges = ',wn(kk), wn(kk+1)
           WRITE(0,*), 'kk = ', kk, ' /', spec_index, ' ::  Edges = ',wn(kk), wn(kk+1)
        ENDIF
        range_min = 100.0D0*range_min                                    ! Convert spectral range to S.I. units
        range_max = 100.0D0*range_max

        anylines = .FALSE.
        
        !c-mm  Block for reading in spectral lines.  Code takes the gas index number and opens up the
        !c-mm     correct single-gas data file (e.g. 01_hit08.par) for each gas.
        num_gases = 0
        DO x=1, max_num_gases
           IF ((gas_index(x) == 1).AND.(variable_gas_mixratio(jj,1) == 0.0D0)) CYCLE       ! If this is water, but mr=0, CYCLE to next gas
           IF (gas_index(x) == 0) then
              CYCLE ! skip blank entries
           ELSE
              num_gases=num_gases+1
           ENDIF
              
           WRITE(gas_label,'(I2)') gas_index(x)
           IF (gas_index(x) < 10) THEN
              hit_data_file=TRIM(ADJUSTL(hit_data_dir))//zero_pad//TRIM(ADJUSTL(gas_label))//'_hit08.par'
           ELSE
              hit_data_file=TRIM(ADJUSTL(hit_data_dir))//TRIM(ADJUSTL(gas_label))//'_hit08.par'
           ENDIF
           write(0,*) 'filename= ',hit_data_file
           OPEN(UNIT=10,FILE=hit_data_file)
           SELECT CASE (gas_index(x))                                                    ! Water and methane have different margin than other gases
              CASE(1)
                 margin(x)=margin_list(5)
              CASE(6)                       
                 margin(x)=margin_list(5)   
              CASE DEFAULT
                 margin(x)=margin_list(1)   !c-mm  This is the line we change to fix the margin at a set value, along with the timing line (7200)
           END SELECT                       !c-mm     and starting_margin_index.  See notes above for margin_list
           DO
!HITRAN2000 I2 , I1 , F12.6 , F6.3  , 1X    ,I3  , F6.3,1X,I3, F5.4,F5.4,   F10.4  , F4.2,F8.6,I3,I3,A9,A9,3I1,3I2
!NAMES      id   iso  cent   [ strength         ]  [ temp1   ] abhw, shbw low_ener  
!HITRAN2K4  I2 , I1 , f12.6 , e10.3 , e10.4 , f5.4, f5.4, f10.4, f4.2
              READ(10,'(I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2)',IOSTAT=instat)      &  ! Read in all necessary data from HITRAN file
                   mol_no, iso_no, line_cen, line_stren, temp1, abhw, sbhw, low_ener, &  
                   abhw_coeff
              line_cen = line_cen*100.0D0                                                ! Convert spectral range to S.I. units
              IF (line_cen < (range_min - margin(x))) CYCLE                              ! If we are lower than our range, CYCLE
              IF (line_cen > (range_max + margin(x))) EXIT                               ! If we are higher than our range, EXIT to next gas
              !c-mm  For speed sake, restrict H2O and CO2 to first isotope only.  This covers 99.7% of water lines and 98.4% of CO2 lines
              !c-mm     A slight error is introduced since we are not counting 100% of the lines, but what'cha gonna do?
              !c-mm  For SO2 and H2S, just do all isotopes since there are many fewer lines.
              IF (iso_no > number_of_active_isotopes(gas_index(x))) CYCLE
              anylines=.TRUE.                                                            ! Lines exist in this spectral band, so TRUE
              raw_data(x,iso_no)%num_rawdat=raw_data(x,iso_no)%num_rawdat+1	         ! Read in relevant data
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,1)=line_cen	 !          |
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,2)=line_stren	 !          |
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,3)=abhw		 !          |
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,4)=low_ener	 !          |
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,5)=abhw_coeff	 !          |
              raw_data(x,iso_no)%rawdat1(raw_data(x,iso_no)%num_rawdat,6)=sbhw		 !          |
              raw_data(x,iso_no)%rawdat2(raw_data(x,iso_no)%num_rawdat)=iso_no		 !          V
           ENDDO
           CLOSE(10)
        ENDDO

        !c-mm get ready to output the data
!c-mm	op_file = TRIM(gk_filename)//'_p'//press_char//'_t'//temp_char//'.dat'
	op_file = TRIM(gk_filename)//'_p'//press_char//'.dat'
	op_margin_file = TRIM(gk_filename)//'_p'//press_char//'_margin.dat'
	OPEN(UNIT=20+jj-1,FILE=op_file)	                                         ! Allow for different files for different mixing
        OPEN(UNIT=(20+jj-1)*10,FILE=op_margin_file)
        !c-mm "i" is the index over pressure
        DO i=p_start_index,p_end_index
           margin_for_new_pressure=.TRUE.
           !c-mm "j" is the index over temperature
!c-mm           DO j=1,t_num
           DO j=t_start_index,t_end_index
              IF (j == 1) write(0,*) 'pressure index= ',i
              temp=layer_temp(j)
              pratio=press(i)
              vstep = (range_max-range_min)
              !c-mm At this point, we know exactly how many lines, so we can allocate exact space to line_dat.  Last 8 is for the 8 CO2 isotopes
              ALLOCATE(line_dat(4,MAXVAL(raw_data%num_rawdat),num_gases,8))

              !c-mm  This block calculates all the line data for all the gases in the atmosphere and determines the narrowest line from which we
              !c-mm     can establish the stepping in wavenumber space.
              DO x=1,num_gases
                 DO y=1,number_of_active_isotopes(gas_index(x))
                    IF (raw_data(x,y)%num_rawdat > 0) THEN
                       !c-mm Call abs_line_data to generate line_dat for each gas in the list, only if there are lines in the chosen band.
                       CALL abs_line_data(line_dat(:,:,x,y), pratio, temp, raw_data(x,y), q_coef(gas_index(x),:,y), &
                                          q_ref(gas_index(x),y), iso_mass(gas_index(x),y),                          &
                                          mixratio_list(x)*iso_fraction(gas_index(x),y))
                       !c-mm  Choose vstep_min to be 20% of the Voigt line halfwidth.  This expression calculates the Voigt halfwidth with great
                       !c-mm     precision.  Taken from Olivero and Longbothum, JQSRT 17(2), 233.
                       vstep_min=0.2D0*MINVAL(0.5D0*(0.5346D0*(2.0D0*line_dat(1,1:raw_data(x,y)%num_rawdat,x,y)*line_dat(4,1:raw_data(x,y)%num_rawdat,x,y))+ &
                                 SQRT(0.2166*4.0D0*(line_dat(1,1:raw_data(x,y)%num_rawdat,x,y)*line_dat(4,1:raw_data(x,y)%num_rawdat,x,y))**2+ &
                                 4.0D0*line_dat(3,1:raw_data(x,y)%num_rawdat,x,y)**2)))
                       IF (vstep_min < vstep) vstep = vstep_min
                    ENDIF
                 ENDDO ! y
              ENDDO ! x

              WRITE(0,*) 'Working on temperature: ',layer_temp(j)
              WRITE(0,*) 'vstep is: ',vstep/100.,'cm-1'
              !c-mm Now check to make sure that vstep allows for at least 1.D4 steps
              IF ( 1.0D+04*vstep > (range_max-range_min)) THEN
                 vstep = (range_max-range_min)/1.00D+04
                 write(0,*), 'Adjusting vstep (1d4) for temp ',layer_temp(j),' to ',vstep/100.,'cm-1'
              ENDIF
              !c-mm  Occasionally will experience an out of memory error if vstep is too small.  This is a problem for the very low frequency bands.
              !c-mm     The value 2.0D-04 is about the limit that many systems can handle considering the number of lines and the band width.
              IF (vstep < 5.0D-03) THEN
                 vstep = 5.0D-03
                 write(0,*) 'Adjusting vstep minimum for temp ',layer_temp(j),' to ',vstep/100.,'cm-1'
              ENDIF
              CALL flush(0)
              !c-mm make room for k-value array -- based on vstep
              kval_dim = CEILING((range_max-range_min)/vstep)                                     ! kval_dim is the size of the kval array 
              ALLOCATE(kvals(kval_dim))
              ALLOCATE(kvals_single_gas(kval_dim))
              kvals=0.0
              !c-mm  k_self_cont and k_frn_cont are only used for water vapor (and CO2, but only if MT_CKD continuum is used)
              !c-mm This loop calls abs_coeff_calc for each gas individually, adding continuum terms where appropriate.
              DO x=1, num_gases
                 IF (margin_for_new_pressure) THEN
                    SELECT CASE (gas_index(x))                                                    ! We are automatically selecting incrementally
                       CASE(1)                                                                    !    smaller margins if it takes too long with
                          margin_index(x)=starting_h2o_margin_index                               !    the current margin.  H2O is initially at
                       CASE DEFAULT                                                               !    25 cm-1, which should be good.  CO2 starts
                          margin_index(x)=starting_margin_index                                   !    at 500 cm-1, which may be too long.  Reset
                    END SELECT                                                                    !    the value for each new pressure index.
                 ENDIF
                 margin(x)=margin_list(margin_index(x))
                 too_long=.TRUE.
                 DO WHILE (too_long)
                    DO y=1,number_of_active_isotopes(gas_index(x))
                       kvals_single_gas=1.D-25                                                    ! Start fresh for each gas/continuum
                       IF (raw_data(x,y)%num_rawdat > 0) THEN
                          CALL cpu_time(start)
                          !c-mm Two separate routines to generating the absorption coefficients.  Both yield basically the same results, although I
                          !c-mm    think the humlik routine is slightly faster.  You can use either one.
!                          CALL abs_coeff_calc(line_dat=line_dat(:,:,x,y), range_min=range_min, range_max=range_max,      & 
!                                              gas=gas_index(x), kvals=kvals_single_gas,vstep=vstep,                      &
!                                              numvals=numvals, numlines=raw_data(x,y)%num_rawdat,temp=temp,maxwid=margin(x), &
!                                              too_long=too_long,failed_wvn=failed_wvn)
                          CALL humlik(line_dat=line_dat(:,:,x,y), range_min=range_min, range_max=range_max,              & 
                                              gas=gas_index(x), kvals=kvals_single_gas,vstep=vstep,                      &
                                              numvals=numvals, numlines=raw_data(x,y)%num_rawdat,temp=temp,maxwid=margin(x), &
                                              too_long=too_long,failed_wvn=failed_wvn,kval_dim=kval_dim)
                          CALL cpu_time(finish)
                          WRITE((20+jj-1)*10,'(a3,f5.1,a9,f5.1,a13,f7.1,a2)') 't= ',temp,' margin= ',margin(x)/100.,' cm-1, time= ',finish-start,' s'
                          CALL flush((20+jj-1)*10)
                          !c-mm  The variable too_long will become true inside abs_coeff_calc/humlik if the cpu_time exceeds some threshold.  It will dump out
                          !c-mm     and return to here where the code will tell you how far along in spectral space the code got before reaching the threshold.
                          IF (too_long) THEN
                             WRITE((20+jj-1)*10,'(a15,f6.3,a24)') 'Only completed ',100.*(failed_wvn-range_min)/(range_max-range_min),'% of band before timeout'
                             CALL flush((20+jj-1)*10)
                             IF (margin_index(x)+1 > 5) THEN
                                WRITE(0,*) 'Cannot get a suitable band margin'
                                STOP
                             ELSE
                                 margin_index(x)=margin_index(x)+1
                                 margin(x)=margin_list(margin_index(x))             !c-mm Automatically adjust next margin if it exceeds threshold...
                             ENDIF
                             EXIT                                                   !c-mm  If we are too long for any of the isotopes, start over, exit y
                          ENDIF
                       ELSE
                          too_long=.FALSE.                                          !c-mm  If there are no lines, we need to force too_long to be false
                       ENDIF
                       !c-mm  Choose how to do the continuum here
                       SELECT CASE (gas_index(x))  ! Both H2O (1) and CO2 (2) have various additions to the line spectrum, so handle specially.
                          CASE (1)
                             !c-mm Note here that we rely on the fact that the parsing block above will not read in any water lines if mixing ratio
                             !c-mm    is set to zero.  That makes num_rawdat=0, preventing execution of this statement for range_min > 1000000 m-1.
                             !c-mm Water vapor continuum data goes out to 10000 cm-1 (1000000 m-1)
                             IF ((raw_data(x,y)%num_rawdat > 0).OR.((range_min < 1000000.).AND.(variable_gas_mixratio(jj,1) > 0.0D0))) THEN  
                                ALLOCATE(k_self_cont(kval_dim))
                                ALLOCATE(k_frn_cont(kval_dim))
                                ALLOCATE(combined_continuum(kval_dim))
                                k_self_cont=0.0D0
                                k_frn_cont=0.0D0
                                combined_continuum=0.0D0
                                !c-mm It seems silly to define these variables just to hold other, perfectly good variables, but if we ever
                                !c-mm    go back to using MT_CKD values for CO2, this will simplify things greatly.
                                ALLOCATE(wn_arr(nS_h2o))
                                ALLOCATE(temp_arr(nT_h2o))
                                ALLOCATE(dim_self_arr(nS_h2o,nT_h2o))
                                ALLOCATE(dim_frn_arr(nS_h2o,nT_h2o))
                                wn_arr=wn_h2o_arr
                                temp_arr=temp_h2o_arr
                                dim_self_arr=dim_h2o_self_arr
                                dim_frn_arr=dim_h2o_frn_arr

                                CALL continuum_calc(range_min=range_min,range_max=range_max,gas=gas_index(x),vstep=vstep, &
                                                    numvals=numvals,temp=temp,pratio=pratio,wn_arr=wn_arr,                &
                                                    temp_arr=temp_arr,dim_self_arr=dim_self_arr,dim_frn_arr=dim_frn_arr,  &
                                                    k_self_cont=k_self_cont,k_frn_cont=k_frn_cont,kval_dim=kval_dim)
                                anylines=.TRUE.
                                !c-mm The continuum term is the combination of self and foreign broadening.  The mixing ratio is what
                                !c-mm    dictates which of the two terms is dominant.
                                combined_continuum = k_self_cont*mixratio_list(x)*iso_fraction(gas_index(x),y) +          & !c-mm  Self-continuum
                                                     k_frn_cont*(1.-mixratio_list(x))*iso_fraction(gas_index(x),y)          !c-mm  Foreign continuum
!c-mm                                kvals = kvals + k_frn_cont + (kvals_single_gas + k_self_cont - k_frn_cont) *            &
!c-mm                                        mixratio_list(x)*iso_fraction(gas_index(x),y)
                                !c-mm  Absorption is scaled by the mixing ratio of the gas 'x' under consideration.  Note that the continuum
                                !c-mm     is being doubly multiplied by the mixing ratio (here and above in combined_continuum), but for
                                !c-mm     different purposes.  This ensures that the continuum term will go to zero as the mixing ratio
                                !c-mm     decreases.
                                kvals = kvals + (kvals_single_gas + combined_continuum)*mixratio_list(x)*iso_fraction(gas_index(x),y)
!c-mm                                kvals = kvals + kvals_single_gas*mixratio_list(x)*iso_fraction(gas_index(x),y)
                                DEALLOCATE(wn_arr)
                                DEALLOCATE(temp_arr)
                                DEALLOCATE(dim_self_arr)
                                DEALLOCATE(dim_frn_arr)
                                DEALLOCATE(k_self_cont)
                                DEALLOCATE(k_frn_cont)
                                DEALLOCATE(combined_continuum)
                             ENDIF
                          CASE (2)
!c-mm Do not do MT_CKD continuum if doing CIA continuum.  Comment this out, but preserve for later use.
!c-mm                          IF ((raw_data(x,y)%num_rawdat > 0).OR.(range_min < 1000000.)) THEN  ! MT_CKD Continuum data goes out to 10000 cm-1
!c-mm                             ALLOCATE(k_self_cont(kval_dim))
!c-mm                             ALLOCATE(k_frn_cont(kval_dim))
!c-mm                             ALLOCATE(wn_arr(nS_co2))
!c-mm                             ALLOCATE(temp_arr(nT_co2))
!c-mm                             ALLOCATE(dim_self_arr(nS_co2,nT_co2))
!c-mm                             wn_arr=wn_co2_arr
!c-mm                             temp_arr=temp_co2_arr
!c-mm                             dim_self_arr=dim_co2_arr
!c-mm                             CALL continuum_calc(range_min=range_min,range_max=range_max,gas=gas_index(x),vstep=vstep, &
!c-mm                                                 numvals=numvals,temp=temp,pratio=pratio,wn_arr=wn_arr,temp_arr=temp_arr,            &
!c-mm                                                 dim_self_arr=dim_self_arr,k_self_cont=k_self_cont,k_frn_cont=k_frn_cont,kval_dim=kval_dim)
!c-mm                             anylines=.TRUE.
!c-mm                             kvals = kvals + (kvals_single_gas + k_self_cont) * mixratio_list(x)*iso_fraction(gas_index(x),y)
!c-mm                             DEALLOCATE(wn_arr)
!c-mm                             DEALLOCATE(temp_arr)
!c-mm                             DEALLOCATE(dim_self_arr)
!c-mm                             DEALLOCATE(k_self_cont)
!c-mm                             DEALLOCATE(k_frn_cont)

                             IF ((raw_data(x,y)%num_rawdat > 0).OR.((range_min <= 25000.).AND.(range_min >= 0.)).OR.((range_max <= 25000.).AND.             &
                                  (range_max >= 0.)).OR.((range_max >= 25000.).AND.(range_min <= 0.)).OR.                 &
                                  ((range_min <= 200000.).AND.(range_min >=70000.)).OR.((range_max <= 200000.).AND.       &
                                  (range_max >= 70000.)).OR.((range_max >= 200000).AND.(range_min <= 70000.))) THEN         ! CIA Continuum data goes from 700-2000 cm-1
                                ALLOCATE(k_self_cont(kval_dim))
                                ALLOCATE(k_frn_cont(kval_dim))
                                k_self_cont=0.0D0
                                k_frn_cont=0.0D0
                                CALL continuum_calc(range_min=range_min,range_max=range_max,gas=gas_index(x),vstep=vstep, &
                                                    numvals=numvals,temp=temp,pratio=pratio,wn_CIA_arr=wn_CIA_arr,        &
                                                    temp_CIA_arr=temp_CIA_arr,dim_CIA_arr=dim_CIA_arr,                    &
                                                    k_self_cont=k_self_cont,k_frn_cont=k_frn_cont,kval_dim=kval_dim)
                                anylines=.TRUE.
                                kvals = kvals + (kvals_single_gas + k_self_cont) * mixratio_list(x)*iso_fraction(gas_index(x),y)
                                DEALLOCATE(k_self_cont)
                                DEALLOCATE(k_frn_cont)
                             ENDIF
                          CASE DEFAULT                             
			     IF (raw_data(x,y)%num_rawdat > 0) THEN
                                anylines=.TRUE.
                                kvals = kvals + kvals_single_gas * mixratio_list(x)*iso_fraction(gas_index(x),y)
			     ENDIF
                       END SELECT
                    ENDDO   ! y
                 ENDDO   ! While
              ENDDO   ! x
              margin_for_new_pressure=.FALSE.
              DEALLOCATE(line_dat)
              DEALLOCATE(kvals_single_gas)
              !c-mm  Here's where we call the binning and sorting routines
              IF (anylines) THEN							  ! If there are no gases of interest in
                 numvals=MIN(numvals,kval_dim)
                 CALL gk_calc(kvals,numvals,gkdat,kk,num_terms)                           !    this band, then no need to sort, as
                 kvals_at_one_p(j,i,:) = gkdat    					  !    abs_coeff_calc has been completely
              ELSE									  !    skipped.  Instead, just bypass
                 kvals_at_one_p(j,i,:) = 1.0D-25					  !    sorting and set kvals_at_one_p to
              ENDIF									  !    1.0D-25 everywhere.
              DEALLOCATE(kvals)							    
              CALL flush(20+jj-1)
	     
           ENDDO    ! j-loop (temperature)
        ENDDO       ! i-loop (pressure)

        DO k=1,num_terms
           DO i=p_start_index,p_end_index
              WRITE(20+jj-1,form1) kvals_at_one_p(:,i,k)
           ENDDO							
        ENDDO								
        CALL flush(20+jj-1)
        CLOSE(20+jj-1)								   !    processing.  Cannot use a namelist variable
        
     ENDDO ! end kk-loop over spectral bands

  ENDDO ! end jj-loop over h2o mixing ratios

END PROGRAM ck_calc


!---------------------------------------------------------------------------------------------------------------------------------
! 2.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to read in basic HITRAN gas data
!c-mm
!c-mm This subroutine reads in the necessary gas data from the file hitran_gas_data
!c-mm List of gases (and number of isotopes) is listed below
!c-mm
!c-mm Gas #1 = H2O  (4)	 Gas #11 = NH3	(2)	Gas #21 = HOCl  (2)	Gas #31 = H2S    (3)
!c-mm Gas #2 = CO2  (8)	 Gas #12 = HNO3	(1)	Gas #22 = N2	(1)	Gas #32 = HCOOH  (1)
!c-mm Gas #3 = O3   (5)	 Gas #13 = OH	(3)	Gas #23 = HCN	(3)	Gas #33 = HO2    (1)
!c-mm Gas #4 = N2O  (5)	 Gas #14 = HF	(1)	Gas #24 = CH3Cl	(2)	Gas #34 = O      (1)
!c-mm Gas #5 = CO   (6)	 Gas #15 = HCl	(2)	Gas #25 = H2O2	(1)	Gas #35 = ClONO2 (2)
!c-mm Gas #6 = CH4  (3)	 Gas #16 = HBr	(2)	Gas #26 = C2H2	(2)	Gas #36 = NO+    (1)
!c-mm Gas #7 = O2   (3)	 Gas #17 = HI	(1)	Gas #27 = C2H6  (1)
!c-mm Gas #8 = NO   (3)	 Gas #18 = ClO	(2)	Gas #28 = PH3   (1)
!c-mm Gas #9 = SO2  (2)	 Gas #19 = OCS	(4)	Gas #29 = COF2  (1)
!c-mm Gas #10 = NO2 (1)	 Gas #20 = H2CO	(3)	Gas #30 = SF6   (1)
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE read_gas_data(q_coef, q_ref, iso_mass, iso_fraction)

  IMPLICIT NONE

  INTEGER, PARAMETER ::                         &
     number_hitran_gases = 36

  INTEGER                                       &
     i, j, k,                                   &
     number_of_isotopes(number_hitran_gases)

  REAL*8, INTENT(  OUT) ::                      &
     q_coef(number_hitran_gases,4,8),           &
     q_ref(number_hitran_gases,8),              &
     iso_mass(number_hitran_gases,8),           &
     iso_fraction(number_hitran_gases,8)

  !-----------------------------------------------------------------------------------------
  DATA number_of_isotopes / 4, 8, 5, 5, 6, 3, 3, 3, 2, 1, 2, 1, 3, 1, 2, 2, 1, 2,  &
                            4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 1, 1, 3, 1, 1, 1, 2, 1 /

  OPEN(UNIT=10,FILE='./ck/hitran_gas_data_isotopes')                

  DO k=1,number_hitran_gases
     READ(10,*) ((q_coef(k,j,i), j=1,4), i=1,number_of_isotopes(k))
  ENDDO
  DO k=1,number_hitran_gases
     READ(10,*) (q_ref(k,i), i=1,number_of_isotopes(k))
  ENDDO
  DO k=1,number_hitran_gases
     READ(10,*) (iso_mass(k,i), i=1,number_of_isotopes(k))
  ENDDO
  DO k=1,number_hitran_gases
     READ(10,*) (iso_fraction(k,i), i=1,number_of_isotopes(k))
  ENDDO
  CLOSE(10)

END SUBROUTINE read_gas_data


!---------------------------------------------------------------------------------------------------------------------------------
! 3.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to read in basic HITRAN data
!c-mm
!c-mm This subroutine is input with the basic absorption line data from HITRAN, as stored in the rawdat1 and rawdat2 data arrays 
!c-mm for the whole band.  It outputs in the line_dat array the necessary spectral information required to calculate absorption 
!c-mm coefficients for the spectral range of the g(k) distribution.  Notice the spectral region has line data up to 25 cm-1 either 
!c-mm side of its limits, so as to allow for overlap effects of absorption lines outside the spectral region.  The partition sums
!c-mm have been extensively validated against the routines written by Gamache.  This polynomial fit is more than adequate.
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE abs_line_data(line_dat, pratio, temp, raw_data, q_coef, q_ref, iso_mass, mixratio)

  IMPLICIT NONE

  TYPE :: hitrandata_onegas
     REAL*8, DIMENSION(80000,6) ::    &
       rawdat1
     INTEGER, DIMENSION(80000)  ::    &
       rawdat2
     INTEGER                    ::    &
       num_rawdat
  END TYPE hitrandata_onegas
  
  TYPE(hitrandata_onegas) ::          &
     raw_data

  INTEGER, PARAMETER ::               &
     num_hitran_gases = 36

  REAL*8, PARAMETER ::                &
     ctwo = 1.438786D-02,             & ! Second radiation constant
     tref = 296.0D0,                  & ! Reference temperature from which data is calculated
     v_light = 3.0D+08,               & ! Speed of light
     rd = 8314.3D0                      ! Universal gas constant

  REAL*8, INTENT(IN   ) ::            &
     temp,                            & ! Temperature at which calculation is done
     pratio,                          & ! The pressure (/1.0D+05, meaning in bars) at which calculation is done
     mixratio,                        &
     iso_mass,                        & ! Mass in kg of one molecule
     q_coef(4),                       & ! Coefficients to calculate the internal partition sum
     q_ref                              ! Internal partition sum at 296 Kelvin

  REAL*8, INTENT(  OUT) ::            &
     line_dat(4,raw_data%num_rawdat)    ! Array of basic line data processed from HITRAN file

  ! Local variables:

  REAL*8                              &
     line_cen,                        & ! Wavenumber of line center (m-1)
     line_stren,                      & ! Line strength (input as cm-1/(molec.cm-2) , output as m/kg)
     abhw,                            & ! Air broadened halfwidth (cm-1/atm from HITRAN, output as m-1/atm)
     sbhw, 			      & ! Self-broadened halfwidth (cm-1/atm from HITRAN, output as m-1/atm)
     low_ener,                        & ! Lower state energy (cm-1 from HITRAN, output as m-1)
     abhw_coeff,                      & ! Coefficient of temperature dependence of air broadened halfwidth (no units)
     gas_amu,                         & ! Mass in kg of one molecule
     qt,                              & ! Internal partition sum at the temperature of the calculation (temp)
     stim_emis_fac,                   & ! Factor governing temperature dependence of line strength due to stimulated emission process
     tips, 			      & ! A factor governing temperature dependence of line strength due to TIPS.
     boltz_fac                          ! A factor governing temperature dependence of line strength due to the Boltzmann population 
                                        !    factors of energy levels
  INTEGER                             &
     i                                  ! Loop variable

  !-----------------------------------------------------------------------------------------
  line_cen = 0.0D0
  DO i = 1, raw_data%num_rawdat
     line_cen = raw_data%rawdat1(i,1)
     line_stren = raw_data%rawdat1(i,2)
     abhw = raw_data%rawdat1(i,3)
     low_ener = raw_data%rawdat1(i,4)
     abhw_coeff = raw_data%rawdat1(i,5)
     sbhw = raw_data%rawdat1(i,6)
     gas_amu = iso_mass*1.667D-27
     line_stren = line_stren/(gas_amu*100.0D0)                    ! Scale values to the appropriate temperature and pressure values
     abhw = abhw*100.0D0                                          ! In some cases HITRAN does not have the energy of the lower 
     sbhw = sbhw*100.0D0				          !    state of the transition, and sets its value to -1.
     low_ener = low_ener*100.0D0                                    
     abhw = ((tref/temp)**abhw_coeff)*(abhw*(pratio*(1.0D0 -    &
          mixratio))+sbhw*(pratio*mixratio))
     qt = q_coef(1) + q_coef(2)*                                & ! Now the internal partition sum at the layer temperature is
          temp + q_coef(3)*temp*temp +                          & !    calculated and stored in qt
          q_coef(4)*temp*temp*temp
     tips = q_ref/qt
     stim_emis_fac = (1.0D0 - dEXP(-ctwo*line_cen/temp))/       &
                       (1.0D0 - dEXP(-ctwo*line_cen/tref))
     IF (low_ener > 0) THEN
        boltz_fac = dEXP(ctwo*low_ener*(temp-tref)/(temp*tref))                   
     ELSE
        boltz_fac = 1.0D0                                         ! Boltzmann population factor =1 as no data available in HITRAN
     ENDIF
     line_stren = line_stren*boltz_fac*stim_emis_fac*tips
     line_dat(1,i) = line_cen                                     ! Write values to the output array
     line_dat(2,i) = line_stren
     line_dat(3,i) = abhw
     line_dat(4,i) = ((2.0D0*rd/iso_mass*temp)**                & ! Coefficient used to determine Doppler line width
                       0.5)/v_light
  ENDDO
END SUBROUTINE abs_line_data


!---------------------------------------------------------------------------------------------------------------------------------
! 4.
!-----------------------------------------------------------------------------------------------------------------------------------
!     A subroutine to return absorption coefficients over a spectral range assuming a Humlicek line shape
!
!     This subroutine calculates the absorption coefficient in the spectral range range_min to range_max, at intervals determined 
!     by vstep.  It uses the Humlicek approximation for the Voigt line shape for each absorption line.  Note that the value of 
!     pressure is no longer required as the line_dat array has already been processed to account for pressure effects.
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE abs_coeff_calc(line_dat, range_min, range_max, gas, kvals, vstep, numvals, numlines, temp, maxwid, too_long, failed_wvn)

  IMPLICIT NONE

  !c-mm Program computes the complex probability function w(z)=EXP(-z**2)*erfc(-i*z) in the upper 
  !c-mm    half-plane z=x+i*y (i.e., for y>=0) Maximum relative error of both parts is < 1e-4.
  !c-mm From Humlicek (1982).

  REAL*8, PARAMETER ::          &
     pi = 3.141592653589793D0

  INTEGER, INTENT(IN   ) ::     &
     numlines,                  & ! Number of lines in total range
     gas                          ! HITRAN number of gas of interest

  REAL*8, INTENT(IN   ) ::      &
     line_dat(:,:),             & ! Array of basic line data processed from HITRAN file
     range_min,                 & ! The min wavenumber of spectral range of g(k) distribution
     range_max,                 & ! The max wavenumber of spectral range of g(k) distribution
     vstep,                     & ! Amount (in m-1) by which the spectrum is stepped over
     temp,                      & ! Temperature used to calculate chi and scale density
     maxwid

  INTEGER, INTENT(  OUT) ::     &
     numvals                      ! Number of values of the absorption coefficient calculated

  REAL*8, INTENT(  OUT) ::      &
     kvals(:)                     ! Array to store absorption coefficient values across the spectral range

  LOGICAL, INTENT(  OUT) ::     &
     too_long                     ! Flag that indicates code is taking too long to run

  REAL*8, INTENT(  OUT) ::      &
     failed_wvn                   ! Calculates wavelength at which calculations exceed time allotted

! Local variables

  REAL*8, DIMENSION(numlines) :: &
     hwd,                        &
     hwl    

  REAL*8                         &
     x,                          &
     y,         	         &
     s,                          &
     vee,                        &
     delta_vees,                 &
     kay, 	                 &
     chi,                        &
     chi_h2o,                    &
     norm_fac,                   &
     f_vvw,                      &
     voigt_strength

  REAL*8                         &
     finish,                     &
     start,                      &
     total_time

  COMPLEX*8                      &
     w4,                         &
     u,                          &
     t

  INTEGER                        &
     a, b, i                    

  !-----------------------------------------------------------------------------------------
  numvals = 0
  vee = range_min
  hwd = line_dat(1,:)*line_dat(4,:)
  hwl = line_dat(3,:)
  a = 1
  b = numlines                                                     ! maxwid corresponds to 700 line halfwidths.  I choose a line
  !G maxwid is a buffer around spectral band width (edge to edge)
  !G   within which to get the wings from lines whose centers are
  !G   not within the band but are close enough that their wings
  !G   contribute.
  !c-mm  maxwid = MIN( MAX(hwl((b-a)/2+a),hwd((b-a)/2+a))*7.0D+03, 2.5D+03 ) !   randomly for estimating the halfwidth.  Choose larger of hwl and hwd.
!c-mm  IF (gas == 1) THEN
!c-mm     maxwid = 2.5D+03    ! H2O continuum from MT_CKD cuts off lines at 25 cm-1.  Set this to avoid double counting from 25-500 cm-1.
!c-mm  ELSE
!c-mm     maxwid = 5.0D+04    ! Set to 500 cm-1 from line center (from Halevy et al., 2009)
!c-mm  ENDIF

  a = 1                                               
  b = 1 

  total_time=0.0
  too_long=.FALSE.
  calc: DO                                                         ! Now the calculation of absorption coefficients begins.
     CALL cpu_time(start)
     IF (vee >= range_max ) EXIT calc
     IF (a > 1) a = a-1                                            ! a and b now define the part of the line_dat array used for the 
     IF (b > 1) b = b-1                                            !    wavenumber at which the absorption coefficient is calculated.
     IF (b == 0) b = 1
     find_a: DO
        IF (line_dat(1,a) >= (vee-maxwid) .OR. a >= numlines) EXIT find_a
        a = a+1
     ENDDO find_a
     find_b: DO
        IF (line_dat(1,b) >= (vee+maxwid) .OR. b >= numlines) EXIT find_b
        b = b+1
     ENDDO find_b
     IF (line_dat(1,b) > (vee+maxwid)) b = b-1
     numvals = numvals + 1
     kay = 1.0D-25
     IF (a <= b) THEN
        DO i = a, b
           delta_vees = vee-line_dat(1,i)
           x = delta_vees/hwd(i)
           y = hwl(i)/hwd(i)
           t = CMPLX(y,-x)
           s = DABS(x) + y
           u = t*t
           IF ( (s < 5.5D0).AND.(y < 0.195D0*DABS(x) - 0.176D0) ) THEN
              w4 = CEXP(u) - t*(36183.31D0 - u*(3321.9905D0 - u*(1540.787D0 - u*(219.0313D0 - u*(35.76683D0 - u*(1.320522D0 -     &
                   u*0.56419D0))))))/(32066.6D0 - u*(24322.84D0 - u*(9022.228D0 - u*(2186.181D0 - u*(364.2191D0 - u*(61.57037D0 - &
                   u*(1.841439D0 - u)))))))
           ELSE IF (s < 5.5d0) THEN
              w4 = (16.4955D0 + t*(20.20933D0 + t*(11.96482D0 + t*(3.778987D0 + t*0.5642236D0))))/ &
                   (16.4955D0 + t*(38.82363D0 + t*(39.27121D0 + t*(21.69274D0 + t*(6.699398D0 + t)))))
           ELSE IF (s < 15.0D0) THEN
              w4 = t*(1.410474D0 + u*0.5641896D0)/(0.75D0 + u*(3.0D0 + u))
           ELSE
              w4 = t*0.5641896D0/(0.5D0 + u)
           ENDIF
           voigt_strength=line_dat(2,i)*DBLE(w4)/(hwd(i)*DSQRT(pi))
           SELECT CASE (gas)
              CASE (1)    ! Has a separate super-Lorentzian chi factor, and uses the Van Vleck Weisskopf lineshape for far wings
                 IF (delta_vees <= 40.*hwd(i)) THEN
                    kay = kay+voigt_strength*norm_fac(vee,line_dat(1,i),temp)*chi_h2o(delta_vees,temp)
                 ELSE
                    kay = kay+line_dat(2,i)*f_vvw(vee,line_dat(1,i),hwl(i))*norm_fac(vee,line_dat(1,i),temp)*chi_h2o(delta_vees,temp)
                 ENDIF
              CASE (2)    ! Has a chi factor, and sticks with the Voigt lineshape
                 kay = kay+voigt_strength*chi(delta_vees,temp)*norm_fac(vee,line_dat(1,i),temp)
              CASE DEFAULT   ! Has no chi factor, but uses Van Vleck Weisskopf lineshape for far wings.
                 IF (delta_vees <= 40.*hwd(i)) THEN
                    kay = kay+voigt_strength*norm_fac(vee,line_dat(1,i),temp)
                 ELSE
                    kay = kay+line_dat(2,i)*f_vvw(vee,line_dat(1,i),hwl(i))*norm_fac(vee,line_dat(1,i),temp)
                 ENDIF
           END SELECT
        ENDDO							                                                      
     ENDIF  ! (a<=b) condition


     kvals(numvals) = kay
     vee = vee + vstep
     CALL cpu_time(finish)
     total_time=total_time+(finish-start)
     IF (total_time > 72000.) THEN    !c-mm  I've set this to *20* hours so that it will be guaranteed to finish for whatever we choose
        too_long=.TRUE.               !c-mm     to be the margin.  This bypasses the 'auto' feature, and should be reverted back to
        failed_wvn=vee-vstep          !c-mm     something more reasonable (3600? 7200?) to use the auto feature again.
        EXIT
     ENDIF
  ENDDO calc
END SUBROUTINE abs_coeff_calc


!---------------------------------------------------------------------------------------------------------------------------------
! 5.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to create precise g(k) distributions
!c-mm
!c-mm This subroutine shall create all the data describing the g(k) distribution.  There are two pieces of data for each g(k) bin: 
!c-mm the k value and the width of the bin (in units of normalised probability).  For each pressure and temperature, one g(k) 
!c-mm distribution has a given number of bins.  The array gkdat will contain this information and pass it back into the main 
!c-mm program for output to file.
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE gk_calc(kvals,numvals,gkdat,kk,num_terms)

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE quicksort(numvals,kvals)
       IMPLICIT NONE
       INTEGER                  &
          numvals
       REAL*8                   &
          kvals(:)
     END SUBROUTINE quicksort
  END INTERFACE

  INTEGER, INTENT(IN   ) ::     &
     numvals,                   & ! Number of values of the absorption coefficient calculated
     kk,                        &
     num_terms

  REAL*8, INTENT(IN   ) ::      &
     kvals(:)                     ! The array of absorption coefficient values 

  REAL*8, INTENT(  OUT) ::      &
     gkdat(num_terms)             ! Array for all g(k) distribution data

  ! Local variables:

  REAL*8 ::                     &
     cum_prob,                  & ! The cumulative probability function's value
     cum_prob_prev,             & ! The cumulative probability function's value at the previous data point
     abscissae(num_terms)         ! The Gaussian abscissae values for doing Gaussian quadrature.  Replaces bin_val

  INTEGER ::                    &
       i                          ! Loop variables

  !-----------------------------------------------------------------------------------------  
  CALL get_gaussian_abscissae(num_terms/2, 0.0D0, 0.95D0, abscissae(1:num_terms/2))
  CALL get_gaussian_abscissae(num_terms/2, 0.95D0, 1.0D0, abscissae(num_terms/2+1:num_terms))

  CALL quicksort(numvals, kvals)     ! kvals is sorted by size, minimum value first, to facilitate construction of g(k).

  DO i = 1, num_terms
     gkdat(i) = kvals(IDNINT(DBLE(numvals)*abscissae(i)))
  ENDDO

END SUBROUTINE gk_calc


!---------------------------------------------------------------------------------------------------------------------------------
! 6.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to calculate the gaussian abscissae for a given number of quadrature points
!c-mm
!c-mm Method:
!c-mm	Originally derived from gauleg.f in Numerical Recipes for F
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE get_gaussian_abscissae(num_terms,x1,x2,abscissae)

  IMPLICIT NONE

  INTEGER, INTENT(IN   ) :: &
     num_terms

  REAL*8, INTENT(IN   )  :: &
     x1,                    &
     x2

  REAL*8, INTENT(  OUT)  :: &     
     abscissae(num_terms)

  ! Local variables

  REAL*8, PARAMETER ::      &
     small = 5.0D-10,       &
     pi = 3.141592653589D0

  INTEGER ::                &
     m, j, k

  REAL*8 ::                 &
     xm,                    &
     xl, 		    &
     z, 		    &
     z1, 		    &
     p1, 		    &
     p2, 		    &
     p3, 		    &
     rn,                    &
     pp

  !-----------------------------------------------------------------------------------------    
  ! Initialize and set initial values
  m = INT((num_terms+1)/2)
  xm = (x2 + x1)/2.0D0
  xl = (x2 - x1)/2.0D0
  rn = DBLE(num_terms)

  DO j = 1, m
     z = DCOS(pi*(DBLE(j)-0.25D0)/(rn + 0.5D0)) 
     z1 = 1.0D6
     DO WHILE (DABS(z-z1) >= small)
        p1 = 1.0D0
        p2 = 0.0D0
        DO k = 1,num_terms
           p3 = p2
           p2 = p1
           p1 = (DBLE(2*k-1)*z*p2 - DBLE(k-1)*p3)/DBLE(k)
        ENDDO
        pp = rn*(z*p1-p2)/(z*z-1.0D0)
        z1 = z
        z = z1-p1/pp
     ENDDO
     abscissae(j) = xm - xl*z
     abscissae(num_terms+1-j) = xm + xl*z
  ENDDO

END SUBROUTINE get_gaussian_abscissae


!---------------------------------------------------------------------------------------------------------------------------------
! 7.
!-----------------------------------------------------------------------------------------------------------------------------------
!+ subroutine to perform a quick sort.
!
! method:
!	the standard quicksort sorting algorithm is implemented.
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE quicksort(numvals,kvals)

  IMPLICIT NONE

  REAL*8, INTENT(INOUT) ::       &
     kvals(:)

  INTEGER, INTENT(IN   ) ::      &
     numvals

  INTEGER, PARAMETER ::          &
     nn = 15,                    &
     nstack = 50

  ! Local variables:

  REAL*8 ::                      &
     a
      
  INTEGER ::                     &
     k,                          &
     i,                          &
     j,                          &
     jstack,                     &
     l,                          &
     r

  INTEGER, DIMENSION(nstack) ::  &
     istack

  !-----------------------------------------------------------------------------------------    
  jstack = 0
  l = 1
  r = numvals
  loop: DO
     IF (r-l < nn) THEN
        DO j = l+1,r
           a = kvals(j)
           do_a: DO i = j-1, l, -1
              IF (kvals(i) <= a) EXIT do_a
              kvals(i+1) = kvals(i)
           ENDDO do_a
           kvals(i+1) = a
        ENDDO
        if (jstack == 0) RETURN
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
     ELSE
        k = (l+r)/2
        CALL swap(kvals(k),kvals(l+1))
        CALL swap_mask(kvals(l),kvals(r),kvals(l)>kvals(r))
        CALL swap_mask(kvals(l+1),kvals(r),kvals(l+1)>kvals(r))
        CALL swap_mask(kvals(l),kvals(l+1),kvals(l)>kvals(l+1))
        i = l+1
        j = r
        a = kvals(l+1)
        loop_1: DO
           loop_2: DO
              i = i+1
              IF (kvals(i) >= a) EXIT loop_2
           ENDDO loop_2
           loop_3: DO
              j = j-1
              IF (kvals(j) <= a) EXIT loop_3
           ENDDO loop_3
           IF (j < i) EXIT loop_1
           CALL swap(kvals(i),kvals(j))
        ENDDO loop_1
        kvals(l+1) = kvals(j)
        kvals(j) = a
        jstack = jstack+2
        IF (jstack > nstack) write(0,*) 'Error, NSTACK is too small'
        IF (r-i+1 >= j-l) THEN
           istack(jstack) = r
           istack(jstack-1) = i
           r = j-1
        ELSE
           istack(jstack)=j-1
           istack(jstack-1)=l
           l = i
        ENDIF
     ENDIF
  ENDDO loop

END SUBROUTINE quicksort


!---------------------------------------------------------------------------------------------------------------------------------
! 8.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to perform a masked swap of two blocks of data
!c-mm
!c-mm Method:
!c-mm	From Numerical Recipes for F90
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE swap(a,b)

  IMPLICIT NONE

  REAL*8, INTENT(INOUT) :: &
     a, b

  ! Local variables

  REAL*8                   &
     swp

  !-----------------------------------------------------------------------------------------    
  swp = a
  a = b
  b = swp

END SUBROUTINE swap


!---------------------------------------------------------------------------------------------------------------------------------
! 9.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to perform a masked swap of two blocks of data
!c-mm
!c-mm Method:
!c-mm	From Numerical Recipes for F90
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE swap_mask(a,b,mask)

  IMPLICIT NONE

  REAL*8, INTENT(INOUT)  :: &
     a, b

  LOGICAL, INTENT(IN   ) :: &
     mask

  ! Local variables

  REAL*8                    &
     swp

  !-----------------------------------------------------------------------------------------    
  IF (mask) THEN
     swp = a
     a = b
     b = swp
  ENDIF

END SUBROUTINE swap_mask

!---------------------------------------------------------------------------------------------------------------------------------
! 10.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Function to calculate a sub-Lorentzian adjustment to line strengths in the far wings
!c-mm
!c-mm Method:
!c-mm	Analytical expression taken from Perrin and Hartmann (1989)
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION chi(delta_vees,temp)

  IMPLICIT NONE

  REAL*8, INTENT(IN   ) :: delta_vees, & !c-mm  Distance [m-1] from line center
                           temp          !c-mm  Temperature


! Local variables

  REAL*8 ::                            &
     chi

  REAL*8, PARAMETER ::                 &
     alpha1 = 0.0888D0,                &
     alpha2 = 0.0D0,                   &
     alpha3 = 0.0232D0,                &
     beta1 = -0.16D0,                  &
     beta2 = 0.0526D0,                 &
     beta3 = 0.0D0,                    &
     eps1 = 0.0041D0,                  &
     eps2 = 0.00152D0,                 &
     eps3 = 0.0D0

  REAL*8 ::                            &
     b1,                               &
     b2,                               &
     b3,                               &
     delta_vees_cm                       !c-mm  Local variable to convert delta_vees from m-1 to cm-1

  !-----------------------------------------------------------------------------------------    


   delta_vees_cm = delta_vees/100.
   b1=alpha1+beta1*DEXP(-eps1*temp)
   b2=alpha2+beta2*DEXP(-eps2*temp)
   b3=alpha3+beta3*DEXP(-eps3*temp)
   IF (ABS(delta_vees_cm) < 3.0) THEN
      chi = 1.0
   ELSE IF ((ABS(delta_vees_cm) >= 3.0).AND.(ABS(delta_vees_cm) < 30.0)) THEN
      chi=DEXP(-b1*(ABS(delta_vees_cm)-3.0))
   ELSE IF ((ABS(delta_vees_cm) >= 30.0).AND.(ABS(delta_vees_cm) < 120.0)) THEN
      chi=DEXP(-b1*(30.0-3.0)-b2*(ABS(delta_vees_cm)-30.0))
   ELSE
      chi=DEXP(-b1*(30.0-3.0)-b2*(120.0-30.0)-b3*(ABS(delta_vees_cm)-120.0))
   ENDIF

END FUNCTION chi


!---------------------------------------------------------------------------------------------------------------------------------
! 11.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Routine to calculate collision-induced absorption (CIA) at the low density limit at temperatures between 200-800 K
!c-mm
!c-mm Method:
!c-mm	The paper (Icarus) describing this computer model, is based on original paper by
!c-mm     M. Gruszka and A. Borysow,"Computer Simulation of the Far Infrared Collision-induced absorption spectra 
!c-mm     of gaseous CO2", Mol. Phys., vol. 93, pp. 1007-1016, 1998.
!c-mm  The actual computer model is derived from the paper:
!c-mm     M. Gruszka and A. Borysow, "Roto-translational collision-induced absorption of CO2 for the atmosphere of Venus
!c-mm     at frequencies from 0 to 250 cm-1, at temperatures from 200 to 800 K", Icarus, 129, 172-177, 1997. 
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE gruszka(temp,wn,abcoef) ! computes the spectrum

   IMPLICIT NONE

   REAL*8, INTENT(IN   ) ::         &
      temp,                         & ! temperature, restricted to the range from 200 K to 800 K
      wn                              ! wavenumber, designed to operate up to 250 cm-1.  No restriction on higher, but degraded accuracy

   REAL*8, INTENT(  OUT) ::         &
      abcoef                          ! absorption coefficient in cm^-1 amagat^-2

! Local variables
   INTEGER, PARAMETER ::            &
      ndeg = 3

   REAL*8, DIMENSION(ndeg) ::       &
      ah,                           &
      bh,                           &
      at,                           &
      bt,                           &
      gam

!c-mm     tau^{L}_{1}:
   DATA ah /0.1586914382D+13,-0.9344296879D+01,0.6943966881D+00/
!c-mm     tau^{L}_{2}:
   DATA bh /0.1285676961D-12,0.9420973263D+01,-0.7855988401D+00/
!c-mm     tau^{H}_{1}:
   DATA at /0.3312598766D-09,0.7285659464D+01,-0.6732642658D+00/
!c-mm     tau^{H}_{2}:
   DATA bt /0.1960966173D+09,-0.6834613750D+01,0.5516825232D+00/
!c-mm     gamma_{1}:
   DATA gam /0.1059151675D+17,-0.1048630307D+02, 0.7321430968D+00/

   REAL*8, PARAMETER ::             &
      w1 = 50.0D0,                  &
      w2 = 100.0D0

   REAL*8 ::                        &
      a1,                           &
      b1,                           &
      a2,                           &
      b2,                           &
      gamma,                        &
      bc1,                          &
      bc2,                          &
      bcbc,                         &
      scon,                         &
      mtot,                         &
      spunit,                       &
      gm0con,                       &
      icm,                          &
      lntemp,                       &
      xk1,                          &
      x1,                           &
      x2

   lntemp=DLOG(temp)
   a1=1.0D0/(ah(1)*DEXP(ah(2)*lntemp+ah(3)*lntemp*lntemp))
   b1=bh(1)*DEXP(bh(2)*lntemp+bh(3)*lntemp*lntemp)
   a2=1.0D0/(at(1)*DEXP(at(2)*lntemp+at(3)*lntemp*lntemp))
   b2=bt(1)*DEXP(bt(2)*lntemp+bt(3)*lntemp*lntemp)
   gamma=gam(1)*DEXP(gam(2)*lntemp+gam(3)*lntemp*lntemp)
   icm=0.1885D0                                             !c-mm   Conversion from wavenumber to frequency  = 2*pi*c  (c in cm/ps)
   x1=b1*DSQRT(a1*a1+icm*icm*wn*wn)
   bc1=DEXP(a1*b1)*a1*xk1(x1)/(a1*a1+icm*icm*wn*wn)
   x2=b2*DSQRT(a2*a2+icm*icm*wn*wn)
   bc2=DEXP(a2*b2)*a2*xk1(x2)/(a2*a2+icm*icm*wn*wn)
!c-mm  This IF-THEN block is used to join the low-frequency and high-frequency fits
   IF (wn < 50.D0) THEN
      bcbc=bc1
   ELSE IF (wn > 100.D0) THEN
      bcbc=bc2
   ELSE
      bcbc=(1.-((wn-50.D0)/50.D0))*bc1+((wn-50.D0)/50.D0)*bc2
   ENDIF
   spunit=1.296917D55               !c-mm   Equals 2*pi*pi/3*n*n/k/c, found in Eq. 1 of Gruszka and Borysow (1997)  n=2.687e19 cm^-3
   gm0con=1.259009D-6               !c-mm   Equals 3ck/pi/pi, found after Eq. 3 in Gruszka and Borysow (1997)
   scon=spunit/temp                 !c-mm   Equals the entirety of the first two terms of Eq. 1 in Gruszka and Borysow (1997)
   mtot=gm0con*temp*gamma*1.0D-56   !c-mm   Equals M_0, found after Eq. 3 in Gruszka and Borysow (1997)
   abcoef=scon*mtot*bcbc*wn*wn      !c-mm   Solution to Eq. 1 of Gruszka and Borysow (1997) in cm^-1/amagat^-2

   RETURN

END SUBROUTINE gruszka

!---------------------------------------------------------------------------------------------------------------------------------
! 12.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Function to calculate modified Bessel function used in Gruszka routine.  Taken from Gruszka & Borysow code.
!c-mm
!c-mm Method:
!c-mm	Expression taken from Abramowitz and Stegun, "Handbook of Mathematical Functions"
!-----------------------------------------------------------------------------------------------------------------------------------

FUNCTION xk1(x)
!c-mm Modified Bessel function K1(x) times x
!c-mm Precision is better than 2.2e-7 everywhere.
!c-mm From Abramowitz and Stegun, p. 379; Tables p. 417

   IMPLICIT NONE

   REAL*8, INTENT(IN   ) ::    &
      x

   REAL*8 ::                   &
      xk1

! Local variables

   REAL*8 ::                   &
      t,                       &
      fi1,                     &
      p,                       &
      y

   IF ((x-2.) <= 0.0) THEN
      t=(x/3.75)**2
      fi1=x*((((((0.00032411*t+0.00301532)*t+0.02658733)*t+0.15084934)*t+0.51498869)*t+0.87890594)*t+.5)
      t=(x/2.)**2
      p=(((((-0.00004686*t-0.00110404)*t-0.01919402)*t-0.18156897)*t-0.67278579)*t+0.15443144)*t+1.
      xk1=x*DLOG(x/2)*fi1+p
   ELSE
      t=2./x
      p=(((((-0.00068245*t+0.00325614)*t-0.00780353)*t+0.01504268)*t-0.03655620)*t+0.23498619)*t+1.25331414
      y=DMIN1(x,330.D0)
      xk1=DSQRT(y)*DEXP(-y)*p
   ENDIF

   RETURN

END FUNCTION xk1

!---------------------------------------------------------------------------------------------------------------------------------
! 13.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Routine for reading Baranov data for CO2 dimer absorption from 1200-1500 cm-1.  This routine only does the read of
!c-mm    the data file CO2_dimer_data, and does it once, before the loop over spectrum is performed.
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE read_baranov_data(wn_arr,temp_arr,dim_arr)

   IMPLICIT NONE

   INTEGER, PARAMETER ::             &
      nS = 3000,                     &   !c-mm  Number of wavenumbers of data available
      nT = 9                             !c-mm  Number of temperatures of data available

   REAL*8, DIMENSION(nS) ::          &
      wn_arr

   REAL*8, DIMENSION(nT) ::          &
      temp_arr

   REAL*8, DIMENSION(nS,nT) ::       &
      dim_arr

   INTEGER ::                        &
      ios

   OPEN(33,FILE='./ck/CO2_dimer_data',FORM='unformatted',STATUS='old',IOSTAT=ios)

   IF (ios /= 0) THEN                  ! file not found
      WRITE(0,*) 'Error!'
      WRITE(0,*) 'CO2 CIA data file could not be found.'
      STOP
   ELSE
      READ(33) wn_arr
      READ(33) temp_arr
      READ(33) dim_arr
   ENDIF
   CLOSE(33)

   RETURN

END SUBROUTINE read_baranov_data

!---------------------------------------------------------------------------------------------------------------------------------
! 14.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Routine for reading water vapor continuum data generated from MT_CKD (AER continuum code).
!c-mm    Reads in both self and foreign continuum and stores them separately.
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE read_h2o_continuum_data(wn_h2o_arr,temp_h2o_arr,dim_h2o_self_arr,dim_h2o_frn_arr)

   IMPLICIT NONE

   INTEGER, PARAMETER ::             &
      nS = 2001,                     &  !c-mm  Number of wavenumbers of data available
      nT = 72                           !c-mm  Number of temperatures of data available

   REAL*8, PARAMETER ::              &
      t_init = 50.0D0,               &  !c-mm  Initial T for constructing continuum (same as for k-values)
      delta_T = 5.D0                    !c-mm  Temperature spacing between code blocks

   REAL*8, DIMENSION(nS) ::          &
      wn_h2o_arr

   REAL*8, DIMENSION(nT) ::          &
      temp_h2o_arr

   REAL*8, DIMENSION(nS,nT) ::       &
      dim_h2o_self_arr,              &
      dim_h2o_frn_arr

   INTEGER ::                        &
      ios,                           &
      i,                             &
      j

   CHARACTER (len=132) ::            &
      junk
   

   OPEN(33,file='./ck/H2O_SELF_COEFF',status='old',iostat=ios)
   OPEN(34,file='./ck/H2O_FRN_COEFF',status='old',iostat=ios)

   IF (ios /= 0) THEN                   ! file not found
      WRITE(0,*) 'Error!'
      WRITE(0,*) 'H2O self continuum data file could not be found.'
      STOP
   ELSE
!c-mm  A bit of a sloppy way to read in the data file.  Data file has form
!c-mm  Temp = XXX
!c-mm     wn(1)      cont_coeff(1)
!c-mm     wn(2)      cont_coeff(2)
!c-mm      .               .
!c-mm      .               .
!c-mm  Temp = YYY
!c-mm     wn(1)      cont_coeff(1)
!c-mm     wn(2)      cont_coeff(2)
!c-mm      .               .
!c-mm      .               .
!c-mm  And so on.  So I'm just reading in (as junk) the temperature value
!c-mm     then hardwiring a value for temp_h2o_arr.  Then values for wn_h2o_arr
!c-mm     and dim_h2o_{self/frn}_arr are read in, but the values for wn_h2o_arr
!c-mm     are just overwritten time and again.  Only the continuum coefficients
!c-mm     are really unique.
      DO i=1,nT
         temp_h2o_arr(i) = t_init+REAL(i-1)*delta_T
         READ(33,'(A)') junk
         READ(34,'(A)') junk
         DO j=1,nS
            READ(33,*) wn_h2o_arr(j),dim_h2o_self_arr(j,i)
            READ(34,*) wn_h2o_arr(j),dim_h2o_frn_arr(j,i)
         ENDDO
      ENDDO
   ENDIF
   CLOSE(33)
   CLOSE(34)

   RETURN

END SUBROUTINE read_h2o_continuum_data


!---------------------------------------------------------------------------------------------------------------------------------
! 15.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Routine for reading CO2 continuum data generated from MT_CKD (AER continuum code).
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE read_co2_continuum_data(wn_co2_arr,temp_co2_arr,dim_co2_arr)

   IMPLICIT NONE

   INTEGER, PARAMETER ::             &
      nS = 5001,                     &  !c-mm  Number of wavenumbers of data available
      nT = 71                           !c-mm  Number of temperatures of data available

   REAL*8, PARAMETER ::              &
      t_init = 50.0D0, &                !c-mm  Initial T for constructing continuum (same as for k-values)
      delta_T = 5.D0                    !c-mm  Temperature spacing between code blocks

   REAL*8, DIMENSION(nS) ::          &
      wn_co2_arr

   REAL*8, DIMENSION(nT) ::          &
      temp_co2_arr

   REAL*8, DIMENSION(nS,nT) ::       &
      dim_co2_arr

   INTEGER ::                        &
      ios,                           &
      i,                             &
      j

   CHARACTER (len=132) ::            &
      junk

   

   OPEN(33,FILE='./ck/CO2_CONT_COEFF',FORM='formatted',STATUS='old',IOSTAT=ios)

   IF (ios /= 0) THEN                   ! file not found
      WRITE(0,*) 'Error!'
      WRITE(0,*) 'CO2 continuum data file could not be found.'
      STOP
   ELSE
!c-mm  A bit of a sloppy way to read in the data file.  Data file has form
!c-mm  Temp = XXX
!c-mm     wn(1)      cont_coeff(1)
!c-mm     wn(2)      cont_coeff(2)
!c-mm      .               .
!c-mm      .               .
!c-mm  Temp = YYY
!c-mm     wn(1)      cont_coeff(1)
!c-mm     wn(2)      cont_coeff(2)
!c-mm      .               .
!c-mm      .               .
!c-mm  And so on.  So I'm just reading in (as junk) the temperature value
!c-mm     then hardwiring a value for temp_co2_arr.  Then values for wn_co2_arr
!c-mm     and dim_co2_arr are read in, but the values for wn_co2_arr
!c-mm     are just overwritten time and again.  Only the continuum coefficients
!c-mm     are really unique.
      DO i=1,nT
         temp_co2_arr(i) = t_init+REAL(i-1)*delta_T
         READ(33,'(A)') junk
         DO j=1,nS
            READ(33,*) wn_co2_arr(j),dim_co2_arr(j,i)
         ENDDO
      ENDDO
   ENDIF
   CLOSE(33)

   RETURN

END SUBROUTINE read_co2_continuum_data

!---------------------------------------------------------------------------------------------------------------------------------
! 16.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Routine to interpolate the Baranov data table to the exact wavelength and temperature we want.
!c-mm
!c-mm Method:
!c-mm	Basic 2-D interpolation scheme.  We already are assured that we are within the wavenumber and
!c-mm      temperature range by the fact that we're here, so no need to check again.
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE bilinear(x_arr,y_arr,nX,nY,f2d_arr,x,y,f)

   IMPLICIT NONE

   INTEGER, INTENT(IN   ) ::                  &
      nX,                                     &
      nY

   REAL*8, INTENT(IN   ) ::                   &
      x,                                      &
      y

   REAL*8, INTENT(IN   ), DIMENSION(nX) ::    &
      x_arr

   REAL*8, INTENT(IN   ), DIMENSION(nY) ::    &
      y_arr

   REAL*8, INTENT(IN   ), DIMENSION(nX,nY) :: &
      f2d_arr

   REAL*8, INTENT(  OUT) ::                   &
      f

! Local variables

   INTEGER ::                                 &
      i,                                      &
      j,                                      &
      a,                                      &
      b

   REAL*8 ::                                  &
      x1,                                     &
      x2,                                     &
      y1,                                     &
      y2,                                     &
      f11,                                    &
      f12,                                    &
      f21,                                    &
      f22,                                    &
      fA,                                     &
      fB
      
   i=1
   DO
      IF (x_arr(i) > x) THEN
         x1=x_arr(i-1)
         x2=x_arr(i)
         a=i-1
         EXIT
      ELSE
         i=i+1
      ENDIF
   ENDDO

!     in the y (temperature) direction 2nd
   j=1
   DO
      IF (y_arr(j) > y) THEN
         y1=y_arr(j-1)
         y2=y_arr(j)
         b=j-1
         EXIT
      ELSE
         j=j+1
      ENDIF
   ENDDO
      
   f11=f2d_arr(a,b)
   f21=f2d_arr(a+1,b)
   f12=f2d_arr(a,b+1)
   f22=f2d_arr(a+1,b+1)
      
!     1st in x-direction
   fA=f11*(x2-x)/(x2-x1)+f21*(x-x1)/(x2-x1)
   fB=f12*(x2-x)/(x2-x1)+f22*(x-x1)/(x2-x1)
      
!     then in y-direction
   f=fA*(y2-y)/(y2-y1)+fB*(y-y1)/(y2-y1)
      
   RETURN

END SUBROUTINE bilinear

!---------------------------------------------------------------------------------------------------------------------------------
! 17.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Function to calculate a normalization factor to account for asymmetry between absorption above and below line center.
!c-mm
!c-mm Method:
!c-mm	Analytical expression taken from Halevy et al. (2009)
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION norm_fac(vee,line_center,temp)

  IMPLICIT NONE

  REAL*8, INTENT(IN   ) ::             &
     vee,                              & !c-mm  Distance [m-1] from line center
     line_center,                      & !c-mm  Line center [m-1]
     temp                                !c-mm  Temperature

! Local variables

  REAL*8 ::                            &
     norm_fac

  REAL*8, PARAMETER ::                 &
     h = 6.626D-34,                    &
     k = 1.38D-23,                     &
     c = 3.0D-8

  norm_fac = (vee*c*TANH(h*c*vee/2./k/temp))/(line_center*c*TANH(h*c*line_center/2./k/temp))

  !-----------------------------------------------------------------------------------------    

END FUNCTION norm_fac


!---------------------------------------------------------------------------------------------------------------------------------
! 18.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Function to calculate a Van-Vleck Weisskopf line shape function for far wings
!c-mm
!c-mm Method:
!c-mm	Analytical expression taken from Halevy et al. (2009)
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION f_vvw(vee,line_center,hwl)

  IMPLICIT NONE

  REAL*8, INTENT(IN   ) ::             &
     vee,                              & !c-mm  Distance [m-1] from line center
     line_center,                      & !c-mm  Line center [m-1]
     hwl                                 !c-mm  Lorentz half width

! Local variables

  REAL*8 ::                            &
     f_vvw

  REAL*8, PARAMETER ::                 &
     pi = 3.14159D0

  f_vvw = (vee/line_center)**2*(((hwl/pi)*(1./((vee-line_center)**2+(hwl)**2)))+((hwl/pi)*(1./((vee+line_center)**2+(hwl)**2))))

  !-----------------------------------------------------------------------------------------    

END FUNCTION f_vvw


!---------------------------------------------------------------------------------------------------------------------------------
! 19.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to calculate all continuum values for CO2 and H2O
!c-mm
!c-mm Method:
!c-mm	Uses data from Baranov et al., Gruszka and Borysow and MT_CKD
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE continuum_calc(range_min, range_max, gas, k_self_cont,k_frn_cont,vstep, numvals, &
                          temp,pratio,wn_arr,temp_arr,dim_self_arr,dim_frn_arr,wn_CIA_arr,temp_CIA_arr,dim_CIA_arr,kval_dim)

  IMPLICIT NONE

  INTEGER, PARAMETER ::         &
     nS = 3000,                 &
     nT = 9,                    &
     nS_h2o = 2001,             &
     nT_h2o = 72,               &
     nS_co2 = 5001,             &
     nT_co2 = 71

  REAL*8, PARAMETER ::          &
     pi = 3.141592653589793D0,  &
     gruszka_llim = 0.0D0,      &  ! The lower limit (m-1) of where CIA processes are calculated from Gruszka model.
     gruszka_ulim = 25000.0D0,  &  ! The upper limit (m-1) of where CIA processes are calculated from Gruszka model.
     baranov_llim = 70000.0D0,  &  ! The lower limit (m-1) of where CIA processes are measured from Baranov data.
     baranov_ulim = 200000.0D0, &  ! The upper limit (m-1) of where CIA processes are measured from Baranov data.
     h2o_cont_llim = 0.0D0,     &  ! The lower limit (m-1) of where h2o continuum processes are calculated from MT_CKD
     h2o_cont_ulim = 2000000.D0,&  ! The upper limit (m-1) of where h2o continuum processes are calculated from MT_CKD
     co2_cont_llim = 0.0D0,     &
     co2_cont_ulim = 1000000.D0,&
     p0 = 1.01325D0,            &  ! Standard pressure (in bars)
     t0 = 273.15D0,             &  ! Standard temperature (in K)
     inv_cm_to_inv_m = 100.0D0, &  ! Conversion from cm-1 to m-1
     Lo = 2.687D25,             &  ! Loschmidt's number (molec/m3 @ STP)
     Na = 6.022D23,             &  ! Avogadro's number (molec/mole)   
     mass_co2 = 0.044D0,        &  ! Stopgap measure to hardwire molecular mass of CO2 (kg/mole)         
     mass_h2o = 0.018D0

  INTEGER, INTENT(IN   ) ::     &
     gas,                       & ! HITRAN number of gas of interest
     kval_dim

  REAL*8, INTENT(IN   ) ::      &
     range_min,                 & ! The min wavenumber of spectral range of g(k) distribution
     range_max,                 & ! The max wavenumber of spectral range of g(k) distribution
     vstep,                     & ! Amount (in m-1) by which the spectrum is stepped over
     temp                         ! Temperature used to calculate chi and scale density

  REAL*8, INTENT(IN   ), OPTIONAL ::                   &
     pratio                       ! Pressure used to scale density.  In bars

  REAL*8, INTENT(IN   ), DIMENSION(nS), OPTIONAL ::    &
     wn_CIA_arr

  REAL*8, INTENT(IN   ), DIMENSION(nT), OPTIONAL ::    &
     temp_CIA_arr

  REAL*8, INTENT(IN   ), DIMENSION(nS,nT), OPTIONAL :: &
     dim_CIA_arr

  REAL*8, INTENT(IN   ), DIMENSION(:,:), OPTIONAL ::   &
     dim_self_arr,                                     &
     dim_frn_arr                

  REAL*8, INTENT(IN   ), DIMENSION(:), OPTIONAL ::     & 
     wn_arr                      

  REAL*8, INTENT(IN   ), DIMENSION(:), OPTIONAL ::     &
     temp_arr                      

  INTEGER, INTENT(  OUT) ::     &
     numvals                      ! Number of values of the absorption coefficient calculated

  REAL*8, INTENT(  OUT) ::      &
       k_self_cont(:),          & ! Array to store self continuum values across the spectral range
       k_frn_cont(:)              ! Array to store foreign continuum values across the spectral range

! Local variables

  REAL*8                         &
     vee,                        &
     kay, 	                 &
     k_gruszka,                  & ! Temporary holding of Gruszka CIA value
     k_baranov,                  & ! Temporary holding of Baranov CIA value
     k_h2o_self,                 & ! Temporary holding of self-continuum value
     k_h2o_frn,                  & ! Temporary holding of foreign continuum value
     k_co2_cont

  !-----------------------------------------------------------------------------------------
  numvals = 0
  vee = range_min
  calc: DO                                                         ! Now the calculation of absorption coefficients begins.
     IF (vee >= range_max ) EXIT calc
     numvals = numvals + 1
     IF (numvals > kval_dim) EXIT calc
     kay = 1.0D-25
     k_gruszka=0.0D0
     k_baranov=0.0D0
     k_h2o_self=0.0D0
     k_h2o_frn=0.0D0

     SELECT CASE (gas)
        CASE (1)
           IF ((vee >= h2o_cont_llim).AND.(vee <= h2o_cont_ulim).AND.(temp >= 50.0).AND.(temp <= 400.0)) THEN
              !c-mm******************************************************************
              !c-mm
              !c-mm  H2O Self Continuum Absorption:  0-10000 cm-1
              !c-mm
              !c-mm******************************************************************
              CALL bilinear(wn_arr,temp_arr,nS_h2o,nT_h2o,dim_self_arr,vee/100.,temp,k_h2o_self)
                 k_h2o_self = k_h2o_self*Na/mass_h2o/1.D4*(pratio/p0)*(t0/temp)                   ! Converts k_h2o_self from cm2/molecule to m2/kg.
                 k_self_cont(numvals)=k_h2o_self                                                  ! Also scales for p and T
              !c-mm******************************************************************
              !c-mm
              !c-mm  H2O Foreign Continuum Absorption:  0-10000 cm-1
              !c-mm
              !c-mm******************************************************************
              CALL bilinear(wn_arr,temp_arr,nS_h2o,nT_h2o,dim_frn_arr,vee/100.,temp,k_h2o_frn)
              k_h2o_frn = k_h2o_frn*Na/mass_h2o/1.D4*(pratio/p0)*(t0/temp)                        ! Converts k_h2o_frn from cm2/molecule to m2/kg
              k_frn_cont(numvals)=k_h2o_frn                                                       ! Also scales for p and T
           ENDIF
        CASE (2)
           IF ((vee >= gruszka_llim).AND.(vee <= gruszka_ulim).AND.(temp >= 200.0).AND.(temp <= 800.0)) THEN
              !c-mm******************************************************************
              !c-mm
              !c-mm  Gruszka and Barysow Collision-Induced Absorption:  0-250 cm-1
              !c-mm
              !c-mm******************************************************************
              CALL gruszka(temp,vee/100.,k_gruszka)                                               ! vee is in m-1, gruszka code wants cm-1, so divide by 100.
              k_gruszka = k_gruszka*inv_cm_to_inv_m*(pratio/p0)*(t0/temp)*Na/Lo/mass_co2          ! Converts k_gruszka from cm-1 amagat-2 to m2/kg
           ENDIF
           IF ((vee >= baranov_llim).AND.(vee <= baranov_ulim).AND.(temp >= 100.0).AND.(temp <= 400.0)) THEN
              !c-mm******************************************************************
              !c-mm
              !c-mm  Baranov Collision-Induced Absorption:  700-2000 cm-1
              !c-mm
              !c-mm******************************************************************
              CALL bilinear(wn_CIA_arr,temp_CIA_arr,nS,nT,dim_CIA_arr,vee/100.,temp,k_baranov)    ! Baranov data is expecting cm-1, so divide vee by 100.
              k_baranov = k_baranov*inv_cm_to_inv_m*(pratio/p0)*(t0/temp)*Na/Lo/mass_co2          ! Converts k_baranov from cm-1 amagat-2 to m2/kg
           ENDIF
!c-mm           IF ((vee >= co2_cont_llim).AND.(vee <= co2_cont_ulim)) THEN
!c-mm              !c-mm******************************************************************
!c-mm              !c-mm
!c-mm              !c-mm  CO2 Continuum Absorption:  0-10000 cm-1
!c-mm              !c-mm
!c-mm              !c-mm******************************************************************
!c-mm              !c-mm  NOTE:  The CO2 continuum data already has an empirical scaling to account for sub-Lorentzian line shape.
!c-mm              CALL bilinear(wn_arr,temp_arr,nS_co2,nT_co2,dim_self_arr,vee/100.,temp,k_co2_cont)  ! vee is in m-1, dim_co2_arr is in cm-1, so divide vee by 100.
!c-mm              k_co2_cont = k_co2_cont*Na/mass_co2/1.D4                                            ! Converts k_co2_continuum from cm2/molecule to m2/kg
!c-mm              k_self_cont(numvals) = k_co2_cont                                                   ! Treat co2 continuum as if it were 'self' continuum
!c-mm              k_frn_cont(numvals) = 0.0D0                                                         ! There is no 'foreign' continuum for CO2 so set to zero.
!c-mm           ENDIF
           k_self_cont(numvals) = k_gruszka + k_baranov
           k_frn_cont(numvals) = 0.0D0
        CASE DEFAULT
              k_self_cont(numvals) = 0.0D0                                                        ! No other gases have a known continuum function
              k_frn_cont(numvals) = 0.0D0
        END SELECT
     vee = vee + vstep
  ENDDO calc

END SUBROUTINE continuum_calc


!---------------------------------------------------------------------------------------------------------------------------------
! 20.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Function to calculate a super-Lorentzian adjustment to line strengths in the far wings for H2O
!c-mm
!c-mm Method:
!c-mm	Analytical expression taken from Clough et al., (1989)
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION chi_h2o(delta_vees,temp)

  IMPLICIT NONE

  REAL*8, INTENT(IN   ) :: delta_vees, & !c-mm  Distance [m-1] from line center
                           temp          !c-mm  Temperature

! Local variables

  REAL*8 ::                            &
     chi_h2o

  REAL*8, PARAMETER ::                 &
     alpha1 = 0.0888D0,                &
     alpha2 = 0.0D0,                   &
     alpha3 = 0.0232D0,                &
     beta1 = -0.16D0,                  &
     beta2 = 0.0526D0,                 &
     beta3 = 0.0D0,                    &
     eps1 = 0.0041D0,                  &
     eps2 = 0.00152D0,                 &
     eps3 = 0.0D0

  REAL*8 ::                            &
     b1,                               &
     b2,                               &
     b3,                               &
     delta_vees_cm                       !c-mm  Local variable to convert delta_vees from m-1 to cm-1

  !-----------------------------------------------------------------------------------------    


   delta_vees_cm = delta_vees/100.
   b1=alpha1+beta1*DEXP(-eps1*temp)
   b2=alpha2+beta2*DEXP(-eps2*temp)
   b3=alpha3+beta3*DEXP(-eps3*temp)
   IF (ABS(delta_vees_cm) < 3.0) THEN
      chi_h2o = 1.0
   ELSE IF ((ABS(delta_vees_cm) >= 3.0).AND.(ABS(delta_vees_cm) < 30.0)) THEN
      chi_h2o=DEXP(-b1*(ABS(delta_vees_cm)-3.0))
   ELSE IF ((ABS(delta_vees_cm) >= 30.0).AND.(ABS(delta_vees_cm) < 120.0)) THEN
      chi_h2o=DEXP(-b1*(30.0-3.0)-b2*(ABS(delta_vees_cm)-30.0))
   ELSE
      chi_h2o=DEXP(-b1*(30.0-3.0)-b2*(120.0-30.0)-b3*(ABS(delta_vees_cm)-120.0))
   ENDIF
   chi_h2o = 1.0      ! Just hardwire this off for now.

END FUNCTION chi_h2o



!---------------------------------------------------------------------------------------------------------------------------------
! 21.
!-----------------------------------------------------------------------------------------------------------------------------------
!c-mm Subroutine to calculate Voigt lineshape
!c-mm
!c-mm Method:
!c-mm	Uses approach from Wells et al., (1999)
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE humlik(line_dat, range_min, range_max, gas, kvals, vstep, numvals, numlines, temp, maxwid, too_long, failed_wvn,kval_dim)

!c-mm To calculate the Faddeeva function with relative error less than 10^(-4).

  IMPLICIT NONE

  REAL*8, PARAMETER ::          &
     pi = 3.141592653589793D0

  INTEGER, INTENT(IN   ) ::     &
     numlines,                  & ! Number of lines in total range
     gas,                       & ! HITRAN number of gas of interest
     kval_dim

  REAL*8, INTENT(IN   ) ::      &
     line_dat(:,:),             & ! Array of basic line data processed from HITRAN file
     range_min,                 & ! The min wavenumber of spectral range of g(k) distribution
     range_max,                 & ! The max wavenumber of spectral range of g(k) distribution
     vstep,                     & ! Amount (in m-1) by which the spectrum is stepped over
     temp,                      & ! Temperature used to calculate chi and scale density
     maxwid

  INTEGER, INTENT(  OUT) ::     &
     numvals                      ! Number of values of the absorption coefficient calculated

  REAL*8, INTENT(  OUT) ::      &
     kvals(:)                     ! Array to store absorption coefficient values across the spectral range

  LOGICAL, INTENT(  OUT) ::     &
     too_long                     ! Flag to indicate that it takes too long to calculate

  REAL*8, INTENT(  OUT) ::      &
     failed_wvn                   ! Carries the wavenumber at which the calculations become 'too long', for diagnostics.

  !Local variables

  REAL*8, PARAMETER ::          &
     rrtpi = 0.56418958D0,      & ! 1/SQRT(pi)
     y0 = 1.5D0,                &
     y0py0 = 2.0*y0,            &
     y0q = y0*y0

  REAL*8, DIMENSION(0:5), SAVE :: &
     c,                           &
     s,                           &
     t

  DATA c / 1.0117281,     -0.75197147,        0.012557727,        &
           0.010022008,   -0.00024206814,     0.00000050084806 /
  DATA s / 1.393237,       0.23115241,       -0.15535147,         &
           0.0062183662,   0.000091908299,   -0.00000062752596 /
  DATA t / 0.31424038,     0.94778839,        1.5976826,          &
           2.2795071,      3.0206370,         3.8897249 /

  INTEGER ::                &
     i,                     &
     j                                                      ! Loop variables
  INTEGER ::                &
     rg1,                   &
     rg2,                   &
     rg3                                                    ! y polynomial flags
  REAL*8 ::                 &
     x,                     &
     y,                     &
     k,                     &
     abx,                   &                               ! ABS(x)
     xq,                    &                               ! x**2
     yq,                    &                               ! y**2
     yrrtpi                                                 ! y/SQRT(pi)
  REAL*8 ::                 &
     xlim0,                 &
     xlim1,                 &
     xlim2,                 &
     xlim3,                 & 
     xlim4                                                  ! |x| on region boundaries
  REAL*8 ::                 &
     a0,                    &
     d0,                    &
     d2,                    &
     e0,                    &
     e2,                    &
     e4,                    &
     h0,                    &
     h2,                    &
     h4,                    &
     h6                                                     ! W4 temporary variables
  REAL*8 ::                 &
     p0,                    &
     p2,                    &
     p4,                    &
     p6,                    &
     p8,                    &
     z0,                    &
     z2,                    &
     z4,                    &
     z6,                    &
     z8    
  REAL*8, DIMENSION(0:5) :: &
     xp,                    &
     xm,                    &
     yp,                    &
     ym,                    &      ! CPF12 temporary values
     mq,                    &
     pq,                    &
     mf,                    &
     pf
  REAL*8 ::                 &
     d,                     &
     yf,                    &
     ypy0,                  &
     ypy0q  

  real*8, dimension(:), allocatable :: hwd, hwl

  REAL*8 ::                      &
     vee,                        &
     delta_vees,                 &
!c-mm     maxwid,                     &
     kay, 	                 &
     chi,                        &
     chi_h2o,                    &
     norm_fac,                   &
     f_vvw,                      &
     voigt_strength

  REAL*8 ::                      &
     start,                      &
     finish,                     &
     total_time

  INTEGER ::                     &
     a, b, perc_threshold

  allocate(hwd(numlines))
  allocate(hwl(numlines))

  numvals = 0
  vee = range_min
  hwd = line_dat(1,:)*line_dat(4,:)
  hwl = line_dat(3,:)
  a = 1
  b = numlines                                                     ! maxwid corresponds to 700 line halfwidths.  I choose a line
  !G maxwid is a buffer around spectral band width (edge to edge)
  !G   within which to get the wings from lines whose centers are
  !G   not within the band but are close enough that their wings
  !G   contribute.
  !c-mm  maxwid = MIN( MAX(hwl((b-a)/2+a),hwd((b-a)/2+a))*7.0D+03, 2.5D+03 ) !   randomly for estimating the halfwidth.  Choose larger of hwl and hwd.
!c-mm  IF (gas == 1) THEN
!c-mm     maxwid = 2.5D+03    ! H2O continuum from MT_CKD cuts off lines at 25 cm-1.  Set this to avoid double counting from 25-500 cm-1.
!c-mm  ELSE
!c-mm     maxwid = 5.0D+04    ! Set to 500 cm-1 from line center (from Halevy et al., 2009)
!c-mm  ENDIF

  a = 1                                               
  b = 1 

  perc_threshold=10
  total_time=0.0
  too_long=.FALSE.
  calc: DO       
      if (100*(vee-range_min)/(range_max-range_min) > perc_threshold) THEN
          write(*,'(A,I3,A)') "humlik: ",perc_threshold,"% complete"
          perc_threshold = perc_threshold + 10
       endif
                                                 ! Now the calculation of absorption coefficients begins.
     CALL cpu_time(start)
     IF (vee >= range_max ) EXIT calc
     IF (a > 1) a = a-1                                            ! a and b now define the part of the line_dat array used for the 
     IF (b > 1) b = b-1                                            !    wavenumber at which the absorption coefficient is calculated.
     IF (b == 0) b = 1
     find_a: DO
        IF (line_dat(1,a) >= (vee-maxwid) .OR. a >= numlines) EXIT find_a
        a = a+1
     ENDDO find_a
     find_b: DO
        IF (line_dat(1,b) >= (vee+maxwid) .OR. b >= numlines) EXIT find_b
        b = b+1
     ENDDO find_b
     IF (line_dat(1,b) > (vee+maxwid)) b = b-1
     numvals = numvals + 1
     IF (numvals > kval_dim) EXIT calc
     kay = 1.0D-25
     IF (a <= b) THEN
        DO i = a, b
           delta_vees = vee-line_dat(1,i)
           x = DSQRT(DLOG(2.0D0))*delta_vees/hwd(i)
           y = DSQRT(DLOG(2.0D0))*hwl(i)/hwd(i)
           yq  = y*y                                                         ! y^2
           yrrtpi = y*rrtpi                                                  ! y/DSQRT(pi)

           IF ( y >= 70.55 ) THEN                                          ! All points
                 xq   = x*x
                 k = yrrtpi / (xq + yq)
           ELSE
              rg1 = 1                                                           ! Set flags
              rg2 = 1
              rg3 = 1

              xlim0 = DSQRT ( 15100.0D0 + y*(40.0 - y*3.6) )                       ! y<70.55
              IF ( y >= 8.425 ) THEN
                 xlim1 = 0.0
              ELSE
                 xlim1 = DSQRT ( 164.0D0 - y*(4.3 + y*1.8) )
              ENDIF
              xlim2 = 6.8 - y
              xlim3 = 2.4*y
              xlim4 = 18.1*y + 1.65
              IF ( y <= 0.000001 ) THEN                                       ! When y<10^-6
                 xlim1 = xlim0                                                    ! avoid W4 algorithm
                 xlim2 = xlim0
              ENDIF

              abx = ABS ( x )                                               ! |x|
              xq  = abx*abx                                                    ! x^2

              IF ( abx >= xlim0 ) THEN                                   ! Region 0 algorithm
                 k = yrrtpi / (xq + yq)
              ELSEIF ( abx >= xlim1 ) THEN                                   ! Humlicek W4 Region 1
                 IF ( rg1 /= 0 ) THEN                                          ! First point in Region 1
                    rg1 = 0
                    a0 = yq + 0.5                                                  ! Region 1 y-dependents
                    d0 = a0*a0
                    d2 = yq + yq - 1.0
                 ENDIF
                 d = rrtpi / (d0 + xq*(d2 + xq))
                 k = d*y *(a0 + xq)
              ELSEIF ( abx > xlim2 ) THEN                                   ! Humlicek W4 Region 2 
                 IF ( rg2 /= 0 ) THEN                                          ! First point in Region 2
                    rg2 = 0
                    h0 =  0.5625 + yq*(4.5 + yq*(10.5 + yq*(6.0 + yq)))            ! Region 2 y-dependents
                    h2 = -4.5 + yq*(9.0 + yq*( 6.0 + yq* 4.0))
                    h4 = 10.5 - yq*(6.0 - yq*  6.0)
                    h6 = -6.0 + yq* 4.0
                    e0 =  1.875 + yq*(8.25 + yq*(5.5 + yq))
                    e2 =  5.25 + yq*(1.0  + yq* 3.0)
                    e4 =  0.75*h6
                 ENDIF
                 d = rrtpi / (h0 + xq*(h2 + xq*(h4 + xq*(h6 + xq))))
                 k = d*y *(e0 + xq*(e2 + xq*(e4 + xq)))
              ELSEIF ( abx < xlim3 ) THEN                                   ! Humlicek W4 Region 3
                 IF ( rg3 /= 0 ) THEN                                          ! First point in Region 3
                    rg3 = 0
                    z0 = 272.1014 + y*(1280.829 + y*(2802.870 + y*(3764.966      &   ! Region 3 y-dependents
                         + y*(3447.629 + y*(2256.981 + y*(1074.409 + y*(369.1989 &
                         + y*(88.26741 + y*(13.39880 + y)))))))))
                    z2 = 211.678 + y*(902.3066 + y*(1758.336 + y*(2037.310       &
                         + y*(1549.675 + y*(793.4273 + y*(266.2987               &
                         + y*(53.59518 + y*5.0)))))))
                    z4 = 78.86585 + y*(308.1852 + y*(497.3014 + y*(479.2576      &
                         + y*(269.2916 + y*(80.39278 + y*10.0)))))
                    z6 = 22.03523 + y*(55.02933 + y*(92.75679 + y*(53.59518      &
                         + y*10.0)))
                    z8 = 1.496460 + y*(13.39880 + y*5.0)
                    p0 = 153.5168 + y*(549.3954 + y*(919.4955 + y*(946.8970      &
                         + y*(662.8097 + y*(328.2151 + y*(115.3772 + y*(27.93941 &
                         + y*(4.264678 + y*0.3183291))))))))
                    p2 = -34.16955 + y*(-1.322256+ y*(124.5975 + y*(189.7730     &
                         + y*(139.4665 + y*(56.81652 + y*(12.79458               &
                         + y*1.2733163))))))
                    p4 = 2.584042 + y*(10.46332 + y*(24.01655 + y*(29.81482      &
                         + y*(12.79568 + y*1.9099744))))
                    p6 = -0.07272979 + y*(0.9377051+ y*(4.266322 + y*1.273316))
                    p8 = 0.0005480304 + y*0.3183291
                 ENDIF
                 d = 1.7724538 / (z0 + xq*(z2 + xq*(z4 + xq*(z6 + xq*(z8+xq)))))
                 k = d*(p0 + xq*(p2 + xq*(p4 + xq*(p6 + xq*p8))))
              ELSE                                                             ! Humlicek CPF12 algorithm
                 ypy0 = y + y0
                 ypy0q = ypy0*ypy0
                 k = 0.0
                 DO j = 0, 5
                    d = x - t(j)
                    mq(j) = d*d
                    mf(j) = 1.0 / (mq(j) + ypy0q)
                    xm(j) = mf(j)*d
                    ym(j) = mf(j)*ypy0
                    d = x + t(j)
                    pq(j) = d*d
                    pf(j) = 1.0 / (pq(j) + ypy0q)
                    xp(j) = pf(j)*d
                    yp(j) = pf(j)*ypy0
                 ENDDO
                 IF ( abx <= xlim4 ) THEN                                      ! Humlicek CPF12 Region I
                    DO j = 0, 5
                       k = k + c(j)*(ym(j)+yp(j)) - s(j)*(xm(j)-xp(j))
                    ENDDO
                 ELSE                                                            ! Humlicek CPF12 Region II
                    yf   = y + y0py0
                    DO j = 0, 5
                       k = k + (c(j)*(mq(j)*mf(j)-y0*ym(j)) + s(j)*yf*xm(j)) /   &
                            (mq(j)+y0q) + (c(j)*(pq(j)*pf(j)-y0*yp(j)) -          &
                            s(j)*yf*xp(j)) / (pq(j)+y0q)
                    ENDDO
                    k = y*k + EXP(-xq)
                 ENDIF
              ENDIF
           ENDIF
           voigt_strength=line_dat(2,i)*DBLE(k)*DSQRT(DLOG(2.0D0))/(hwd(i)*DSQRT(pi))
           SELECT CASE (gas)
              CASE (1)    ! Has a separate super-Lorentzian chi factor, and uses the Van Vleck Weisskopf lineshape for far wings
                 IF (delta_vees <= 40.*hwd(i)) THEN
                    kay = kay+voigt_strength*norm_fac(vee,line_dat(1,i),temp)*chi_h2o(delta_vees,temp)
                 ELSE
                    kay = kay+line_dat(2,i)*f_vvw(vee,line_dat(1,i),hwl(i))*norm_fac(vee,line_dat(1,i),temp)*chi_h2o(delta_vees,temp)
                 ENDIF
              CASE (2)    ! Has a chi factor, and sticks with the Voigt lineshape
                 kay = kay+voigt_strength*chi(delta_vees,temp)*norm_fac(vee,line_dat(1,i),temp)
              CASE DEFAULT   ! Has no chi factor, but uses Van Vleck Weisskopf lineshape for far wings.
                 IF (delta_vees <= 40.*hwd(i)) THEN
                    kay = kay+voigt_strength*norm_fac(vee,line_dat(1,i),temp)
                 ELSE
                    kay = kay+line_dat(2,i)*f_vvw(vee,line_dat(1,i),hwl(i))*norm_fac(vee,line_dat(1,i),temp)
                 ENDIF
           END SELECT
        ENDDO
     ENDIF  ! (a<=b) condition
     kvals(numvals) = kay
     vee = vee + vstep
     CALL cpu_time(finish)
     total_time=total_time+(finish-start)
     IF (total_time > 7200.) THEN
        too_long=.TRUE.
        failed_wvn=vee-vstep
        EXIT
     ENDIF
  ENDDO calc
  deallocate(hwd)
  deallocate(hwl)
END SUBROUTINE humlik
