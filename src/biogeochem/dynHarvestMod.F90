module dynHarvestMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the harvest data, as well as the state updates that happen as a
  ! result of harvest.
  !
  !Polly Buotte, University of Idaho, added code to allow prescribed bark beetle mortality 11/15/2016
  !slevis of SLevis Consulting LLC added code to allow prognostic bark beetle tree mortality 12/2017
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils              , only : endrun
  use dynFileMod              , only : dyn_file_type
  use dynVarTimeUninterpMod   , only : dyn_var_time_uninterp_type
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use CNVegStateType          , only : cnveg_state_type  ! for beetle modules
  use SoilBiogeochemStateType , only : soilbiogeochem_state_type
  use pftconMod               , only : pftcon
  use clm_varcon              , only : grlnd
  use GridcellType            , only : grc  ! for beetle modules
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  use clm_varctl              , only : iulog  ! slevis added for debugging
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynHarvest_init    ! initialize data structures for harvest information
  public :: dynHarvest_interp  ! get harvest data for current time step, if needed
  public :: CNHarvest          ! harvest mortality routine for CN code
  public :: dynProgBB          ! calc. tree mortality from beetles (prognostic)
  !PBuotte: begin prescribe BB
  public :: dynPrescBB_init    ! initialize data structures for prescribed beetle mortality information
  public :: dynPrescBB_interp  ! get beetle mortality data for current time step, if needed
  public :: CNPrescBB          ! beetle mortality routine for CN code
  !PBuotte: end prescribed BB
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNHarvestPftToColumn   ! gather patch-level harvest fluxes to the column level
  !PBuotte: begin prescribed BB
  private :: CNBBmortPftToColumn   ! gather patch-level beetle mortality fluxes to the column level
  !PBuotte: end prescribed BB
  !
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dynHarvest_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing harvest data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dynHarvest_file ! information for the file containing harvest data

  ! Define the underlying harvest variables
  integer, parameter :: num_harvest_inst = 5
  character(len=64), parameter :: harvest_varnames(num_harvest_inst) = &
       [character(len=64) :: 'HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']
  
  type(dyn_var_time_uninterp_type) :: harvest_inst(num_harvest_inst)   ! value of each harvest variable

  real(r8) , allocatable   :: harvest(:) ! harvest rates
  logical                  :: do_harvest ! whether we're in a period when we should do harvest

  !PBuotte: begin prescribed BB
  type(dyn_file_type), target :: dynBeetle_file ! information for the file containing beetle data
  ! Define the underlying beetle mortality variables
  integer, parameter :: num_prescbb_inst = 5
  character(len=64), parameter :: prescbb_varnames(num_prescbb_inst) = &
       [character(len=64) :: 'BEETLE_VH1', 'BEETLE_VH2', 'BEETLE_SH1','BEETLE_SH2', 'BEETLE_SH3']
  type(dyn_var_time_uninterp_type) :: prescbb_inst(num_prescbb_inst)   ! value of each beetle variable
  real(r8) , allocatable   :: beetle_mort_rates(:)  ! grid cell level variable
  logical                  :: do_presc_bb ! whether we're in a period when we should do beetle mortality
  !PBuotte end prescribed BB; slevis renamed presc_bb to beetle_mort_rates

  character(len=*), parameter :: string_not_set = "not_set"  ! string to initialize with to indicate string wasn't set
  character(len=64)        :: harvest_units = string_not_set ! units from harvest variables 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynHarvest_init(bounds, harvest_filename)
    !
    ! !DESCRIPTION:
    ! Initialize data structures for harvest information.
    ! This should be called once, during model initialization.
    ! 
    ! !USES:
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: harvest_filename ! name of file containing harvest information
    !
    ! !LOCAL VARIABLES:
    integer :: varnum     ! counter for harvest variables
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    character(len=64) :: units = string_not_set
    
    character(len=*), parameter :: subname = 'dynHarvest_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    allocate(harvest(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for harvest'//errMsg(sourcefile, __LINE__))
    end if

    ! Get the year from the START of the timestep for consistency with other dyn file
    ! stuff (though it shouldn't actually matter for harvest, given the way the
    ! once-per-year harvest is applied).
    dynHarvest_file = dyn_file_type(harvest_filename, YEAR_POSITION_START_OF_TIMESTEP)
    
    ! Get initial harvest data
    num_points = (bounds%endg - bounds%begg + 1)
    do varnum = 1, num_harvest_inst
       harvest_inst(varnum) = dyn_var_time_uninterp_type( &
            dyn_file=dynHarvest_file, varname=harvest_varnames(varnum), &
            dim1name=grlnd, conversion_factor=1.0_r8, &
            do_check_sums_equal_1=.false., data_shape=[num_points])
       call harvest_inst(varnum)%get_att("units",units)
       if ( trim(units) == string_not_set ) then
          units = "unitless"
       else if ( trim(units) == "unitless" ) then

       else if ( trim(units) /= "gC/m2/yr" ) then
          call endrun(msg=' bad units read in from file='//trim(units)//errMsg(sourcefile, __LINE__))
       end if
       if ( varnum > 1 .and. trim(units) /= trim(harvest_units) )then
          call endrun(msg=' harvest units are inconsitent on file ='// &
               trim(harvest_filename)//errMsg(sourcefile, __LINE__))
       end if
       harvest_units = units
       units = string_not_set
    end do
    
  end subroutine dynHarvest_init

  !PBuotte: begin prescribed BB------------------------------------------------
  subroutine dynPrescBB_init(bounds, beetle_filename)
    !
    ! !DESCRIPTION:
    ! Initialize data structures for beetle information.
    ! This should be called once, during model initialization.
    !
    ! !USES:
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_END_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: beetle_filename  ! name of file containing beetle information
    !
    ! !LOCAL VARIABLES:
    integer :: varnum     ! counter for beetle variables
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code

    character(len=*), parameter :: subname = 'dynPrescBB_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    allocate(beetle_mort_rates(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for beetle mortality'//errMsg(__FILE__,__LINE__))
    end if

    ! Get the year from the END of the timestep for compatibility with how things were
    ! done before, even though that doesn't necessarily make the most sense conceptually.
    dynBeetle_file = dyn_file_type(beetle_filename,YEAR_POSITION_END_OF_TIMESTEP)

    ! Get initial beetle mortality data
    num_points = (bounds%endg - bounds%begg + 1)
    do varnum = 1, num_prescbb_inst
       prescbb_inst(varnum) = dyn_var_time_uninterp_type( &
            dyn_file=dynBeetle_file, varname=prescbb_varnames(varnum), &
            dim1name=grlnd, conversion_factor=1.0_r8, &
            do_check_sums_equal_1=.false., data_shape=[num_points])
    end do

  end subroutine dynPrescBB_init
  !PBuotte: end prescribed BB---------------------------------------------

  !-----------------------------------------------------------------------
  subroutine dynHarvest_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get harvest data for model time, when needed.
    !
    ! Note that harvest data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the harvest rate is held constant through the year.  This is
    ! consistent with the treatment of changing PFT weights, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.
    !
    ! !USES:
    use dynTimeInfoMod , only : time_info_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for harvest variables
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable

    character(len=*), parameter :: subname = 'dynHarvest_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynHarvest_file%time_info%set_current_year()

    ! Get total harvest for this time step
    harvest(bounds%begg:bounds%endg) = 0._r8

    if (dynHarvest_file%time_info%is_before_time_series()) then
       ! Turn off harvest before the start of the harvest time series
       do_harvest = .false.
    else
       ! Note that do_harvest stays true even past the end of the time series. This
       ! means that harvest rates will be maintained at the rate given in the last
       ! year of the file for all years past the end of this specified time series.
       do_harvest = .true.
       allocate(this_data(bounds%begg:bounds%endg))
       do varnum = 1, num_harvest_inst
          call harvest_inst(varnum)%get_current_data(this_data)
          harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + &
               this_data(bounds%begg:bounds%endg)
       end do
       deallocate(this_data)
    end if

  end subroutine dynHarvest_interp

!PBuotte: begin prescribed BB
  subroutine dynPrescBB_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get beetle mortality data for model time, when needed.
    !
    ! Beetle mortality is applied all at once on Aug 1
    ! !USES:
    use dynTimeInfoMod , only : time_info_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for beetle variables
    real(r8), allocatable :: this_data(:) ! data for a single beetle variable

    character(len=*), parameter :: subname = 'dynPrescBB_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynBeetle_file%time_info%set_current_year()

    ! Get total beetle mortality for this time step
    beetle_mort_rates(bounds%begg:bounds%endg) = 0._r8

    if (dynBeetle_file%time_info%is_before_time_series()) then
       ! Turn off beetle mortality before the start of the beetle time series
       do_presc_bb = .false.
    else
       ! Note that do_presc_bb stays true even past the end of the time series.
       ! This means that mortality rates will be maintained at the rate given in the
       ! last year of the file for all years past the end of this specified time series.
       do_presc_bb = .true.
       allocate(this_data(bounds%begg:bounds%endg))
       do varnum = 1, num_prescbb_inst
          call prescbb_inst(varnum)%get_current_data(this_data)
          beetle_mort_rates(bounds%begg:bounds%endg) = beetle_mort_rates(bounds%begg:bounds%endg) +&
               this_data(bounds%begg:bounds%endg)
       end do
       deallocate(this_data)
    end if

  end subroutine dynPrescBB_interp
  !PBuotte: end prescribed BB---------------------------------------------

! -----------------------------------------------------------------------------

  subroutine dynProgBB(bounds, temperature_inst, atm2lnd_inst, cnveg_state_inst, cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Simulate beetle mortality according to Jeff Hicke's MMOBB (Mechanistic
    ! Model of Outbreaking Bark Beetles) version 0.8 received 2017/12/2.
    !
    ! Beetle mortality is applied all at once on Aug 1
    !
    ! The code consists of two main loops:
    ! 1) Where beetles "live" in the focal cell; out-dispersal is determined
    ! 2) Where beetles disperse to in-dispersing cells
    ! Loops and commands determining the correspondence between out- and
    ! in-dispersing beetles are located between the main loops.

    ! !HISTORY:
    ! Sam Levis (slevis) introduced MMOBB to the CLM using as template
    ! the prescribed beetle model introduced to the CLM by Polly Buotte.
    ! Jeff Hicke wrote MMOBB in R and slevis translated to Fortran.

    ! Differences from the R code as implemented in the CLM:
    !
    ! 1) The latest version of MMOBB-MPB in R used an inverse distance
    !    weighting (IDW) scheme, with the weights from 1 at 0 km declining to
    !    0 at 10 km. For FMEC's 4-km2 grid we have agreed that nearest neighbor
    !    is sufficient, and that we will not implement IDW at this time.


    ! !NOTES:
    ! - Restart test 2017/10/27: PASSES
    ! - Test for another pe count 2017/10/26: FAILS  w  call random_number
    !                                         PASSES wo call random_number
    ! - TODO when there's time :-)
    ! 1) Here ave_num_trees_perha = 1113._r8, while stocking = 1000._r8
    !    in CNVegStructUpdateMod.F90 for height calculation.
    !    Recommendation: make consistent.
    ! 2) Consider changing B_od and other beetle variables to integer
    ! 3) Functionalize to avoid repetition of codes
    !    E.g. code from the out-dispersal loop is duplicated in the
    !         in-dispersal loop
    ! ----------------------------------------------------
    ! !USES:
    use spmdMod         , only : MPI_REAL8, MPI_SUM, mpicom
    use clm_time_manager, only : get_curr_date
    use pftconMod       , only : lodgepole
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use decompMod       , only : ldecomp, get_proc_global
    use TemperatureType , only : temperature_type
    use atm2lndType     , only : atm2lnd_type
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds  ! proc-level bounds
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_inst
    type(temperature_type)      , intent(in)    :: temperature_inst
    type(cnveg_carbonstate_type), intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_state_type)      , intent(inout) :: cnveg_state_inst

    ! !LOCAL VARIABLES:
    real(r8), pointer :: neighbors_count(:)  ! complete grid cell array of count
    real(r8), pointer :: B_id_glob(:)  ! complete grid cell array of B_id
    real(r8), pointer :: B_od_long(:)  ! B_od array for all grid cells
    real(r8), pointer :: B_od_glob(:)  ! complete grid cell array of B_od
    real(r8), pointer :: lodgepole_wtgcell_long(:)  ! same pair for...
    real(r8), pointer :: lodgepole_wtgcell_glob(:)  ! ...lodgepole_wtgcell
    integer :: ng          ! total number of grid cells
    integer :: p, g        ! patch, gridcell indices
    integer :: g_id, g_od  ! gridcell indices in/out-dispersing cells
    integer :: ier         ! error code
    integer :: yr          ! year
    integer :: mon         ! month
    integer :: day         ! day
    integer :: tod         ! seconds

    real(r8) :: Tsprsum     ! average temperature April to July (C)
    real(r8) :: Tfall       ! average temperature September to November (C)
    real(r8) :: Tmin        ! min monthly average temperature Dec to Feb (C)
    real(r8) :: PPTdrought  ! cumulative precipitation since Oct 1 (mm)

    real(r8) :: ag_tree_carbon_perha  ! pre-beetle aboveground carbon (Mg/ha)
    real(r8) :: killed_carbon_perha   ! killed aboveground carbon (Mg/ha)
    real(r8) :: new_ag_tree_carbon_perha  ! post beetle aboveground carbon (Mg/ha)
    real(r8) :: ag_c_1    ! equal to min_ag_carbon_perha
    real(r8) :: shape_g   ! begin: used in calculating the gamma distribution
    real(r8) :: shape_slope
    real(r8) :: shape_offset
    real(r8) :: scale_g
    real(r8) :: scale_slope
    real(r8) :: scale_offset
    real(r8) :: x         ! end: used in calculating the gamma distribution

    real(r8) :: B_endemic_local  ! (beetles)
    real(r8) :: BPT_stressed
    real(r8) :: BPT_a_u     ! (beetles)
    real(r8) :: B_e         ! number of emerging (beetles)
    real(r8) :: B_a         ! number of attacking (beetles)
    real(r8) :: B_a_u       ! (beetles) uncoordinated
    real(r8) :: B_a_c       ! (beetles)   coordinated
    real(r8) :: B_k         ! number that killed trees (beetles)
    real(r8) :: B_t         ! number (beetles)
    real(r8) :: B_od_1      ! number of out-dispersing (beetles) part 1
    real(r8) :: B_od_2      ! number of out-dispersing (beetles) part 2
    real(r8) :: B_od_loop2  ! number of out-dispersing (beetles) 2nd do p loop
    real(r8) :: T_k_max     ! max# trees killed by beetles (trees)
    real(r8) :: T_k_u       ! (trees) killed by uncoordinated beetles
    real(r8) :: T_k_c       ! (trees) killed by   coordinated beetles

    real(r8) :: P_Tfall      ! probability of beetle survival given Tfall
    real(r8) :: P_Twin       ! probability of beetle survival given Tmin
    real(r8) :: P_Tsprsum    ! probability of beetle survival given Tsprsum
    real(r8) :: P_attack_success  ! probability of successful attack
    real(r8) :: P_avail  ! probability of beetle survival given susceptibility of host trees from stand structure (carbon)
    real(r8) :: P_host_suscept_drought  ! probability of beetle survival given susceptibility of host trees from drought
    real(r8) :: P_unacc_beetle_survival_rate   ! probability of unaccounted-for
    real(r8) :: P_unacc_beetle_mortality_rate  ! survival and mortality

    real(r8) :: num_female_eggs_last_year
    real(r8) :: od_calc_max_beetle_pop
    real(r8) :: fraction_outdispersing

    ! !PARAMETERS
    real(r8), parameter :: min_probability = 0.01_r8  ! threshold for various factors
    real(r8), parameter :: wet_min_probability = min_probability
    ! assume that the number of trees in a stand is fixed, and that the trees
    ! grow, which is represented by C
    ! calculate average C in a tree in a stand that is fully stocked with
    ! mature, highly suitable lodgepole pines
    real(r8), parameter :: ave_num_trees_perha = 1113._r8  ! number of trees in stand that is fully stocked with mature, highly suitable lodgepole pines (trees/ha)
    real(r8), parameter :: lower_bound_ag_carbon_perha = 3._r8  ! aboveground carbon (Mg/ha)
    ! min fraction of BPT_healthy required to kill one tree (I can't
    ! find reference to support this)
    real(r8), parameter :: min_frac_BPT_healthy = 0.1_r8  ! require some beetles to kill a tree even in highly susceptible conditions determined by drought and stand structure
    real(r8), parameter :: num_beetles_per_barkarea_kill_lodgepole = 6.04_r8  ! females per ft2; equal to 65/m2; from Raffa and Berryman, Ecol. Monog., 1983; Reid 1963; and Raffa email sent 11/14/17
    real(r8), parameter :: bark_surf_area_average_lodgepole = 91.2_r8  ! ft2; Lynch dissertation, p 92, from DBH = 12.5" = 31.75 cm
    real(r8), parameter :: BPT_healthy = &
     num_beetles_per_barkarea_kill_lodgepole * bark_surf_area_average_lodgepole  ! BPT_healthy = 6.04 * 91.2 = 551, avg number of beetles required to kill 1 tree given drought, stand conditions
    real(r8), parameter :: BPT_min = min_frac_BPT_healthy * BPT_healthy  ! a constant, as coded now; BPT_min = 0.1 * 551 = 55 beetles per tree, min number of beetles required to kill 1 tree given drought, stand conditions
    real(r8), parameter :: overshoot_fraction = 0.24_r8  ! number of beetles actually fill up a tree, which is greater than the required to kill it; from Raffa and Berryman 1983 (according to email from Ken on 11/15/17)
    real(r8), parameter :: B_endemic = 1.e4_r8  ! (beetles/km2)
! From Safranyik and Carroll, Chapter in CFS MPB book, page 39 (60 eggs/female)
! Reid 1962:  up to 75 eggs/female in natural settings, p 610 and 612
! 2/3s of eggs are female, from Reid 1958, cited in Reid 1962, Can. Ent., p 610
    real(r8), parameter :: num_female_eggs_per_female = 60._r8 * 2._r8 / 3._r8
    real(r8), parameter :: od_calc_max_num_trees = 100._r8  ! basis for x value of 2nd point
    real(r8), parameter :: od_calc_max_od_fraction = 0.1_r8  ! y value of 2nd point
    real(r8), parameter :: m = 1.4_r8      ! used in P_attack_success
    real(r8), parameter :: b = -1.6474_r8  ! used in P_attack_success
    real(r8) :: point_1_x  ! parameters defining lines (begin)
    real(r8) :: point_2_x  ! all these have been hardwired directly where...
    real(r8) :: point_1_y  ! ...they get used
    real(r8) :: point_2_y
    real(r8) :: point_3_x
    real(r8) :: point_4_x
    real(r8) :: point_3_y
    real(r8) :: point_4_y
    real(r8) :: slope
    real(r8) :: offset     ! parameters defining lines (end)
    real(r8), parameter :: PPTdrought_x1 = 150._r8
    real(r8), parameter :: PPTdrought_x2 = 500._r8
    real(r8), parameter :: shape_1 = 0.2_r8   ! used in calculating the gamma
    real(r8), parameter :: scale_1 = 0.01_r8  ! distribution
    real(r8), parameter :: ag_c_2 = 60._r8
    real(r8), parameter :: shape_2 = 1.50_r8
    real(r8), parameter :: scale_2 = 0.15_r8

    character(len=*), parameter :: subname = 'dynProgBB'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    associate(                                                                &
         ! Input:  [integer (:)]  pft vegetation type
         ivt                 =>    patch%itype                              , &
         ! Input:  [real(r8) (:)] (gC/m2) aboveground carbon
         abovegroundc        =>    cnveg_carbonstate_inst%abovegroundc_patch, &
         ! Input:  [real(r8) (:)] (K) September to November average temperature
         tbotav_son          =>    atm2lnd_inst%tbotav_son                  , &
         ! Input:  [real(r8) (:)] (K) April to July average temperature
         tbotav_amjj         =>    atm2lnd_inst%tbotav_amjj                 , &
         ! Input:  [real(r8) (:)] (K) December average of daily min temperature
         tminav_dec          =>    temperature_inst%tminav_dec              , &
         ! Input:  [real(r8) (:)] (K) January average of daily min temperature
         tminav_jan          =>    temperature_inst%tminav_jan              , &
         ! Input:  [real(r8) (:)] (K) February average of daily min temperature
         tminav_feb          =>    temperature_inst%tminav_feb              , &
         ! Input:  [real(r8) (:)] (mm) water year precipitation (Oct 1 - Sep 30)
         prec_water_year     =>    atm2lnd_inst%prec_water_year             , &
         ! Output history only: [real(r8) (:)] (beetles) in-dispersing beetles
         B_id                =>    cnveg_state_inst%B_id                    , &
         ! Output history only: [real(r8) (:)] (beetles) out-dispersing beetles
         B_od                =>    cnveg_state_inst%B_od                    , &
         ! Output history only: [real(r8) (:)] (beetles) needed to kill one tree
         BPT                 =>    cnveg_state_inst%BPT                     , &
         ! Output history only: [real(r8) (:)] (trees) available trees
         T_a                 =>    cnveg_state_inst%T_a                     , &
         ! Output history only: [real(r8) (:)] (trees) killed by beetles
         T_k                 =>    cnveg_state_inst%T_k                     , &
         ! Output history only: [real(r8) (:)] (trees) total (max?) number of
         ! trees that could be present; corresponds to T_t in manuscript
         T_total             =>    cnveg_state_inst%T_total                 , &
         ! Output history only: [real(r8) (:)] (Mg C/ha) min amount of carbon
         ! permitted, associated with surviving dominants and subdominants,
         ! saplings, shrubs, herbaceous vegetation
         min_ag_carbon_perha =>    cnveg_state_inst%min_ag_carbon_perha     , &
         ! Output history only: [real(r8) (:)] (carbon per year) beetle induced mortality rate
         beetle_mort_rates_patch => cnveg_state_inst%beetle_mort_rates_patch, &
         ! Input/Output: [real(r8) (:)] (beetles) egg-laying females last year
         num_egglayingfemales_last_year => cnveg_state_inst%num_egglayingfemales_last_year                                                                   , &
         ! Input/Output: [real(r8) (:)] (unitless) random numb. generator seed
         beetle_seed         => cnveg_state_inst%beetle_seed)

    ! Risky to hardwire month, day, and tod values separately here and
    ! separately in subroutines CNPrescBB and UpdateAccVars.
    ! Recommendation: set in one place.
    call get_curr_date(yr, mon, day, tod)
    if (mon == 8 .and. day == 1 .and. tod == 1800) then

       ! Initialize for the mpi_allreduce located between the out-dispersal and
       ! in-dispersal loops
       call get_proc_global(ng=ng)
       ! Variables to gather from all PEs, while between the out-dispersal and
       ! in-dispersal loops
       allocate(B_id_glob(ng))
       allocate(B_od_long(ng))
       allocate(B_od_glob(ng))
       allocate(lodgepole_wtgcell_long(ng))
       allocate(lodgepole_wtgcell_glob(ng))
       ! Initialize to 0 so as to MPI_SUM zeros in all PEs but one per grid cell
       B_od_long(:) = 0._r8
       lodgepole_wtgcell_long(:) = 0._r8
       ! Initialize to 1e36 to help make errors stand out
       B_od_glob(:) = 1.e36_r8
       lodgepole_wtgcell_glob(:) = 1.e36_r8
       ! Initialize to 0 to sum correctly
       B_id_glob(:) = 0._r8

       do p = bounds%begp,bounds%endp  ! first p loop: out-dispersal
          g = patch%gridcell(p)

          ! Variables calculated by the CLM
          ! Aboveground carbon
          ! Change units gC/m2 to Mg C/ha, ie to 1e6 gC/1e4 m2
          ag_tree_carbon_perha = abovegroundc(p) * 1.e-6_r8 * 1.e4_r8

          ! Calculate beetle mortality for lodgepole pine only
          if (patch%wtgcell(p) > 0._r8 .and. ivt(p) == lodgepole .and. ag_tree_carbon_perha > 0._r8) then

             ! Variables calculated by the CLM
             ! Temperature variables: change units from K to degC
             ! Tfall & Tmin = -273.15 degC in year 1 when using arbitrary
             ! initial conditions
             Tsprsum   = tbotav_amjj(g) - SHR_CONST_TKFRZ
             Tfall = tbotav_son(g) - SHR_CONST_TKFRZ
             Tmin = min(tminav_dec(p), tminav_jan(p), &
                        tminav_feb(p)) - SHR_CONST_TKFRZ
             ! Cumulative precipitation since Oct 1
             PPTdrought = prec_water_year(g)

             ! ----------------------------------------------------
             ! MMOBB: Mechanistic Model of Outbreaking Bark Beetles
             !        Starts here
             ! ----------------------------------------------------

             ! from mean of Pfeifer et al., GCB, 2011, Figure 2a (constant in
             ! space and time)
             ! Hicke and Jenkins, FEM, 2008, Figure 2c (=1000) and
             ! Briggs et al., For Sci, 2015 (=1050)
             T_total(p) = ave_num_trees_perha * grc%area(g) * 100._r8 * patch%wtgcell(p)  ! (trees) where area is in km2 and there are 100 ha/km2; constant in time and varies only due to varying gridcell area and pft weight in the gridcell

             ! --------------------------------------
             ! STEP 0. Unaccounted-for mortality/survival
             ! Scale gamma distribution parameters by AG C
             ! More beetle mortality in low C stands, less in high C
             !    designed to allow stands to recover more and to
             !    increase eruptive dynamics/explosive population growth
             min_ag_carbon_perha(p) = min(ag_tree_carbon_perha, &
                                       lower_bound_ag_carbon_perha)
             ag_c_1 = min_ag_carbon_perha(p)
             shape_slope = (shape_2 - shape_1) / (ag_c_2 - ag_c_1)
             scale_slope = (scale_2 - scale_1) / (ag_c_2 - ag_c_1)
             shape_offset = shape_2 - shape_slope * ag_c_2
             scale_offset = scale_2 - scale_slope * ag_c_2
             if (ag_tree_carbon_perha < ag_c_1) then
                shape_g = shape_1
                scale_g = scale_1
             else if (ag_tree_carbon_perha > ag_c_2) then
                shape_g = shape_2
                scale_g = scale_2
             else
                shape_g = shape_offset + &
                          shape_slope * ag_tree_carbon_perha
                scale_g = scale_offset + &
                          scale_slope * ag_tree_carbon_perha
             end if
             ! slevis: What I found out about rgamma at...
             ! http://ugrad.stat.ubc.ca/R/library/base/html/GammaDist.html
             ! rgamma generates random deviates
             ! The Gamma distribution with parameters shape = a and scale = s
             ! has density f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
             ! for a > 0 and s > 0

             ! In Jeff Hicke's R code I found this comment and equation:
             ! results in most values near 0, so this is survival
             ! P_unacc_beetle_survival_rate = max(0.01, min(1, rgamma(1, shape, scale=scale)))
             ! I replaced it with the next few lines, where:
             ! gamma is the gamma function
             ! call random_number(x) returns x values from 0 to 1
             ! call random_seed results in bfb restarts in point simulation

             ! beetle_seed(p,1) == huge(1) at the beginning of a simulation with
             ! no restart data. I offset the initial value to avoid all p
             ! starting the simulation with the same seed.
             if (beetle_seed(p,1) == huge(1)) then
                beetle_seed(p,:) = huge(1) + p
             end if
             call random_seed(put=beetle_seed(p,:))  ! from previous p
             call random_number(x)
             call random_seed(get=beetle_seed(p,:))  ! for next p
!            x = 0.9_r8  ! com. out call random* to complete the PE test

             P_unacc_beetle_survival_rate = max(min_probability, min(1._r8, &
                1._r8 / (scale_g**shape_g * gamma(shape_g)) *               &
                x**(shape_g - 1._r8) * exp(-1._r8 * x / scale_g) ))
             P_unacc_beetle_mortality_rate = 1._r8 - P_unacc_beetle_survival_rate

             ! --------------------------------------
             ! STEP 1. Compute number of beetles required to kill a tree in this
             ! cell in this year. Varies as a function of drought alone.

             ! a) Host susceptibility due to drought
             !    (0 = healthy, 1 = stressed)
             ! PPTdrought will be compared against these parameters
             ! in function 2 below
             if (PPTdrought < PPTdrought_x1) then
                P_host_suscept_drought = 1._r8
             else if ( (PPTdrought >= PPTdrought_x1) .and. &
                       (PPTdrought < PPTdrought_x2) ) then
                P_host_suscept_drought = max(wet_min_probability,  &
                           1._r8 - (PPTdrought - PPTdrought_x1) *  &
                                   (1._r8 - wet_min_probability) / &
                                   (PPTdrought_x2 - PPTdrought_x1))
             else
                P_host_suscept_drought = wet_min_probability
             end if

             ! b) Number of beetles required to kill one tree in this cell in
             !    this year decreases as tree drought stress increases
             BPT_stressed = BPT_healthy * (1. - P_host_suscept_drought)
             BPT(p) = max(BPT_min, BPT_stressed * (1._r8 + overshoot_fraction))

             ! --------------------------------------
             ! STEP 2. Beetle development.
             ! Compute number of emerging beetles this year.

             ! a) Number of female eggs
             B_endemic_local = B_endemic * grc%area(g) * patch%wtgcell(p)
             num_egglayingfemales_last_year(p) = &
              max(B_endemic_local, num_egglayingfemales_last_year(p))
             num_female_eggs_last_year = num_egglayingfemales_last_year(p) * &
                                         num_female_eggs_per_female

             ! b) Climate effects on beetles ranging in phase of development
             ! from eggs (previous fall) through to emerging adults (current
             ! summer)

             ! i) Effects of fall temperature
             ! Three piecewise linear functions
             point_1_x = -2._r8
             point_2_x = 2._r8
             point_1_y = 0._r8
             point_2_y = 1._r8

             point_3_x = 5._r8
             point_4_x = 18._r8
             point_3_y = 1._r8
             point_4_y = 0._r8

             if (Tfall < point_1_x) then
                P_Tfall = min_probability
             else if ( (Tfall >= point_1_x) .and. (Tfall < point_2_x) ) then
                slope = (point_2_y - point_1_y) / (point_2_x - point_1_x)
                offset = point_1_y - slope * point_1_x
                P_Tfall = max(min_probability, min(1._r8, slope * Tfall + offset))
             else if ( (Tfall >= point_2_x) .and. (Tfall < point_3_x) ) then
                P_Tfall = 1._r8
             else if ( (Tfall >= point_3_x) .and. (Tfall < point_4_x) ) then
                slope = (point_4_y - point_3_y) / (point_4_x - point_3_x)
                offset = point_3_y - slope * point_3_x
                P_Tfall = max(min_probability, min(1._r8, slope * Tfall + offset))
             else
                P_Tfall = min_probability
             end if

             ! ii) Effects of winter temperature
             ! One linear function between Points 1 and 2
             ! y = 0.01 for x < Point 1 x
             ! y = 1 for x > Point 2 x
             point_1_x = -18._r8
             point_2_x = -12._r8
             point_1_y = 0._r8
             point_2_y = 1._r8

             if (Tmin < point_1_x) then
                P_Twin = min_probability
             else if ( (Tmin >= point_1_x) .and. (Tmin < point_2_x) ) then
                slope = (point_2_y - point_1_y) / (point_2_x - point_1_x)
                offset = point_1_y - slope * point_1_x
                P_Twin = max(min_probability, min(1._r8, slope * Tmin + offset))
             else if ( (Tmin >= point_2_x) ) then
                P_Twin = 1._r8
             end if

             ! iii) Effects of spring-summer temperature
             point_1_x = 6._r8
             point_2_x = 11._r8
             point_1_y = 0._r8
             point_2_y = 1._r8

             point_3_x = 12._r8
             point_4_x = 24._r8
             point_3_y = 1._r8
             point_4_y = 0._r8

             if (Tsprsum < point_1_x) then
                P_Tsprsum = min_probability
             else if ( (Tsprsum >= point_1_x) .and. (Tsprsum < point_2_x) ) then
                slope = (point_2_y - point_1_y) / (point_2_x - point_1_x)
                offset = point_1_y - slope * point_1_x
                P_Tsprsum = max(min_probability, min(1._r8, slope * Tsprsum + offset))
             else if ( (Tsprsum >= point_2_x) .and. (Tsprsum < point_3_x) ) then
                P_Tsprsum = 1._r8
             else if ( (Tsprsum >= point_3_x) .and. (Tsprsum < point_4_x) ) then
                slope = (point_4_y - point_3_y) / (point_4_x - point_3_x)
                offset = point_3_y - slope * point_3_x
                P_Tsprsum = max(min_probability, min(1._r8, slope * Tsprsum + offset))
             else
                P_Tsprsum = min_probability
             end if

             ! iv) Apply climate effects on emerging beetles
             B_e = num_female_eggs_last_year * P_Tfall * P_Twin * P_Tsprsum * &
                   P_unacc_beetle_survival_rate

             ! --------------------------------------
             ! STEP 3. Compute number of attacking beetles this year from
             ! focal cell and from neighborhood. Also part 1 of the number of
             ! out-dispersing beetles
             ! Fraction of out-dispersing beetles increases linearly with the
             ! beetle population

             ! compute parameters for first out-dispersal, which represents
             ! out-dispersal for some fraction of emerging beetles
             ! include density dependency (fraction that disperses
             ! increases with population density)
             ! define a linear relationship specified by two points:
             ! 0,0 defines one point (0 fraction out-dispersing when
             ! beetle population is 0)
             ! second point:  the fraction of the beetle population
             ! that out-disperses (specified by below parmeter) occurs
             ! when the population >= that needed to kill 100 trees
             ! x value of second point based on current stand conditions
             od_calc_max_beetle_pop = od_calc_max_num_trees * BPT(p)
             fraction_outdispersing = min(od_calc_max_od_fraction,  &
                                          od_calc_max_od_fraction / &
                                          od_calc_max_beetle_pop * B_e)

             ! number of out-dispersing beetles (part 1)
             B_od_1 = floor(B_e * fraction_outdispersing)

             ! number of attacking beetles
             B_a = B_e - B_od_1

             ! --------------------------------------
             ! STEP 4. Compute number of available trees in cell based on
             ! tree/stand structure.
             point_1_x = ag_c_1
             point_2_x = 50._r8
             point_1_y = 0._r8
             point_2_y = 1._r8

             if (ag_tree_carbon_perha < point_1_x) then
                P_avail = 0._r8
             else if (ag_tree_carbon_perha >= point_1_x .and. &
                      ag_tree_carbon_perha < point_2_x) then
                slope = (point_2_y - point_1_y) / (point_2_x - point_1_x)
                offset = point_2_y - slope * point_2_x
                P_avail = max(0._r8, &
                              min(1._r8, slope * ag_tree_carbon_perha + offset))
             else
                P_avail = 1._r8
             end if

             ! Number of available trees
             T_a(p) = T_total(p) * P_avail

             ! --------------------------------------
             ! From here to the end of the do p loop:
             !    Code is duplicated in the in-dispersal loop below
             ! STEP 5. Prepare to compute number of killed trees
             ! a) Max killed trees
             T_k_max = min(T_a(p), T_total(p))  ! seems redundant; T_a enough
             ! b) Attack success probability
             P_attack_success = 1._r8 / &
                               (1._r8 + exp(-(m * log10(B_a / 100._r8) + b)))

             ! --------------------------------------
             ! STEP 6. Compute number of attacked trees
             ! use P_attack_success to decide how many beetles are uncoordinated
             ! versus coordinated

             ! beetles associated with perfectly uncoordinated attack:
             ! beetles dispersed equally among T_a
             B_a_u = B_a * (1._r8 - P_attack_success)  ! uncoordinated

             ! beetles associated with perfectly coordinated attack:
             ! beetles only attack the min number of trees needed to kill them
             B_a_c = B_a * P_attack_success  ! coordinated

             ! how many trees are killed by uncoordinated attacks
             BPT_a_u = B_a_u / T_a(p)  ! equally divided among all trees
             if (BPT_a_u > BPT(p)) then
                T_k_u = min(T_a(p), T_k_max)  ! seems redundant; T_a enough
             else
                T_k_u = 0._r8
             end if

             ! how many trees are killed by coordinated attacks
             T_k_c = floor(B_a_c / BPT(p))  ! the remainder of beetles, those that cannot kill one tree ( B_a_c* %% BPT), die
             T_k_c = min(T_k_max, T_k_c)

             ! combine trees killed by uncoordinated and by coordinated attacks
             T_k(p) = min(T_a(p), T_k_u + T_k_c)
             B_k = T_k(p) * BPT(p)
             B_t = B_k
!            write(iulog,*) 'B_a, B_a_u, B_a_c, BPT_a_u, B_t = '  ! slevis diag
!            write(iulog,*) B_a, B_a_u, B_a_c, BPT_a_u, B_t  ! slevis diag

             ! --------------------------------------
             ! STEP 7. Compute number of out-dispersing beetles > 0 if max
             ! #killed trees is reached and there is a resulting excess of
             ! beetles; else = 0
             B_od_2 = floor(max(0._r8, (B_a - B_t)))

             ! are beetle populations lower than endemic level?
             ! require at least B_endemic beetles to be present
             ! slevis: Require B_t >= B_endemic to simulate T_k > 0
             if (B_t < B_endemic_local) then
                B_e = B_endemic_local
                B_a = B_endemic_local
                B_k = 0._r8
                B_t = B_endemic_local
                T_k(p) = 0._r8
                B_od_1 = 0._r8
                B_od_2 = 0._r8
             end if

             B_od(p) = B_od_1 + B_od_2

             ! --------------------------------------
             ! STEP 8. Compute amount of carbon in killed trees so as to
             ! ensure a minimum amount of aboveground carbon
             ! --------------------------------------
             if (T_a(p) > 0._r8) then
                killed_carbon_perha = ag_tree_carbon_perha * T_k(p) / T_a(p)
!               write(iulog,*) 'Tsprsum, Tfall, t_dec, t_jan, t_feb, prec_water_year, g, p =', Tsprsum, Tfall, tminav_dec(p), tminav_jan(p), tminav_feb(p), prec_water_year(g), g, p  ! slevis diag
!               write(iulog,*) 'T_total, ave_num_trees_perha, area, wtgcell, g, p ='  ! slevis diag
!               write(iulog,*) T_total(p), ave_num_trees_perha, grc%area(g), patch%wtgcell(p), g, p  ! slevis diag
!               write(iulog,*) 'x, P_unacc_beetle_survival_rate, yr, p =', x, P_unacc_beetle_survival_rate, yr, p ! slevis diag
!               write(iulog,*) 'BPT, BPT_min, BPT_stressed, BPT_healthy, P_host_suscept_drought, p =', BPT(p), BPT_min, BPT_stressed, BPT_healthy, P_host_suscept_drought, p  ! slevis diag
!               write(iulog,*) 'num_egglayingfemales_last_year coming in, p =', num_egglayingfemales_last_year(p), p  ! slevis diag
!               write(iulog,*) 'P_Tfall, P_Twin, P_Tsprsum, p =', P_Tfall, P_Twin, P_Tsprsum, p  ! slevis diag
!               write(iulog,*) 'B_a, B_e, B_od_1, p =', B_a, B_e, B_od_1, p  ! slevis diag
!               write(iulog,*) 'T_k, T_a, T_k_u, T_k_c, B_a, P_attack_success, BPT, P_avail, p = '  ! slevis diag
!               write(iulog,*) T_k(p), T_a(p), T_k_u, T_k_c, B_a, P_attack_success, BPT(p), P_avail, p  ! slevis diag
             else
                killed_carbon_perha = 0._r8
             end if
             new_ag_tree_carbon_perha = ag_tree_carbon_perha - killed_carbon_perha
             ! ensure a minimum amount of AG C
             if (new_ag_tree_carbon_perha < min_ag_carbon_perha(p)) then
                ! recalculate important variables
                new_ag_tree_carbon_perha = min_ag_carbon_perha(p)
                killed_carbon_perha = ag_tree_carbon_perha - new_ag_tree_carbon_perha
                T_k(p) = T_a(p) * killed_carbon_perha / ag_tree_carbon_perha
!               write(iulog,*) 'T_k, T_a, killed_carbon_perha, ag_tree_carbon_perha, p ='  ! slevis diag
!               write(iulog,*) T_k(p), T_a(p), killed_carbon_perha, ag_tree_carbon_perha, p  ! slevis diag

                B_k = T_k(p) * BPT(p)  ! number of beetles that killed trees
                B_t = B_endemic_local  ! collapse population to endemic levels
                B_od_2 = floor(max(0._r8, B_a - B_t))
                B_od(p) = B_od_1 + B_od_2
             end if

             if (T_a(p) > 0._r8) then
                ! Pass to the CLM to simulate effect of beetles on the
                ! carbon and nitrogen cycles. Introduced by slevis to mimic
                ! the implementation of the prescribed beetle model.
                ! Cap at 0.99 to ensure regrowth of this pft.
                beetle_mort_rates_patch(p) = min(0.99_r8, T_k(p) / T_a(p))
!               write(iulog,*) 'B_t, beetle_mort_rates_patch, T_k, T_a, g, p leaving 1st do loop ='  ! slevis diag
!               write(iulog,*)  B_t, beetle_mort_rates_patch(p), T_k(p), T_a(p), g, p  ! slevis diag
             end if

             ! Saved in CLM restart files
             num_egglayingfemales_last_year(p) = B_t

             ! Saving for mpi_allreduce coming up between the do p loops
             B_od_long(g) = B_od(p)  ! B_od goes to history
             lodgepole_wtgcell_long(g) = patch%wtgcell(p)

          end if  ! lodgepole with wtgcell > 0 and  aboveground biomass > 0
       end do  ! first p loop: out-dispersal

       ! ------------------------------------------------
       ! Between the out-dispersal and in-dispersal loops
       ! all processors share all the relevant data needed to come up with B_id
       ! Could the same be done in if (iam == 0) followed by call mpi_bcast?
       ! ------------------------------------------------

       ! Gather B_od_glob from B_od_long, ultimately
       !                  from B_od
       ! Gather lodgepole_wtgcell_glob from lodgepole_wtgcell_long, ultimately
       !                               from patch%wtgcell

       call mpi_allreduce(B_od_long, B_od_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)
       call mpi_allreduce(lodgepole_wtgcell_long, lodgepole_wtgcell_glob, ng, &
                          MPI_REAL8, MPI_SUM, mpicom, ier)

       ! From out-dispersing beetles (B_od_glob) get in-dispersing beetles
       ! (B_id) by finding nearest neighbors with the help of ixy and jxy

       allocate(neighbors_count(ng))
       neighbors_count = 0._r8  ! initialize counter vector

       ! Loop to find all out-dispersing beetles' neighbors and
       ! sum each out-dispersing cell's neighbors' lodgepole pine weights so
       ! that neighbors_count is the weighted sum of grid cells with lodgepole
       ! pine > 0
       do g_od = 1, ng
          if (B_od_glob(g_od) > 0._r8) then
             do g_id = 1, ng
                ! identify neighbors with the ixy, jxy indices of grid cells
                if (ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1) then
                   neighbors_count(g_od) = neighbors_count(g_od) + &
                                           lodgepole_wtgcell_glob(g_id)
                end if  ! find surrounding neighbors
             end do  ! g_id loop
          end if  ! B_od > 0
          ! B_id gets beetles from its neighbors
          ! Disperse beetles from out-dispersing cells to in-dispersing cells
          ! Weight out-dispersing beetles by lodgepole pine weights divided by
          ! the neighbors_count sum
          if (neighbors_count(g_od) > 0._r8) then
             do g_id = 1, ng
                if (ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) - 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) + 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id)     .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id) + 1 .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1 .or. &

                    ldecomp%ixy(g_od) == ldecomp%ixy(g_id)     .and.  &
                    ldecomp%jxy(g_od) == ldecomp%jxy(g_id) - 1) then
                   B_id_glob(g_id) = B_id_glob(g_id) + &
                                     lodgepole_wtgcell_glob(g_id) * &
                                     B_od_glob(g_od) / neighbors_count(g_od)
                end if  ! find surrounding neighbors
             end do  ! g_id
          end if  ! neighbors_count > 0
       end do  ! g_od loop

       do p = bounds%begp,bounds%endp  ! second p loop: in-dispersal
          g = patch%gridcell(p)

          B_id(p) = B_id_glob(g)

          ! Calculate beetle mortality for lodgepole pine only
          ! B_id > 0 establishes that wtgcell > 0 and ag_tree_carbon > 0
          ! as a result of earlier if statements
          if (B_id(p) > 0._r8 .and. ivt(p) == lodgepole) then

             ! number of available trees remaining from out-dispersal loop
             T_a(p) = max(0._r8, T_a(p) - T_k(p))
             ! Aboveground carbon calculated by the CLM
             ! Change units gC/m2 to Mg C/ha, ie to 1e6 gC/1e4 m2
             ag_tree_carbon_perha = abovegroundc(p) * 1.e-6_r8 * 1.e4_r8
             ! Aboveground tree carbon remaining from out-dispersal loop
             ag_tree_carbon_perha = ag_tree_carbon_perha * &
                                    (1._r8 - beetle_mort_rates_patch(p))

             ! Warnings
             if (B_od_glob(g) /= B_od(p)) then
                write(iulog,*) 'Warning: Expected these two variables to be equal: B_od_glob, B_od, g, p =', B_od_glob(g), B_od(p), g, p
             end if

             if (lodgepole_wtgcell_glob(g) /= patch%wtgcell(p)) then
                write(iulog,*) 'Warning: Expected these two variables to be equal: lodgepole_wtgcell_glob, patch%wtgcell, g, p =', lodgepole_wtgcell_glob(g), patch%wtgcell(p), g, p
             end if

             ! -----------------------------------------------------------
             ! In this do p loop, the next 6 lines replace STEPs 2-4
             ! of the out-dispersal loop
             ! -----------------------------------------------------------
             ! #in-dispersing beetles from neighboring grid cells
             B_a = B_id(p)
             B_e = 0._r8
             B_od_1 = 0._r8

             ! --------------------------------------
             ! From here to the end of this do p loop:
             !    Code is duplicated from the out-dispersal loop above
             !    REPEAT as in the first do p loop STEP 5 TO END OF LOOP
             ! --------------------------------------
             ! STEP 5. Compute number of killed trees
             ! a) Max killed
             T_k_max = min(T_a(p), T_total(p))  ! seems redundant
             ! b) Attack success probability
             P_attack_success = 1._r8 / &
                               (1._r8 + exp(-(m * log10(B_a / 100._r8) + b)))

             ! --------------------------------------
             ! STEP 6. Compute number of attacked trees
             ! use P_attack_success to decide how many beetles are uncoord.
             ! versus coordinated

             ! beetles associated with perfectly uncoordinated attack:
             ! beetles dispersed equally among T_a
             B_a_u = B_a * (1._r8 - P_attack_success)  ! uncoordinated

             ! beetles associated with perfectly coordinated attack:
             ! beetles only attack min number of trees needed to kill them
             B_a_c = B_a * P_attack_success  ! coordinated

             ! how many trees are killed by uncoordinated attacks
             BPT_a_u = B_a_u / T_a(p)  ! equally divided among all trees
             if (BPT_a_u > BPT(p)) then
                T_k_u = min(T_a(p), T_k_max)
             else
                T_k_u = 0._r8
             end if

             ! how many trees are killed by coordinated attacks
             T_k_c = floor(B_a_c / BPT(p))  ! the remainder of beetles, those that cannot kill one tree ( B_a_c* %% BPT), die
             T_k_c = min(T_k_max, T_k_c)

             ! combine trees killed by uncoord. and by coordinated attacks
             T_k(p) = min(T_a(p), T_k_u + T_k_c)
             B_k = T_k(p) * BPT(p)
             B_t = B_k

             ! --------------------------------------
             ! STEP 7. Compute number of out-dispersing beetles > 0 if max
             ! #killed trees is reached and there is a resulting excess of
             ! beetles; else = 0
             B_od_2 = floor(max(0._r8, (B_a - B_t)))

             ! are beetle populations lower than endemic level?
             ! require at least B_endemic beetles to be present
             ! slevis: Require B_t >= B_endemic to simulate T_k > 0
             B_endemic_local = B_endemic * grc%area(g) * patch%wtgcell(p)
             if (B_t < B_endemic_local) then
                B_e = B_endemic_local  ! set but not used after this line
                B_a = B_endemic_local  ! set but not used after this line
                B_k = 0._r8
                B_t = B_endemic_local
                T_k(p) = 0._r8
                B_od_1 = 0._r8
                B_od_2 = 0._r8
             end if

             B_od_loop2 = B_od_1 + B_od_2

             ! --------------------------------------
             ! STEP 8. Compute amount of carbon in killed trees so as to
             ! ensure a minimum amount of aboveground carbon
             ! --------------------------------------
             if (T_a(p) > 0._r8) then
                killed_carbon_perha = ag_tree_carbon_perha * T_k(p) / T_a(p)
!               write(iulog,*) 'LOOP2: B_a, p =', B_a, p  ! slevis diag
!               write(iulog,*) 'LOOP2: T_k, T_a, T_k_u, T_k_c, B_a, P_attack_success, BPT, p = '  ! slevis diag
!               write(iulog,*) T_k(p), T_a(p), T_k_u, T_k_c, B_a, P_attack_success, BPT(p), p  ! slevis diag
             else
                killed_carbon_perha = 0._r8
             end if
             new_ag_tree_carbon_perha = ag_tree_carbon_perha - killed_carbon_perha
             ! ensure a minimum amount of AG C
             if (new_ag_tree_carbon_perha < min_ag_carbon_perha(p)) then
                ! recalculate important variables
                new_ag_tree_carbon_perha = min_ag_carbon_perha(p)
                killed_carbon_perha = ag_tree_carbon_perha - new_ag_tree_carbon_perha
                T_k(p) = T_a(p) * killed_carbon_perha / ag_tree_carbon_perha
!               write(iulog,*) 'LOOP2: T_k, T_a, killed_carbon_perha, ag_tree_carbon_perha, p ='  ! slevis diag
!               write(iulog,*) T_k(p), T_a(p), killed_carbon_perha, ag_tree_carbon_perha, p  ! slevis diag

                B_k = T_k(p) * BPT(p)  ! number of beetles that killed trees
                B_t = B_endemic_local  ! collapse population to endemic levels
                B_od_2 = floor(max(0._r8, B_a - B_t))
                B_od_loop2 = B_od_1 + B_od_2
             end if

             if (T_a(p) > 0._r8) then
                ! Pass to the CLM to simulate effect of beetles on the
                ! carbon and nitrogen cycles. Introduced by slevis to mimic
                ! the implementation of the prescribed beetle model.
                ! Cap at 0.99 to ensure regrowth of this pft.
                ! Second time through, ie for in-dispersing beetles, I
                ! assume we must add the mortality rates calculated in the
                ! out-dispersal loop
!               write(iulog,*) 'LOOP2: B_t, beetle_mort_rates_patch, T_k, T_a, g, p ='  ! slevis diag
!               write(iulog,*)  B_t, beetle_mort_rates_patch(p), T_k(p), T_a(p), g, p  ! slevis diag
                beetle_mort_rates_patch(p) = min(0.99_r8, beetle_mort_rates_patch(p) + T_k(p) / T_a(p))
             end if

             ! Saved in CLM restart files
             num_egglayingfemales_last_year(p) = B_t

          end if  ! lodgepole and in-dispersing beetles
       end do  ! second p loop: in-dispersal

       deallocate(neighbors_count)
       deallocate(B_id_glob)
       deallocate(B_od_long)
       deallocate(B_od_glob)
       deallocate(lodgepole_wtgcell_long)
       deallocate(lodgepole_wtgcell_glob)

       ! When absolutely satisfied that prognostic beetles work correctly,
       ! please delete these 10 lines:
       ! Hardwired values resulting in bfb same answers as running the case
       ! /glade/p/work/slevis/FMEC/point_runs/presc_bb_lpp1
       ! Uncomment so as to overwrite earlier calculation of the same variable
       !if (yr == 1902) then
       !   beetle_mort_rates_patch(bounds%begp:bounds%endp) = 0.99_r8
       !else
       !   beetle_mort_rates_patch(bounds%begp:bounds%endp) = 0.00_r8
       !end if

    else  ! not Aug 1st's 1800 time step
       ! This works because MMOBB (above) will run only on Aug 1st each yr
       ! Initialize variables that go to history
       B_id(bounds%begp:bounds%endp) = 0._r8
       B_od(bounds%begp:bounds%endp) = 0._r8
       beetle_mort_rates_patch(bounds%begp:bounds%endp) = 0._r8
       ! Initialize other vector variables
       T_k(bounds%begp:bounds%endp) = 0._r8
       T_a(bounds%begp:bounds%endp) = 0._r8
       BPT(bounds%begp:bounds%endp) = 0._r8
    end if  ! Aug 1st's 1800 time step

    end associate

  end subroutine dynProgBB  ! slevis: prognostic beetle module


  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftconMod       , only : noveg, nbrdlf_evr_shrub
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_step_size_real, is_beg_curr_year
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p                         ! patch index
    integer :: g                         ! gridcell index
    integer :: fp                        ! patch filter index
    real(r8):: thistreec                 ! carbon in this tree for calculating harvest fraction (gC/m2)
    real(r8):: cm                        ! rate for carbon harvest mortality (gC/m2/yr)
    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
    real(r8):: dtime                     ! model time step (s)
    !-----------------------------------------------------------------------

    associate(& 
         ivt                                 =>    patch%itype                                                      , & ! Input:  [integer (:)]  pft vegetation type                                
         
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                            , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
         livestemc                           =>    cnveg_carbonstate_inst%livestemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
         deadstemc                           =>    cnveg_carbonstate_inst%deadstemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
         livecrootc                          =>    cnveg_carbonstate_inst%livecrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
         deadcrootc                          =>    cnveg_carbonstate_inst%deadcrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
         xsmrpool                            =>    cnveg_carbonstate_inst%xsmrpool_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
         
         leafn                               =>    cnveg_nitrogenstate_inst%leafn_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
         frootn                              =>    cnveg_nitrogenstate_inst%frootn_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
         livestemn                           =>    cnveg_nitrogenstate_inst%livestemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
         deadstemn                           =>    cnveg_nitrogenstate_inst%deadstemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
         livecrootn                          =>    cnveg_nitrogenstate_inst%livecrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
         deadcrootn                          =>    cnveg_nitrogenstate_inst%deadcrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
         
         hrv_leafc_to_litter                 =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_to_litter                =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_to_litter             =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
         wood_harvestc                       =>    cnveg_carbonflux_inst%wood_harvestc_patch                      , & ! Output: [real(r8) (:)]
         hrv_livecrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         
         hrv_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_to_litter                =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         wood_harvestn                       =>    cnveg_nitrogenflux_inst%wood_harvestn_patch                    , & ! Output: [real(r8) (:)]
         hrv_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_retransn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch      & ! Output: [real(r8) (:)]                                                    
         )

      dtime = get_step_size_real()

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         g = patch%gridcell(p)

         ! If this is a tree pft, then
         ! get the annual harvest "mortality" rate (am) from harvest array
         ! and convert to rate per second
!P.Buotte changed if statement for FMEC tree PFT slots
         if (ivt(p) > noveg .and. ivt(p) < 16) then

            if (do_harvest) then
               if (harvest_units == "gC/m2/yr") then
                  thistreec = leafc(p) + frootc(p) + livestemc(p) + deadstemc(p) + livecrootc(p) + deadcrootc(p) + xsmrpool(p)
                  cm = harvest(g)
                  if (thistreec > 0.0_r8) then
                     am = min(0.98_r8,cm/thistreec)    ! Only harvest up to 98% so regrowth is possible PJL
                  else
                     am = 0._r8
                  end if
               else
                  am = harvest(g)
               end if

               ! Apply all harvest at the start of the year
               if (is_beg_curr_year()) then
                  m  = am/dtime
               else
                  m = 0._r8
               end if
            else
               m = 0._r8
            end if

            ! patch-level harvest carbon fluxes
            ! displayed pools
            hrv_leafc_to_litter(p)               = leafc(p)               * m
            hrv_frootc_to_litter(p)              = frootc(p)              * m
            hrv_livestemc_to_litter(p)           = livestemc(p)           * m
            wood_harvestc(p)                     = deadstemc(p)           * m
            hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
            hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
            hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

            ! storage pools
            hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
            hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
            hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
            hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
            hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
            hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
            hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

            ! transfer pools
            hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
            hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
            hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
            hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
            hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
            hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
            hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

            ! patch-level harvest mortality nitrogen fluxes
            ! displayed pools
            hrv_leafn_to_litter(p)               = leafn(p)               * m
            hrv_frootn_to_litter(p)              = frootn(p)              * m
            hrv_livestemn_to_litter(p)           = livestemn(p)           * m
            wood_harvestn(p)                     = deadstemn(p)           * m
            hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
            hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
            hrv_retransn_to_litter(p)            = retransn(p)            * m

            ! storage pools
            hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
            hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
            hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
            hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
            hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
            hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

            ! transfer pools
            hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
            hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
            hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
            hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
            hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
            hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m

         end if  ! end tree block

      end do ! end of pft loop

      ! gather all patch-level litterfall fluxes from harvest to the column
      ! for litter C and N inputs

      call CNHarvestPftToColumn(num_soilc, filter_soilc, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

    end associate 

  end subroutine CNHarvest

  !PBuotte: begin prescribed BB----------------------------------------------------------
  subroutine CNPrescBB (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, cnveg_state_inst)
    !
    ! !DESCRIPTION:
    ! Prescribed bark beetle mortality routine for coupled carbon-nitrogen code (CN)
    ! Different from harvest in that beetle mortality is applied on one radiation
    ! timestep, not throughout the year.
    ! Mortality occurs on Aug 1, but transfers occur on all radiation timesteps
    !
    ! !USES:
    use pftconMod       , only : noveg, lodgepole !PBuotte: not referencing shrub pft slot to determine if tree PFT. Hardwired last tree slot.
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_days_per_year
    use clm_time_manager, only : get_curr_date
    use dynSubgridControlMod, only : get_do_prog_bb
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_state_type)          , intent(in)    :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p                         ! patch index
    integer :: g                         ! gridcell index
    integer :: fp                        ! patch filter index
    real(r8):: am                        ! rate for fractional beetle mortality (1/yr)
    real(r8):: m                         ! rate for fractional beetle mortality (1/time step)
    real(r8):: days_per_year             ! days per year
    integer :: yr                        ! year
    integer :: mon                       ! month
    integer :: day                       ! day
    integer :: tod                       ! seconds
    integer :: offset                    ! offset

    !-----------------------------------------------------------------------

    associate(&
         ivt                                 =>    patch%itype                                                      , & ! Input:  [integer (:)]  pft vegetation type
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C  
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                            , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C
         livestemc                           =>    cnveg_carbonstate_inst%livestemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C
         deadstemc                           =>    cnveg_carbonstate_inst%deadstemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C
         livecrootc                          =>    cnveg_carbonstate_inst%livecrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C
         deadcrootc                          =>    cnveg_carbonstate_inst%deadcrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C

         xsmrpool                            =>    cnveg_carbonstate_inst%xsmrpool_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer
         leafsnag1c                          =>    cnveg_carbonstate_inst%leafsnag1c_patch                        , &
         leafsnag2c                          =>    cnveg_carbonstate_inst%leafsnag2c_patch                        , &
         leafsnag3c                          =>    cnveg_carbonstate_inst%leafsnag3c_patch                        , &
         snag1c                              =>    cnveg_carbonstate_inst%snag1c_patch                            , &
         snag2c                              =>    cnveg_carbonstate_inst%snag2c_patch                            , &
         snag3c                              =>    cnveg_carbonstate_inst%snag3c_patch                            , &
         snag4c                              =>    cnveg_carbonstate_inst%snag4c_patch                            , &
         snag5c                              =>    cnveg_carbonstate_inst%snag5c_patch                            , &
         snag6c                              =>    cnveg_carbonstate_inst%snag6c_patch                            , &

         leafn                               =>    cnveg_nitrogenstate_inst%leafn_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N  
         frootn                              =>    cnveg_nitrogenstate_inst%frootn_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N
         livestemn                           =>    cnveg_nitrogenstate_inst%livestemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N
         deadstemn                           =>    cnveg_nitrogenstate_inst%deadstemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N
         livecrootn                          =>    cnveg_nitrogenstate_inst%livecrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N
         deadcrootn                          =>    cnveg_nitrogenstate_inst%deadcrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer
         leafsnag1n                          =>    cnveg_nitrogenstate_inst%leafsnag1n_patch                       , &
         leafsnag2n                          =>    cnveg_nitrogenstate_inst%leafsnag2n_patch                       , &
         leafsnag3n                          =>    cnveg_nitrogenstate_inst%leafsnag3n_patch                       , &
         snag1n                              =>    cnveg_nitrogenstate_inst%snag1n_patch                           , &
         snag2n                              =>    cnveg_nitrogenstate_inst%snag2n_patch                           , &
         snag3n                              =>    cnveg_nitrogenstate_inst%snag3n_patch                           , &
         snag4n                              =>    cnveg_nitrogenstate_inst%snag4n_patch                           , &
         snag5n                              =>    cnveg_nitrogenstate_inst%snag5n_patch                           , &
         snag6n                              =>    cnveg_nitrogenstate_inst%snag6n_patch                           , &

         bb_frootc_to_litter                =>    cnveg_carbonflux_inst%bb_frootc_to_litter_patch                , & ! Output: [real(r8) (:)]                   
         bb_leafc_to_leafsnag1              =>    cnveg_carbonflux_inst%bb_leafc_to_leafsnag1c_patch              , & ! Output: [real(r8) (:)]                  

         bb_livestemc_to_snag1              =>    cnveg_carbonflux_inst%bb_livestemc_to_snag1c_patch              , & ! Output: [real(r8) (:)]                  
         bb_deadstemc_to_snag1              =>    cnveg_carbonflux_inst%bb_deadstemc_to_snag1c_patch              , & ! Output: [real(r8) (:)]
         bb_deadstemc_to_prod10c            =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod10c_patch           , & ! Output: [real(r8) (:)]                    
         bb_deadstemc_to_prod100c           =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod100c_patch          , & ! Output: [real(r8) (:)]                    
         bb_liverootc_to_litter            =>    cnveg_carbonflux_inst%bb_liverootc_to_litter_patch           , & ! Output: [real(r8) (:)]                      
         bb_deadrootc_to_litter            =>    cnveg_carbonflux_inst%bb_deadrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                      
         bb_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%bb_xsmrpool_to_atm_patch                , & ! Output: [real(r8) (:)]                    
         leafsnag1c_to_leafsnag2c           =>    cnveg_carbonflux_inst%leafsnag1c_to_leafsnag2c_patch          , &
         leafsnag2c_to_leafsnag3c           =>    cnveg_carbonflux_inst%leafsnag2c_to_leafsnag3c_patch          , &
         leafsnag3c_to_litter               =>    cnveg_carbonflux_inst%leafsnag3c_to_litter_patch              , &
         snag1c_to_snag2c                   =>    cnveg_carbonflux_inst%snag1c_to_snag2c_patch                  , &
         snag2c_to_snag3c                   =>    cnveg_carbonflux_inst%snag2c_to_snag3c_patch                  , &
         snag3c_to_snag4c                   =>    cnveg_carbonflux_inst%snag3c_to_snag4c_patch                  , &
         snag4c_to_snag5c                   =>    cnveg_carbonflux_inst%snag4c_to_snag5c_patch                  , &
         snag5c_to_snag6c                   =>    cnveg_carbonflux_inst%snag5c_to_snag6c_patch                  , &
         snag6c_to_litter                   =>    cnveg_carbonflux_inst%snag6c_to_litter_patch                  , &

         bb_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%bb_frootc_storage_to_litter_patch       , & ! Output: [real(r8) (:)]                    
         bb_leafc_storage_to_leafsnag1      =>    cnveg_carbonflux_inst%bb_leafc_storage_to_leafsnag1c_patch     , & ! Output: [real(r8) (:)]
         bb_livestemc_storage_to_snag1      =>    cnveg_carbonflux_inst%bb_livestemc_storage_to_snag1c_patch     , & ! Output: [real(r8) (:)]
         bb_deadstemc_storage_to_snag1      =>    cnveg_carbonflux_inst%bb_deadstemc_storage_to_snag1c_patch     , & ! Output: [real(r8) (:)]                   
         bb_liverootc_storage_to_litter    =>    cnveg_carbonflux_inst%bb_liverootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                      
         bb_deadrootc_storage_to_litter    =>    cnveg_carbonflux_inst%bb_deadrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                      
         bb_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%bb_gresp_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                    

         bb_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%bb_frootc_xfer_to_litter_patch          , & ! Output: [real(r8) (:)]                    
         bb_leafc_xfer_to_leafsnag1         =>    cnveg_carbonflux_inst%bb_leafc_xfer_to_leafsnag1c_patch        , & ! Output: [real(r8) (:)]                   
         bb_livestemc_xfer_to_snag1         =>    cnveg_carbonflux_inst%bb_livestemc_xfer_to_snag1c_patch        , & ! Output: [real(r8) (:)]                   
         bb_deadstemc_xfer_to_snag1         =>    cnveg_carbonflux_inst%bb_deadstemc_xfer_to_snag1c_patch        , & ! Output: [real(r8) (:)]                   
         bb_liverootc_xfer_to_litter       =>    cnveg_carbonflux_inst%bb_liverootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                      
         bb_deadrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%bb_deadrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                      
         bb_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%bb_gresp_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                    
         leafsnag1n_to_leafsnag2n           =>    cnveg_nitrogenflux_inst%leafsnag1n_to_leafsnag2n_patch          , &
         leafsnag2n_to_leafsnag3n           =>    cnveg_nitrogenflux_inst%leafsnag2n_to_leafsnag3n_patch          , &
         leafsnag3n_to_litter               =>    cnveg_nitrogenflux_inst%leafsnag3n_to_litter_patch              , &
         snag1n_to_snag2n                   =>    cnveg_nitrogenflux_inst%snag1n_to_snag2n_patch                  , &
         snag2n_to_snag3n                   =>    cnveg_nitrogenflux_inst%snag2n_to_snag3n_patch                  , &
         snag3n_to_snag4n                   =>    cnveg_nitrogenflux_inst%snag3n_to_snag4n_patch                  , &
         snag4n_to_snag5n                   =>    cnveg_nitrogenflux_inst%snag4n_to_snag5n_patch                  , &
         snag5n_to_snag6n                   =>    cnveg_nitrogenflux_inst%snag5n_to_snag6n_patch                  , &
         snag6n_to_litter                   =>    cnveg_nitrogenflux_inst%snag6n_to_litter_patch                  , &

         bb_frootn_to_litter                =>    cnveg_nitrogenflux_inst%bb_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                    
         bb_leafn_to_leafsnag1              =>    cnveg_nitrogenflux_inst%bb_leafn_to_leafsnag1n_patch             , & ! Output: [real(r8) (:)]                 
         bb_livestemn_to_snag1              =>    cnveg_nitrogenflux_inst%bb_livestemn_to_snag1n_patch           , & ! Output: [real(r8) (:)]
         bb_deadstemn_to_snag1              =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_snag1n_patch           , & ! Output: [real(r8) (:)]
         bb_deadstemn_to_prod10n            =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod10n_patch         , & ! Output: [real(r8) (:)]                    
         bb_deadstemn_to_prod100n           =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod100n_patch        , & ! Output: [real(r8) (:)]                    
         bb_liverootn_to_litter            =>    cnveg_nitrogenflux_inst%bb_liverootn_to_litter_patch         , & ! Output: [real(r8) (:)]                      
         bb_deadrootn_to_litter            =>    cnveg_nitrogenflux_inst%bb_deadrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                      
         bb_retransn_to_litter              =>    cnveg_nitrogenflux_inst%bb_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                    

         bb_leafn_storage_to_leafsnag1      =>    cnveg_nitrogenflux_inst%bb_leafn_storage_to_leafsnag1n_patch   , & ! Output: [real(r8) (:)]                   
         bb_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%bb_frootn_storage_to_litter_patch     , & ! Output: [real(r8) (:)]                    
         bb_livestemn_storage_to_snag1      =>    cnveg_nitrogenflux_inst%bb_livestemn_storage_to_snag1n_patch   , & ! Output: [real(r8) (:)]                   
         bb_deadstemn_storage_to_snag1      =>    cnveg_nitrogenflux_inst%bb_deadstemn_storage_to_snag1n_patch   , & ! Output: [real(r8) (:)]                   
         bb_liverootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%bb_liverootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                      
         bb_deadrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%bb_deadrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                      
         bb_leafn_xfer_to_leafsnag1         =>    cnveg_nitrogenflux_inst%bb_leafn_xfer_to_leafsnag1n_patch      , & ! Output: [real(r8) (:)]                   
         bb_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%bb_frootn_xfer_to_litter_patch        , & ! Output: [real(r8) (:)]                    
         bb_livestemn_xfer_to_snag1         =>    cnveg_nitrogenflux_inst%bb_livestemn_xfer_to_snag1n_patch      , & ! Output: [real(r8) (:)]                   
         bb_deadstemn_xfer_to_snag1         =>    cnveg_nitrogenflux_inst%bb_deadstemn_xfer_to_snag1n_patch      , & ! Output: [real(r8) (:)]                   
         bb_liverootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%bb_liverootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                      
         bb_deadrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%bb_deadrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                      

         ! Input: [real(r8) (:)] (unitless) beetle induced mortality rate
         beetle_mort_rates_patch => cnveg_state_inst%beetle_mort_rates_patch &
         )

      days_per_year = get_days_per_year()
      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         g = patch%gridcell(p)

         ! If lodgepole pine, do snag pool transfer and beetle mortality
         if (ivt(p) == lodgepole) then
            if (do_presc_bb .or. get_do_prog_bb()) then
               ! Transfer between yearly snag pools on Aug 1
               ! This happens before current year mortality
               ! Transfer total pool at once. Divide by length of time step
               ! here, then multiply back in state update routine
               call get_curr_date(yr, mon, day, tod)
               if(mon == 8 .and. day == 1 .and. tod == 1800) then
                  snag1c_to_snag2c(p)       = snag1c(p) * 1._r8 / 1800
                  snag2c_to_snag3c(p)       = snag2c(p) * 1._r8 / 1800
                  snag3c_to_snag4c(p)       = snag3c(p) * 1._r8 / 1800
                  snag4c_to_snag5c(p)       = snag4c(p) * 1._r8 / 1800
                  snag5c_to_snag6c(p)       = snag5c(p) * 1._r8 / 1800
                  leafsnag1c_to_leafsnag2c(p)   = leafsnag1c(p) * 1._r8 / 1800
                  leafsnag2c_to_leafsnag3c(p)   = leafsnag2c(p) * 1._r8 / 1800
                  snag1n_to_snag2n(p)       = snag1n(p) * 1._r8 / 1800
                  snag2n_to_snag3n(p)       = snag2n(p) * 1._r8 / 1800
                  snag3n_to_snag4n(p)       = snag3n(p) * 1._r8 / 1800
                  snag4n_to_snag5n(p)       = snag4n(p) * 1._r8 / 1800
                  snag5n_to_snag6n(p)       = snag5n(p) * 1._r8 / 1800
                  leafsnag1n_to_leafsnag2n(p)   = leafsnag1n(p) * 1._r8 / 1800
                  leafsnag2n_to_leafsnag3n(p)   = leafsnag2n(p) * 1._r8 / 1800
                else
                  snag1c_to_snag2c(p)       = 0._r8
                  snag2c_to_snag3c(p)       = 0._r8
                  snag3c_to_snag4c(p)       = 0._r8
                  snag4c_to_snag5c(p)       = 0._r8
                  snag5c_to_snag6c(p)       = 0._r8
                  leafsnag1c_to_leafsnag2c(p)   = 0._r8
                  leafsnag2c_to_leafsnag3c(p)   = 0._r8
                  snag1n_to_snag2n(p)       = 0._r8
                  snag2n_to_snag3n(p)       = 0._r8
                  snag3n_to_snag4n(p)       = 0._r8
                  snag4n_to_snag5n(p)       = 0._r8
                  snag5n_to_snag6n(p)       = 0._r8
                  leafsnag1n_to_leafsnag2n(p)   = 0._r8
                  leafsnag2n_to_leafsnag3n(p)   = 0._r8
                end if ! date check

              ! apply decay rate for transfer to litter pool
              ! 10-year half life is: ln(2)/10 = 0.0693147_r8
              ! transfer same amount every timestep, not all at once as
              ! transfer between snag pools
              snag6c_to_litter(p) = snag6c(p) * 0.0693147_r8 / (days_per_year * secspday)
              snag6n_to_litter(p) = snag6n(p) * 0.0693147_r8 / (days_per_year * secspday)
              ! apply decay rate for transfer to litter pool
              ! 1-year half life is: ln(2)/1  = 0.693147_r8
              ! transfer same amount every timestep, not all at once as
              ! transfer between leafsnag pools
              leafsnag3c_to_litter(p) = leafsnag3c(p) * 0.693147_r8 / (days_per_year * secspday)
              leafsnag3n_to_litter(p) = leafsnag3n(p) * 0.693147_r8 / (days_per_year * secspday)

              ! Apply beetle mortality after running the prog. beetle code.
              ! Get the annual beetle mortality rate (am) from beetle array
              ! and convert to rate per length of time step
              if(mon == 8 .and. day == 1 .and. tod == 1800) then
                 ! Prescribed beetle module works with grid cell level variable
                 ! Prognostic beetle module works with pft level variable
                 if (do_presc_bb) then
                    am = beetle_mort_rates(g)
                 else
                    am = beetle_mort_rates_patch(p)
!                   if (beetle_mort_rates_patch(p) > 0._r8) write(iulog,*) 'Outside beetle model now: beetle_mort_rates_patch, p =', beetle_mort_rates_patch(p), p  ! slevis diag
                 end if
                 ! Apply beetle mortality all at once. Divide by seconds/tstep,
                 ! then multiply back in state update routine.
                 m  = am / 1800
              else
                 m = 0._r8
              end if ! second date check
             end if ! do_presc_bb or do_prog_bb check

            ! patch-level beetle mortality carbon fluxes
            ! displayed pools
            bb_leafc_to_leafsnag1(p)            = leafc(p)               * m
            bb_frootc_to_litter(p)              = frootc(p)              * m
            bb_livestemc_to_snag1(p)            = livestemc(p)           * m
            bb_deadstemc_to_snag1(p)            = deadstemc(p)           * m
!commenting this out for now because I don't know why pftcon%pprodharv10 is greater than 0
!will need to fix before we can allocate harvest to beelte kill
! slevis recommnedation: replace 1800s throughout with dt
!            bb_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * pftcon%pprodharv10(ivt(p))
!            bb_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * (1.0_r8 - pftcon%pprodharv10(ivt(p)))
            bb_liverootc_to_litter(p)          = livecrootc(p)          * m
            bb_deadrootc_to_litter(p)          = deadcrootc(p)          * m
            bb_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

            ! storage pools
            bb_leafc_storage_to_leafsnag1(p)    = leafc_storage(p)       * m
            bb_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
            bb_livestemc_storage_to_snag1(p)    = livestemc_storage(p)   * m
            bb_deadstemc_storage_to_snag1(p)    = deadstemc_storage(p)   * m
            bb_liverootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
            bb_deadrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
            bb_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

            ! transfer pools
            bb_leafc_xfer_to_leafsnag1(p)       = leafc_xfer(p)          * m
            bb_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
            bb_livestemc_xfer_to_snag1(p)       = livestemc_xfer(p)      * m
            bb_deadstemc_xfer_to_snag1(p)       = deadstemc_xfer(p)      * m
            bb_liverootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
            bb_deadrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
            bb_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

            ! patch-level beetle mortality nitrogen fluxes
            ! displayed pools
            bb_leafn_to_leafsnag1(p)            = leafn(p)               * m
            bb_frootn_to_litter(p)              = frootn(p)              * m
            bb_livestemn_to_snag1(p)            = livestemn(p)           * m
            bb_deadstemn_to_snag1(p)            = deadstemn(p)           * m
!            bb_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * pftcon%pprodharv10(ivt(p))
!            bb_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * (1.0_r8 - pftcon%pprodharv10(ivt(p)))
            bb_liverootn_to_litter(p)          = livecrootn(p)          * m
            bb_deadrootn_to_litter(p)          = deadcrootn(p)          * m
            bb_retransn_to_litter(p)            = retransn(p)            * m

            ! storage pools
            bb_leafn_storage_to_leafsnag1(p)    = leafn_storage(p)       * m
            bb_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
            bb_livestemn_storage_to_snag1(p)    = livestemn_storage(p)   * m
            bb_deadstemn_storage_to_snag1(p)    = deadstemn_storage(p)   * m
            bb_liverootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
            bb_deadrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

            ! transfer pools
            bb_leafn_xfer_to_leafsnag1(p)       = leafn_xfer(p)          * m
            bb_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
            bb_livestemn_xfer_to_snag1(p)       = livestemn_xfer(p)      * m
            bb_deadstemn_xfer_to_snag1(p)       = deadstemn_xfer(p)      * m
            bb_liverootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
            bb_deadrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         end if  ! end tree block

      end do ! end of pft loop
      ! gather all patch-level litterfall fluxes from beetle mortality to the column
      ! for litter C and N inputs

      call CNBBMortPftToColumn(num_soilc, filter_soilc, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

    end associate

  end subroutine CNPrescBB
 ! PBuotte: end Prescribed BB

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilc, filter_soilc, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all patch-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
   type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
   type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                                   & 
        ivt                              =>    patch%itype                                                      , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_flab                          =>    pftcon%lf_flab                                                 , & ! Input:  leaf litter labile fraction                       
        lf_fcel                          =>    pftcon%lf_fcel                                                 , & ! Input:  leaf litter cellulose fraction                    
        lf_flig                          =>    pftcon%lf_flig                                                 , & ! Input:  leaf litter lignin fraction                       
        fr_flab                          =>    pftcon%fr_flab                                                 , & ! Input:  fine root litter labile fraction                  
        fr_fcel                          =>    pftcon%fr_fcel                                                 , & ! Input:  fine root litter cellulose fraction               
        fr_flig                          =>    pftcon%fr_flig                                                 , & ! Input:  fine root litter lignin fraction                  
        
        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        hrv_leafc_to_litter              =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_to_litter             =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_to_litter          =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
        pwood_harvestc                   =>    cnveg_carbonflux_inst%wood_harvestc_patch                      , & ! Input:  [real(r8) (:)   ]
        hrv_livecrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        cwood_harvestc                   =>    cnveg_carbonflux_inst%wood_harvestc_col                        , & ! InOut:  [real(r8) (:)   ]
        harvest_c_to_litr_met_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_met_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
        harvest_c_to_litr_cel_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_cel_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
        harvest_c_to_litr_lig_c          =>    cnveg_carbonflux_inst%harvest_c_to_litr_lig_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
        harvest_c_to_cwdc                =>    cnveg_carbonflux_inst%harvest_c_to_cwdc_col                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
        
        hrv_leafn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_to_litter          =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        pwood_harvestn                   =>    cnveg_nitrogenflux_inst%wood_harvestn_patch                    , & ! Input:  [real(r8) (:)   ]
        hrv_livecrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransn_to_litter           =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_storage_to_litter      =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_xfer_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        cwood_harvestn                   =>    cnveg_nitrogenflux_inst%wood_harvestn_col                      , & ! InOut:  [real(r8) (:)   ]
        harvest_n_to_litr_met_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_met_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
        harvest_n_to_litr_cel_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_cel_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
        harvest_n_to_litr_lig_n          =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_lig_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
        harvest_n_to_cwdn                =>    cnveg_nitrogenflux_inst%harvest_n_to_cwdn_col                    & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        )

     do j = 1, nlevdecomp
        do pi = 1,maxpatch_pft
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col%npatches(c)) then
                 p = col%patchi(c) + pi - 1

                 if (patch%active(p)) then

                    ! leaf harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                         hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                         hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood harvest mortality carbon fluxes
                    harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                         hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j) 
                    harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                         hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                         hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

                    ! storage harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                         hrv_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                    ! transfer harvest mortality carbon fluxes
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                         hrv_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                         hrv_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                         hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                         hrv_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                    ! leaf harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                         hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root litter nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                         hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood harvest mortality nitrogen fluxes
                    harvest_n_to_cwdn(c,j)  = harvest_n_to_cwdn(c,j)  + &
                         hrv_livestemn_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                         hrv_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                         hrv_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! retranslocated N pool harvest mortality fluxes
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                    ! storage harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                         hrv_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                         hrv_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_livestemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_deadstemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! transfer harvest mortality nitrogen fluxes
                    harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                         hrv_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                         hrv_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_livestemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                         hrv_deadstemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                         hrv_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)

                 end if
              end if

           end do

        end do
     end do
   
     do pi = 1,maxpatch_pft
        do fc = 1,num_soilc
           c = filter_soilc(fc)

           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1

              if (patch%active(p)) then
                 ! wood harvest mortality carbon fluxes to product pools
                 cwood_harvestc(c)  = cwood_harvestc(c)  + &
                      pwood_harvestc(p)  * wtcol(p)

                 ! wood harvest mortality nitrogen fluxes to product pools
                 cwood_harvestn(c)  = cwood_harvestn(c)  + &
                      pwood_harvestn(p)  * wtcol(p)
              end if
           end if

        end do

     end do

   end associate 

 end subroutine CNHarvestPftToColumn

 ! PBuotte: begin Prescribed BB
 !------------------------------------------------------------------------------------
 subroutine CNBBMortPftToColumn (num_soilc, filter_soilc, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! called at the end of CNPrescBB gather all patch-level beetle mortality litter fluxes
   ! to the column level and assign them to the three litter pools and one CWD pool
   !
   ! !USES:
   use clm_varpar , only : maxpatch_pft, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
   type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
   type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                                   &
        ivt                              =>    patch%itype                                                      , & ! Input:  [integer  (:)   ]  pft vegetation type
        wtcol                            =>    patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)

        lf_flab                          =>    pftcon%lf_flab                                                 , & ! Input:  leaf litter labile fraction         
        lf_fcel                          =>    pftcon%lf_fcel                                                 , & ! Input:  leaf litter cellulose fraction      
        lf_flig                          =>    pftcon%lf_flig                                                 , & ! Input:  leaf litter lignin fraction         
        fr_flab                          =>    pftcon%fr_flab                                                 , & ! Input:  fine root litter labile fraction    
        fr_fcel                          =>    pftcon%fr_fcel                                                 , & ! Input:  fine root litter cellulose fraction 
        fr_flig                          =>    pftcon%fr_flig                                                 , & ! Input:  fine root litter lignin fraction    

        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems

        leafsnag3c_to_litter           =>    cnveg_carbonflux_inst%leafsnag3c_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                       
        snag6c_to_litter             =>    cnveg_carbonflux_inst%snag6c_to_litter_patch                  , & ! Input: [real(r8) (:)   ]
        bb_frootc_to_litter             =>    cnveg_carbonflux_inst%bb_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemc_to_snag1           =>    cnveg_carbonflux_inst%bb_livestemc_to_snag1c_patch             , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemc_to_snag1           =>    cnveg_carbonflux_inst%bb_deadstemc_to_snag1c_patch             , & ! Input: [real(r8) (:)   ]
        pbb_deadstemc_to_prod10c        =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod10c_patch           , & ! Input:  [real(r8) (:)   ]                     
        pbb_deadstemc_to_prod100c       =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod100c_patch          , & ! Input:  [real(r8) (:)   ]                     
        bb_liverootc_to_litter         =>    cnveg_carbonflux_inst%bb_liverootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootc_to_litter         =>    cnveg_carbonflux_inst%bb_deadrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                       
        bb_leafc_storage_to_leafsnag1   =>    cnveg_carbonflux_inst%bb_leafc_storage_to_leafsnag1c_patch     , & ! Input:  [real(r8) (:)   ]                    
        bb_frootc_storage_to_litter     =>    cnveg_carbonflux_inst%bb_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemc_storage_to_snag1   =>    cnveg_carbonflux_inst%bb_livestemc_storage_to_snag1c_patch     , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemc_storage_to_snag1   =>    cnveg_carbonflux_inst%bb_deadstemc_storage_to_snag1c_patch     , & ! Input:  [real(r8) (:)   ]                    
        bb_liverootc_storage_to_litter =>    cnveg_carbonflux_inst%bb_liverootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootc_storage_to_litter =>    cnveg_carbonflux_inst%bb_deadrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                       
        bb_gresp_storage_to_litter      =>    cnveg_carbonflux_inst%bb_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                     
        bb_leafc_xfer_to_leafsnag1      =>    cnveg_carbonflux_inst%bb_leafc_xfer_to_leafsnag1c_patch        , & ! Input:  [real(r8) (:)   ]                    
        bb_frootc_xfer_to_litter        =>    cnveg_carbonflux_inst%bb_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemc_xfer_to_snag1      =>    cnveg_carbonflux_inst%bb_livestemc_xfer_to_snag1c_patch        , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemc_xfer_to_snag1      =>    cnveg_carbonflux_inst%bb_deadstemc_xfer_to_snag1c_patch        , & ! Input:  [real(r8) (:)   ]                    
        bb_liverootc_xfer_to_litter    =>    cnveg_carbonflux_inst%bb_liverootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%bb_deadrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                       
        bb_gresp_xfer_to_litter         =>    cnveg_carbonflux_inst%bb_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                     
        cbb_deadstemc_to_prod10c        =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod10c_col             , & ! InOut:  [real(r8) (:)   ]                     
        cbb_deadstemc_to_prod100c       =>    cnveg_carbonflux_inst%bb_deadstemc_to_prod100c_col            , & ! InOut:  [real(r8) (:)   ]                     
        beetle_c_to_litr_met_c          =>    cnveg_carbonflux_inst%beetle_c_to_litr_met_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with beetles to litter metabolic pool (gC/m3/s)
        beetle_c_to_litr_cel_c          =>    cnveg_carbonflux_inst%beetle_c_to_litr_cel_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with beetles to litter cellulose pool (gC/m3/s)
        beetle_c_to_litr_lig_c          =>    cnveg_carbonflux_inst%beetle_c_to_litr_lig_c_col              , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with beetles to litter lignin pool (gC/m3/s)
        beetle_c_to_cwdc                =>    cnveg_carbonflux_inst%beetle_c_to_cwdc_col                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with beetles to CWD pool (gC/m3/s)

       ! Nitrogen
        leafsnag3n_to_litter           =>    cnveg_nitrogenflux_inst%leafsnag3n_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                       
        snag6n_to_litter                =>    cnveg_nitrogenflux_inst%snag6n_to_litter_patch               , & ! Input [real(r8) (:) ]
        bb_frootn_to_litter             =>    cnveg_nitrogenflux_inst%bb_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemn_to_snag1           =>    cnveg_nitrogenflux_inst%bb_livestemn_to_snag1n_patch           , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemn_to_snag1           =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_snag1n_patch           , & ! Input:  [real(r8) (:)   ]
        pbb_deadstemn_to_prod10n        =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod10n_patch         , & ! Input:  [real(r8) (:)   ]                     
        pbb_deadstemn_to_prod100n       =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod100n_patch        , & ! Input:  [real(r8) (:)   ]                     
        bb_liverootn_to_litter         =>    cnveg_nitrogenflux_inst%bb_liverootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootn_to_litter         =>    cnveg_nitrogenflux_inst%bb_deadrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                       
        bb_retransn_to_litter           =>    cnveg_nitrogenflux_inst%bb_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                     
        bb_leafn_storage_to_leafsnag1   =>    cnveg_nitrogenflux_inst%bb_leafn_storage_to_leafsnag1n_patch   , & ! Input:  [real(r8) (:)   ]                    
        bb_frootn_storage_to_litter     =>    cnveg_nitrogenflux_inst%bb_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemn_storage_to_snag1   =>    cnveg_nitrogenflux_inst%bb_livestemn_storage_to_snag1n_patch   , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemn_storage_to_snag1   =>    cnveg_nitrogenflux_inst%bb_deadstemn_storage_to_snag1n_patch   , & ! Input:  [real(r8) (:)   ]                    
        bb_liverootn_storage_to_litter =>    cnveg_nitrogenflux_inst%bb_liverootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%bb_deadrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                       
        bb_leafn_xfer_to_leafsnag1      =>    cnveg_nitrogenflux_inst%bb_leafn_xfer_to_leafsnag1n_patch      , & ! Input:  [real(r8) (:)   ]                    
        bb_frootn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%bb_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                     
        bb_livestemn_xfer_to_snag1      =>    cnveg_nitrogenflux_inst%bb_livestemn_xfer_to_snag1n_patch      , & ! Input:  [real(r8) (:)   ]                    
        bb_deadstemn_xfer_to_snag1      =>    cnveg_nitrogenflux_inst%bb_deadstemn_xfer_to_snag1n_patch      , & ! Input:  [real(r8) (:)   ]                    
        bb_liverootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%bb_liverootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                       
        bb_deadrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%bb_deadrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                       
        cbb_deadstemn_to_prod10n        =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod10n_col           , & ! InOut:  [real(r8) (:)   ]                     
        cbb_deadstemn_to_prod100n       =>    cnveg_nitrogenflux_inst%bb_deadstemn_to_prod100n_col          , & ! InOut:  [real(r8) (:)   ]                     
        beetle_n_to_litr_met_n          =>    cnveg_nitrogenflux_inst%beetle_n_to_litr_met_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
        beetle_n_to_litr_cel_n          =>    cnveg_nitrogenflux_inst%beetle_n_to_litr_cel_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
        beetle_n_to_litr_lig_n          =>    cnveg_nitrogenflux_inst%beetle_n_to_litr_lig_n_col            , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
        beetle_n_to_cwdn                =>    cnveg_nitrogenflux_inst%beetle_n_to_cwdn_col                    & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        )

     do j = 1, nlevdecomp
        do pi = 1,maxpatch_pft
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col%npatches(c)) then
                 p = col%patchi(c) + pi - 1

                 if (patch%active(p)) then

                    ! leaf beetle mortality carbon fluxes
                    ! leaves go to leafsnag1 pool first. Leaf fall to litter begins on 3rd year after mortality,
                    beetle_c_to_litr_met_c(c,j) = beetle_c_to_litr_met_c(c,j) + &
                         leafsnag3c_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    beetle_c_to_litr_cel_c(c,j) = beetle_c_to_litr_cel_c(c,j) + &
                         leafsnag3c_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    beetle_c_to_litr_lig_c(c,j) = beetle_c_to_litr_lig_c(c,j) + &
                         leafsnag3c_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root beetle mortality carbon fluxes to litter
                    beetle_c_to_litr_met_c(c,j) = beetle_c_to_litr_met_c(c,j) + &
                         bb_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    beetle_c_to_litr_cel_c(c,j) = beetle_c_to_litr_cel_c(c,j) + &
                         bb_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    beetle_c_to_litr_lig_c(c,j) = beetle_c_to_litr_lig_c(c,j) + &
                         bb_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood beetle mortality carbon fluxes to cwd
                    ! Roots and transfer from last snag pool go here. Stems go to snag1 first
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_liverootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_deadrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         snag6c_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! storage beetle mortality carbon fluxes
                    ! PBuotte: Not sure this is necessary because for evergreen PFT,
                    ! storage and transfer pools are allocated to display at
                    ! beginning of each timestep, I think, so storage and transfer pools
                    ! should be empty. Beetle mortality does not occur in
                    ! non-evergreen PFTs. Keeping in case I'm wrong.

                    ! c storage fluxes to litter
                    beetle_c_to_litr_met_c(c,j)     = beetle_c_to_litr_met_c(c,j)     + &
                         bb_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    beetle_c_to_litr_met_c(c,j)      = beetle_c_to_litr_met_c(c,j)      + &
                         bb_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                    ! c storage fluxes to CWD
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_liverootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_deadrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! c transfer fluxes to litter
                    beetle_c_to_litr_met_c(c,j)     = beetle_c_to_litr_met_c(c,j)     + &
                         bb_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                    beetle_c_to_litr_met_c(c,j)      = beetle_c_to_litr_met_c(c,j)      + &
                         bb_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)

                    ! c transfer fluxes to CWD
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_liverootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_c_to_cwdc(c,j) = beetle_c_to_cwdc(c,j) + &
                         bb_deadrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)


                    ! leaf beetle mortality nitrogen fluxes
                    beetle_n_to_litr_met_n(c,j) = beetle_n_to_litr_met_n(c,j) + &
                         leafsnag3n_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    beetle_n_to_litr_cel_n(c,j) = beetle_n_to_litr_cel_n(c,j) + &
                         leafsnag3n_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                    beetle_n_to_litr_lig_n(c,j) = beetle_n_to_litr_lig_n(c,j) + &
                         leafsnag3n_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                    ! fine root litter nitrogen fluxes
                    beetle_n_to_litr_met_n(c,j) = beetle_n_to_litr_met_n(c,j) + &
                         bb_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    beetle_n_to_litr_cel_n(c,j) = beetle_n_to_litr_cel_n(c,j) + &
                         bb_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                    beetle_n_to_litr_lig_n(c,j) = beetle_n_to_litr_lig_n(c,j) + &
                         bb_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                    ! wood beetle mortality nitrogen fluxes
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_liverootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_deadrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         snag6n_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! retranslocated N pool beetle mortality fluxes
                    beetle_n_to_litr_met_n(c,j) = beetle_n_to_litr_met_n(c,j) + &
                         bb_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                    ! n storage fluxes to litter
                    beetle_n_to_litr_met_n(c,j)     = beetle_n_to_litr_met_n(c,j)     + &
                         bb_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)

                    ! n storage fluxes to CWD
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_liverootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_deadrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! n transfer fluxes to litter
                    beetle_n_to_litr_met_n(c,j)     = beetle_n_to_litr_met_n(c,j)     + &
                         bb_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)

                    ! n transfer fluxes to CWD
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_liverootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    beetle_n_to_cwdn(c,j) = beetle_n_to_cwdn(c,j) + &
                         bb_deadrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)


                 end if
              end if

           end do

        end do
     end do

     do pi = 1,maxpatch_pft
        do fc = 1,num_soilc
           c = filter_soilc(fc)

           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1

              if (patch%active(p)) then

                 ! wood harvest mortality carbon fluxes to product pools
                 cbb_deadstemc_to_prod10c(c)  = cbb_deadstemc_to_prod10c(c)  + &
                      pbb_deadstemc_to_prod10c(p)  * wtcol(p)
                 cbb_deadstemc_to_prod100c(c)  = cbb_deadstemc_to_prod100c(c)  + &
                      pbb_deadstemc_to_prod100c(p)  * wtcol(p)

                 ! wood harvest mortality nitrogen fluxes to product pools
                 cbb_deadstemn_to_prod10n(c)  = cbb_deadstemn_to_prod10n(c)  + &
                      pbb_deadstemn_to_prod10n(p)  * wtcol(p)
                 cbb_deadstemn_to_prod100n(c)  = cbb_deadstemn_to_prod100n(c)  + &
                      pbb_deadstemn_to_prod100n(p)  * wtcol(p)

              end if
           end if

        end do

     end do

   end associate
  end subroutine CNBBMortPftToColumn
! PBuotte: end Prescribed BB

end module dynHarvestMod
