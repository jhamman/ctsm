module lnd_comp_esmf

  !----------------------------------------------------------------------------
  This is the ESMF cap for CTSM
  !----------------------------------------------------------------------------

  use ESMF
  use shr_kind_mod          , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod           , only : shr_sys_abort
  use shr_log_mod           , only : shr_log_Unit
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel, shr_file_getUnit
  use shr_orb_mod           , only : shr_orb_decl
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date

  use lnd_import_export     , only : advertise_fields, realize_fields
  use lnd_import_export     , only : import_fields, export_fields
  use spmdMod               , only : masterproc, mpicom, spmd_init
  use decompMod             , only : bounds_type, ldecomp, get_proc_bounds
  use domainMod             , only : ldomain
  use controlMod            , only : control_setNL
  use clm_varorb            , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl            , only : inst_index, inst_suffix, inst_name
  use clm_varctl            , only : single_column, clm_varctl_set, iulog
  use clm_varctl            , only : nsrStartup, nsrContinue, nsrBranch
  use clm_time_manager      , only : set_timemgr_init, advance_timestep
  use clm_time_manager      , only : set_nextsw_cday, update_rad_dtime
  use clm_time_manager      , only : get_nstep, get_step_size
  use clm_time_manager      , only : get_curr_date, get_curr_calday
  use clm_initializeMod     , only : initialize1, initialize2
  use clm_driver            , only : clm_drv
  use perf_mod              , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  ! Module routines
  public  :: land_setvm
  public  :: land_register

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer, parameter         :: dbug = 1
  character(*),parameter     :: modName =  "(lnd_comp_esmf)"
  character(*),parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine land_setvm(comp, rc)
    type(ESMF_GridComp) :: comp
    integer, intent(out) :: rc

    ! Initialize return code
    rc = ESMF_SUCCESS

  end subroutine


  subroutine land_register(comp, rc)
    type(ESMF_GridComp) :: comp
    integer, intent(out) :: rc

    ! Initialize return code
    rc = ESMF_SUCCESS

    print *, "Land Register starting"

    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=land_init, &
      rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=land_run, &
      rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=land_final, &
      rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    print *, "Registered Initialize, Run, and Finalize routines"
    print *, "Land Register returning"
    
  end subroutine


  subroutine land_init(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: comp
    type(ESMF_State) :: importState, exportState
    type(ESMF_Clock) :: clock
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_ArraySpec)  :: arrayspec
    type(ESMF_DistGrid)   :: distgrid
    type(ESMF_Grid)       :: grid
    type(ESMF_Field)      :: field
    type(ESMF_VM)         :: vm
    integer               :: petCount, localPet
    integer               :: ierr
    integer               :: compid      ! component id

    ! Initialize return codes
    rc = ESMF_SUCCESS
    ierr = MPI_SUCCESS

    ! Determine petCount
    call ESMF_GridCompGet(comp, vm=vm, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    
    print *, "Land Init starting, localPet =", localPet
    
    !----------------------------------------------------------------------------
    ! initialize CTSM MPI stuff
    !----------------------------------------------------------------------------
    call mpi_comm_dup(lmpicom, mpicom, ierr)  ! duplicate mpicom
    ! check error?

    ! TODO: get compid from comp
    compid = 1

    call spmd_init(mpicom, compid)  ! MPI initialization


    ! Create the source Field and add it to the export State
    call ESMF_ArraySpecSet(arrayspec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    !
    grid = make_grid_sph(45,120,8.,1.5,0.,-90., msx=200., mex=360., msy=-90., mey=90., &
          maskvalue=4, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    !
    field = ESMF_FieldCreate(arrayspec=arrayspec, grid=grid, &
      indexflag=ESMF_INDEX_GLOBAL, name="F_lnd", rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_StateAdd(exportState, (/field/), rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    call ESMF_StatePrint(exportState, rc=rc)
   
    print *, "land Init returning"

  end subroutine land_init


  !-------------------------------------------------------------------------
!   !  The Run routine where data is computed.
!   !
 
  subroutine land_run(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: comp
    type(ESMF_State) :: importState, exportState
    type(ESMF_Clock) :: clock
    integer, intent(out) :: rc

    ! Local variables
    real(ESMF_KIND_R8)    :: pi, kx, ky
    type(ESMF_Field)      :: field
    real(ESMF_KIND_R8), pointer :: farrayPtr(:,:)   ! matching F90 array pointer
    integer               :: i, j, elb(2), eub(2)
    
    ! Initialize return code
    rc = ESMF_SUCCESS

    print *, "Land Run starting"

    pi = 3.14159d0

    ! Get the source Field from the export State
    call ESMF_StateGet(exportState, "F_lnd", field, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    ! Gain access to actual data via F90 array pointer
    call ESMF_FieldGet(field, localDe=0, farrayPtr=farrayPtr, &
      exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    ! Fill source Field with data
    kx = 2.*pi/(eub(1)-elb(1))
    ky = 2.*pi/(eub(2)-elb(2))
    do i = elb(1), eub(1)
      do j = elb(2), eub(2)
        farrayPtr(i,j) = cos(kx*i)*sin(ky*j)
      enddo
    enddo
 
    print *, "Land Run returning"

  end subroutine land_run


!-------------------------------------------------------------------------
!   !  The Finalization routine where things are deleted and cleaned up.
!   !
 
  subroutine land_final(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: comp
    type(ESMF_State) :: importState, exportState
    type(ESMF_Clock) :: clock
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_Grid) :: grid
    type(ESMF_Field) :: field
    
    ! Initialize return code
    rc = ESMF_SUCCESS

    print *, "Land Final starting"

    call ESMF_StateGet(exportState, "F_lnd", field, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_FieldGet(field, grid=grid, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_FieldDestroy(field, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out
    call ESMF_GridDestroy(grid, rc=rc)
    if (rc/=ESMF_SUCCESS) return ! bail out

    print *, "Land Final returning"

  end subroutine land_final


end module lnd_comp_esmf
