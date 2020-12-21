! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
! IVT caculation from:
! equation 8a and 8b from Zhu and Newell (1998)
! https://journals.ametsoc.org/doi/pdf/10.1175/1520-0493%281998%29126%3C0725%3AAPAFMF%3E2.0.CO%3B2
!
! Ming Ge: reference isobaric_diagnostics.F
! Oct 2018

module ivt

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_LOG_ERR, MPAS_LOG_CRIT
    use mpas_kind_types, only : RKIND
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_clock_type), pointer :: clock

    public :: ivt_setup, &
              ivt_compute, &

    private


    contains


    !-----------------------------------------------------------------------
    !  routine ivt_setup
    !
    !> \brief Initialize the diagnostic
    !> \author Ming Ge
    !> \date   Oct. 2018
    !> \details
    !>  Initialize the IVT diagnostic. Copy subroutine isobaric_diagnostics_setup 
    !
    !-----------------------------------------------------------------------
    subroutine ivt_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        ! Perform initialization, memory allocation, etc.

        call mpas_pool_get_subpool(all_pools, 'mesh',  mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag',  diag)

        clock => simulation_clock
   
    end subroutine ivt_setup


    !-----------------------------------------------------------------------
    !  routine ivt_compute
    !
    !> \brief Compute diagnostic before model output is written
    !> \author Ming Ge
    !> \date   Oct. 2018
    !> \details
    !>  Compute IVT diagnostic before model output is written
    !
    !-----------------------------------------------------------------------
    subroutine ivt_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

        ! Called immediately before diagnostics will be written
        ! Compute the diagnostic
       
       integer :: iCell, k
       real (kind=RKIND), parameter :: g = 9.81
       real (kind=RKIND) :: dp, qv_avg, uzonal_avg, umeridional_avg
       real (kind=RKIND), dimension(:),   pointer :: ivt_uzonal, ivt_umeridional
       real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p
       real (kind=RKIND), dimension(:,:), pointer :: uzonal, umeridional, qv
       real (kind=RKIND), dimension(:,:,:),pointer:: scalars

       integer, pointer :: nCells, nVertLevels
       integer, pointer :: index_qv
       
       logical :: need_ivt_umeridional, need_ivt_uzonal, need_any_diags

       need_any_diags       = .false.
       need_ivt_uzonal      = .false.
       need_ivt_umeridional = .false.

       need_ivt_uzonal = MPAS_field_will_be_written('ivt_uzonal')
       need_any_diags = need_any_diags .or. need_ivt_uzonal
       need_ivt_umeridional = MPAS_field_will_be_written('ivt_umeridional')
       need_any_diags = need_any_diags .or. need_ivt_umeridional
       
       if (need_any_diags) then
         call mpas_pool_get_dimension(mesh,  'nCells',       nCells)
         call mpas_pool_get_dimension(mesh,  'nVertLevels',  nVertLevels)
         call mpas_pool_get_dimension(state, 'index_qv',     index_qv) 

         call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
         call mpas_pool_get_array(diag, 'uReconstructZonal',      uzonal) 
         call mpas_pool_get_array(diag, 'ivt_uzonal',             ivt_uzonal)
         call mpas_pool_get_array(diag, 'ivt_umeridional',        ivt_umeridional)
         call mpas_pool_get_array(diag, 'pressure_base',          pressure_b)
         call mpas_pool_get_array(diag, 'pressure_p',             pressure_p)

         call mpas_pool_get_array(state,'scalars', scalars, 1)

         ivt_uzonal(:)      = 0.0
         ivt_umeridional(:) = 0.0

         do iCell = 1, nCells
           do k = 2, nVertLevels 
             uzonal_avg      = 0.5*(uzonal     (k,iCell) + uzonal     (k-1,iCell))
             umeridional_avg = 0.5*(umeridional(k,iCell) + umeridional(k-1,iCell)) 
             qv_avg          = 0.5*(scalars(index_qv,k,iCell) +  scalars(index_qv, k-1, iCell))
             dp = ((pressure_p(k,iCell) + pressure_b(k,iCell))-(pressure_p(k-1,iCell) + pressure_b(k-1,iCell)))
             ivt_uzonal(iCell)      = ivt_uzonal     (iCell) + qv_avg*uzonal_avg*dp/g
             ivt_umeridional(iCell) = ivt_umeridional(iCell) + qv_avg*umeridional_avg*dp/g
           end do
         end do

      end if
    end subroutine ivt_compute

end module ivt
