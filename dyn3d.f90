!------------------------------------------------------------------------
!  Program     :            dyn3d
!------------------------------------------------------------------------
!  Purpose      : The main routine
!
!  Details      ：
!
!  Input        :
!
!  Input/output :
!
!  Output       :
!
!  Remarks      : This program is written based on the Matlab version by
!                 Prof. Eldredge
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2017 Aug
!------------------------------------------------------------------------
PROGRAM dyn3d

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_init_system
    USE module_write_structure
    USE module_artic_rhs_3d

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  ARGUMENT
    !--------------------------------------------------------------------
    REAL(dp),DIMENSION(:),ALLOCATABLE         :: y_init

    ! add_body, add_joint and assemble them
    CALL config_3d_hinged

    ! initialize system
    ALLOCATE(y_init(2*system%np))
    CALL init_system(y_init)
    WRITE(*,*) y_init
    DEALLOCATE(y_init)

    ! construct rhs for the first timestep
!    system%soln%t(1) = system%params%dt
!    CALL artic_rhs_3d(system%soln%t(1),y_init,system%soln%y)

    ! write data
    CALL write_structure

END PROGRAM dyn3d