!------------------------------------------------------------------------
!  Subroutine     :          write_Matlab_plot
!------------------------------------------------------------------------
!  Purpose      : This subroutine writes the position of every body verts_i
!                 at regularized time, for Matlab plot.
!
!  Details      ：
!
!  Input        :
!  Input/output :
!
!  Output       :
!
!  Remarks      : One line of the output .dat file is referring to one timestep:
!                 First column: current time
!                                body(1)%vert(1)%x, body(1)%vert(1)%y, body(1)%vert(1)%z
!                                body(1)%vert(2)%x, body(1)%vert(2)%y
!                                ... until vert(4)
!                                body(2)%vert(1)%x, body(2)%vert(1)%y
!                                body(2)%vert(2)%x, body(2)%vert(2)%y
!                                ... until body(n)
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!  Ruizhi Yang, 2018 Feb
!------------------------------------------------------------------------

SUBROUTINE write_Matlab_plot(i)

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                 :: i

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    INTEGER                                 :: j,k,Mfile_idx

    !--------------------------------------------------------------------
    !  Write verts_i to file
    !--------------------------------------------------------------------

    Mfile_idx = system%Mfile_idx

    ! current time in the solution
    WRITE(Mfile_idx,'(F9.5)', ADVANCE="NO") system%soln%t(i)

    DO j = 1, system%nbody
        DO k = 1, body_system(j)%nverts
            WRITE(Mfile_idx,'(F12.7)', ADVANCE="NO") body_system(j)%verts_i(k,1) ! x
            WRITE(Mfile_idx,'(F12.7)', ADVANCE="NO") body_system(j)%verts_i(k,2) ! y
            WRITE(Mfile_idx,'(F12.7)', ADVANCE="NO") body_system(j)%verts_i(k,3) ! z
        END DO
    END DO
    WRITE(Mfile_idx,'(/)')

END SUBROUTINE write_Matlab_plot