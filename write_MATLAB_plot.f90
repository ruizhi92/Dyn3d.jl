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
!  Ruizhi Yang, 2018 Dec
!------------------------------------------------------------------------

SUBROUTINE write_Matlab_plot

    !--------------------------------------------------------------------
    !  MODULE
    !--------------------------------------------------------------------
    USE module_constants
    USE module_data_type
    USE module_trans_matrix

IMPLICIT NONE

    !--------------------------------------------------------------------
    !  Local variables
    !--------------------------------------------------------------------
    CHARACTER(LEN = max_char)               :: outfile,outfolder,fullname
    INTEGER                                 :: i,j,k
    REAL(dp),DIMENSION(6,6)                 :: Xi_to_b
    REAL(dp),DIMENSION(6,1)                 :: q_temp

    !--------------------------------------------------------------------
    !  Set write to file
    !--------------------------------------------------------------------
    ! note here that we are building in the cmake build folder
    outfolder = '../Matlab_plot'

    !--------------------------------------------------------------------
    !  Write verts_i
    !--------------------------------------------------------------------
    outfile = 'verts.dat'
    fullname = TRIM(outfolder)//'/'//TRIM(outfile)
    OPEN(2018,file = fullname)

    DO i = 1, system%params%nstep + 1

        ! current time in the solution
        WRITE(2018,'(F9.5)', ADVANCE="NO") system%soln%t(i)

        ! apply the solution, only q is needed
        DO j = 1, system%nbody
            body_system(j)%q(:,1) = system%soln%y(i,6*(j-1)+1:6*j)
        END DO

        ! calculate verts_i based on solution at the current time
        DO j = 1, system%nbody
            ! Xb_to_i
            CALL trans_matrix(body_system(j)%q(4:6,1),body_system(j)%q(1:3,1), &
                              Xi_to_b,body_system(j)%Xb_to_i)

            ! x_0
            body_system(j)%x_0 = body_system(j)%q(4:6,1)

            ! verts_i
            DO k = 1, body_system(j)%nverts
                q_temp(1:3,1) = 0.0_dp
                q_temp(4:6,1) = body_system(j)%verts(k,:)
                q_temp = MATMUL(body_system(j)%Xb_to_i,q_temp)
                body_system(j)%verts_i(k,:) = q_temp(4:6,1) + body_system(j)%x_0
            END DO
        END DO

        ! write verts_i to file
        DO j = 1, system%nbody
            DO k = 1, body_system(j)%nverts
                WRITE(2018,'(F9.5)', ADVANCE="NO") body_system(j)%verts_i(k,1) ! x
                WRITE(2018,'(F9.5)', ADVANCE="NO") body_system(j)%verts_i(k,2) ! y
                WRITE(2018,'(F9.5)', ADVANCE="NO") body_system(j)%verts_i(k,3) ! z
            END DO
        END DO

        WRITE(2018,'(/)')

    END DO

    CLOSE(2018)

    WRITE(*,*) 'Matlab plot info output done.'

END SUBROUTINE write_Matlab_plot