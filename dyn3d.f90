PROGRAM dyn3d

    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_trans_matrix

IMPLICIT NONE

    REAL,DIMENSION(6,6)           :: X,Xinv,rot,tr
    REAL,DIMENSION(3)           :: r,theta
    INTEGER                    :: max_iterations
    REAL                       :: tol


    r = (/1.0, 1.0, 1.0/)
    theta = (/ 0.0, 0.0, 0.0/)

    !CALL trans_matrix(r,theta,X)
    CALL trans_matrix(r,theta,X,Xinv,rot,tr)
    CALL write_matrix(X,6)


END PROGRAM dyn3d