PROGRAM dyn3d

    USE module_basic_matrix_operations
    USE module_constants
    USE module_data_type
    USE module_trans_matrix

IMPLICIT NONE

    REAL,ALLOCATABLE           :: X(:,:)!,Xinv(:,:),rot(:,:),tr(:,:)
    REAL,ALLOCATABLE           :: r(:),theta(:)
    INTEGER                    :: max_iterations
    REAL                       :: tol


    ALLOCATE(X(6,6))
    !ALLOCATE(Xinv(6,6))
    !ALLOCATE(rot(6,6))
    !ALLOCATE(tr(6,6))
    ALLOCATE(r(3))
    ALLOCATE(theta(3))

    r = (/1.0, 1.0, 1.0/)
    theta = (/ 0.0, 0.0, 0.0/)

    CALL trans_matrix(r,theta,X)
    !CALL trans_matrix(r,theta,X,Xinv,rot,tr)
    CALL write_matrix(X)


END PROGRAM dyn3d