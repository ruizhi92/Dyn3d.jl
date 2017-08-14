!------------------------------------------------------------------------
!  Module	    :            module_data_type
!------------------------------------------------------------------------
!  Purpose      :
!
!  Details      ：
!
!  Input        :
!
!  Input/output :
!
!  Output       :
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------
!  whirl vortex-based immersed boundary library
!  SOFIA Laboratory
!  University of California, Los Angeles
!  Los Angeles, California 90095  USA
!------------------------------------------------------------------------

MODULE module_data_type

    USE module_constants

IMPLICIT NONE

    TYPE dof
    !----------------------------- TYPE dof -----------------------------
    ! the total number dof of the system is number of joint*6, because each joint
    ! has 6 degree of freedoms
    ! Provide a type to describe property of each degree of freedom
    !                 {dof_id,dof_type,stiffness,damping}
    !                or
    !                 {dof_id,dof_type,motion_type,motion_params}
    !                depending on whether dof_type is 'passive' or 'active'.
    !                     dof_id:    ID of the degree of freedom
    !                                whose motion is to be specified. This
    !                                ID corresponds to the overall degree.
    !                The other entries for a 'passive' degree of freedom are
    !                     stiff: Stiffness of a spring associated with the
    !                                degree of freedom.
    !                     damp:   Damping coefficient of a spring associated
    !                                with the degree of freedom.
    !                 or for an 'active' degree of freedom
    !                     motion_type: Can be 'ramp','oscillatory',
    !                                        'velocity','hold'
    !                     motion_params: An array of parameters specifying the
    !                                    behavior of the specified motion.
    !                                    1. For 'ramp', this comes as a list
    !                                     [a,amp1,t_hold1,rate1,amp2,t_hold2,...]
    !                                    where a is the smoothing value,
    !                                    each pair ampj,t_holdj specifies
    !                                    to hold for at a certain value ampj
    !                                    for duration t_holdj, and each of
    !                                    these hold pairs is joined by a rate
    !                                    that specifies the duration of the ramp.
    !                                    2. For 'oscillatory',the list is
    !                                     [amp,freq,phase]
    !                                    3. For 'velocity',the list is
    !                                      [velocity]
    !                                    which specifies the constant velocity
    !                                    at which to move that degree.
    !                                    4. For 'hold', the list is a single value
    !                                     [amp]
    !                                    at which to hold the degree of freedom.
    !                  If no sub-array is given for a degree of freedom, then
    !                  the default is 'passive' with no stiffness/damping,
    !                  except for cases in which body1=0, in which case
    !                  default is 'active' with a hold at zero.
        INTEGER                         :: dof_id
        CHARACTER(LEN = max_char)       :: dof_type
        REAL,ALLOCATABLE                :: stiff(:),damp(:)
        CHARACTER(LEN = max_char)       :: motion_type
        REAL,ALLOCATABLE                :: motion_params(:)
    END TYPE

END MODULE module_data_type