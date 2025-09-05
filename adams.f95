module adams_m
    implicit none
    real, parameter :: PI = 3.141592653589793
    real, parameter :: TOLERANCE = 1.0e-12
contains

function quat_mult(q1, q2) result(ans)
    implicit none
    real, intent(in) :: q1(4), q2(4)
    real :: ans(4)
    real :: e01, ex1, ey1, ez1
    real :: e02, ex2, ey2, ez2
    e01 = q1(1) 
    ex1 = q1(2) 
    ey1 = q1(3)
    ez1 = q1(4)
    e02 = q2(1) 
    ex2 = q2(2) 
    ey2 = q2(3)
    ez2 = q2(4)
    ans(1) = e01*e02 - ex1*ex2 - ey1*ey2 -ez1*ez2
    ans(2) = e01*ex2 + ex1*e02 + ey1*ez2 -ez1*ey2
    ans(3) = e01*ey2 - ex1*ez2 + ey1*e02 + ez1*ex2
    ans(4) = e01*ez2 + ex1*ey2 - ey1*ex2 + ez1*e02
end function quat_mult

function quat_base_to_dependent(vec, quat) result(v2)
    ! Eq. (1.5.4 in the book which uses the quaternion product twice to get Eq. 1.5.7)
    implicit none
    real, intent(in) :: vec(3), quat(4)
    real :: v2(3)
    real :: t0, tx, ty, tz
    real :: e0, ex, ey, ez
    real :: v0, vx, vy, vz
    v0 = 0
    vx = vec(1)
    vy = vec(2)
    vz = vec(3)
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    t0 = -vx*ex - vy* ey - vz*ez 
    tx =  vx*e0 + vy* ez - vz* ey 
    ty =  -vx* ez + vy* e0 + vz* ex 
    tz = vx* ey - vy* ex + vz* e0
    v2(1) = e0*tx - ex*t0 - ey*tz + ez*ty 
    v2(2) = e0*ty + ex*tz - ey*t0 - ez*tx 
    v2(3) = e0*tz - ex*ty + ey*tx - ez*t0
end function quat_base_to_dependent

function quat_dependent_to_base(v2, quat) result(v1)
    implicit none
    real, intent(in) :: v2(3), quat(4)
    real :: v1(3)
    real :: t0, tx, ty, tz
    real :: e0, ex, ey, ez
    real :: v0, vx, vy, vz
    v0 = 0
    vx = v2(1)
    vy = v2(2)
    vz = v2(3)
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    t0 = vx*ex + vy* ey + vz*ez 
    tx =  vx*e0 - vy* ez + vz* ey 
    ty =  vx* ez + vy* e0 - vz* ex 
    tz = vx* ey - vy* ex + vz* e0
    v1(1) = e0 * tx + ex *t0 + ey* tz - ez *ty 
    v1(2) = e0 * ty - ex *tz + ey* t0 + ez *tx 
    v1(3) = e0 *tz + ex *ty - ey *tx + ez *t0
end function quat_dependent_to_base

subroutine quat_norm(quat)
    ! normalizes a quaternion such that it becomes a unit quaternion. 
    ! Done the EXACT way a unit vector is the vector divided by its magnitude
    implicit none
    real, intent(inout) :: quat(4)
    real :: magnitude 
    magnitude = sqrt(quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2)
    quat = quat/magnitude
end subroutine quat_norm

function euler_to_quat(euler) result(quat)
    ! Eq. 1.6.2 in the book. 
    real, intent(in) :: euler(3) 
    real :: quat(4)
    real :: phi, theta, psi 
    real :: cos_phi, cos_theta, cos_psi
    real :: sin_phi, sin_theta, sin_psi
    real :: cos_phi_cos_theta, sin_phi_sin_theta, sin_phi_cos_theta, cos_phi_sin_theta
    real :: hphi, htheta, hpsi
    phi = euler(1)
    theta = euler(2)
    psi = euler(3)
    hphi = 0.5*phi
    htheta = 0.5*theta
    hpsi = 0.5*psi
    cos_phi = cos(hphi)
    cos_theta = cos(htheta)
    cos_psi = cos(hpsi)
    sin_phi = sin(hphi)
    sin_theta = sin(htheta)
    sin_psi = sin(hpsi)
    cos_phi_cos_theta = cos_phi*cos_theta
    sin_phi_sin_theta = sin_phi*sin_theta
    sin_phi_cos_theta = sin_phi*cos_theta
    cos_phi_sin_theta = cos_phi*sin_theta
    quat(1) = cos_phi_cos_theta*cos_psi + sin_phi_sin_theta*sin_psi
    quat(2) = sin_phi_cos_theta*cos_psi - cos_phi_sin_theta*sin_psi
    quat(3) = cos_phi_sin_theta*cos_psi + sin_phi_cos_theta*sin_psi
    quat(4) = cos_phi_cos_theta*sin_psi - sin_phi_sin_theta*cos_psi
end function euler_to_quat

function quat_to_euler(quat) result(euler)
    ! Eq. 1.6.3 in the book. 
    real, intent(in) :: quat(4) 
    real :: euler(3)
    real :: cos_pi_4, arcsin_quat_cos_pi_4
    real :: e0, ex, ey, ez
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    cos_pi_4 = ex/cos(PI*0.25)
    arcsin_quat_cos_pi_4 = asin(cos_pi_4)
    if (e0*ey-ex*ez == 0.5) then
        euler(1) = 2*arcsin_quat_cos_pi_4
        euler(2) = 0.5*PI
        euler(3) = 0.0
    else if (e0*ey-ex*ez == 0.5) then
        euler(1) = 2*arcsin_quat_cos_pi_4
        euler(2) = -0.5*PI
        euler(3) = 0.0
    else 
        euler(1) = atan2(2*(e0*ex + ey*ez), (e0**2 + ez**2 -ex**2 -ey**2))
        euler(2) = asin(2*(e0*ey-ex*ez))
        euler(3) = atan2(2*(e0*ez + ex*ey), (e0**2 + ex**2 -ey**2 -ez**2))
    end if
end function quat_to_euler

end module adams_m
