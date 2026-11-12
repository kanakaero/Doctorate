! Variable Declaration Module

module var_dec
    implicit none

    ! Define Double Precision
    integer, parameter :: dp = kind(1.d0)

    ! Number of Variables
    integer, parameter :: nvar  = 3
    integer, parameter :: nprim = 4

    ! Initial Condition, Flux Scheme, Time Integrator
    integer :: icase, fcase, tcase

    ! Define Domain and Grid Parameters
    integer, parameter :: Lx = 2
    integer, parameter :: N = 40

    real(dp) :: xstart = -1.0_dp
    real(dp) :: xend   =  1.0_dp

    ! Define Gas Parameters
    real(dp), parameter :: gamma = 1.4_dp
    real(dp), parameter :: Rgas  = 287.0_dp

    ! Define Time Parameters
    real(dp) :: dx, dt
    real(dp) :: t, tend = 0.002_dp

    ! CFL Number
    real(dp), parameter :: CFL = 0.2_dp

    ! Define Grid Arrays
    real(dp), allocatable :: x(:), xc(:)

    ! Conservative Variables
    ! Qc(1,:) = Density
    ! Qc(2,:) = Momentum
    ! Qc(3,:) = Total Energy
    real(dp), allocatable :: Qc(:,:), Qc_new(:,:)

    ! Primitive Variables
    ! Qp(1,:) = Density
    ! Qp(2,:) = Velocity
    ! Qp(3,:) = Pressure
    ! Qp(4,:) = Temperature
    real(dp), allocatable :: Qp(:,:)

    ! Flux Arrays
    real(dp), allocatable :: flux(:,:)

    ! Residual Arrays
    real(dp), allocatable :: res(:,:)

    ! Output Arrays
    real(dp), allocatable :: Qc_out_01(:,:), Qc_out_05(:,:)
    real(dp), allocatable :: Qp_out_01(:,:), Qp_out_05(:,:)

    ! Exact Solution Arrays
    real(dp), allocatable :: Qp_exact(:,:)

    ! Helper Variables
    integer(kind=8) :: i, iface
    logical :: case_failed

end module var_dec