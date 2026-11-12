! Variable Declaration Module

module var_dec
    implicit none
    
    ! Define Double Precision
    integer, parameter :: dp = kind(1.d0)

    ! Initial Condition Case
    integer :: icase, fcase

    ! Define Domain and Grid Parameters
    integer, parameter :: Lx = 2
    integer, parameter :: N = 21
    real(dp) :: xstart = -1.0_dp, xend = 1.0_dp

    ! Define Time Parameters
    real(dp) :: dx, dt = 0.01
    real(dp) :: t, tend = 0.5_dp

    ! Define Arrays
    real(dp), allocatable :: x(:), xc(:)

    real(dp), allocatable :: u(:), u_new(:)
    real(dp), allocatable :: u_out_01(:,:), u_out_05(:,:)

    real(dp), allocatable :: flux(:)
    real(dp), allocatable :: res(:)
    real(dp), allocatable :: r_out_01(:,:), r_out_05(:,:)
    
    ! Helper Variables
    integer(kind=8) :: i, iface
end module var_dec