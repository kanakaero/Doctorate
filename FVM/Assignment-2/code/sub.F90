! Subroutine Definitions

subroutine ic(ic_case)
    use var_dec
    implicit none

    integer, intent(in) :: ic_case
    select case (ic_case)
        ! Case 1: u = 1 for x <= 0 and u = 0 for x > 0
        case (1)
            do i = 2, N
                if (xc(i) <= 0.0_dp) then
                    u(i) = 1.0_dp
                else
                    u(i) = 0.0_dp
                end if
            end do

        ! Case 2: u = 0 for x <= 0 and u = 1 for x > 0
        case (2)
            do i = 2, N
                if (xc(i) <= 0.0_dp) then
                    u(i) = 0.0_dp
                else
                    u(i) = 1.0_dp
                end if
            end do

        ! Case 3: u = -1 for x <= 0 and u =  1 for x > 0
        case (3)
            do i = 2, N
                if (xc(i) <= 0.0_dp) then
                    u(i) = -1.0_dp
                else
                    u(i) = 1.0_dp
                end if
            end do

        ! Case 4: u =  1 for x <= 0 and u = -1 for x > 0
        case (4)
            do i = 2, N
                if (xc(i) <= 0.0_dp) then
                    u(i) = 1.0_dp
                else
                    u(i) = -1.0_dp
                end if
            end do

        ! Case 5:
        ! u = 1       for x <= -0.5 and 
        ! u = 0       for x >=  0.5
        ! linear      for -0.5 < x < 0.5
        case (5)
            do i = 2, N
                if (xc(i) <= -0.5_dp) then
                    u(i) = 1.0_dp
                else if (xc(i) >= 0.5_dp) then
                    u(i) = 0.0_dp
                else
                    u(i) = 0.5_dp - xc(i)
                end if
            end do
    end select
end subroutine ic

subroutine cflux(uL, uR, scheme, F)

    use var_dec
    implicit none

    real(dp), intent(in) :: uL, uR
    integer, intent(in) :: scheme
    real(dp), intent(out) :: F

    real(dp) :: fL, fR
    real(dp) :: s

    ! fL and fR correspond to f- and f+ respectively

    select case (scheme)

    ! Flux Splitting
    case (1)

        if (uL > 0.0_dp) then; fL = 0.5_dp * uL**2
            else; fL = 0.0_dp
        end if

        if (uR < 0.0_dp) then; fR = 0.5_dp * uR**2
            else; fR = 0.0_dp
        end if

        F = fL + fR

    ! Flux Difference Splitting
    case (2)

        fL = 0.5_dp * uL**2
        fR = 0.5_dp * uR**2

        s = 0.5_dp * (uL + uR)
        F = 0.5_dp * (fL + fR - (abs(s) * (uR - uL)))

    ! Godunov (Exact)
    case (3)

        if (uL <= uR) then

            ! Rarefaction
            if (uL >= 0.0_dp) then; F = 0.5_dp * uL**2
                else if (uR <= 0.0_dp) then; F = 0.5_dp * uR**2
                else; F = 0.0_dp ! Rarefaction crosses zero
            end if

        else

            ! Shock
            if (uL + uR >= 0.0_dp) then; F = 0.5_dp * uL**2
                else; F = 0.5_dp * uR**2
            end if

        end if

    end select

end subroutine cflux