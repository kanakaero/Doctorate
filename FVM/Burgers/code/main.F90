! Main Program

program burgers
    use var_dec
    implicit none

    character(len=100) :: filename
    integer :: nstep

    ! Allocate Arrays
    allocate(x(N), xc(N), u(N+1), u_new(N+1), flux(N), res(N))
    allocate(u_out_01(N,3), u_out_05(N,3))
    allocate(r_out_01(N,3), r_out_05(N,3))

    ! Initialize Grid
    dx = Lx / real(N - 1, dp)
    do i = 1, N; x(i) = xstart + (i - 1) * dx; end do
    do i = 2, N; xc(i) = 0.5 * (x(i) + x(i-1)); end do

    do icase = 1, 5
        do fcase = 1, 3

            ! Assign Initial Condition
            call ic(icase)

            t = 0.0_dp
            nstep = 0

            do while (t < tend)
                nstep = nstep + 1

                ! Boundary Condition
                u(1) = u(2); u(N+1) = u(N)

                ! Compute Flux
                do iface = 1, N; i = iface; call cflux(u(i), u(i+1), fcase, flux(i)); end do

                ! Compute Residual
                do i = 2, N; res(i) = flux(i) - flux(i-1); end do

                ! Write t = 0.1
                if (nstep == 10) then
                    do i = 2, N
                        u_out_01(i,fcase) = u(i)
                        r_out_01(i,fcase) = res(i)
                    end do
                end if

                ! Write t = 0.5
                if (nstep == 50) then
                    do i = 2, N
                        u_out_05(i,fcase) = u(i)
                        r_out_05(i,fcase) = res(i)
                    end do
                end if

                ! Update Solution
                do i = 2, N; u_new(i) = u(i) - dt/dx * res(i); end do

                ! Copy updated solution
                do i = 2, N; u(i) = u_new(i); end do

                ! Advance time
                t = t + dt
            end do
        end do

        ! Write results for the associated IC case
        write(filename,'("results-",I0,".dat")') icase
        open(unit=20, file=filename, status="replace", action="write")

        ! t = 0.1
        write(20,'(A)') '# t = 0.1'
        write(20,'(A)') '# x  u_FVS  R_FVS  u_FDS  R_FDS  u_Godunov  R_Godunov'

        do i = 2, N
            write(20,'(ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3)') &
                xc(i), &
                u_out_01(i,1), r_out_01(i,1), &
                u_out_01(i,2), r_out_01(i,2), &
                u_out_01(i,3), r_out_01(i,3)
        end do

        write(20,*)
        write(20,*)

        ! t = 0.5
        write(20,'(A)') '# t = 0.5'
        write(20,'(A)') '# x  u_FVS  R_FVS  u_FDS  R_FDS  u_Godunov  R_Godunov'

        do i = 2, N
            write(20,'(ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3,1X,ES24.15E3)') &
                xc(i), &
                u_out_05(i,1), r_out_05(i,1), &
                u_out_05(i,2), r_out_05(i,2), &
                u_out_05(i,3), r_out_05(i,3)
        end do

        close(20)
    end do

    ! Deallocate Arrays
    deallocate(x, xc, u, u_new, u_out_01, u_out_05, r_out_01, r_out_05, flux, res)
end program burgers