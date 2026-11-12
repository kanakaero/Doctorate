! Main program for 1D Euler equations

program euler

    use var_dec
    use sub
    use exact
    implicit none

    integer :: nstep
    integer :: output_count
    character(len=100) :: filename
    real(dp) :: QpL(4), QpR(4)
    real(dp) :: x0, next_output_time
    integer :: exact_unit
    character(len=100) :: exact_filename
    character(len=100) :: residual_filename
    integer :: residual_unit
    real(dp), parameter :: output_dt = 0.0002_dp
    integer :: k

    allocate(x(N+1), xc(N))
    allocate(Qc(nvar,N+2), Qc_new(nvar,N+2))
    allocate(Qp(nprim,N))
    allocate(flux(nvar,N+1))
    allocate(res(nvar,N))
    allocate(Qp_exact(nprim,N))

    dx = (xend - xstart) / real(N,dp)
    do i = 1, N+1; x(i) = xstart + (i-1)*dx; end do
    do i = 1, N; xc(i) = 0.5_dp*(x(i) + x(i+1)); end do

    do icase = 1, 3

        call ic(icase)
        call get_riemann_states(icase, QpL, QpR)

        x0 = 0.0_dp

        write(exact_filename,'("exact-IC",I0,".dat")') icase
        open(newunit=exact_unit, file=trim(exact_filename), status="replace", action="write")
        write(exact_unit,'(A)') '# x rho u p T'

        do k = 0, 10
            call exact_riemann(xc, real(k,dp)*output_dt, x0, QpL, QpR, gamma, Qp_exact)
            call write_exact_solution(exact_unit)
        end do

        close(exact_unit)

        do fcase = 1, 5
            do tcase = 1, 2

                case_failed = .false.

                call ic(icase)
                call prim_to_cons(Qp,Qc)

                t = 0.0_dp
                nstep = 0
                output_count = 0

                write(filename,'("results-IC",I0,"-F",I0,"-T",I0,".dat")') icase, fcase, tcase
                open(unit=20, file=trim(filename), status="replace", action="write")
                write(20,'(A)') '# x rho u p T'

                write(residual_filename,'("residual-IC",I0,"-F",I0,"-T",I0,".dat")') icase, fcase, tcase
                open(newunit=residual_unit, file=trim(residual_filename), status="replace", action="write")
                write(residual_unit,'(A)') '# t max_rho max_momentum max_energy'

                call cons_to_prim(Qc,Qp)
                call write_solution(20)

                next_output_time = 0.0002_dp

                ! Time Stepping Loop
                do while (t < tend)

                    nstep = nstep + 1

                    call compute_dt(Qc,dt)

                    if (case_failed) then
                        write(*,'(A,I0,A,I0,A,I0)') &
                            'Skipping IC=',icase,' F=',fcase,' T=',tcase
                        exit
                    end if

                    dt = min(dt, next_output_time - t)

                    call compute_rhs(Qc,res)
                    call write_residual(residual_unit, t, res)
                    
                    call time_integrator(Qc,Qc_new,dt,tcase)

                    if (case_failed) then
                        write(*,'(A,I0,A,I0,A,I0)') &
                            'Skipping IC=',icase,' F=',fcase,' T=',tcase
                        exit
                    end if

                    Qc = Qc_new
                    t = t + dt

                    call cons_to_prim(Qc,Qp)

                    if (t >= next_output_time - 1.0e-10_dp) then
                        call write_solution(20)
                        next_output_time = next_output_time + 0.0002_dp
                    end if

                end do

                close(20)
                close(residual_unit)
            end do
        end do
    end do

    deallocate(x,xc,Qc,Qc_new,Qp,Qp_exact,flux,res)

end program euler