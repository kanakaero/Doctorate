! Subroutines Module

module sub
    use var_dec
    implicit none

contains

    subroutine ic(which)
        integer, intent(in) :: which
        real(dp) :: rhoL,uL,pL,TL,rhoR,uR,pR,TR,M,aL,pRatio,rhoRatio
        integer :: j

        Qp=0.0_dp
        select case(which)
        case(1)
            rhoL=1.0_dp; uL=0.0_dp; TL=500.0_dp
            rhoR=2.5_dp; uR=0.0_dp; TR=200.0_dp
            pL=rhoL*Rgas*TL; pR=rhoR*Rgas*TR
        case(2)
            pL=100.0_dp; TL=401.8_dp; M=3.0_dp
            rhoL=pL/(Rgas*TL); aL=sqrt(gamma*Rgas*TL); uL=M*aL
            pRatio=1.0_dp+2.0_dp*gamma/(gamma+1.0_dp)*(M*M-1.0_dp)
            rhoRatio=(gamma+1.0_dp)*M*M/((gamma-1.0_dp)*M*M+2.0_dp)
            pR=pL*pRatio; rhoR=rhoL*rhoRatio
            TR=TL*pRatio/rhoRatio; uR=uL/rhoRatio
        case(3)
            rhoL=2.5_dp; uL=0.0_dp; TL=500.0_dp
            rhoR=1.0_dp; uR=0.0_dp; TR=250.0_dp
            pL=rhoL*Rgas*TL; pR=rhoR*Rgas*TR
        case default
            error stop "ic: invalid case"
        end select
        do j=1,N
            if (xc(j)<=0.0_dp) then
                Qp(:,j)=[rhoL, uL, pL, TL]
            else
                Qp(:,j)=[rhoR, uR, pR, TR]
            end if
        end do
    end subroutine ic

    subroutine prim_to_cons(Qp,Qc)
        real(dp),intent(in)::Qp(:,:)
        real(dp),intent(out)::Qc(:,:)
        integer::j

        Qc=0.0_dp
        do j=1,N
            Qc(:,j+1)=primitive_to_conservative(Qp(:,j))
        end do
    end subroutine prim_to_cons

    function primitive_to_conservative(Qp) result(Qc)
        real(dp),intent(in)::Qp(4)
        real(dp)::Qc(3),rho,u,p,T

        rho=Qp(1);u=Qp(2);p=Qp(3);T=Qp(4)
        Qc=[rho,rho*u,p/(gamma-1.0_dp)+0.5_dp*rho*u*u]
    end function primitive_to_conservative

    subroutine cons_to_prim(Qc,Qp)
        real(dp),intent(in)::Qc(:,:)
        real(dp),intent(out)::Qp(:,:)
        integer::j
        real(dp)::rho,u,p

        do j=1,N
            rho=Qc(1,j+1)
            if (rho<=1e-12_dp) error stop "cons_to_prim: non-positive density"
            u=Qc(2,j+1)/rho
            p=(gamma-1.0_dp)*(Qc(3,j+1)-0.5_dp*Qc(2,j+1)*u)
            if (p<=1.0e-10_dp) then
                write(*,'(A,I0,A,ES14.6,A,ES14.6)') "bad cell ",j," rho=",rho," p=",p
                case_failed=.true.
                return
            end if
            Qp(:,j)=[rho, u, p, p/(rho*Rgas)]
        end do
    end subroutine

    subroutine boundary_conditions(Qc)
        real(dp),intent(inout)::Qc(:,:)

        ! Transmissive (zero-gradient) boundary conditions for ghost cells
        Qc(:,1)=Qc(:,2)
        Qc(:,N+2)=Qc(:,N+1)
    end subroutine

    subroutine compute_dt(Qc,dtloc)
        real(dp),intent(in)::Qc(:,:)
        real(dp),intent(out)::dtloc
        integer::j
        real(dp)::rho,u,p,a,smax

        smax=0.0_dp
        do j=2,N+1
            rho=Qc(1,j);u=Qc(2,j)/rho
            p=(gamma-1.0_dp)*(Qc(3,j)-0.5_dp*Qc(2,j)*u)
            if (rho <= 1.0e-12_dp .or. p <= 1.0e-10_dp) then
                case_failed = .true.
                dtloc = 0.0_dp
                return
            end if
            a=sqrt(gamma*p/rho)
            smax=max(smax,abs(u)+a)
        end do
        dtloc=CFL*dx/smax
    end subroutine

    subroutine compute_rhs(Qc,rhs)
        real(dp),intent(in)::Qc(:,:)
        real(dp),intent(out)::rhs(:,:)
        real(dp)::Ub(3,N+2),F(3,N+1)
        integer::j

        Ub=Qc
        call boundary_conditions(Ub)
        do j=1,N+1
            F(:,j)=numerical_flux(Ub(:,j),Ub(:,j+1),fcase)
            if(case_failed) then
                rhs=0.0_dp
                return
            end if
        end do
        do j=1,N
            rhs(:,j)=-(F(:,j+1)-F(:,j))/dx
        end do
    end subroutine

    function numerical_flux(QcLin,QcRin,scheme) result(Fout)
        real(dp),intent(in)::QcLin(3),QcRin(3)
        integer,intent(in)::scheme
        real(dp)::Fout(3),Qpleft(4),Qpright(4)
        real(dp)::HL,HR,rhol,uloc,pl,rhor,uright,pr,al,ar
        
        Qpleft=to_primitive(QcLin);Qpright=to_primitive(QcRin)
        if (case_failed) then
            Fout = 0.0_dp
            return
        end if
        rhol=Qpleft(1);uloc=Qpleft(2);pl=Qpleft(3)
        rhor=Qpright(1);uright=Qpright(2);pr=Qpright(3)
        al=sqrt(gamma*pl/rhol);ar=sqrt(gamma*pr/rhor)
        HL=(QcLin(3)+pl)/rhol;HR=(QcRin(3)+pr)/rhor

        select case(scheme)
            case(1)
                Fout=sw_flux(QcLin,QcRin)
            case(2)
                Fout=vanleer_flux(Qpleft,Qpright,al,ar,.false.)
            case(3)
                Fout=vanleer_flux(Qpleft,Qpright,al,ar,.true.)
            case(4)
                Fout=ausm_flux(Qpleft,Qpright,al,ar)
            case(5)
                Fout=roe_flux(QcLin,QcRin,Qpleft,Qpright,HL,HR)
        end select
    end function

    function to_primitive(Qc) result(Qp)
        real(dp),intent(in)::Qc(3)
        real(dp)::Qp(4),rho,u,p

        rho=Qc(1)

        if(rho<=1.0e-12_dp) then
            case_failed=.true.
            Qp=0.0_dp
            return
        end if

        u=Qc(2)/rho
        p=(gamma-1.0_dp)*(Qc(3)-0.5_dp*Qc(2)*u)

        if(p<=1.0e-10_dp) then
            case_failed=.true.
            Qp=0.0_dp
            return
        end if

        Qp=[rho, u, p, p/(rho*Rgas)]
    end function

    function sw_flux(QcL,QcR) result(F)
        real(dp),intent(in)::QcL(3),QcR(3)
        real(dp)::F(3),QpL(4),QpR(4),FLp(3),FLm(3),FRp(3),FRm(3)

        QpL=to_primitive(QcL);QpR=to_primitive(QcR)
        FLp=sw_split(QcL,QpL,.true.);FLm=sw_split(QcL,QpL,.false.)
        FRp=sw_split(QcR,QpR,.true.);FRm=sw_split(QcR,QpR,.false.)
        F=FLp+FRm
    end function

    function sw_split(Qc,Qp,positive) result(Fs)
        real(dp),intent(in)::Qc(3),Qp(3)
        logical,intent(in)::positive
        real(dp)::Fs(3),rho,u,p,a,H,lam(3),Rmat(3,3),Lmat(3,3),coef(3),lp(3)
        integer::k
        
        rho=Qp(1);u=Qp(2);p=Qp(3);a=sqrt(gamma*p/rho)
        H=(Qc(3)+p)/rho

        Rmat(:,1)=[1.0_dp,u,0.5_dp*u*u]
        Rmat(:,2)=[1.0_dp,u+a,H+u*a]
        Rmat(:,3)=[1.0_dp,u-a,H-u*a]

        call invert3(Rmat,Lmat)
        coef=matmul(Lmat,Qc)
        lam=[u,u+a,u-a]

        do k=1,3
            if(positive) then
                lp(k)=0.5_dp*(lam(k)+abs(lam(k)))
            else
                lp(k)=0.5_dp*(lam(k)-abs(lam(k)))
            end if
        end do

        Fs=matmul(Rmat,lp*coef)
    end function

    subroutine invert3(A,B)
        real(dp),intent(in)::A(3,3)
        real(dp),intent(out)::B(3,3)
        real(dp)::det

        det=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
           -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
           +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

        if(abs(det)<1.0e-14_dp) error stop "invert3: singular matrix"

        B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
        B(1,2)=(A(1,3)*A(3,2)-A(1,2)*A(3,3))/det
        B(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
        B(2,1)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
        B(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
        B(2,3)=(A(1,3)*A(2,1)-A(1,1)*A(2,3))/det
        B(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
        B(3,2)=(A(1,2)*A(3,1)-A(1,1)*A(3,2))/det
        B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
    end subroutine

    function vanleer_flux(Qpleft,Qpright,aL,aR,common) result(Fout)
        real(dp),intent(in)::Qpleft(4),Qpright(4),aL,aR
        logical,intent(in)::common
        real(dp)::Fout(3),alocal,arocal,ML,MR

        if(common) then
            alocal=max(aL,aR);arocal=alocal
        else
            alocal=aL;arocal=aR
        end if

        ML=Qpleft(2)/alocal;MR=Qpright(2)/arocal
        Fout=vl_split(Qpleft,alocal,ML,.true.)+vl_split(Qpright,arocal,MR,.false.)
    end function

    function vl_split(Qp,a,M,plus) result(Fp)
        real(dp),intent(in)::Qp(4),a,M
        logical,intent(in)::plus
        real(dp)::Fp(3),rho,u,p,fac

        rho=Qp(1);u=Qp(2);p=Qp(3)

        if(plus) then
            if(M>=1.0_dp) then
                Fp=[rho*u,rho*u*u+p,u*(rho*(0.5_dp*u*u)+gamma*p/(gamma-1.0_dp))]
            else if(M<=-1.0_dp) then
                Fp=0.0_dp
            else
                fac=0.25_dp*rho*a*(M+1.0_dp)**2
                Fp=fac*[1.0_dp, 2.0_dp*a/gamma*(M+(gamma-1.0_dp)/2.0_dp), &
                    2.0_dp*a*a/(gamma*gamma-1.0_dp)*(0.5_dp*(gamma-1.0_dp)*M*M+M+1.0_dp/(gamma+1.0_dp))]
            end if
        else
            if(M<=-1.0_dp) then
                Fp=[rho*u,rho*u*u+p,u*(rho*(0.5_dp*u*u)+gamma*p/(gamma-1.0_dp))]
            else if(M>=1.0_dp) then
                Fp=0.0_dp
            else
                fac=-0.25_dp*rho*a*(M-1.0_dp)**2
                Fp=fac*[1.0_dp, 2.0_dp*a/gamma*(M-(gamma-1.0_dp)/2.0_dp), &
                    2.0_dp*a*a/(gamma*gamma-1.0_dp)*(0.5_dp*(gamma-1.0_dp)*M*M-M+1.0_dp/(gamma+1.0_dp))]
            end if
        end if
    end function

    function ausm_flux(Qpleft,Qpright,aL,aR) result(Fout)
        real(dp),intent(in)::Qpleft(4),Qpright(4),aL,aR
        real(dp)::Fout(3)
        real(dp)::rhoL,uL,pL,HL, rhoR,uR,pR,HR
        real(dp)::a_half,ML,MR,Mhalf,phalf
        real(dp)::massflux,uface,Hface

        rhoL=Qpleft(1);uL=Qpleft(2);pL=Qpleft(3)
        rhoR=Qpright(1);uR=Qpright(2);pR=Qpright(3)

        HL=0.5_dp*uL*uL + gamma/(gamma-1.0_dp)*pL/rhoL
        HR=0.5_dp*uR*uR + gamma/(gamma-1.0_dp)*pR/rhoR

        a_half=max(aL,aR)
        ML=uL/a_half; MR=uR/a_half

        Mhalf=mach_plus(ML)+mach_minus(MR)
        phalf=pressure_plus(ML)*pL + pressure_minus(MR)*pR

        if(Mhalf>=0.0_dp) then
            massflux=a_half*Mhalf*rhoL
            uface=uL; Hface=HL
        else
            massflux=a_half*Mhalf*rhoR
            uface=uR; Hface=HR
        end if

        Fout(1)=massflux; Fout(2)=massflux*uface+phalf; Fout(3)=massflux*Hface
    end function ausm_flux

    function mach_plus(M) result(v)
        real(dp),intent(in)::M
        real(dp)::v
        
        if(abs(M)>=1.0_dp) then
            v=0.5_dp*(M+abs(M))
        else
            v=0.25_dp*(M+1.0_dp)**2
        end if
    end function

    function mach_minus(M) result(v)
        real(dp),intent(in)::M
        real(dp)::v

        if(abs(M)>=1.0_dp) then
            v=0.5_dp*(M-abs(M))
        else
            v=-0.25_dp*(M-1.0_dp)**2
        end if
    end function

    function pressure_plus(M) result(v)
        real(dp),intent(in)::M
        real(dp)::v

        if(abs(M)>=1.0_dp) then
            v=0.5_dp*(1.0_dp+sign(1.0_dp,M))
        else
            v=0.25_dp*(M+1.0_dp)**2*(2.0_dp-M)
        end if
    end function

    function pressure_minus(M) result(v)
        real(dp),intent(in)::M
        real(dp)::v

        if(abs(M)>=1.0_dp) then
            v=0.5_dp*(1.0_dp-sign(1.0_dp,M))
        else
            v=0.25_dp*(M-1.0_dp)**2*(2.0_dp+M)
        end if
    end function

    function roe_flux(QcL,QcR,QpL,QpR,HL,HR) result(F)
        real(dp),intent(in)::QcL(3),QcR(3),QpL(4),QpR(4),HL,HR
        real(dp)::F(3),ubar,Hbar,rbar,abar,Rmat(3,3),alpha(3),lam(3)
        real(dp)::s1,s2,s3
        real(dp)::drho,duvel,dpres

        rbar=sqrt(QpL(1)*QpR(1))
        ubar=(sqrt(QpL(1))*QpL(2)+sqrt(QpR(1))*QpR(2))/(sqrt(QpL(1))+sqrt(QpR(1)))
        Hbar=(sqrt(QpL(1))*HL+sqrt(QpR(1))*HR)/(sqrt(QpL(1))+sqrt(QpR(1)))
        abar=sqrt(max(0.0_dp,(gamma-1.0_dp)*(Hbar-0.5_dp*ubar*ubar)))
    
        Rmat(:,1)=[1.0_dp,ubar,0.5_dp*ubar*ubar]
        Rmat(:,2)=[1.0_dp,ubar+abar,Hbar+ubar*abar]
        Rmat(:,3)=[1.0_dp,ubar-abar,Hbar-ubar*abar]

        drho  = QpR(1)-QpL(1)
        duvel = QpR(2)-QpL(2)
        dpres = QpR(3)-QpL(3)

        alpha(1)=drho-dpres/(abar*abar)
        alpha(2)=(dpres+rbar*abar*duvel)/(2.0_dp*abar*abar)
        alpha(3)=(dpres-rbar*abar*duvel)/(2.0_dp*abar*abar)

        lam=[ubar,ubar+abar,ubar-abar]
        s1=abs(lam(1));s2=abs(lam(2));s3=abs(lam(3))
        F=0.5_dp*(physical_flux(QcL)+physical_flux(QcR)) &
          -0.5_dp*(Rmat(:,1)*s1*alpha(1)+Rmat(:,2)*s2*alpha(2)+Rmat(:,3)*s3*alpha(3))
    end function

    function physical_flux(Qc) result(f)
        real(dp),intent(in)::Qc(3)
        real(dp)::f(3),rho,u,p

        rho=Qc(1);u=Qc(2)/rho
        p=(gamma-1.0_dp)*(Qc(3)-0.5_dp*Qc(2)*u)
        f=[Qc(2),Qc(2)*u+p,u*(Qc(3)+p)]
    end function

    subroutine time_integrator(Qc,Qcout,dtloc,which)
        real(dp),intent(in)::Qc(:,:),dtloc
        real(dp),intent(out)::Qcout(:,:)
        integer,intent(in)::which
        real(dp)::k1(3,N),k2(3,N),k3(3,N),k4(3,N),Us(3,N+2)
        
        select case(which)
            case(1)
                call compute_rhs(Qc,k1)
                if(case_failed) then
                    Qcout=Qc
                    return
                end if
                Qcout=Qc
                Qcout(:,2:N+1)=Qc(:,2:N+1)+dtloc*k1
                call boundary_conditions(Qcout)
            case(2)
                call compute_rhs(Qc,k1)
                if(case_failed) then
                    Qcout=Qc
                    return
                end if
                Us=Qc;Us(:,2:N+1)=Qc(:,2:N+1)+0.5_dp*dtloc*k1
                call boundary_conditions(Us)
                call compute_rhs(Us,k2)
                if(case_failed) then
                    Qcout=Qc
                    return
                end if
                Us=Qc;Us(:,2:N+1)=Qc(:,2:N+1)+0.5_dp*dtloc*k2
                call boundary_conditions(Us)
                call compute_rhs(Us,k3)
                if(case_failed) then
                    Qcout=Qc
                    return
                end if
                Us=Qc;Us(:,2:N+1)=Qc(:,2:N+1)+dtloc*k3
                call boundary_conditions(Us)
                call compute_rhs(Us,k4)
                if(case_failed) then
                    Qcout=Qc
                    return
                end if
                Qcout=Qc
                Qcout(:,2:N+1)=Qc(:,2:N+1)+dtloc/6.0_dp*(k1+2.0_dp*k2+2.0_dp*k3+k4)
                call boundary_conditions(Qcout)
        end select
    end subroutine

    subroutine write_solution(unitno)
        integer,intent(in)::unitno
        integer::j
        real(dp)::rho,u,p,T

        do j=1,N
            rho=Qp(1,j);u=Qp(2,j);p=Qp(3,j);T=Qp(4,j)
            write(unitno,'(5(ES22.13,1X))') xc(j),rho,u,p,T
        end do

        write(unitno,*)
    end subroutine

    subroutine write_exact_solution(unitno)
        integer,intent(in)::unitno
        integer::j
        real(dp)::rho,u,p,T

        do j=1,N
            rho=Qp_exact(1,j);u=Qp_exact(2,j);p=Qp_exact(3,j);T=Qp_exact(4,j)
            write(unitno,'(5(ES22.13,1X))') xc(j),rho,u,p,T
        end do

        write(unitno,*)
    end subroutine

    subroutine write_residual(unit,t,res)
        integer,intent(in) :: unit
        real(dp),intent(in) :: t,res(:,:)
        real(dp) :: max_res(3)
        integer :: k

        do k = 1, 3
            max_res(k) = maxval(abs(res(k,:)))
        end do

        write(unit,'(4ES20.10)') t, max_res(1), max_res(2), max_res(3)
    end subroutine

    subroutine get_riemann_states(which, QpL, QpR)
        integer, intent(in) :: which
        real(dp), intent(out) :: QpL(4), QpR(4)

        call ic(which)

        QpL = Qp(:,1)
        QpR = Qp(:,N)
    end subroutine get_riemann_states

end module sub