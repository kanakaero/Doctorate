! Exact solution Module

module exact
    use var_dec, only: dp, Rgas
    implicit none
    private
    public :: exact_riemann

contains

    pure real(dp) function pf(p,rho,pk,a,g) result(f)
        real(dp),intent(in)::p,rho,pk,a,g
        real(dp)::coef,B

        if(p>pk)then
            coef=2.0_dp/((g+1.0_dp)*rho); B=(g-1.0_dp)*pk/(g+1.0_dp)
            f=(p-pk)*sqrt(coef/(p+B))
        else
            f=2.0_dp*a/(g-1.0_dp)*((p/pk)**((g-1.0_dp)/(2.0_dp*g))-1.0_dp)
        endif
    end function

    subroutine starp(rL,uL,pL,aL,rR,uR,pR,aR,g,ps)
        real(dp),intent(in)::rL,uL,pL,aL,rR,uR,pR,aR,g
        real(dp),intent(out)::ps
        real(dp)::lo,hi,mid,f
        integer::i

        lo=0.0_dp; hi=max(pL,pR)

        do while(pf(hi,rL,pL,aL,g)+pf(hi,rR,pR,aR,g)+uR-uL<0.0_dp)
            hi=2.0_dp*hi
        enddo

        do i=1,100
            mid=0.5_dp*(lo+hi)
            f=pf(mid,rL,pL,aL,g)+pf(mid,rR,pR,aR,g)+uR-uL
            if(abs(f)<1.0e-12_dp)exit
            if(f>0.0_dp)then
                hi=mid
            else
                lo=mid
            endif
        enddo

        ps=mid
    end subroutine

    subroutine exact_riemann(x,t,x0,QL,QR,g,Q)
        real(dp),intent(in)::x(:),t,x0,QL(4),QR(4),g
        real(dp),intent(out)::Q(:,:)
        real(dp)::rL,uL,pL,aL,rR,uR,pR,aR,ps,us,rsL,rsR
        real(dp)::asL,asR,SL,SR,SHL,STL,SHR,STR,xi,r,u,p,a,rg
        integer::i

        rL=QL(1); uL=QL(2); pL=QL(3); aL=sqrt(g*pL/rL)
        rR=QR(1); uR=QR(2); pR=QR(3); aR=sqrt(g*pR/rR)
        rg=Rgas
        if(t<=0.0_dp)then
            do i=1,size(x)
                Q(:,i)=merge(QL,QR,x(i)<x0)
            enddo
            return
        endif
        call starp(rL,uL,pL,aL,rR,uR,pR,aR,g,ps)
        us=0.5_dp*(uL+uR+pf(ps,rR,pR,aR,g)-pf(ps,rL,pL,aL,g))
        if(ps>pL)then
            rsL=rL*(ps/pL+(g-1.0_dp)/(g+1.0_dp))/ &
                ((g-1.0_dp)/(g+1.0_dp)*ps/pL+1.0_dp)
            SL=uL-aL*sqrt((g+1.0_dp)*ps/(2.0_dp*g*pL)+(g-1.0_dp)/(2.0_dp*g))
        else
            rsL=rL*(ps/pL)**(1.0_dp/g); asL=aL*(ps/pL)**((g-1.0_dp)/(2.0_dp*g))
            SHL=uL-aL; STL=us-asL
        endif
        if(ps>pR)then
            rsR=rR*(ps/pR+(g-1.0_dp)/(g+1.0_dp))/ &
                ((g-1.0_dp)/(g+1.0_dp)*ps/pR+1.0_dp)
            SR=uR+aR*sqrt((g+1.0_dp)*ps/(2.0_dp*g*pR)+(g-1.0_dp)/(2.0_dp*g))
        else
            rsR=rR*(ps/pR)**(1.0_dp/g); asR=aR*(ps/pR)**((g-1.0_dp)/(2.0_dp*g))
            SHR=uR+aR; STR=us+asR
        endif
        do i=1,size(x)
            xi=(x(i)-x0)/t
            if(xi<=us)then
                if(ps>pL)then
                    if(xi<=SL)then
                        r=rL; u=uL; p=pL
                    else
                        r=rsL; u=us; p=ps
                    endif
                else
                    if(xi<=SHL)then
                        r=rL; u=uL; p=pL
                    else if(xi>=STL)then
                        r=rsL; u=us; p=ps
                    else
                        u=2.0_dp/(g+1.0_dp)*(aL+0.5_dp*(g-1.0_dp)*uL+xi)
                        a=2.0_dp/(g+1.0_dp)*(aL+0.5_dp*(g-1.0_dp)*(uL-xi))
                        r=rL*(a/aL)**(2.0_dp/(g-1.0_dp))
                        p=pL*(a/aL)**(2.0_dp*g/(g-1.0_dp))
                    endif
                endif
            else
                if(ps>pR)then
                    if(xi>=SR)then
                        r=rR; u=uR; p=pR
                    else
                        r=rsR; u=us; p=ps
                    endif
                else
                    if(xi>=SHR)then
                        r=rR; u=uR; p=pR
                    else if(xi<=STR)then
                        r=rsR; u=us; p=ps
                    else
                        u=2.0_dp/(g+1.0_dp)*(-aR+0.5_dp*(g-1.0_dp)*uR+xi)
                        a=2.0_dp/(g+1.0_dp)*(aR-0.5_dp*(g-1.0_dp)*(uR-xi))
                        r=rR*(a/aR)**(2.0_dp/(g-1.0_dp))
                        p=pR*(a/aR)**(2.0_dp*g/(g-1.0_dp))
                    endif
                endif
            endif
            Q(:,i)=[r,u,p,p/(r*rg)]
        enddo
    end subroutine
end module exact