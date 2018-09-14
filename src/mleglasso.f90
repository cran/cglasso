!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: August 13, 2018
!
subroutine mleglasso(p,S,nrho,maxit,thr,Sgm,Tht,R2,nit,conv,trace)
implicit none
integer :: p,nrho,maxit,nit(nrho),conv,trace
double precision :: S(p,p),thr,Sgm(p,p,nrho),Tht(p,p,nrho),R2(nrho)
! internal variables
integer :: i,j,k,ncomp,Ck(p),pk(p)
double precision, parameter :: big = huge(1.d0)
double precision :: Sgm_k(p,p),Tht_k(p,p),rho_gl(p,p),R2num,R2den
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing R2den             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
R2den = 0.d0
do j = 1, p
    do i = j + 1, p
        R2den = R2den + S(i, j)**2
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Sgm_k = Sgm(:, :, 1)
Tht_k = Tht(:, :, 1)
do k = 1, nrho
    rho_gl = 0.d0
    do j = 1, p
        do i = j + 1, p
            if(abs(Tht(i, j, k)).eq.0.d0) then
                rho_gl(i, j) = big
                rho_gl(j, i) = big
            end if
        end do
    end do
    ncomp = 0
    Ck = 0
    pk = 0
    if(trace.eq.2) call mleglasso_trace2(k)
    call glassosub(p,S,0,rho_gl,maxit,thr,Sgm_k,Tht_k,ncomp,Ck,pk,nit(k),conv,trace)
    if(trace.eq.1) call mleglasso_trace1(k, nit(k))
    if(conv.ne.0) then
        nrho = k - 1
        exit
    end if
    Sgm(:, :, k) = Sgm_k
    Tht(:, :, k) = Tht_k
    R2num = 0.d0
    do j = 1, p
        do i = j, p
            R2num = R2num + (S(i, j) - Sgm_k(i, j))**2
        end do
    end do
    R2(k) = 1.d0 - R2num / R2den
end do
end subroutine mleglasso
