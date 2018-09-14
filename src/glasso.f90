!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 27, 2018
!
subroutine glasso(p,S,w,pendiag,nrho,rhoratio,rho,maxR2,maxit,thr,Sgm,Tht,Adj,df,R2,ncomp,Ck,pk,nit,conv,trace)
implicit none
integer :: p,pendiag,nrho,maxit,Adj(p,p,nrho),df(nrho),ncomp(nrho),Ck(p,nrho),pk(p,nrho),nit(nrho),conv,trace
double precision :: S(p,p),w(p,p),rhoratio,rho(nrho),maxR2,thr,Sgm(p,p,nrho),Tht(p,p,nrho),R2(nrho)
! internal variables
integer :: i,j,k
double precision :: Sgm_k(p,p),Tht_k(p,p),maxrho,minrho,ratio,rho_gl(p,p),R2num,R2den
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing maxrho            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
maxrho = 0.d0
do j = 1, p
    do i = j + 1, p
        if(w(1, 1).lt.0.d0) then
            maxrho = max(maxrho, abs(S(i, j)))
        else
            if(w(i, j).gt.0.d0) then
                maxrho = max(maxrho, abs(S(i, j) / w(i, j)))
            end if
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing R2den             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
R2den = 0.d0
do j = 1, p
    if(pendiag.eq.1) R2den = R2den + (maxrho * w(j, j))**2
    do i = j + 1, p
        R2den = R2den + S(i, j)**2
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing rho-values        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(rho(1).lt.0.d0) then
    rho(1) = maxrho
    minrho = rhoratio * maxrho
    ratio = exp((log(minrho) - log(maxrho)) / (nrho - 1.d0))
    do k = 2, nrho
        rho(k) = ratio * rho(k - 1)
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inizializing Sgm_k and Tht_k            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Sgm_k = 0.d0
Tht_k = 0.d0
do j = 1, p
    Sgm_k(j, j) = S(j, j)
    if(pendiag.eq.1) then
        if(w(j, j).lt.0.d0) then
            Sgm_k(j, j) = Sgm_k(j, j) + maxrho
        else
            Sgm_k(j, j) = Sgm_k(j, j) + maxrho * w(j, j)
        end if
    end if
    Tht_k(j, j) = 1.d0 / Sgm_k(j, j)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k = 1, nrho
    if(trace.eq.2) call glasso_trace2_1(k, rho(k))
    rho_gl = rho(k)
    if(w(1, 1).ge.0.d0) rho_gl = rho_gl * w
    call glassosub(p,S,pendiag,rho_gl,maxit,thr,Sgm_k,Tht_k,ncomp(k),Ck(:,k),pk(:,k),nit(k),conv,trace)
    if(trace.eq.1) call glasso_trace1(k, rho(k), nit(k))
    if(conv.ne.0) then
        nrho = k - 1
        exit
    end if
    Sgm(:, :, k) = Sgm_k
    Tht(:, :, k) = Tht_k
    do j = 1, p
        do i = j + 1, p
            if(abs(Tht_k(i, j)).gt.0.d0) then
                df(k) = df(k) + 1
                Adj(i, j, k) = 1
                Adj(j, i, k) = 1
            end if
        end do
    end do
    R2num = 0.d0
    do j = 1, p
        do i = j, p
            R2num = R2num + (S(i, j) - Sgm_k(i, j))**2
        end do
    end do
    R2(k) = 1.d0 - R2num / R2den
    if(R2(k).ge.maxR2) then
        nrho = k
        exit
    end if
end do
end subroutine glasso




































