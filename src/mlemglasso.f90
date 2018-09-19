!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: September 17, 2018
!
subroutine mlemglasso(n,p,X,R,startmis,nrho,maxit_em,thr_em,maxit_bcd,thr_bcd,Xipt,S,mu,Sgm,Tht,R2,nit,conv,subrout,trace)
implicit none
integer :: n,p,R(n,0:p),startmis,nrho,maxit_em,maxit_bcd,nit(nrho,2),conv,subrout,trace
double precision :: X(n,p),thr_em,thr_bcd
double precision :: Xipt(n,p,nrho),S(p,p,nrho),mu(p,nrho),Sgm(p,p,nrho),Tht(p,p,nrho),R2(nrho)
! internal variables
integer :: i,j,k,ii,nnit,ncomp,Ck(p),pk(p)
double precision, parameter :: big = huge(1.d0)
double precision :: T1o(p),T2o(p,p),T1(p),T2(p,p),mu_o(p),Sgm_o(p,p),Tht_o(p,p),S_k(p,p),mu_k(p),Sgm_k(p,p),Tht_k(p,p)
double precision :: Xipt_k(n,p),rho_gl(p,p),dmu,dTht,R2num,R2den
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing the observed statistics T1o and T2o !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T1o = 0.d0
do j = 1, p
    do i = 1, n
        if(R(i, j).eq.0) then
            T1o(j) = T1o(j) + X(i, j)
        end if
    end do
end do
T2o = 0.d0
do i = 1, p
    do j = i, p
        do k = 1, n
            if((R(k, i).eq.0).and.(R(k, j).eq.0)) then
                T2o(i, j) = T2o(i, j) + X(k, i) * X(k, j)
            end if
        end do
        T2o(j, i) = T2o(i, j)
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Xipt_k = X
mu_k = mu(:, 1)
Sgm_k = Sgm(:, :, 1)
Tht_k = Tht(:, :, 1)
do k = 1, nrho
    if(trace.eq.2) call mlemglasso_trace2_1(k)
    rho_gl = 0.d0
    do j = 1, p
        do i = j + 1, p
            if(abs(Tht(i, j, k)).eq.0.d0) then
                rho_gl(i, j) = big
                rho_gl(j, i) = big
            end if
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! starting EM algorithm !
    !!!!!!!!!!!!!!!!!!!!!!!!!
    mu_o = mu_k
    Sgm_o = Sgm_k
    Tht_o = Tht_k
    do ii = 1, maxit_em
        if(trace.eq.2) call cglasso_trace2_2(ii)
        call impute(startmis, n, p, X, R, T1o, T2o, mu_o, Tht_o, Xipt_k, T1, T2, conv)
        if(conv.ne.0) then
            subrout = 1
            exit
        end if
        mu_k = T1 / n
        do j = 1, p
            S_k(:, j) = T2(:, j) / n - mu_k * mu_k(j)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! inizializing the matrix Sgm_k !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Sgm_k = S_k
        do j = 1, p
            do i = j + 1, p
                if(abs(Tht_o(i, j)).gt.0.d0) then
                    Sgm_k(i, j) = Sgm_k(i, j) + rho_gl(i, j) * sign(1.d0, Tht_o(i, j))
                    Sgm_k(j, i) = Sgm_k(i, j)
                end if
            end do
        end do
        if(trace.eq.2) call cglasso_trace2_3()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! solving the glasso problem  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nnit = 0
        ncomp = 0
        Ck = 0
        pk = 0
        if(trace.eq.2) call mlemglasso_trace2_4()
        call glassosub(p, S_k, 0, rho_gl, maxit_bcd, thr_bcd, Sgm_k, Tht_k, ncomp, Ck, pk, nnit, conv, trace)
        if(conv.ne.0) then
            subrout = 2
            exit
        end if
        if(trace.eq.2) call cglasso_trace2_5()
        nit(k, 1) = ii
        nit(k, 2) = nit(k, 2) + nnit
        Sgm_k = (Sgm_k + transpose(Sgm_k)) / 2.d0
        Tht_k = (Tht_k + transpose(Tht_k)) / 2.d0
        dmu = maxval(abs(mu_o - mu_k))
        dTht = 0.d0
        do j = 1, p
            do i = 1, j
                dTht = max(dTht, abs(Tht_o(i, j) - Tht_k(i, j)))
            end do
        end do
        if(trace.eq.2) call cglasso_trace2_6(thr_em, dmu, dTht)
        if(max(dmu, dTht).le.thr_em) then
            exit
        else
            mu_o = mu_k
            Tht_o = Tht_k
        end if
    end do
    if(ii.ge.maxit_em) then
        conv = 1
        subrout = 3
    end if
    if(conv.eq.0) then
        if(trace.eq.1) call mlemglasso_trace1(k, nit(k, 1), nit(k, 2))
        Xipt(:, :, k) = Xipt_k
        S(:, :, k) = S_k
        mu(:, k) = mu_k
        Sgm(:, :, k) = Sgm_k
        Tht(:, :, k) = Tht_k
        R2num = 0.d0
        R2den = 0.d0
        do j = 1, p
            R2num = R2num + (S_k(j, j) - Sgm_k(j, j))**2
            do i = j + 1, p
                R2num = R2num + (S_k(i, j) - Sgm_k(i, j))**2
                R2den = R2den + S_k(i, j)**2
            end do
        end do
        R2(k) = 1.d0 - R2num / R2den
    else
        nrho = k - 1
        exit
    end if
end do
end subroutine mlemglasso
