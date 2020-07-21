!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: September 13, 2018
!
! DESCRIPTION
!
! INPUT
! n = sample size                                                       (integer)
! p = number of variables                                               (integer)
! X = n x p dimensional matrix                                          (double)
! R = n x (0:p) dimensional matrix used to encode the                   (integer)
!       patterns of non-observed values
! startmis = the starting row of the missing values                     (integer)
! w = p x p dimensional matrix of positive weights                      (double)
! pendiag = integer specifying if the diagonal elements                 (integer)
!           of the concentration matrix are penalized
!           '0' unpenalized
!           '1' penalized
! xm = starting values of the marginal expected values                  (double)
! nrho = number of tuning values                                        (double)
! rhoratio = smallest value for the tuning parameter                    (double)
!            is defined as fraction of maxrho, i.e.
!            minrho = rhoratio * maxrho
! rho = (optional) A user supplied rho sequence of length               (double)
!       nrho
! maxR2 = upper limit for the R2 statistic                              (double)
! maxit_em = maximum number of iterations for the EM algorithm          (integer)
! thr_em = convergence threshold for EM algorithm                       (double)
! maxit_bcd = maximum number of iterations for the BCD algorithm        (integer)
! thr_bcd = convergence threshold for BCD algorithm                     (double)
! trace = integer for printing out information on video                 (integer)
!
! OUTPUT
! Xipt = n x p x nrho dimensional array;                                (double)
!           X(:, :, k) is the matrix X with censored values
!           imputed using the kth cglasso model
! S = n x p x nrho dimensional array;                                   (double)
! mu = p x k dimensional matrix;                                        (double)
! Sgm = n x p x nrho dimensional array;                                 (double)
! Tht = n x p x nrho dimensional array;                                 (double)
! Adj = n x p x nrho dimensional array;                                 (integer)
! df = the number of nonzero partial corr for each value of rho         (integer)
! R2 = the value of the R2 statistic for each value of rho              (double)
! ncomp = number of connected components                                (integer)
! Ck = matrix used to identify the connected components                 (integer)
! pk = p-dimensiona vector used to encode the connected components      (integer)
! nit = number of steps untile convergence is met                       (integer)
! conv = integer used to encode the convergence of the algorithm        (integer)
!           '-1' error in memory allocation
!            '0' convergence is met
!            '1' maximum number of iterations has been exceeded
!            '2' error in computing tmean and tvar
!            '3' Tht_mm inversion failed
! subrout = the integer used to encode the subroutine                   (integer)
!           where the error occurred:
!           '0' convergence is met
!           '1' Fit_MarginalDistributions
!           '2' update
!           '3' glasso
!           '4' cglasso
subroutine mglasso_algo1(n,p,X,R,startmis,w,pendiag,nrho,rhoratio,rho,maxR2,maxit_em,thr_em,maxit_bcd,thr_bcd,&
Xipt,S,mu,Sgm,Tht,Adj,df,R2,ncomp,Ck,pk,nit,conv,subrout,trace)
implicit none
integer :: n,p,R(n,0:p),startmis,pendiag,nrho,maxit_em,maxit_bcd,Adj(p,p,nrho),df(nrho),ncomp(nrho),Ck(p,nrho),pk(p,nrho)
integer :: nit(nrho,2),conv,subrout,trace
double precision :: X(n,p),w(p,p),rhoratio,rho(nrho),maxR2,thr_em,thr_bcd
double precision :: Xipt(n,p,nrho),S(p,p,nrho),mu(p,nrho),Sgm(p,p,nrho),Tht(p,p,nrho),R2(nrho)
! internal variables
integer :: no,i,j,k,ii,nnit
double precision :: minrho,maxrho,T1o(p),T2o(p,p),xm(p),Xipt_k(n,p),S_k(p,p),mu_k(p),Sgm_k(p,p),Tht_k(p,p)
double precision :: ratio,T1(p),T2(p,p),mu_o(p),Sgm_o(p,p),Tht_o(p,p),dmu,dTht,rho_gl(p,p),R2num,R2den
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing the observed statistics T1o and T2o !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T1o = 0.d0
do j = 1, p
    no = 0
    do i = 1, n
        if(R(i, j).eq.0) then
            no = no + 1
            T1o(j) = T1o(j) + X(i, j)
        end if
    end do
    xm(j) = T1o(j) / no
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing Xipt_k                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Xipt_k = X
do j = 1, p
    do i = 1, n
        if(R(i, j).eq.1) Xipt_k(i, j) = xm(j)
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing the Sbar matrix and maxrho   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
maxrho = 0.d0
do j = 1, p
    S_k(j, j) = dot_product(Xipt_k(:, j), Xipt_k(:, j)) / n - xm(j) * xm(j)
    do i = j + 1, p
        S_k(i, j) = dot_product(Xipt_k(:, i), Xipt_k(:, j)) / n - xm(i) * xm(j)
        S_k(j, i) = S_k(i, j)
        if(w(1, 1).lt.0.d0) then
            maxrho = max(maxrho, abs(S_k(i, j)))
        else
            if(w(i, j).gt.0.d0) then
                maxrho = max(maxrho, abs(S_k(i, j) / w(i, j)))
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
        R2den = R2den + S_k(i, j)**2
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing rho-values                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(rho(1).lt.0.d0) then
    rho(1) = maxrho
    minrho = rhoratio * maxrho
    ratio = exp((log(minrho) - log(maxrho)) / (nrho - 1.d0))
    do k = 2, nrho
        rho(k) = ratio * rho(k - 1)
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inizializing mu_k, Sgm_k and Tht_k      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mu_k = xm
Tht_k = 0.d0
Sgm_k = 0.d0
do j = 1, p
    Sgm_k(j, j) = S_k(j, j)
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
    if(trace.eq.2) call mglasso_trace2_1(k, rho(k))
    if(rho(k).lt.maxrho) then
        mu_o = mu_k
        Sgm_o = Sgm_k
        Tht_o = Tht_k
        rho_gl = rho(k)
        if(w(1, 1).ge.0.d0) rho_gl = rho_gl * w
        do ii = 1, maxit_em
            if(trace.eq.2) call cglasso_trace2_2(ii)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! updating summary statistics T1 and T2 imputing missing values !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call impute(startmis, n, p, X, R, T1o, T2o, mu_o, Tht_o, Xipt_k, T1, T2, conv)
            if(conv.ne.0) then
                subrout = 1
                exit
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! updating 'mu_k' and 'S_k'!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            if(trace.eq.2) call cglasso_trace2_4()
            nnit = 0
            call glassosub(p,S_k,pendiag,rho_gl,maxit_bcd,thr_bcd,Sgm_k,Tht_k,ncomp(k),Ck(:,k),pk(:,k),nnit,conv,trace)
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
    else
        rho_gl = rho(k)
        if(w(1, 1).ge.0.d0) rho_gl = rho_gl * w
        call Find_ConnectedComp(p, S_k, rho_gl, ncomp(k), Ck(:,k), pk(:,k))
        if(trace.eq.2) then
            call cglasso_trace2_2(1)
            call cglasso_trace2_3()
            call cglasso_trace2_5()
        end if
    end if
    if(conv.eq.0) then
        if(trace.eq.1) call mglasso_trace1(k, rho(k), nit(k, 1), nit(k, 2))
        Xipt(:, :, k) = Xipt_k
        S(:, :, k) = S_k
        mu(:, k) = mu_k
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
            R2num = R2num + (S_k(j, j) - Sgm_k(j, j))**2
            do i = j + 1, p
                R2num = R2num + (S_k(i, j) - Sgm_k(i, j))**2
            end do
        end do
        R2(k) = 1.d0 - R2num / R2den
        if(R2(k).ge.maxR2) then
            nrho = k
            exit
        end if
    else
        nrho = k - 1
        exit
    end if
end do
end subroutine mglasso_algo1
