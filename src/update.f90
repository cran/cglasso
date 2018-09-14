!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 2.0.0
! Data: August 18, 2018
!
! DESCRIPTION
! 'update' is used to compute the conditional expected vaules of the sufficients statistics using
! the approximation proposed in Guo et al. (2015). See also Section 3.2 in Augugliaro et al. (2018).
!
! REFERENCES
! Augugliaro L., Abbruzzo A., and Vinciotti V. (2018) l1-Penalized Censored Gaussian Graphical Model.
!   Biostatistics (to appear)
! Guo, J., Levina, E., Michailidis, G. and Zhu, J. (2015). Graphical models for ordinal data.
!   Journal of Computational and Graphical Statistics 24(1), 183–204.
!
! INPUT
! starmis =                                                                 (integer)
! n = sample size                                                           (integer)
! p = number of variables                                                   (integer)
! X = n x p dimensional matrix                                              (double)
! lo = p dimensional vector of left-censoring values                        (double)
! up = p dimensional vector of right-censoring values                       (double)
! T1o = observed sufficient statistic for mu                                (double)
! T2o = observed sufficient statistic for Sgm                               (double)
! R = (0:n) x (-1:p) matrix used to encode the censoring patterns:          (integer)
!       (see 'setup.f90')
! mu = current estimate of the expected value                               (double)
! Tht = current estimate of the precision matrix                            (double)
!
! Output
! Xipt = the matrix X with missing values imputed using                     (double)
!           the kth cglasso model
! T1 = conditional expected value of the statistic T1                       (double)
! T2 = conditional expected value of the statistic T2 computed using the a  (double)
!       pproximation proposed in Guo et al. (2015).
! conv = integer used to encode errors                                      (intger)
!       '-1' error in memory allocation
!        '2' error in computing tmean and tvar
!        '3' Tht_mm inversion failed
subroutine update(startmis,n,p,X,lo,up,R,T1o,T2o,mu,Tht,tol,Xipt,T1,T2,conv)
implicit none
integer :: startmis,n, p,R(0:n,-1:p),conv
double precision :: X(n,p),lo(p),up(p),T1o(p),T2o(p,p),mu(p),Tht(p,p),tol,Xipt(n,p),T1(p),T2(p,p)
! internal variables
integer :: i,j,j1,j2,j3,no,nlc,nrc,nmis
integer, dimension (:), allocatable :: j_ob,j_lc,j_rc,j_mis
double precision, external :: dnorm, pnorm
double precision :: z,ratio
double precision, dimension (:), allocatable :: zeta,mu_m_given_o,tmean,tvar
double precision, dimension (:,:), allocatable :: Sgm_m_given_o
T1 = T1o
T2 = T2o
i = startmis
nlc = R(i, -1)
nrc = R(i, 0)
nmis = nlc + nrc
no = p - nmis
allocate(j_ob(1:no),j_mis(1:nmis),zeta(1:nmis),mu_m_given_o(1:nmis),tmean(1:nmis),tvar(1:nmis),stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
allocate(Sgm_m_given_o(1:nmis,1:nmis),stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
if(nlc.ne.0) then
    allocate(j_lc(1:nlc),stat=conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end if
if(nrc.ne.0) then
    allocate(j_rc(1:nrc),stat=conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end if
j1 = 0
j2 = 0
j3 = 0
do j = 1, p
    if(R(i, j).eq.0) then
        j1 = j1 + 1
        j_ob(j1) = j
    end if
    if(R(i, j).eq.-1) then
        j2 = j2 + 1
        j_lc(j2) = j
    end if
    if(R(i, j).eq.1) then
        j3 = j3 + 1
        j_rc(j3) = j
    end if
end do
j1 = 0
do j2 = 1, nlc
    j1 = j1 + 1
    j_mis(j1) = j_lc(j2)
end do
do j3 = 1, nrc
    j1 = j1 + 1
    j_mis(j1) = j_rc(j3)
end do
call inv(nmis, Tht(j_mis, j_mis), Sgm_m_given_o, conv)
if(conv.ne.0) then
    conv = 3
    return
end if
do
    zeta = matmul(Tht(j_mis, j_ob), X(i, j_ob) - mu(j_ob))
    mu_m_given_o = mu(j_mis) - matmul(Sgm_m_given_o, zeta)
    j1 = 0
    j2 = 0
    do j = 1, nmis
        if(j.le.nlc) then
            j1 = j1 + 1
            if(R(i, j_lc(j1)).ne.0) then
                z = (lo(j_lc(j1)) - mu_m_given_o(j)) / sqrt(Sgm_m_given_o(j, j))
                if(dnorm(z).gt.tol * pnorm(z, .false.)) then
                    ratio = dnorm(z) / pnorm(z, .true.)
                    tmean(j) = mu_m_given_o(j) - sqrt(Sgm_m_given_o(j, j)) * ratio
                    tvar(j) = Sgm_m_given_o(j, j) * (1.d0 - z * ratio - ratio**2)
                else
                    tmean(j) = mu_m_given_o(j)
                    tvar(j) = Sgm_m_given_o(j, j)
!                    R(i, j_lc(j1)) = 0
                end if
            else
                tmean(j) = mu_m_given_o(j)
                tvar(j) = Sgm_m_given_o(j, j)
            end if
        else
            j2 = j2 + 1
            if(R(i, j_rc(j2)).ne.0) then
                z = (up(j_rc(j2)) - mu_m_given_o(j)) / sqrt(Sgm_m_given_o(j, j))
                if(dnorm(z).gt.tol * pnorm(z, .true.)) then
                    ratio = dnorm(z) / pnorm(z, .false.)
                    tmean(j) = mu_m_given_o(j) + sqrt(Sgm_m_given_o(j, j)) * ratio
                    tvar(j) = Sgm_m_given_o(j, j) * (1.d0 + z * ratio - ratio**2)
                else
                    tmean(j) = mu_m_given_o(j)
                    tvar(j) = Sgm_m_given_o(j, j)
!                    R(i, j_rc(j2)) = 0
                end if
            else
                tmean(j) = mu_m_given_o(j)
                tvar(j) = Sgm_m_given_o(j, j)
            end if
        end if
    end do
    Xipt(i, j_mis) = tmean
    T1(j_mis) = T1(j_mis) + tmean
    do j = 1, nmis
        T2(j_mis, j_mis(j)) = T2(j_mis, j_mis(j)) + tmean * tmean(j)
        T2(j_mis(j), j_mis(j)) = T2(j_mis(j), j_mis(j)) + tvar(j)
        T2(j_ob, j_mis(j)) = T2(j_ob, j_mis(j)) + X(i, j_ob) * tmean(j)
    end do
    T2(j_mis, j_ob) = transpose(T2(j_ob, j_mis))
    i = i + 1
    if(i.gt.n) exit
    if(R(i, 0).ne.0.or.R(i, -1).ne.0) then
        deallocate(j_ob,j_mis,zeta,mu_m_given_o,tmean,tvar,Sgm_m_given_o,stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        if(nlc.ne.0) then
            deallocate(j_lc,stat=conv)
            if(conv.ne.0) then
                conv = -1
                return
            end if
        end if
        if(nrc.ne.0) then
            deallocate(j_rc,stat=conv)
            if(conv.ne.0) then
                conv = -1
                return
            end if
        end if
        nlc = R(i, -1)
        nrc = R(i, 0)
        nmis = nlc + nrc
        no = p - nmis
        allocate(j_ob(1:no),j_mis(1:nmis),zeta(1:nmis),mu_m_given_o(1:nmis),tmean(1:nmis),tvar(1:nmis),stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        allocate(Sgm_m_given_o(1:nmis,1:nmis), stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        if(nlc.ne.0) then
            allocate(j_lc(1:nlc),stat=conv)
            if(conv.ne.0) then
                conv = -1
                return
            end if
        end if
        if(nrc.ne.0) then
            allocate(j_rc(1:nrc),stat=conv)
            if(conv.ne.0) then
                conv = -1
                return
            end if
        end if
        j1 = 0
        j2 = 0
        j3 = 0
        do j = 1, p
            if(R(i, j).eq.0) then
                j1 = j1 + 1
                j_ob(j1) = j
            end if
            if(R(i, j).eq.-1) then
                j2 = j2 + 1
                j_lc(j2) = j
            end if
            if(R(i, j).eq.1) then
                j3 = j3 + 1
                j_rc(j3) = j
            end if
        end do
        j1 = 0
        do j2 = 1, nlc
            j1 = j1 + 1
            j_mis(j1) = j_lc(j2)
        end do
        do j3 = 1, nrc
            j1 = j1 + 1
            j_mis(j1) = j_rc(j3)
        end do
        call inv(nmis, Tht(j_mis, j_mis), Sgm_m_given_o, conv)
        if(conv.ne.0) then
            conv = 3
            return
        end if
    end if
end do
deallocate(j_ob,j_mis,zeta,mu_m_given_o,tmean,tvar,Sgm_m_given_o,stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
if(nlc.ne.0) then
    deallocate(j_lc,stat=conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end if
if(nrc.ne.0) then
    deallocate(j_rc,stat=conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end if
end subroutine update














