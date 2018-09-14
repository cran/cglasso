!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 16, 2018
!
! Description
! 'Fit_MarginalDistributions' is used to estimates the parameters of the marginal normal distribution
! under the assumption that a part of boservations are right-censored
!
! Input
! n = sample size                                                                       (integer)
! p = number of variables                                                               (integer)
! X = n x p dimensional matrix                                                          (double)
! lo = p-dimensional vector of left-censoring values                                    (double)
! up = p-dimensional vector of right-censoring values                                   (double)
! R = (0:n) x (-1:p) dimensional matrix used to save the information about              (integer)
!       the censored values (see 'setup.f90')
! T1o = observed sufficient statistic                                                   (double)
! T2o = observed sufficient statistic                                                   (double)
! nstp = number of steps for the EM algorithm                                           (integer)
! eps = threshold value used for the convergence of the EM algorithm                    (double)
! tol = value use to compute the conditional expected value                             (double)
!
! Output
! xm = p - dimensinal vector of marginal expected values                                (double)
! vm = p - dimensinal vector of marginal variances                                      (double)
! Xipt = X matrix with imputed censored values                                          (double)
! conv = integer used to encode the convergence of the algorithm:                       (integer)
!           '-1' error in memory allocation
!            '0' convergence is met
!            '1' maximum number of steps is reached
!            '2' error in computing tmean and tvar
subroutine Fit_MarginalDistributions(n,p,X,lo,up,R,T1o,T2o,nstp,eps,tol,xm,vm,Xipt,conv)
integer :: n,p,R(0:n,-1:p),nstp,conv
double precision :: X(n,p),lo(p),up(p),T1o(p),T2o(p,p),eps,tol,xm(p),vm(p),Xipt(n,p)
! internal variables
integer :: i,j,k,no,ncr,ncl
double precision :: z,ratio,tmean_r,tmean_l,tvar,xm_o,xm_n,vm_o,vm_n,T1b,T2b,dm,dv
double precision, external :: dnorm, pnorm
double precision, dimension(:), allocatable :: xobs
Xipt = X
do j = 1, p
    no = count(R(1:n, j).eq.0)
    ncr = count(R(1:n, j).eq.1)
    ncl = n - no - ncr
    if(no.ne.n) then
        allocate(xobs(1:no), stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        k = 0
        do i = 1, n
            if(R(i, j).eq.0) then
                k = k + 1
                xobs(k) = X(i, j)
            end if
        end do
        xm_o = xm(j)
        vm_o = vm(j)
        xm_n = 0.d0
        vm_n = 0.d0
        tmean_r = 0.d0
        tmean_l = 0.d0
        do i = 1, nstp
            T1b = T1o(j)
            T2b = T2o(j, j)
            if(ncr.ne.0) then
                z = (up(j) - xm_o) / sqrt(vm_o)
                if(dnorm(z).gt.tol * pnorm(z, .true.)) then
                    ratio = dnorm(z) / pnorm(z, .false.)
                    tmean_r = xm_o + sqrt(vm_o) * ratio
                    tvar = vm_o * (1.d0 + z * ratio - ratio**2)
                    T1b = T1b + ncr * tmean_r
                    T2b = T2b + ncr * (tvar + tmean_r**2)
                else
                    R(0, j) = 0
                    xm_n = xm(j)
                    vm_n = vm(j)
                    tmean_l = xm_n
                    tmean_r = xm_n
                    exit
                end if
            end if
            if(ncl.ne.0) then
                z = (lo(j) - xm_o) / sqrt(vm_o)
                if(dnorm(z).gt.tol * pnorm(z, .false.)) then
                    ratio = dnorm(z) / pnorm(z, .true.)
                    tmean_l = xm_o - sqrt(vm_o) * ratio
                    tvar = vm_o * (1.d0 - z * ratio - ratio**2)
                    T1b = T1b + ncl * tmean_l
                    T2b = T2b + ncl * (tvar + tmean_l**2)
                else
                    R(0, j) = 0
                    xm_n = xm(j)
                    vm_n = vm(j)
                    tmean_l = xm_n
                    tmean_r = xm_n
                    exit
                end if
            end if
            xm_n = T1b / n
            vm_n = T2b / n - xm_n**2
            dm = abs(xm_n - xm_o)
            dv = abs(vm_n - vm_o)
            if(max(dm, dv).le.eps) exit
            xm_o = xm_n
            vm_o = vm_n
            if(xm_o.le.lo(j).or.xm_o.ge.up(j)) then
                conv = 2
                return
            end if
        end do
        if(i.ge.nstp) then
            conv = 1
            return
        end if
        xm(j) = xm_n
        vm(j) = vm_n
        do i = 1, n
            if(R(i, j).eq.-1) Xipt(i, j) = tmean_l
            if(R(i, j).eq.1)  Xipt(i, j) = tmean_r
        end do
        deallocate(xobs, stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
    end if
end do
end subroutine Fit_MarginalDistributions
