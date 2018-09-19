!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: September 13, 2018
!
! INPUT
! starmis = integer identifying the row where are starting the missing      (integer)
!           values
! n = sample size                                                           (integer)
! p = number of variables                                                   (integer)
! X = n x p dimensional matrix                                              (double)
! T1o = observed sufficient statistic for mu                                (double)
! T2o = observed sufficient statistic for Sgm                               (double)
! R = n x (0:p) matrix used to encode the censoring patterns:               (integer)
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
subroutine impute(startmis,n,p,X,R,T1o,T2o,mu,Tht,Xipt,T1,T2,conv)
implicit none
integer :: startmis,n, p,R(n,0:p),conv
double precision :: X(n,p),T1o(p),T2o(p,p),mu(p),Tht(p,p),Xipt(n,p),T1(p),T2(p,p)
! internal variables
integer :: i,j,jo,jm,no,nmis
integer, dimension (:), allocatable :: j_ob,j_mis
double precision, dimension (:), allocatable :: zeta,mu_m_given_o
double precision, dimension (:,:), allocatable :: Sgm_m_given_o
T1 = T1o
T2 = T2o
i = startmis
nmis = R(i, 0)
no = p - nmis
allocate(j_ob(1:no),j_mis(1:nmis),zeta(1:nmis),mu_m_given_o(1:nmis),stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
allocate(Sgm_m_given_o(1:nmis,1:nmis),stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
jo = 0
jm = 0
do j = 1, p
    if(R(i, j).eq.0) then
        jo = jo + 1
        j_ob(jo) = j
    else
        jm = jm + 1
        j_mis(jm) = j
    end if
end do
call inv(nmis, Tht(j_mis, j_mis), Sgm_m_given_o, conv)
if(conv.ne.0) then
    conv = 3
    return
end if
do
    zeta = matmul(Tht(j_mis, j_ob), X(i, j_ob) - mu(j_ob))
    mu_m_given_o = mu(j_mis) - matmul(Sgm_m_given_o, zeta)
    Xipt(i, j_mis) = mu_m_given_o
    T1(j_mis) = T1(j_mis) + mu_m_given_o
    do j = 1, nmis
        T2(j_mis, j_mis(j)) = T2(j_mis, j_mis(j)) + mu_m_given_o * mu_m_given_o(j)
        T2(j_ob, j_mis(j)) = T2(j_ob, j_mis(j)) + X(i, j_ob) * mu_m_given_o(j)
    end do
    T2(j_mis, j_ob) = transpose(T2(j_ob, j_mis))
    i = i + 1
    if(i.gt.n) exit
    if(R(i, 0).ne.0) then
        deallocate(j_ob,j_mis,zeta,mu_m_given_o,Sgm_m_given_o,stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        nmis = R(i, 0)
        no = p - nmis
        allocate(j_ob(1:no),j_mis(1:nmis),zeta(1:nmis),mu_m_given_o(1:nmis),Sgm_m_given_o(1:nmis,1:nmis),stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        jo = 0
        jm = 0
        do j = 1, p
            if(R(i, j).eq.0) then
                jo = jo + 1
                j_ob(jo) = j
            else
                jm = jm + 1
                j_mis(jm) = j
            end if
        end do
        call inv(nmis, Tht(j_mis, j_mis), Sgm_m_given_o, conv)
        if(conv.ne.0) then
            conv = 3
            return
        end if

    end if
end do
deallocate(j_ob,j_mis,zeta,mu_m_given_o,Sgm_m_given_o,stat=conv)
if(conv.ne.0) then
    conv = -1
    return
end if
end subroutine impute
