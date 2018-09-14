!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 19, 2018
!
! DESCRIPTION
! 'setup' is used to the matrices X and R
!
! INPUT
! n = sample size                                                                                   (integer)
! p = number of variables                                                                           (integer)
! X = n x p dimensional matrix                                                                      (double)
! lo = p dimensional vector of left censored values                                                 (double)
! up = p dimensional vector of right censored values                                                (double)
!
! OUTPUT
! X = in output X is ordered to identify the patterns of non-observed values                        (double)
! R = (0:n) x (-1:p) dimensional matrix used to save the information about                          (integer)
!       the censored values:
!       'R(i, j) = -1' (i = startmis...n and j = 1...p) means that x_{ij} is left censored
!       'R(i, j) =  0' (i = startmis...n and j = 1...p) means that x_{ij} is a missing value
!       'R(i, j) = +1' (i = startmis...n and j = 1...p) means that x_{ij} is right censored
!       'R(0, j) = 1' means that the jth variable should be treat as censored
!       'R(0, j) = 0' means that the jth variable should be treat as missing
!       'R(i, -1)' is the number of left censored values in x_i
!       'R(i, 0)'  is the number of right censored values in x_i
! startmis = the starting row of the censored values                                                (integer)
subroutine setup(n, p, X, lo, up, R, startmis)
implicit none
integer :: n, p, R(0:n, -1:p), startmis
double precision :: X(n, p), lo(p), up(p)
! internal variables
integer :: i,j,jj, tmp1(p)
double precision :: tmp2(p)
R = 0
do j = 1, p
    where (X(:, j).ge.up(j)) R(1:n, j) = 1
    where (X(:, j).le.lo(j)) R(1:n, j) = -1
end do
jj = -1
do
    jj = jj + 1
    if(jj.le.n) R(jj, 0) = 1
    if(jj.ge.n) exit
    j = jj
    do i = j + 1, n
        if(all(R(j, 1:p).eq.R(i, 1:p))) then
            jj = jj + 1
            tmp1 = R(jj, 1:p)
            R(jj, 1:p) = R(i, 1:p)
            R(i, 1:p) = tmp1
            tmp2 = X(jj, :)
            X(jj, :) = X(i, :)
            X(i, :) = tmp2
        end if
    end do
end do
R(0,  1:p) = 1
R(0, -1:0) = 0
do i = 1, n
    if(R(i, 0).eq.1) then
        R(i, 0) = count(R(i, 1:p).eq.1)
        R(i, -1) = count(R(i, 1:p).eq.-1)
    end if
end do
startmis = 0
do i = 1, n
    if(R(i, -1).ne.0 .or. R(i, 0).ne.0) then
        startmis = i
        exit
    end if
end do
end subroutine setup
