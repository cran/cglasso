function dnorm(x)
implicit none
double precision :: x, z, dnorm
double precision, parameter :: max_dbl = huge(1.d0)
double precision, parameter :: m_1_sqrt_2pi = 0.398942280401432677939946059934d0
z = x
z = abs(z)
if(z >= 2.d0 * sqrt(max_dbl)) then
    dnorm = 0.d0
end if
dnorm = m_1_sqrt_2pi * exp(-0.5d0 * z * z )
end
