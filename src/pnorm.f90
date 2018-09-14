function pnorm(x, lower)
implicit none

real ( kind = 8 ), parameter, dimension (5) :: a = (/ &
    2.2352520354606839287D00, &
    1.6102823106855587881D02, &
    1.0676894854603709582D03, &
    1.8154981253343561249D04, &
    6.5682337918207449113D-2 /)

real ( kind = 8 ), parameter, dimension (4) :: b = (/ &
    4.7202581904688241870D01, &
    9.7609855173777669322D02, &
    1.0260932208618978205D04, &
    4.5507789335026729956D04 /)

real ( kind = 8 ), parameter, dimension (9) :: c = (/ &
    3.9894151208813466764D-1, &
    8.8831497943883759412D00, &
    9.3506656132177855979D01, &
    5.9727027639480026226D02, &
    2.4945375852903726711D03, &
    6.8481904505362823326D03, &
    1.1602651437647350124D04, &
    9.8427148383839780218D03, &
    1.0765576773720192317D-8 /)

real ( kind = 8 ), parameter, dimension (8) :: d = (/ &
    2.2266688044328115691D01, &
    2.3538790178262499861D02, &
    1.5193775994075548050D03, &
    6.4855582982667607550D03, &
    1.8615571640885098091D04, &
    3.4900952721145977266D04, &
    3.8912003286093271411D04, &
    1.9685429676859990727D04 /)

real ( kind = 8 ), parameter, dimension (6) :: p = (/ &
    2.1589853405795699D-1,  &
    1.274011611602473639D-1,&
    2.2235277870649807D-2,  &
    1.421619193227893466D-3,&
    2.9112874951168792D-5,  &
    2.307344176494017303D-2 /)

real ( kind = 8 ), parameter, dimension (5) :: q = (/ &
    1.28426009614491121D00, &
    4.68238212480865118D-1, &
    6.59881378689285515D-2, &
    3.78239633202758244D-3, &
    7.29751555083966205D-5 /)

integer :: i
logical :: lower, upper, check
double precision :: x, xden, xnum, temp, del, xsq, y, cum, ccum, pnorm
double precision, parameter :: m_1_sqrt_2pi = 0.398942280401432677939946059934d0
double precision, parameter :: m_sqrt_32 = 5.656854249492380195206754896838d0
double precision, parameter :: min_dbl = tiny(1.d0)
double precision, parameter :: eps_dbl = epsilon(1.d0) * 0.5d0

! prevedere riga nel caso in cui isnan(x) = TRUE

upper = .not.lower

y = abs(x)
if(y <= 0.67448975) then
    if(y > eps_dbl) then
        xsq = x * x
        xnum = a(5) * xsq
        xden = xsq
        do i = 1, 3
            xnum = (xnum + a(i)) * xsq
            xden = (xden + b(i)) * xsq
        end do
    else
        xnum = 0.0d0
        xden = 0.0d0
    end if
    temp = x * (xnum + a(4)) / (xden + b(4));
    cum = 0.5d0 + temp
    ccum = 0.5d0 - temp
else
    if(y <= m_sqrt_32) then
        xnum = c(9) * y
        xden = y
        do i = 1, 7
            xnum = (xnum + c(i)) * y
            xden = (xden + d(i)) * y
        end do
        temp = (xnum + c(8)) / (xden + d(8))
        xsq = aint(y * 16.d0) / 16.d0
        del = (y - xsq) * (y + xsq)
        cum = exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp
        ccum = 1.d0 - cum
        if(x > 0.d0) then
            temp = cum
            cum = ccum
            ccum = temp
        end if
    else
        check = (lower.and.(-37.5193 < x).and.(x < 8.2924)).or.(upper.and.(-8.2924 < x).and.(x < 37.5193))
        if(check) then
            xsq = 1.0d0 / (x * x)
            xnum = p(6) * xsq
            xden = xsq
            do i = 1, 4
                xnum = (xnum + p(i)) * xsq
                xden = (xden + q(i)) * xsq
            end do
            temp = xsq * (xnum + p(5)) / (xden + q(5))
            temp = (m_1_sqrt_2pi - temp) / y
            xsq = aint(x * 16.d0) / 16.d0
            del = (x - xsq) * (x + xsq)
            cum = exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp
            ccum = 1.d0 - cum
            if(x > 0.d0) then
                temp = cum
                cum = ccum
                ccum = temp
            end if
        else
            if(x > 0.0d0) then
                cum = 1.0d0
                ccum = 0.0d0
            else
                cum = 0.0d0
                ccum = 1.0d0
            end if
        end if
    end if
end if

if (cum < min_dbl) cum = 0.0d0
if (ccum < min_dbl) ccum = 0.0d0

if(lower) then
    pnorm = cum
else
    pnorm = ccum
end if
end
