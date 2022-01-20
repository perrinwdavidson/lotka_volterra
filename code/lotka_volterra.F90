!***********************************************************************
! lotka-volterra model
!-----------------------------------------------------------------------
! perrin w. davidson
! 19.01.21
! perrinwdavidson@gmail.com
!***********************************************************************
! lk_model
!-----------------------------------------------------------------------
! purpose: reads in initial values and parameters and integrates a set
! differential equations forward in time.
!
! inputs:
! x0 - initial prey value
! y0 - initial predator value
! alpha - prey birth rate
! beta - prey death rate from interaction with predators
! delta - predator growth rate from feeding on prey
! gamma - predator death rate
!
! outputs:
! x - integrated prey population
! y - integrated predator population
!***********************************************************************
! Define program -------------------------------------------------------
PROGRAM lotka_volterra

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Inputs:
    DOUBLE PRECISION, DIMENSION(6) :: args = (/ 0.5, 0.5, 2./3., 4./3., 1., 1. /)  ! [x0, y0, alpha, beta, delta, gamma]
    DOUBLE PRECISION, DIMENSION(4) :: argchange = (/ 0.5, 0.75, 1., 1.5 /)
    INTEGER, PARAMETER :: ichange = 2  ! index for arg that is going to be changed

    ! Model parameters:
    DOUBLE PRECISION, PARAMETER :: dt = 10E-5
    INTEGER, PARAMETER :: t = 10E4
    INTEGER, PARAMETER :: int_type = 2

    ! Internal scalars:
    DOUBLE PRECISION :: x1
    DOUBLE PRECISION :: y1
    DOUBLE PRECISION :: dx
    DOUBLE PRECISION :: dy
    DOUBLE PRECISION :: dxdt
    DOUBLE PRECISION :: dydt
    INTEGER, PARAMETER :: numchange = SIZE(argchange)

    ! Internal indices:
    INTEGER :: i
    INTEGER :: j

    ! Paths:
    CHARACTER(len=48) :: basePath
    CHARACTER(len=15) :: directoryPath
    CHARACTER(len=67) :: outputPath

    ! Outputs:
    DOUBLE PRECISION, DIMENSION(t, numchange) :: x
    DOUBLE PRECISION, DIMENSION(t, numchange) :: y

    ! Title card -------------------------------------------------------
    WRITE(6, 800)
800 FORMAT(1x, '------------------------', /, &
           1x, '  Lotka-Volterra Model  ', /, &
           1x, '      for MATH 273      ', /, &
           1x, '------------------------')

    ! Set paths --------------------------------------------------------
    ! Set base and directory paths:
    basePath = '/Users/perrindavidson/Research/uchicago/current/'
    directoryPath = 'lotka_volterra/'

    ! Set output path:
    outputPath = basePath//directoryPath//'run/'

    ! Integrate model --------------------------------------------------
    ! Open parameter:
    OPEN(UNIT=3, FILE=outputPath//'args.txt', STATUS='unknown')

    ! Loop through all variable change:
    DO j = 1, numchange, 1

        ! Set variable for interaction:
        args(ichange) = argchange(j)

        ! Initialize data:
        x1 = args(1)
        y1 = args(2)

        ! Store data:
        x(1, j) = args(1)
        y(1, j) = args(2)

        ! Loop through all time:
        DO i = 2, t, 1

            ! Integrate - Euler:
            IF (int_type .EQ. 1) THEN

                ! Integrate:
                CALL EULER(x1, y1, args, dt, x1, y1)

            ! Integrate - Runge-Kutta:
            ELSEIF (int_type .EQ. 2) THEN

                ! Integrate:
                CALL RK4(x1, y1, args, dt, x1, y1)

            ENDIF

            ! Store:
            x(i, j) = x1
            y(i, j) = y1

        ENDDO

        ! Write parameters:
        WRITE(3, *) args

    ENDDO

    ! Write out --------------------------------------------------------
    ! Open:
    OPEN(UNIT=1, FILE=outputPath//'lk_prey.txt', STATUS='unknown')
    OPEN(UNIT=2, FILE=outputPath//'lk_pred.txt', STATUS='unknown')

    ! Write phase space:
    DO i = 1, t, 1

        WRITE(1, *) x(i, :)
        WRITE(2, *) y(i, :)

    ENDDO

    ! Close:
    CLOSE(1)
    CLOSE(2)
    CLOSE(3)

    ! Ending card ------------------------------------------------------
    WRITE(6, 802)
802 FORMAT(1x, ' Done with calculations ', /, &
           1x, '------------------------')

END PROGRAM

!***********************************************************************
! lk
!-----------------------------------------------------------------------
! purpose: calculate differential in L-K equations.
!
! inputs:
! x - prey value
! y - predator value
! args = [alpha, beta, delta, gamma] - arguments as follows...
!   alpha - prey birth rate
!   beta - prey death rate from interaction with predators
!   delta - predator growth rate from feeding on prey
!   gamma - predator death rate
!
! outputs:
! dxdt - integrated prey population
! dydt - integrated predator population
!***********************************************************************
! Define program -------------------------------------------------------
SUBROUTINE LK(x, y, args, dxdt, dydt)

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Inputs:
    DOUBLE PRECISION, INTENT(IN) :: x
    DOUBLE PRECISION, INTENT(IN) :: y
    DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: args  ! [x0, y0, alpha, beta, delta, gamma]

    ! Outputs:
    DOUBLE PRECISION, INTENT(OUT) :: dxdt
    DOUBLE PRECISION, INTENT(OUT) :: dydt

    ! Calculate current values -----------------------------------------
    dxdt = (args(3) * x) - (args(4) * x * y)
    dydt = (args(5) * x * y) - (args(6) * y)

ENDSUBROUTINE LK

!***********************************************************************
! euler
!-----------------------------------------------------------------------
! purpose: integrate with euler method.
!
! inputs:
! x0 - prey value
! y0 - predator value
! args = [alpha, beta, delta, gamma] - arguments as follows...
!   alpha - prey birth rate
!   beta - prey death rate from interaction with predators
!   delta - predator growth rate from feeding on prey
!   gamma - predator death rate
! dt - time step
!
! outputs:
! x1 - integrated prey population
! y1 - integrated predator population
!***********************************************************************
! Define program -------------------------------------------------------
SUBROUTINE EULER(x0, y0, args, dt, x1, y1)

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Inputs:
    DOUBLE PRECISION, INTENT(IN) :: x0
    DOUBLE PRECISION, INTENT(IN) :: y0
    DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: args
    DOUBLE PRECISION, INTENT(IN) :: dt

    ! Internal scalars:
    DOUBLE PRECISION :: dxdt
    DOUBLE PRECISION :: dydt
    DOUBLE PRECISION :: dx
    DOUBLE PRECISION :: dy

    ! Outputs:
    DOUBLE PRECISION, INTENT(OUT) :: x1
    DOUBLE PRECISION, INTENT(OUT) :: y1

    ! Calculate current values -----------------------------------------
    ! Model:
    CALL LK(x0, y0, args, dxdt, dydt)

    ! Euler step:
    dx = dxdt * dt
    dy = dydt * dt

    ! Integrate:
    x1 = x0 + dx
    y1 = y0 + dy

ENDSUBROUTINE EULER

!***********************************************************************
! rk4
!-----------------------------------------------------------------------
! purpose: integrate with fourth order runge-kutta method.
!
! inputs:
! x - prey value
! y - predator value
! args = [alpha, beta, delta, gamma] - arguments as follows...
!   alpha - prey birth rate
!   beta - prey death rate from interaction with predators
!   delta - predator growth rate from feeding on prey
!   gamma - predator death rate
! dt - time step
!
! outputs:
! dxdt - integrated prey population
! dydt - integrated predator population
!***********************************************************************
! Define program -------------------------------------------------------
SUBROUTINE RK4(x0, y0, args, dt, x1, y1)

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Inputs:
    DOUBLE PRECISION, INTENT(IN) :: x0
    DOUBLE PRECISION, INTENT(IN) :: y0
    DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: args
    DOUBLE PRECISION, INTENT(IN) :: dt

    ! Internal scalars:
    DOUBLE PRECISION :: dxdt
    DOUBLE PRECISION :: dydt
    DOUBLE PRECISION :: k1x
    DOUBLE PRECISION :: k1y
    DOUBLE PRECISION :: k2x
    DOUBLE PRECISION :: k2y
    DOUBLE PRECISION :: k3x
    DOUBLE PRECISION :: k3y
    DOUBLE PRECISION :: k4x
    DOUBLE PRECISION :: k4y

    ! Outputs:
    DOUBLE PRECISION, INTENT(OUT) :: x1
    DOUBLE PRECISION, INTENT(OUT) :: y1

    ! Calculate current values -----------------------------------------
    ! Calculate coefficients - k1:
    CALL LK(x0, y0, args, dxdt, dydt)
    k1x = dxdt * dt
    k1y = dydt * dt

    ! Calculate coefficients - k2:
    CALL LK(x0 + (k1x / 2), y0 + (k1y / 2), args, dxdt, dydt)
    k2x = dxdt * dt
    k2y = dydt * dt

    ! Calculate coefficients - k3:
    CALL LK(x0 + (k2x / 2), y0 + (k2y / 2), args, dxdt, dydt)
    k3x = dxdt * dt
    k3y = dydt * dt

    ! Calculate coefficients - k4:
    CALL LK(x0 + k3x, y0 + k3y, args, dxdt, dydt)
    k4x = dxdt * dt
    k4y = dydt * dt

    ! Estimate values:
    x1 = x0 + ((k1x + (2 * k2x) + (2 * k3x) + k4x) / 6)
    y1 = y0 + ((k1y + (2 * k2y) + (2 * k3y) + k4y) / 6)

ENDSUBROUTINE RK4

!***********************************************************************
