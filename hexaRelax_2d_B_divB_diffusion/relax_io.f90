MODULE io

  USE var_global
  USE grid

  IMPLICIT NONE

CONTAINS

  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n

    WRITE (filename, FMT = '(a, i5.5)')  output_file, n

    OPEN(UNIT = 42, &
         FILE = FILENAME, &
         FORM = 'UNFORMATTED', &
         STATUS = 'UNKNOWN')
    WRITE (42) t
    WRITE (42) bbx
    WRITE (42) bby
    WRITE (42) bbz
    CLOSE (42)

  END SUBROUTINE writedata

  SUBROUTINE output_diag

    REAL(num), DIMENSION(nx, 2:nz) :: b_c, bm_c ! bb and bbm at cell centre
    REAL(num), DIMENSION(nx) :: exbase, eybase, bxbase, bybase, poynting
    REAL(num), DIMENSION(nx+1, nz+1) :: vv ! Velocity squared
    REAL(num), DIMENSION(nx, 2:nz) :: diss
    REAL(num), DIMENSION(nx+1, nz) :: bx_err
    REAL(num), DIMENSION(nx, nz) :: by_err
    REAL(num), DIMENSION(nx, nz+1) :: bz_err
    REAL(num), DIMENSION(nx, 2:nz) :: err_array
    REAL(num), DIMENSION(nx, nz) :: divb, divb_norm, bbm_c
    REAL(num), DIMENSION(nx) :: divb_base
    REAL(num) :: mag_eng, sz, q, err, min_divb, max_divb, mean_divb
    REAL(num) :: max_divb_norm, mean_divb_norm, divb_flux, divb_diss

    ! Calculate magnetic energy
    b_c = 0.25_num * (bb(1:nx  , 2:nz  ) + bb(2:nx+1, 2:nz  ) + &
                      bb(1:nx  , 3:nz+1) + bb(2:nx+1, 3:nz+1))
    bm_c = 0.25_num * (bbm(1:nx  , 2:nz  ) + bbm(2:nx+1, 2:nz  ) + &
                       bbm(1:nx  , 3:nz+1) + bbm(2:nx+1, 3:nz+1))
    mag_eng = 0.5_num * SUM(b_c) * delx * delz

    !poynting flux (Ex * By - Ey * Bx)
    bxbase = 0.25_num * (bbx(1:nx, 1) + bbx(2:nx+1, 1) + &
                         bbx(1:nx, 2) + bbx(2:nx+1, 2))
    bybase = 0.25_num * (bby(1:nx, 1) + bby(1:nx, 1) + &
                         bby(1:nx, 2) + bby(1:nx, 2))
    exbase = eex(:, 2)
    eybase = 0.5_num * (eey(1:nx, 2) + eey(2:nx+1, 2))
    poynting = (exbase * bybase - eybase * bxbase) * delx
    sz = SUM(poynting)

    ! Dissipation Q = B^2 * nu * v^2
    vv = vx * vx + vy * vy + vz * vz
    diss = 0.25_num * (vv(1:nx  , 2:nz  ) + vv(2:nx+1, 2:nz  ) + &
                       vv(1:nx  , 3:nz+1) + vv(2:nx+1, 3:nz+1))
    diss = diss / frc
    diss = diss * bm_c
    q = SUM(diss) * delx * delz

    ! Analytic solutions vs. numerical solution
    bx_err = (bbx(:, 1:nz) - bx_a) ** 2.0_num
    by_err = (bby(1:nx, 1:nz) - by_a) ** 2.0_num
    bz_err = (bbz(1:nx, :) - bz_a) ** 2.0_num
    err_array = 0.5_num * (bx_err(1:nx, 2:nz) + bx_err(2:nx+1, 2:nz)) &
              + by_err(:, 2:nz) &
              + 0.5_num * (bz_err(:, 2:nz) + bz_err(:, 3:nz+1))
    err = SUM(err_array) * delx * delz

    ! divb
    divb = (bbx(2:nx+1, 1:nz) - bbx(1:nx, 1:nz)) / delx &
         + (bbz(1:nx, 2:nz+1) - bbz(1:nx, 1:nz)) / delz
    min_divb = MINVAL(ABS(divb))
    max_divb = MAXVAL(ABS(divb))
    mean_divb = SUM(ABS(divb)) / REAL(nx * nz, num)

    ! divb norm
    bbm_c = 0.25_num * (bbm(1:nx  , 1:nz  ) + bbm(2:nx+1, 1:nz  ) + &
                        bbm(1:nx  , 2:nz+1) + bbm(2:nx+1, 2:nz+1))
    bbm_c = SQRT(bbm_c)
    divb_norm = divb * delx / (bbm_c * 6.0_num)
    max_divb_norm = MAXVAL(ABS(divb_norm))
    mean_divb_norm = SUM(ABS(divb_norm)) / REAL(nx * nz, num)

    ! divb flux
    divb_base = 0.5_num * (divb(:, 1) + divb(:, 2))
    divb_flux = etad * SUM(divb_base * bbz(:, 2)) * delx

    ! divB dissipation
    divb_diss = etad * SUM(divB(:, 2:nz) * divB(:, 2:nz)) * delx * delz

    ! Current-field angle
    

    WRITE(50, *) t, &
                 dt, &
                 mag_eng, &
                 sz, &
                 q, &
                 err, &
                 min_divb, &
                 max_divb, &
                 mean_divb, &
                 max_divb_norm, &
                 mean_divb_norm, &
                 divb_flux, &
                 divb_diss, &
                 iter_no

  END SUBROUTINE output_diag

END MODULE io
