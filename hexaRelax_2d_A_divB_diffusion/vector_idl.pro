nx = 256
nz = 256

time = 1.0d0
bbx = dblarr(nx+1, nz+2)
bby = dblarr(nx+2, nz+2)
bbz = dblarr(nx+2, nz+1)

OPENR, 10, "run1/relax_00001", /F77_UNFORMATTED

READU, 10, time
PRINT, 10, time
READU, 10, bbx
READU, 10, bby
READU, 10, bbz

CLOSE, 10

END
