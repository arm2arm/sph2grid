function  smartlog, data


Ndims= size(data,/DIMENSIONS)
;if(ndims lt 1) then return ,0
vollog=reform(data, n_elements(data))

iro=where(vollog gt 0, COMPLEMENT=niro)

vollog(iro)=alog10(vollog(iro))
minvol=min(vollog(iro))
if(n_elements(niro) gt 1)then begin
    vollog(niro)=minvol
endif
vollog=reform(vollog,Ndims)
return, vollog
end

device,retain=2,decomposed=0

file='../src/LGR2Mpc_sfrh_143_grid'
file='../_sigdiv_RZA3D_1024_500_s2_zinit120_011_grid'
;file='C:\arm2arm\DATA\MODEL7\MODELS\MODEL7\RUNG2\SNAPS\snap_gal_0.25_sfr_030_grid'
openr, 1, file
type=0l
GRID=0l
XC=fltarr(4)
loadct, 4
readu, 1, type
readu, 1, XC
readu, 1, GRID
print, type, XC, GRID
data=fltarr(GRID, GRID, GRID)
readu, 1, data
close, 1
RC=XC[3]
;stop
RCobj=1
scR=30.0
rangx=[XC[0]-scR, XC[0]+scR]
rangy=[XC[1]-scR, XC[1]+scR]
scR=1.0
ngR=min([RCobj*scR/(RC*2.0)*GRID, GRID/2.0-1])
Print, "slice with is :", 2*RCobj*scR, " Kpc/h (to control it change scR is in units RCobj)"
slice=smartlog(total(data(*,*,GRID/2-ngR:GRID/2+ngR),3, /double))
xv=findgen(GRID)/GRID*2*RC-RC+XC[0]
yv=findgen(GRID)/GRID*2*RC-RC+XC[1]
zv=findgen(GRID)/GRID*2*RC-RC+XC[2]


;contour, slice,xv,yv,$
;  /fill, nlevels=255, /xs, /ys, /iso,$;xrange=rangx, yrange=rangy,$
;MIN_VALUE=-12, MAX_VALUE=-0.6
window, xsize=512, ysize=512
loadct, 0
contour, slice,xv,yv, /fill, nlevels=255, MIN_VALUE=-1.5, /xs, /ys, /iso, /nodata
loadct, 13
contour, slice,xv,yv, /fill, nlevels=255, /overplot
loadct, 0
contour, congrid(slice, grid, grid, /interp),xv,yv,$
  nlevels=5,/overplot
;Tvcircle,  RCobj, XC[0], XC[1],  255, /data
write_png, 'c:/arm2arm/progs/sph2grid/isol_gal.png', tvrd(true=1)
end
