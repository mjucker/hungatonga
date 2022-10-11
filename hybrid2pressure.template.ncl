; ****************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ****************************************************

begin

; *********************
; Reading Data
; *********************
diro = ""  ; Directory
fin  = "INFILE"     ; Input file
fout = "OUTFILE"    ; Output file
vname = "VARNAME"

data1 = addfile(diro + fin,"r")
data2 = addfile("GRIDFILE","r")
data3 = addfile(diro + "DATA3FILE","r")

VAR1 = data1->$vname$
VAR1@missing_value= -999
VAR1@_FillValue = -999
P0 = 0.01*data2->P0       ; reference pressure [mb]
PS = data3->PS
hyam = data2->hyam  
hybm = data2->hybm
lat = data1->lat    
lon = data1->lon
time = data1->time
ts = data3->TS
phis = data3->PHIS

; desired output levels [hPa/mb] 
;pnew = (/  10, 20, 30, 50, 70,100,150,200,250 \
;           , 300,400,500,600,700,800,850,900,925,950,1000 /)  
pnew = (/ 1, 2, 3, 5, 7, \
    10, 15, 20, 30, 40, 50, 60, 75, 93, \
    115, 140, 170, 200, 250, 300, 350, 400, 450, \
    500, 550, 600, 700, 750, 800, 850, 900, 950 /)


;**************************************************
; Interpolating from CAM3 hybrid to pressure levels
;**************************************************
; type of interpolation: 1 = linear, 2 = log, 3 = loglog
; is extrapolation desired:   extrap = True
; the 7th argument is not used, and so is set to 1.

VARnew = vinth2p_ecmwf(VAR1,hyam,hybm,pnew,PS,1,P0,1,True,0,ts,phis)


; ***********************************************
; Saving variable into netcdf file 
; ***********************************************

system("/bin/rm -f " + diro + fout)    ; remove if exists
ncdf     = addfile( diro + fout ,"c")  ; open output netCDF file

; make time and UNLIMITED dimension      ; recommended  for most applications
filedimdef(ncdf,"time",-1,True) 

; output variables directly
VARnew!0    = "time"
VARnew!1    = "lev"
VARnew!2    = "lat"                 ; assign named dimensions
VARnew!3    = "lon"
VARnew&time =  time                     
VARnew&lev  =  pnew                     
VARnew&lat  =  lat                     
VARnew&lon  =  lon 

VARnew@long_name = "Geopotential Height"
VARnew@units     = "m"

ncdf->$vname$ = VARnew


end

