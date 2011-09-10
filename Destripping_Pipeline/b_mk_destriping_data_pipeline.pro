pro b_mk_destriping_data_pipeline, dirl

;***********************************************************************
;NAME: b_mk_destriping_data_pipeline
;***********************************************************************
;PURPOSE: Save Level 2 data into fits format for Andrea's destrippping code
;***********************************************************************
;AUTHOR: Brian Williams - CalPoly - December 2010
;***********************************************************************
;CALLING: b_mk_destriping_data_pipeline
;***********************************************************************
;INPUTS:
;Level 2 data from Data2 
;***********************************************************************
;OPTIONAL INPUTS:
;
;***********************************************************************
;OUTPUT:
;2 fits files one with data the other with pointing in the form of Paralactic angle
;        and pixel number for nside equals 256
;***********************************************************************
;EXAMPLE:
;  b_mk_destriping_data_pipeline,['20080826']
;***********************************************************************
;This code uses tools from:
;The IDL Astronomy Users Library. http://idlastro.gsfc.nasa.gov/
;The HEALPix software package. http://healpix.jpl.nasa.gov/
;***********************************************************************

compile_opt idl2

print,systime()
;!quiet=1
Latitude=37.583375
Longitude=118.237095

RTD=180./!pi
dtr=!pi/180.
miss = -1.63750e30

ndirs=n_elements(dirl)
channum = [1,2,3,6]

nside=256.0
npix=12.0*nside*nside

j=0

l2path='/smbmount/data2/BMachine/Data/Level2New/'
dspath_short='/smbmount/data2/BMachine/Data/Destripe_Data/'

for i=0, ndirs-1 do begin
   dspath=strcompress(dspath_short+dirl[i]+'/',/remove_all)
   file_mkdir, dspath
  for j=0,3 do begin
    dat=dirl[i]
    year=float(strmid(dat,0,4))
    month=float(strmid(dat,4,2))
    day=float(strmid(dat,6,2))
    ch=channum[j]
    lvl2d=b_read_lvl2_chan(dirl[i],channel=ch)
    timed=lvl2d.time
    az=lvl2d.azimuth
    el=lvl2d.elevation
    gooddata=where(lvl2d.flag eq 0.0)
    
    Temp=reform(lvl2d.data[0,*])*lvl2d.cals[0]
    
    Q=reform(lvl2d.data[1,*]-mean(lvl2d.data[1,gooddata]))*lvl2d.cals[1]
    U=reform(lvl2d.data[2,*]-mean(lvl2d.data[2,gooddata]))*lvl2d.cals[1]
   
    para_angle= get_parallactic_angle(month, day, year, timed, az, el, latitude, longitude, altitude)
    lst=sidetime(month,day,year,timed,longitude)
    Azel2radec,az,el,ra,dec,latitude,lst
    ra=ra mod(24)
    rad=ra*360./24.
    phi = rad *( (2*!DPI)/360. )
    theta = (90. - dec) * (!DPI / 180.)
    ang2pix_ring,nside,theta,phi,ipring
    output_file_data=strcompress(dspath+'ch'+string(ch)+'.fits',/remove_all)
    output_file_pointing=strcompress(dspath+'ch'+string(ch)+'_pix.fits',/remove_all)
    baddata=where(lvl2d.flag gt 0.0)
    ipring[baddata]=npix
    
    
      fxhmake,header,/initialize,/extend,/date
      fxwrite,output_file_data,header
      fxbhmake,header,n_elements(timed),'Data'

      fxbaddcol,1,header,timed[0],'Time', 'Time', TUNIT    = 'Fractional Hours'
      fxbaddcol,2,header,Temp[0],'Temp', 'Temp', TUNIT    = 'Kelvin'
      fxbaddcol,3,header,Q[0],'Q', 'Q', TUNIT    = 'Kelvin'
      fxbaddcol,4,header,U[0],'U', 'U', TUNIT    = 'Kelvin'
     
;      (4) Write extension header to FITS file
      fxbcreate,unit,output_file_data,header

;      (5) Use FXBWRITM to write all data to the extension in a single call
      fxbwritm,unit,['Time','Temp','Q','U'],timed, Temp, Q, U

      FXBFINISH, UNIT
     
;************************************ Write pointing fits**************************************
      fxhmake,headerp,/initialize,/extend,/date
      fxwrite,output_file_pointing,headerp
      fxbhmake,headerp,n_elements(timed),'Data'
      
      fxbaddcol,1,headerp,timed[0],'Time', 'Time', TUNIT    = 'Fractional Hours'
      fxbaddcol,2,headerp,ipring[0],'Pixel', 'Pixel Number'
      fxbaddcol,3,headerp,para_angle[0],'Alpha','Paralactic Angle', TUNIT    = 'Degrees'
;      (4) Write extension header to FITS file
      fxbcreate,unitp,output_file_pointing,headerp
;      (5) Use FXBWRITM to write all data to the extension in a single call
      fxbwritm,unitp,['Time','Pixel','Alpha'],timed,ipring,Para_angle
      FXBFINISH, UNITp

;ftab_help, output_file
  endfor
endfor
print,systime()
end

