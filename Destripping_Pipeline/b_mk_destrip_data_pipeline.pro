pro b_mk_destrip_data_pipeline,dirlist=dirlist

;writehistory='First try 4-21-2009 by Brian Williams'
;writehistory='Second try 8-9-2009 by Brian Williams'
;writehistory='Third  try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
;writehistory='Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
  
; Sample Call
;     b_mk_destrip_data_Pipeline, dirlist=['20080826']

;compile_opt idl2
starttime=systime()

if not(keyword_set(dirlist)) then dirlist=['20080826']

;White Mountain 
Latitude=37.583375
Longitude=118.237095
altitude=3700

RTD=180./!pi
dtr=!pi/180.

;********* Radish directories
;lvl1path = '/smbmount/data2/BMachine/Data/Level1New/'
;lvl2path = '/smbmount/data2/BMachine/Data/Level2New/'
;restore, '/smbmount/data2/BMachine/Data/Calibration/CalCorrectionStructure.sav'

;********** Laptop directories
lvl1path='/Users/brianw/B_Machine/Level1_new/'
lvl2path='/Users/brianw/B_Machine/LevelDestrip/'
restore, '/Users/brianw/B_Machine/calibration/CalCorrectionStructure.sav'
ndir=n_elements(dirlist)

hr=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']

for i=0,ndir-1 do begin

	file_mkdir,strcompress(lvl2path+dirlist[i],/remove_all)
	  dat=dirlist[i]
    year=float(strmid(dat,0,4))
    month=float(strmid(dat,4,2))
    day=float(strmid(dat,6,2))

	;******************************* Getting Calibrations******************************
	calnum= where(strmid(calcorrections.dates[*],24,8) eq strcompress(dirlist[i]),count)
	if count gt 0 then begin
		calsACT=calcorrections.cals[*,calnum]
		polcor=calcorrections.cals[*,calnum]/calcorrections.Cals0807[*]
		CalsACP=[-1.6572*polcor[0],-1.4822*polcor[1],-2.790*polcor[2],-1.4164*polcor[3]]
	endif else begin
		CalsACT=[0.9203,0.8319,1.378,0.7329]
		CalsACP=[-1.6572,-1.4822,-2.790,-1.4164]
	endelse
	;********************************************************************************************

	;************************************************Get pointing offsets*******************************
	pointoff=strcompress('/smbmount/data2/BMachine/Data/pointing_offsets/'+dirlist[i]+'.txt',/remove_all)
	;pointoff=strcompress('/Users/brianw/B_Machine/pointing_offsets/'+dirlist[i]+'.txt',/remove_all)
	stuff=file_info(pointoff)
	if stuff.exists then begin
		pointingoffs=findgen(2,16)
		pointingoff=findgen(2,4)
 		openr,2,pointoff
		readf,2,pointingoffs
		close,2
		pointingoff[0,*]=[pointingoffs[0,1],pointingoffs[0,2],pointingoffs[0,3],pointingoffs[0,6]]
		pointingoff[1,*]=[pointingoffs[1,1],pointingoffs[1,2],pointingoffs[1,3],pointingoffs[1,6]]
	endif else begin
		pointingoff=findgen(2,4)
		pointingoff[0,*]=[0.0,0.0,0.0,0.0]
		pointingoff[1,*]=[0.0,0.0,0.0,0.0]
	endelse
	;****************************************************************************************************


    for hour=0,23 do begin

       ;lvl1filelist=file_search(lvl1path+dirlist[i]+'\'+hr[hour]+'*.fts')
       lvl1filelist=file_search(lvl1path+dirlist[i]+'/'+hr[hour]+'*.fts')
       if (n_elements(lvl1filelist) gt 1 ) then begin
	;***************************************************Data Cuts******************************************************
			lvl1data=b_read_lvl1_fits_pipeline(files=lvl1filelist,keep=[0,1,2,3,4,5,6])
			lvl1data=B_pointing_interp(lvl1data)
			flags=lvl2_Data_selection(lvl1data)
			time=reform(lvl1data.time)

			Status=lvl1data.status


;*************************************************Create data variables for chan 1 fits******************************

		data_chan1=reform(lvl1data.demoddata[1,*,*])
		T_chan1=reform(data_chan1[0,*])*calsACT[0]
    U_chan1=reform(data_chan1[1,*])*calsACP[0]
    Q_chan1=reform(data_chan1[2,*])*calsACP[0]

		Elevation_chan1=lvl1data.elevation[*]*rtd+ pointingoff[1,0]
		Azimuth_chan1=lvl1data.Azimuth[*]*rtd + pointingoff[0,0]/cos(Elevation_chan1*dtr)
		azimuth_chan1=azimuth_chan1 mod 360.0
		aztest=where(azimuth_chan1 lt 0,countaz)
			if countaz gt 0 then begin
				azimuth_chan1[aztest]=azimuth_chan1[aztest]+360.0
			endif
			
      para_angle=azimuth_chan1[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,azimuth_chan1,elevation_chan1,ra_chan1,dec_chan1,latitude,lst
      ra_chan1=ra_chan1 mod(24)
      rad_chan1=ra_chan1*360./24.
      phi_chan1 = rad_chan1 *( (2*!DPI)/360. )          ;Right Acsension
      theta_chan1 = (90. - dec_chan1) * (!DPI / 180.)   ;Declination
      Psi_chan1= get_parallactic_angle(month, day, year, time, azimuth_chan1, elevation_chan1, latitude, longitude, altitude)
	
      Revs_chan1=lvl1data.revolution
      time_chan1=lvl1data.time

;********************************************************************************************************************


 ;*************************************************Create data variables for chan 2 fits******************************
    data_chan2=reform(lvl1data.demoddata[2,*,*])
    T_chan2=reform(data_chan2[0,*])*calsACT[1]
    U_chan2=reform(data_chan2[1,*])*calsACP[1]
    Q_chan2=reform(data_chan2[2,*])*calsACP[1]

    Elevation_chan2=lvl1data.elevation[*]*rtd+ pointingoff[1,0]
    Azimuth_chan2=lvl1data.Azimuth[*]*rtd + pointingoff[0,0]/cos(Elevation_chan2*dtr)
    azimuth_chan2=azimuth_chan2 mod 360.0
    aztest=where(azimuth_chan2 lt 0,countaz)
      if countaz gt 0 then begin
        azimuth_chan2[aztest]=azimuth_chan2[aztest]+360.0
      endif
      
      para_angle=azimuth_chan2[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,azimuth_chan2,elevation_chan2,ra_chan2,dec_chan2,latitude,lst
      ra_chan2=ra_chan2 mod(24)
      rad_chan2=ra_chan2*360./24.
      phi_chan2 = rad_chan2 *( (2*!DPI)/360. )          ;Right Acsension
      theta_chan2 = (90. - dec_chan2) * (!DPI / 180.)   ;Declination
      Psi_chan2= get_parallactic_angle(month, day, year, time, azimuth_chan2, elevation_chan2, latitude, longitude, altitude)
  
      Revs_chan2=lvl1data.revolution
      time_chan2=lvl1data.time
;********************************************************************************************************************


;*************************************************Create data variables for chan 3 fits******************************
    data_chan3=reform(lvl1data.demoddata[3,*,*])
    T_chan3=reform(data_chan3[0,*])*calsACT[2]
    U_chan3=reform(data_chan3[1,*])*calsACP[2]
    Q_chan3=reform(data_chan3[2,*])*calsACP[2]

    Elevation_chan3=lvl1data.elevation[*]*rtd+ pointingoff[1,0]
    Azimuth_chan3=lvl1data.Azimuth[*]*rtd + pointingoff[0,0]/cos(Elevation_chan3*dtr)
    azimuth_chan3=azimuth_chan3 mod 360.0
    aztest=where(azimuth_chan3 lt 0,countaz)
      if countaz gt 0 then begin
        azimuth_chan3[aztest]=azimuth_chan3[aztest]+360.0
      endif
      
      para_angle=azimuth_chan3[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,azimuth_chan3,elevation_chan3,ra_chan3,dec_chan3,latitude,lst
      ra_chan3=ra_chan3 mod(24)
      rad_chan3=ra_chan3*360./24.
      phi_chan3 = rad_chan2 *( (2*!DPI)/360. )          ;Right Acsension
      theta_chan3 = (90. - dec_chan3) * (!DPI / 180.)   ;Declination
      Psi_chan3= get_parallactic_angle(month, day, year, time, azimuth_chan3, elevation_chan3, latitude, longitude, altitude)
  
      Revs_chan3=lvl1data.revolution
      time_chan3=lvl1data.time
;********************************************************************************************************************

;*************************************************Create data variables for chan 6 fits******************************
    data_chan6=reform(lvl1data.demoddata[6,*,*])
    T_chan6=reform(data_chan6[0,*])*calsACT[3]
    U_chan6=reform(data_chan6[1,*])*calsACP[3]
    Q_chan6=reform(data_chan6[2,*])*calsACP[3]

    Elevation_chan6=lvl1data.elevation[*]*rtd+ pointingoff[1,0]
    Azimuth_chan6=lvl1data.Azimuth[*]*rtd + pointingoff[0,0]/cos(Elevation_chan6*dtr)
    azimuth_chan6=azimuth_chan6 mod 360.0
    aztest=where(azimuth_chan6 lt 0,countaz)
      if countaz gt 0 then begin
        azimuth_chan6[aztest]=azimuth_chan6[aztest]+360.0
      endif
      
      para_angle=azimuth_chan6[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,azimuth_chan6,elevation_chan6,ra_chan6,dec_chan6,latitude,lst
      ra_chan6=ra_chan6 mod(24)
      rad_chan6=ra_chan6*360./24.
      phi_chan6 = rad_chan6 *( (2*!DPI)/360. )          ;Right Acsension
      theta_chan6 = (90. - dec_chan6) * (!DPI / 180.)   ;Declination
      Psi_chan6= get_parallactic_angle(month, day, year, time, azimuth_chan6, elevation_chan6, latitude, longitude, altitude)
  
      Revs_chan6=lvl1data.revolution
      time_chan6=lvl1data.time
;********************************************************************************************************************


;*******************************************Write Fits file for chan 1***********************************
		unit=2
		;l2file_chan1=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan1.fits',/Remove_all)
		l2file_chan1=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan1.fits',/Remove_all)
         	fxhmake,header,/initialize,/extend,/date
        	fxwrite,l2file_chan1,header
         	fxbhmake,header,n_elements(azimuth_chan1),'Destriping','Destriping'
         	fxbaddcol,1,header,time_chan1[0],'Time'
         	fxbaddcol,2,header,Phi_chan1[0],'Phi'         ;RA
         	fxbaddcol,3,header,Theta_chan1[0],'Theta'     ;Dec
        	fxbaddcol,4,header,Psi_chan1[0],'Psi'             ;Paralactic angle
        	fxbaddcol,5,header,flags.chan1[0],'Flag'
         	fxbaddcol,6,header,  T_chan1[0],'Temp'
         	fxbaddcol,7,header,  Q_chan1[0],'Q'
         	fxbaddcol,8,header,  U_chan1[0],'U'
         	
         	sxaddpar,header, '          ','          '
         	sxaddpar,header,'Band', '38-45Ghz'
         	sxaddpar,header,'DataUnits', 'mK'
         	sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
         	sxaddpar,header,'Pointunits', 'Degrees'
         	sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'History','Demodutlated data for Andrea de-striping code'
          sxaddpar,header,'FileDate', dirlist[i]

       		fxbcreate,unit,l2file_chan1,header
         	fxbwritm,unit,['Time','Phi','Theta','Psi','Flag', 'Temp','Q','U'],Time_chan1,Phi_chan1,Theta_chan1,Psi_chan1,flags.chan1,T_chan1,Q_chan1,U_chan1
         	;fxbwritm,unit,['Data'],data_chan1
        	fxbfinish,unit
  ;********************************************************************************************************************

;*******************************************Write Fits file for chan 2***********************************
    ;l2file_chan2=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan2.fits',/Remove_all)
    l2file_chan2=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan2.fits',/Remove_all)
          fxhmake,header,/initialize,/extend,/date
          fxwrite,l2file_chan2,header
          fxbhmake,header,n_elements(azimuth_chan2),'Destriping','Destriping'
          fxbaddcol,1,header,time_chan2[0],'Time'
          fxbaddcol,2,header,Phi_chan2[0],'Phi'         ;RA
          fxbaddcol,3,header,Theta_chan2[0],'Theta'     ;Dec
          fxbaddcol,4,header,Psi_chan2[0],'Psi'             ;Paralactic angle
          fxbaddcol,5,header,flags.chan2[0],'Flag'
          fxbaddcol,6,header,  T_chan2[0],'Temp'
          fxbaddcol,7,header,  Q_chan2[0],'Q'
          fxbaddcol,8,header,  U_chan2[0],'U'
          
          sxaddpar,header, '          ','          '
          sxaddpar,header,'Band', '38-45Ghz'
          sxaddpar,header,'DataUnits', 'mK'
          sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
          sxaddpar,header,'Pointunits', 'Degrees'
          sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'History','Demodutlated data for Andrea de-striping code'
          sxaddpar,header,'FileDate', dirlist[i]

          fxbcreate,unit,l2file_chan2,header
          fxbwritm,unit,['Time','Phi','Theta','Psi','Flag', 'Temp','Q','U'],Time_chan2,Phi_chan2,Theta_chan2,Psi_chan2,flags.chan2,T_chan2,Q_chan2,U_chan2
          ;fxbwritm,unit,['Data'],data_chan2
          fxbfinish,unit
  ;********************************************************************************************************************

 ;*******************************************Write Fits file for chan 3***********************************
    ;l2file_chan3=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan3.fits',/Remove_all)
    l2file_chan3=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan3.fits',/Remove_all)
          fxhmake,header,/initialize,/extend,/date
          fxwrite,l2file_chan3,header
          fxbhmake,header,n_elements(azimuth_chan3),'Destriping','Destriping'
          fxbaddcol,1,header,time_chan3[0],'Time'
          fxbaddcol,2,header,Phi_chan3[0],'Phi'         ;RA
          fxbaddcol,3,header,Theta_chan3[0],'Theta'     ;Dec
          fxbaddcol,4,header,Psi_chan3[0],'Psi'             ;Paralactic angle
          fxbaddcol,5,header,flags.chan3[0],'Flag'
          fxbaddcol,6,header,  T_chan3[0],'Temp'
          fxbaddcol,7,header,  Q_chan3[0],'Q'
          fxbaddcol,8,header,  U_chan3[0],'U'
          
          sxaddpar,header, '          ','          '
          sxaddpar,header,'Band', '38-45Ghz'
          sxaddpar,header,'DataUnits', 'mK'
          sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
          sxaddpar,header,'Pointunits', 'Degrees'
          sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'History','Demodutlated data for Andrea de-striping code'
          sxaddpar,header,'FileDate', dirlist[i]

          fxbcreate,unit,l2file_chan3,header
          fxbwritm,unit,['Time','Phi','Theta','Psi','Flag', 'Temp','Q','U'],Time_chan3,Phi_chan3,Theta_chan3,Psi_chan3,flags.chan3,T_chan3,Q_chan3,U_chan3
          ;fxbwritm,unit,['Data'],data_chan3
          fxbfinish,unit
  ;********************************************************************************************************************


   ;*******************************************Write Fits file for chan1***********************************
    ;l2file_chan6=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan6.fits',/Remove_all)
    l2file_chan6=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan6.fits',/Remove_all)
          fxhmake,header,/initialize,/extend,/date
          fxwrite,l2file_chan6,header
          fxbhmake,header,n_elements(azimuth_chan6),'Destriping','Destriping'
          fxbaddcol,1,header,time_chan6[0],'Time'
          fxbaddcol,2,header,Phi_chan6[0],'Phi'         ;RA
          fxbaddcol,3,header,Theta_chan6[0],'Theta'     ;Dec
          fxbaddcol,4,header,Psi_chan6[0],'Psi'             ;Paralactic angle
          fxbaddcol,5,header,flags.chan6[0],'Flag'
          fxbaddcol,6,header,  T_chan6[0,0],'Temp'
          fxbaddcol,7,header,  Q_chan6[0,0],'Q'
          fxbaddcol,8,header,  U_chan6[0,0],'U'
          
          sxaddpar,header, '          ','          '
          sxaddpar,header,'Band', '38-45Ghz'
          sxaddpar,header,'DataUnits', 'mK'
          sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
          sxaddpar,header,'Pointunits', 'Degrees'
          sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'History','Demodutlated data for Andrea de-striping code'
          sxaddpar,header,'FileDate', dirlist[i]

          fxbcreate,unit,l2file_chan6,header
          fxbwritm,unit,['Time','Phi','Theta','Psi','Flag', 'Temp','Q','U'],Time_chan6,Phi_chan6,Theta_chan6,Psi_chan6,flags.chan6,T_chan6,Q_chan6,U_chan6
          ;fxbwritm,unit,['Data'],data_chan6
          fxbfinish,unit
  ;********************************************************************************************************************

       endif 	 ;end if for filedata check
    endfor   ; end of hour loop
endfor    ; end of directory loop

print,'Starttime is ' + starttime
print,systime()

end



