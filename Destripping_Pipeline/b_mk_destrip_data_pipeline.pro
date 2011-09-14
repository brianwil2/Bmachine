pro b_mk_destrip_data_pipeline,dirlist

;writehistory='First try 4-21-2009 by Brian Williams'
;writehistory='Second try 8-9-2009 by Brian Williams'
;writehistory='Third  try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
;writehistory='Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
  
; Sample Call
;     b_mk_lvl2_data_Pipeline_c, ['20080826']

;compile_opt idl2
starttime=systime()

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
lvl1path='/Users/brianw/B_Machine/Level1/'
lvl2path='/Users/brianw/B_Machine/LevelDestrip/'
restore, '/Users/brianw/B_Machine/calibration/CalCorrectionStructure.sav'
ndir=n_elements(dirlist)

hr=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']

for i=0,ndir-1 do begin

	file_mkdir,strcompress(lvl2path+dirlist[i],/remove_all)
	  dat=dirl[i]
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


;*************************************************Create data variables for chan1 fits******************************

		data_chan1=reform(lvl1data.demoddata[1,*,*])

		Elevation_chan1=lvl1data.elevation[*]*rtd+ pointingoff[1,0]
		Azimuth_chan1=lvl1data.Azimuth[*]*rtd + pointingoff[0,0]/cos(Elevation_chan1*dtr)
		azimuth_chan1=azimuth_chan1 mod 360.0
		aztest=where(azimuth_chan1 lt 0,countaz)
			if countaz gt 0 then begin
				azimuth_chan1[aztest]=azimuth_chan1[aztest]+360.0
			endif
			
      para_angle=azimuth_chan1[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,azimuth_chan1,elevations_chan1,ra_chan1,dec_chan1,latitude,lst
      ra=ra mod(24)
      rad=ra*360./24.
      phi = rad *( (2*!DPI)/360. )          ;Right Acsension
      theta = (90. - dec) * (!DPI / 180.)   ;Declination
      para_angle= get_parallactic_angle(month, day, year, time, azimuth_chan1, elevation_chan1, latitude, longitude, altitude)
	
      Revs_chan1=lvl1data.revolution
      time_chan1=lvl1data.time

;********************************************************************************************************************


 ;*************************************************Create data variables for chan2 fits******************************

		data_chan2=reform(lvl1data.demoddata[2,*,*])

		Elevation_chan2=lvl1data.elevation[*]*rtd + pointingoff[1,1]
		Azimuth_chan2=lvl1data.Azimuth[*]*rtd + pointingoff[0,1]/cos(Elevation_chan2*dtr)
		azimuth_chan2=azimuth_chan2 mod 360.0
		aztest=where(azimuth_chan2 lt 0,countaz)
			if countaz gt 0 then begin
				azimuth_chan2[aztest]=azimuth_chan2[aztest]+360.0
			endif
		Revs_chan2=lvl1data.revolution
		time_chan2=lvl1data.time


;********************************************************************************************************************


 ;*************************************************Create data variables for chan3 fits******************************

        	data_chan3=reform(lvl1data.demoddata[3,*,*])

		Elevation_chan3=lvl1data.elevation[*]*rtd + pointingoff[1,2]
		Azimuth_chan3=lvl1data.Azimuth[*]*rtd + pointingoff[0,2]/cos(Elevation_chan3*dtr)
		azimuth_chan3=azimuth_chan3 mod 360.0
		aztest=where(azimuth_chan1 lt 0,countaz)
			if countaz gt 0 then begin
				azimuth_chan3[aztest]=azimuth_chan3[aztest]+360.0
			endif
		Revs_chan3=lvl1data.revolution
		time_chan3=lvl1data.time

;********************************************************************************************************************

 ;*************************************************Create data variables for chan6 fits******************************

		data_chan6=reform(lvl1data.demoddata[6,*,*])

		Elevation_chan6=(lvl1data.elevation[*] *rtd)+ pointingoff[1,3]
		Azimuth_chan6=lvl1data.Azimuth[*]*rtd + (pointingoff[0,3])/cos(Elevation_chan6*dtr)
		azimuth_chan6=azimuth_chan6 mod 360.0
		aztest=where(azimuth_chan6 lt 0,countaz)
			if countaz gt 0 then begin
				azimuth_chan6[aztest]=azimuth_chan6[aztest]+360.0
			endif
		Revs_chan6=lvl1data.revolution
		time_chan6=lvl1data.time

;********************************************************************************************************************


;*******************************************Write Fits file for chan1***********************************
		unit=2
		;l2file_chan1=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan1.fits',/Remove_all)
		l2file_chan1=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan1.fits',/Remove_all)
         	fxhmake,header,/initialize,/extend,/date
        	fxwrite,l2file_chan1,header
         	fxbhmake,header,n_elements(azimuth_chan1),'CutData','CutData'
         	fxbaddcol,1,header,time_chan1[0],'Time'
         	fxbaddcol,2,header,azimuth_chan1[0],'Azimuth'
         	fxbaddcol,3,header,elevation_chan1[0],'Elevation'
        	fxbaddcol,4,header,revs_chan1[0],'Revolution'
        	fxbaddcol,5,header,flags.chan1[0],'Flag'
        	fxbaddcol,6,header, status[0],'Status'
         	fxbaddcol,7,header,  data_chan1[*,0],'Data'

         	sxaddpar,header, '          ','          '
         	sxaddpar,header,'Band', '38-45Ghz'
         	sxaddpar,header,'DataUnits', 'Volts'
         	sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
         	sxaddpar,header,'Pointunits', 'Degrees'
         	sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'CalT',string( calsacT[0])
          sxaddpar,header,'CalP', string(calsacP[0])
          sxaddpar,header,'History','Third try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
          sxaddpar,header,'History','Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
          sxaddpar,header,'FileDate', dirlist[i]


       		fxbcreate,unit,l2file_chan1,header
         	fxbwritm,unit,['Time','Azimuth','Elevation','Revolution','Status','Flag'],Time_chan1,Azimuth_chan1,Elevation_chan1,revs_chan1,status,flags.chan1
         	fxbwritm,unit,['Data'],data_chan1
        	fxbfinish,unit
          ;********************************************************************************************************************

	;*******************************************Write Fits file for chan2***********************************
		;l2file_chan2=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan2.fits',/Remove_all)
		l2file_chan2=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan2.fits',/Remove_all)
         	fxhmake,header,/initialize,/extend,/date
        	fxwrite,l2file_chan2,header
         	fxbhmake,header,n_elements(azimuth_chan2),'CutData','CutData'
         	fxbaddcol,1,header,time_chan2[0],'Time'
         	fxbaddcol,2,header,azimuth_chan2[0],'Azimuth'
         	fxbaddcol,3,header,elevation_chan2[0],'Elevation'
        	fxbaddcol,4,header,revs_chan2[0],'Revolution'
        	fxbaddcol,5,header,flags.chan2[0],'Flag'
        	fxbaddcol,6,header,status[0],'Status'
          fxbaddcol,7,header,  data_chan2[*,0],'Data'

          sxaddpar,header, '          ','          '
         	sxaddpar,header,'Band', '37-44Ghz'
         	sxaddpar,header,'DataUnits', 'Kelvin'
         	sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
         	sxaddpar,header,'Pointunits', 'Degrees'
         	sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'CalT',string( calsacT[1])
          sxaddpar,header,'CalP', string(calsacP[1])
          sxaddpar,header,'History','Third try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
          sxaddpar,header,'History','Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
          sxaddpar,header,'FileDate', dirlist[i]

       		fxbcreate,unit,l2file_chan2,header
         	fxbwritm,unit,['Time','Azimuth','Elevation','Revolution','Status','Flag'],Time_chan2,Azimuth_chan2,Elevation_chan2,revs_chan2,Status,flags.chan2
         	fxbwritm,unit,['Data'],data_chan2
        	fxbfinish,unit
          ;********************************************************************************************************************
        ;*******************************************Write Fits file for chan3***********************************
		;l2file_chan3=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan3.fits',/Remove_all)
		l2file_chan3=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan3.fits',/Remove_all)
         	fxhmake,header,/initialize,/extend,/date
        	fxwrite,l2file_chan3,header
         	fxbhmake,header,n_elements(azimuth_chan3),'CutData','CutData'
         	fxbaddcol,1,header,time_chan3[0],'Time'
         	fxbaddcol,2,header,azimuth_chan3[0],'Azimuth'
         	fxbaddcol,3,header,elevation_chan3[0],'Elevation'
        	fxbaddcol,4,header,revs_chan3[0],'Revolution'
        	fxbaddcol,5,header,flags.chan2[0],'Flag'
        	fxbaddcol,6,header,status[0],'Status'
         	fxbaddcol,7,header,  data_chan3[*,0],'Data'

          sxaddpar,header, '          ','          '
         	sxaddpar,header,'Bands', ' 37-44Ghz'
         	sxaddpar,header,'DataUnits', 'Volts'
         	sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
         	sxaddpar,header,'Pointunits', 'Degrees'
         	sxaddpar,header,'Choptype', 'SquareWave'
         	sxaddpar,header,'CalT',string(calsacT[2])  ;Kelvin/Volt
          sxaddpar,header,'CalP', string(calsacP[2])  ;Kelvin/Volt
          sxaddpar,header,'History','Third try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
          sxaddpar,header,'History','Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
          sxaddpar,header,'FileDate', dirlist[i]


       		fxbcreate,unit,l2file_chan3,header
         	fxbwritm,unit,['Time','Azimuth','Elevation','Revolution','Status','Flag'],Time_chan3,Azimuth_chan3,Elevation_chan3,revs_chan3,Status,flags.chan3
         	fxbwritm,unit,['Data'],data_chan3
        	fxbfinish,unit
          ;********************************************************************************************************************


          ;*******************************************Write Fits file for chan6***********************************
		;l2file_chan6=strcompress(lvl2path+dirlist[i]+'\'+dirlist[i]+'_'+hr[hour]+'_chan6.fits',/Remove_all)
		l2file_chan6=strcompress(lvl2path+dirlist[i]+'/'+dirlist[i]+'_'+hr[hour]+'_chan6.fits',/Remove_all)
         	fxhmake,header,/initialize,/extend,/date
        	fxwrite,l2file_chan6,header
         	fxbhmake,header,n_elements(azimuth_chan6),'CutData','CutData'
         	fxbaddcol,1,header,time_chan6[0],'Time'
         	fxbaddcol,2,header,azimuth_chan6[0],'Azimuth'
         	fxbaddcol,3,header,elevation_chan6[0],'Elevation'
        	fxbaddcol,4,header,revs_chan6[0],'Revolution'
        	fxbaddcol,5,header,flags.chan6[0],'Flag'
        	fxbaddcol,6,header,status[0],'Status'
          fxbaddcol,7,header,data_chan6[*,0],'Data'

          sxaddpar,header, '          ','          '
         	sxaddpar,header,'Band', '38-45Ghz'
         	sxaddpar,header,'DataUnits', 'Volts'
         	sxaddpar,header,'Timeunits', 'FractionalHoursOfDay'
         	sxaddpar,header,'Pointunits', 'Degrees'
         	sxaddpar,header,'Choptype', 'SquareWave'
          sxaddpar,header,'CalT',string( calsacT[3])
          sxaddpar,header,'CalP', string(calsacP[3])
          sxaddpar,header,'History','Third try 9-20-2009 changing structure of files including upto date pointing corrections and calibrations by Brian Williams'
          sxaddpar,header,'History','Fourth try 9-10-2011 corrected Level1 files and removed spaces for Destriping codes by Brian Williams'
          sxaddpar,header,'FileDate', dirlist[i]

       		fxbcreate,unit,l2file_chan6,header
         	fxbwritm,unit,['Time','Azimuth','Elevation','Revolution','Status','Flag'],Time_chan6,Azimuth_chan6,Elevation_chan6,revs_chan6,Status,flags.chan6
         	fxbwritm,unit,['Data'],data_chan6
        	fxbfinish,unit
          ;********************************************************************************************************************
      
       endif 	 ;end if for filedata check
    endfor   ; end of hour loop
endfor    ; end of directory loop

print,'Starttime is ' + starttime
print,systime()



para_angle=az[*]*0.0
      lst=sidetime(month,day,year,time,longitude)
      Azel2radec,az,elev,ra,dec,latitude,lst
      ra=ra mod(24)
      rad=ra*360./24.
      phi = rad *( (2*!DPI)/360. )
      theta = (90. - dec) * (!DPI / 180.)

end



