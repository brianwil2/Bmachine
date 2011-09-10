pro b_map_tqu_lvl2_pipeline_d,dirl,filebase=filebase,nside=nside,dchp=dchp,achp=achp,qprime=qprime,dcfilter=dcfilter

;Naive mapmaking tool, uses leve21 fits files from dirl directory list
;start off assuming map of all 16 chan, I,Q,U in one big healpix matrix
; returns outmap as array outmap[number of channels,{I,Q,U},{map, error per pixel},Values of pixels)
; returns nobsmap as array nobsmap[number of channels, number of samples in a given pixel]
;test version uses _f filtering, forces first 8 channels for QU and next 8 for I
;set ignorecal keyword to skip filtering with the status word

;SB!  Latitude=34.425804
;SB!  longitude=119.714189

print,systime()
!quiet=1
Latitude=37.583375
Longitude=118.237095

RTD=180./!pi
dtr=!pi/180.
miss = -1.63750e30



if not(keyword_set(achp)) then achp=.015
if not(keyword_set(nside)) then nside=256
if not(keyword_set(filebase)) then filebase='NoBase'

npix=nside2npix(nside)
ndir=n_elements(dirl)
hr=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']

;note that dirl should be without the path, so we can parse the date

chans=[1,2,3,6]


;mappath_new='/project/projectdirs/cmb/cofe/bmachine/Level2_maps/'
;lvl2path='/project/projectdirs/cmb/cofe/bmachine/Level2/'

mappath_new='Z:\BMachine\Data\Level2_maps\'
lvl2path='Z:\BMachine\Data\Level2\'

for d=0,ndir-1 do begin ;cycle thru selected days
	outmap=fltarr(npix,3,2)
	summap=fltarr(3,npix)
	sumsqmap=fltarr(3,npix)
	nobsmap=fltarr(npix)

	mappath=strcompress(mappath_new+dirl[d]+'/')
	file_mkdir,mappath
	dat=dirl[d]
  	for c=0,3 do begin

  		lvl2d=b_read_lvl2_chan(dat,channel=chans[c])
  		demod=lvl2d.data
		year=float(strmid(dat,0,4))
   		month=float(strmid(dat,4,2))
		day=float(strmid(dat,6,2))

		if (keyword_set(dcfilter)) then begin
		    if not(keyword_set(dchp)) then dchp=.1
      		     demod[0,*]=highpass_f(demod[0,*],37.4,dchp)
      		endif else begin
      		     demod[0,*]=demod[0,*]-mean(demod[0,*])
		     demod[1,*]=demod[1,*]-mean(demod[1,*])
		     demod[2,*]=demod[2,*]-mean(demod[2,*])
      		endelse

      		demod[1,*]=highpass_f(demod[1,*],37.4,achp)
      		demod[2,*]=highpass_f(demod[2,*],37.4,achp)

		flag=Map_data_selection(lvl2d,channel=chans[c])
		gooddata=where(flag eq 0.0,counts)

		print, counts
		print, counts/n_elements(flag)

		dgood=fltarr(3,counts)
		dgood[0,*]=demod[0,gooddata]
      		dgood[1,*]=demod[1,gooddata]
      		dgood[2,*]=demod[2,gooddata]
		timed=lvl2d.time[gooddata]
		az=lvl2d.azimuth[gooddata]
		el=lvl2d.elevation[gooddata]


           	if not(keyword_set(qprime)) then begin
			num=floor(n_elements(timed)/2)
			aza=az[0:num]  &  azb=az[num+1:*]
			ela=el[0:num]  &  elb=el[num+1:*]
			timeda=timed[0:num]  &  timedb=timed[num+1:*]
      para_angle_a= get_parallactic_angle(month, day, year, timeda, aza, ela, latitude, longitude, altitude)
			para_angle_b= get_parallactic_angle(month, day, year, timedb, azb, elb, latitude, longitude, altitude)
			para_angle=[para_angle_a,para_angle_b]
                	QUsky=get_sky_Q_and_U(dgood[1,*],dgood[2,*],para_angle)
               		dgood[1,*]=QUsky[0,*]
                	dgood[2,*]=QUsky[1,*]
           	endif


            	lst=sidetime(month,day,year,timed,longitude)
            	Azel2radec,az,el,ra,dec,latitude,lst
            	ra=ra mod(24)
            	rad=ra*360./24.
            	phi = rad *( (2*!DPI)/360. )
            	theta = (90. - dec) * (!DPI / 180.)
            	ang2pix_ring,nside,theta,phi,ipring
		print,'The channel is '+string(chans[c])
		help, uniq(ipring)

           	for r=0l,counts-1 do begin
            		nobsmap(ipring(r))=nobsmap(ipring(r))+1
           		for j=0,2 do begin 						;cycle thru I,Q,U
             			summap[j,ipring[r]]=summap[j,ipring[r]]+dgood[j,r]
            			sumsqmap[j,ipring[r]]=sumsqmap[j,ipring[r]]+[dgood[j,r]]^2
             		endfor
            	endfor




    		for s=0,2 do begin
      			 h=where(nobsmap gt 0.0,complement=nothit)
      			 outmap(h,s,0)=summap(s,h)/nobsmap(h)
      			 outmap(h,s,1)=(sumsqmap(s,h)/nobsmap(h)-summap(s,h)*summap(s,h))/sqrt(nobsmap(h))
      			 outmap(nothit,s,*)=miss
			 nobsmap[nothit]=miss
    		endfor


			mapfilename=strcompress(mappath+'map_'+filebase+'_'+string(chans[c])+'.fits')
    			nobsfilename=strcompress(mappath+'nobs_'+filebase+'_'+string(chans[c])+'.fits')
    			write_tqu,mapfilename,reform(outmap[*,*,*]),coordsys='C',/ring
    			e=create_struct('HDR',' ', 'N_OBS',REFORM(nobsmap[*]))
    			write_fits_sb,nobsfilename,0,e,/ring,coord='C'

	endfor		;channel
	Print,systime()
endfor                 ;directory

print, systime()
end
