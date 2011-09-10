function b_read_lvl2_chan,day,channel=channel

if not(keyword_set(channel)) then channel =6



;lvl2path = 'Z:\Bmachine\data\level2\'
lvl2path = '/smbmount/data2/BMachine/Data/Level2New/'
l2files=file_search(strcompress(lvl2path+day+'/'+'*_chan'+string(channel)+'.fits',/remove_all))
num=n_elements(l2files)

for i=0,num-1 do begin
	print,l2files[i],'     '  ,i
	read_fits_s,l2files[i],ext,prim

	if i eq 0 then begin
		datab=reform(prim.data,3,n_elements(prim.time))
		azb=prim.azimuth
		elb=prim.elevation
		revb=prim.revolution
		timeb=prim.time
		statusb=prim.status
		flagb=prim.flag
	endif else begin
		datab=[[datab],[reform(prim.data,3,n_elements(prim.time))]]
		azb=[azb,prim.azimuth]
		elb=[elb,prim.elevation]
		timeb=[timeb,prim.time]
		statusb=[statusb,prim.status]
		revb=[revb,prim.revolution]
		flagb=[flagb,prim.flag]
	endelse
endfor



nump=n_elements(prim.time)

tp=strpos(prim.hdr[30],'/')
d=execute(strmid(prim.hdr[30],0,tp-1))

pp=strpos(prim.hdr[31],'/')
d=execute(strmid(prim.hdr[31],0,pp-1))
cals=[calt,calp]

alldata=create_struct('Time',timeb,'Azimuth',azb,'Elevation',elb,'status',statusb,'flag',flagb,'cals',cals,'data',datab)
return, alldata
end
