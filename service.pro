pro service
  for i=0,1 do begin
    print,'ooo',i
    ;i+=2
  endfor
  if -1 then print,'ooolk'
  ;print,norm(540.184,62.7,198)
end
;+
; :光谱角:
;
;-
function SAM,dataset,core
  data=dataset
  s=size(data,/dimensions)
  colum=s[0]&line=s[1]&bands=s[2]
  ;x*S
  data=TRANSPOSE(data,[2,0,1])
  up=fltarr(bands,colum,line)
  core_core=fltarr(colum,line)
  data_data=fltarr(colum,line)
  for i=0,line-1 do begin
    coredata=core#TRANSPOSE(intarr(colum)+1)
    core_core[*,i]=sqrt(total(1L*coredata*coredata,1))
    data_data[*,i]=sqrt(total(1L*data[*,*,i]*data[*,*,i],1))
    up[*,*,i]=coredata*data[*,*,i];;this is the virtual step -> eclipse 56
  endfor
  S=total(up,1)/(core_core*data_data)
  return,ACOS(S)
end

;+
; :计算欧式距离:
;
;-
function calEuclideanDis,data,Center
  if(size(data,/N_DIMENSIONS) eq 3) then nb=(size(data))[3] else nb=1
  EuDis=0.0
  ;print,'nb:',nb,'center:',center[0]
  for bandI=0,nb-1 do begin
    EuDis += ((data[*,*,bandI]-Center[bandI])*1L)^2
  endfor
  ;print,"eudis[0,0]:",eudis[0,0]
  return, EuDis
end

;+
; :标准正态分布:
;
;-
function st_norm,u
  x=abs(u)/sqrt(2)
  T=[0.0705230784,0.0422820123,0.0092705272,$
    0.0001520143,0.0002765672,0.0000430638]
  s=x*0
  foreach a,T,i do begin
    s+=a*(x^(i+1))
  endforeach
  E=1-DOUBLE((1+s)^(-16))
  cc=(u lt 0)+(-1)*(u ge 0)
  p=0.5-0.5*cc*E
  return, float(p)
end

;+
; :正态分布:
;
;-
function norm,a,sigma,x
  u=(x-a)/sigma
  return, st_norm(u)
end

;给定概率分布
function setFloate,data,hardcore,extent,std
  print,format='($,a)','*'
  s=size(data,/dimensions)
  colum=s[0]&line=s[1]&bands=s[2]
  preDtb = fltarr(colum,line,bands)
  assureDtb = intarr(colum,line)+1
  ;首先根据core曲线先给定分类概率
  U=SAM(data,hardcore)
  ;其次按照extend补齐概率
  for i=0,bands-1 do begin
    Z=data[*,*,i]
    exist=(Z lt extent[i,0]) * (Z gt extent[i,1])
    ;print,'exist:',total(exist),mean(exist),stddev(exist)
    assureDtb*=exist
    Z*=exist
    p=norm(hardcore[i],std[i],Z)
    ;print,'p:',total(p),mean(p),stddev(p)
    preDtb[*,*,i]=p*exist
    ;print,'ppreDtb:',i,total(preDtb[*,*,i]),mean(preDtb[*,*,i]),stddev(preDtb[*,*,i])
  endfor
  preDtb=min(preDtb,dimension=3)
  ;print,'ppreDtb:',total(preDtb),mean(preDtb),stddev(preDtb)
  reU=preDtb*assureDtb
  All=[[[U]],[[reU]]]
  All=max(All,dimension=3)
  return,All
end


;保存模型********************
function sv_IMG,map_info,path,data
  compile_opt IDL2,hidden
  ENVI,/RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT

;  envi_open_file, path, r_fid=fid
;  IF (fid EQ -1) THEN return,error

  ; envi_file_query,fid,ns=ns,data_type=dt,nl=nl
  dims=size(data)
  ns=dims[1]
  nl=dims[2]
  if size(data,/N_DIMENSIONS) eq 3 then begin
    nb=dims[3]
    dt=dims[4]
  endif else begin
    nb =1
    dt=dims[3]
  endelse

  r_Pos = STRPOS(path,'.tif')
  IF r_Pos[0] NE -1 THEN BEGIN
    fileName= STRMID(path,0,r_Pos)
    out_name = fileName+'_grad.img'
  ENDIF else begin
    out_name = path
  Endelse
  
  openw,lun,out_name,/get
  writeu,lun,data
  free_lun,lun

  ;map_info = envi_get_map_info(fid = fid)
  ENVI_SETUP_HEAD, fname=out_name, $
    ns=ns, nl=nl, nb=nb, $
    interleave=0, data_type=dt, $
    offset=0, /write,$
    MAP_INFO = map_info
  ;ENVI_BATCH_EXIT
  print,'save done:',out_name
end