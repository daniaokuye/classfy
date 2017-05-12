;+
; :Description:
;
;-
pro subarea,hardCore
  compile_opt idl2
  path = 'F:\Data\test\L028.tif'
  path = 'F:\Data\test\L028_Index.tif'
  outputfile = 'D:\360Downloads\test2.tif'
  outputfile = 'D:\360Downloads\5index\test2.tif'
  subNum=6.0;分块的在行列上的次数

  envi_open_file, path, r_fid=fid
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  ;map_info = envi_get_map_info(fid = fid)

  classed=intarr(ns,nl);用来装载分类后的结果
  colSZ=ceil(ns/subNum)
  lineSZ=ceil(nl/subNum)

  ;setbackCore,hardCore,std,extent,available
  DELfid=[];\\\\\\\\\\\\
  for L=0,subNum-1 do begin
    star_l=lineSZ*L
    end_l=lineSZ*(L+1)-1
    if(end_l ge nl) then end_l=nl-1
    for S=0,subNum-1 do begin
      star_s=colSZ*S
      end_s=colSZ*(S+1)-1
      if(end_s ge ns) then end_s=ns-1
      dim_new=[-1,star_s,end_s,star_l,end_l]
      data=[]
      for i=0,nb-1 do begin
        new = envi_get_data(fid = fid,dims = dim_new,pos = i)
        data = [[[data]],[[new]]]
      endfor
      print,format='($,a)','*:'
      sa_readDBF,data,hardCore,r_fid;,result
      ;////////////
      envi_file_query,r_fid,dims=dim;,nb=nb
      ;dim=[0,0,col,0,line]
      dima=[-1,0,end_s-star_s,0,end_l-star_l];/////////////
      if ARRAY_EQUAL(dim,dima) ne 1 then alert=DIALOG_MESSAGE('dims ne dim',TITLE='ERROR',/ERROR)
      result = envi_get_data(fid = r_fid,dims =dim,pos = 0)
      DELfid=[DELfid,r_Fid]
      ;////////////
      ;classfy,data,hardCore,std,extent,result
      classed[star_s:end_s,star_l:end_l]=result
    endfor
  endfor
  ;obtain meaadata of a classification
  ;只为获取数据头
  b=ENVIFIDToRaster(r_fid)
  meta=b.METADATA
  ;
  ;  data=[]
  ;  for i=0,nb-1 do begin
  ;    new = envi_get_data(fid = fid,dims = dim_new,pos = i)
  ;    data = [[[data]],[[new]]]
  ;  endfor
  ;  sa_readDBF,fid,hardCore,result
  ;sv=sv_IMG(map_info,outputfile,classed)
  a=ENVIFIDToRaster(fid)
  newRaster = ENVIRaster(classed, URI=outputfile,SPATIALREF=a.SPATIALREF,METADATA=meta)
  newRaster.Save
  ok=deleteFID(DELfid)

end


;+
; :Description:按照avaliable来重新组合影像
;
;-
function combineClass,file,ava
  compile_opt idl2
  file = 'D:\360Downloads\5index\test2.tif'
  envi_open_file, file, r_fid=fid
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  data = byte(envi_get_data(fid = fid,dims =dims,pos = 0))
  for j=0,N_ELEMENTS(ava)-1 do begin
    key=ava[j]
    print,key
    for i=0,N_ELEMENTS(key)-1 do begin
      data=(data eq key[i])*byte(j+101) + (data ne key[i])*data
    endfor
  endfor
  a=ENVIFIDToRaster(fid)
  fileName= STRMID(file,0,STRPOS(file,'.tif'))
  outputfile= fileName+'_combine.tif'
  newRaster = ENVIRaster(data-100, URI=outputfile,SPATIALREF=a.SPATIALREF)
  newRaster.Save
  return, 1
end


