;分类后处理
pro postClass
  compile_opt idl2
  path = 'D:\360Downloads\test2_grad.img'
  path = 'D:\360Downloads\5index\test2_grad.img'
  out_name = 'D:\360Downloads\test2_gradClass.tif'
  out_name = 'D:\360Downloads\5index\test2_gradClass.tif'
  envi_open_file, path, r_fid=fid
  ;envi_open_file,out_name , r_fid=r_fid
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  map_info = envi_get_map_info(fid = fid)
  ;;;;
  setbackCore,hardCore,std,extent,available
  min_pop=1500;seg时最小的群落个数
  e = ENVI()

  a=list()
  number=1
  for i=0,7 do begin;7是类别的个数，参考generateCore
    av=where(available[*,i] gt 0,count)
    a.Add,indgen(count-1)+number
    number+=(count-1)
    ;print,number
  endfor
  ;;;;;
  b=a.Count()
  SN=max(a[b-1])
  hardCore=indgen(SN)+1
  ;将输入转成分类图，后面优化的时候放入分类中去
  envi_doit, 'class_doit', fid=fid, pos=indgen(nb),$
    dims=dims, r_fid=r_fid, $
    out_bname='min', method=1, out_name=out_name, $;
    mean=TRANSPOSE(hardCore), class_names=indgen(SN+1), $
    lookup= bytarr(3,SN+1), in_memory=0
  print,'1st：','大类'
  S_fid=intarr(N_ELEMENTS(a));保存seg的各个返回fid，注意有新生成的文件
  ;  object=bytarr(ns,nl)
  for i=0,N_ELEMENTS(a)-1 do begin
    tempFile = e.GetTemporaryFilename()
    ENVI_DOIT, 'ENVI_SEGMENT_DOIT', ALL_NEIGHBORS=1,$
      CLASS_PTR=a[i], DIMS=dims, FID=r_fid,$
      MIN_POPULATION=min_pop, $
      OUT_NAME=tempFile, POS=0, R_FID=temp
    S_fid[i]=temp
  endfor
  print,'2nd：','遗留点'
  svImgList=list();遗漏点的临时文件，亚类数目
  print,'now'
  for i=0,N_ELEMENTS(a)-1 do begin
    for j=0,N_ELEMENTS(a[i])-1 do begin
      ret = reservedData(fid,S_fid[i],a[i,j],e)
      svImgList.add,ret
    endfor
  endfor
  print,svImgList
  ok=deleteFID(S_fid);第一步的文件删除，大类数目

  ;the third step
  print,'3th:','同类转换'
  water=[1,12,9,10,11,13];2,8,
  nonwater=[3,2,8,4,5,6,7]
  water=[0,3,4,5,6];5index
  nonwater=[1,2];5index
  x=[]
  for i=0,N_ELEMENTS(water)-1 do begin
    x=[x,a[water[i]]]
  endfor
  water = x
  x=[]
  for i=0,N_ELEMENTS(nonwater)-1  do begin
    x=[x,a[nonwater[i]]]
  endfor
  nonwater = x
  
  v=getbody(svImgList,water,waterbody);water
  print,'max of water:',max(waterbody),mean(waterbody)
  ;第二步的文件通过getbody（）删除，亚类数目，不包括最后一大类
  waterbody *=2
  v=getbody(svImgList,nonwater,waterbody);nonwater
  print,'max of nonwater:',max(waterbody),mean(waterbody)
  waterbody +=SN
  waterbody[where(waterbody eq SN)]=0


  print,'4th:','遗留点加入'
  data = envi_get_data(fid = fid,dims = dims,pos = 0)
  data[where(waterbody gt 0)]=0
  data += waterbody
  SN=max(data)
  waterbody=[]
  outputfile = e.GetTemporaryFilename()
  print,'data of 4th output:',outputfile
  newRaster = ENVIRaster(data, URI=outputfile)
  newRaster.Save
  Mfid = ENVIRasterToFID(newRaster);新生成临时文件，4th——1，mfid
  data=[]
  ;sv=sv_IMG(map_info,outputfile,data);check the statics


  hardCore=indgen(SN)+1
  ;将输出转成分类图——新生成临时文件，4th——2，a_fid
  tempFile = e.GetTemporaryFilename()
  print,'class with output:',tempFile
  envi_doit, 'class_doit', fid=Mfid, pos=0,$
    dims=dims, r_fid=a_fid, $
    out_bname='min', method=1, out_name=tempFile, $;
    mean=TRANSPOSE(hardCore), class_names=indgen(SN+1), $
    lookup= bytarr(3,SN+1), in_memory=0

  ;调整大类组合
  for i=0,N_ELEMENTS(a)-1 do begin
    if where(a[i,0] eq water) ne -1 then begin
      a[i]=[a[i],SN];大号是water
    endif else if where(a[i,0] eq nonwater) ne -1 then begin
      a[i]=[a[i],SN-1]
    endif
  endfor
  print,'now of a',a

  N_fid=intarr(N_ELEMENTS(a));保存【new 大类 】的各个返回fid，注意有新生成的文件
  for i=0,N_ELEMENTS(a)-1 do begin
    tempFile = e.GetTemporaryFilename()
    print,'seg of ',i,'th output:',tempFile
    ENVI_DOIT, 'ENVI_SEGMENT_DOIT', ALL_NEIGHBORS=1,$
      CLASS_PTR=a[i], DIMS=dims, FID=a_fid,$
      MIN_POPULATION=min_pop, $
      OUT_NAME=tempFile, POS=0, R_FID=temp
    N_fid[i]=temp
  endfor
  ;delete temporate images,4th
  ;ok=deleteFID([Mfid,a_fid])
  print,'5fh'
  data=bytarr(ns,nl)
  for i=0,N_ELEMENTS(N_fid)-1 do begin
    temp = envi_get_data(fid = N_fid[i],dims = dims,pos = 0)
    data += temp*(i+1)
    data[where(data gt i+1)]=i+1
    temp=[]
  endfor
  ok=deleteFID(N_fid);删掉4th
  sv=sv_IMG(map_info,out_name,data)
  print,'od'
end

;+
; get reserved points or objects
;
;-
function reservedData,fid,S_fid,k,e;k用来控制分类中对应的类别
  compile_opt idl2
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  map_info = envi_get_map_info(fid = fid)
  ;  water=[1,2,8,12,9,10,11,13]
  ;  nonwater=[3,4,5,6,7]
  tempData = envi_get_data(fid = S_fid,dims = dims,pos = 0)
  tempdata = tempdata gt 0
  data = envi_get_data(fid = fid,dims = dims,pos = 0)
  ;test=bytarr(ns,nl)
  ;for j=0,N_ELEMENTS(a[i])do begin
  data = data eq k;该类别对应的像元位置
  ;print,'q',size(data)
  data = data-(tempData gt 0);这个矩阵应该只有（0，1，-1->255）,1代表遗漏点
  data[ where(data eq 255) ]=0
  ;data = byte(data)
  ;data = (data le 1)*data
  ;print,'q2',size(data)
  ;  if where(water eq 2 ne -1 then begin
  ;    object+=reserve gt 0
  ;  endif else if where(nonwater eq 2) ne -1 then begin
  ;    object+=(reserve gt 0)*2
  ;  endif
  ;  test=[]
  ;  reserve=[]
  ;endfor
  ;ok=deleteFID(temp)

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  outputfile = e.GetTemporaryFilename()
  sv=sv_IMG(map_info,outputfile,data)
  tempData=[]
  data=[]
  return,outputfile
end

;+
; :Description:
;
;-
function getbody,listImg,index,body
  compile_opt idl2
  fidList=list();最后删除这些文件
  for i=0, N_ELEMENTS(index)-1 do begin
    if i gt N_ELEMENTS(listImg)-1 then continue
    envi_open_file, listImg[index[i]], r_fid=fid
    fidList.add,fid
    envi_file_query,fid,dims=dims,DATA_TYPE=dt;nl=nl,ns=ns,nb=nb,
    ;map_info = envi_get_map_info(fid = fid)
    data = envi_get_data(fid = fid,dims = dims,pos = 0)
    ;print,'maxb',max(data)

    print,'maxb2:',max(data),' dt:',dt,' i:',index[i]
    ;print,'maxb__:',N_ELEMENTS(body)
    if N_ELEMENTS(body) gt 0 then body+=data else body=data
    print,'maxb3:',max(body);,size(body)
    ;    raster = ENVIFIDToRaster(fid)
    ;    raster.close
    data=[]
  endfor
  ok=deleteFID(fidList)
  ;print,'os'
  return, 1
end
