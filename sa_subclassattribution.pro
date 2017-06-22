;亚类整合
pro sa_subClassAttribution
  ;1.获得大类代号，并统计亚类内相关统计
  maskFile='F:\Data\IsoCopy\tpCode.tif'
  path = 'F:\Data\test\L028_Index.tif'
  classPath = 'D:\360Downloads\5index\test2_gradClass.tif'
  envi_open_file, path, r_fid=zzfid
  envi_open_file, maskFile, r_fid=fid
  envi_open_file, classPath, r_fid=c_fid
  newFid=obClass(fid,c_fid);zzfid,
  logFile = 'D:\test\log.txt'
  openW,lun,logFile,/GET_LUN
  for i=0,N_ELEMENTS(newfid)-1 do begin
    PRINTF,lun,format='(/,i-,a)',i+1,'is the value in calculating now'
    IsSaveLog = staticClass(newfid[i],zzfid,lun)
  endfor
  FREE_LUN,lun

  newFid=rt()
  allCluster=coreA(logFile)
  p=generClass(allCluster,newFid)

  ;2.统计相关关系
  ;  inout='D:\test\allInside.txt'
  ;  fs=fileread(logFile);
  ;  hardcore=make_array(5,n_elements(fs),/DOUBLE)
  ;  for n=0L,n_elements(fs)-1 DO BEGIN
  ;    rbool=StringToDoubleArray(fs[n],dataArray,Count,Pos);
  ;    hardcore[*,n] = dataArray
  ;  endfor
  ;  fs=[]

  print,'o'



end

;+
; :Description:adjust the relative core
;注意一下，是不是每次都是从2开始，检测错误
;-
function coreA,logFile
  compile_opt idl2
  meanC=pickColumn(logFile,'mean');meanC是从0开始的序数
  h=HISTOGRAM(meanc[0,*],min=min(meanc[0,*]));h为频数
  binsize=( max(meanc[0,*])-min(meanc[0,*]) )/ (N_ELEMENTS(h)-1)

  reSave= 'D:\test\clusters.txt'
  openW,lun,reSave,/GET_LUN

  allCluster=list()
  for i =0,N_ELEMENTS(h)-1 do begin;特别注意，有不足22个的。
    x=i*binsize+min(meanc[0,*]);对应的大类序号
    print,'now is :',x
    start=total(h[0:i])-h[i]
    hardcore=FLTARR(N_ELEMENTS(meanC[2:-1,0]),max(h))
    ;调整到相应的位置，无值的地方被填充为0，填满总类数目
    hardcore[*,meanC[1,start:start+h[i]-1]-1]=meanC[2:-1,start:start+h[i]-1]
    ;有一个无效值0,也可能是最大值
    hardcore=hardcore[*,1:-1];注意序数还是从0 开始哦。
    sa_matrixCore,realloc,hardcore,SIDmatrix;sidmatrix序数也还是从0 开始哦
    hardcore=[]
    clusters=cluster(SIDmatrix)
    allCluster.add,clusters
    line=[]
    PRINTF,lun,x
    foreach am, clusters do begin
      f=strcompress('('+String(N_ELEMENTS(am))+'(i,:,","))',/REMOVE_ALL)
      ;print,'ffffffff:::',f
      PRINTF,lun,format=f,am
    endforeach

    ;这儿是对核心的一些归属相似的计算
    ;inner action
    ;target class action
    ;a[i]

    ;none target class action

    ;outer action
  endfor
  FREE_LUN,lun
  print,'w'
  return, allCluster
end



;+
; :obtain all the given class:
;
;-
function obClass,fid,c_fid
  compile_opt idl2
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  envi_file_query,c_fid,NUM_CLASSES =SN,CLASS_NAMES =clsname
  x=where(clsname eq 0 or clsname eq 'Unclassified',count)
  if count then SN-=1

  ;obtain meaadata of a classification
  ;只为获取数据头
  e=ENVI()
  tempFile = e.GetTemporaryFilename();
  envi_doit, 'class_doit', fid=c_fid, pos=0,$
    dims=dims, r_fid=r_fid, $;m_fid=Mfid,
    out_bname='min', method=1, out_name=tempFile, $;
    mean=TRANSPOSE(INDGEN(SN)+1), class_names=sindgen(SN+1), $
    lookup= byte(randomu(1,[3,SN+1])*255), in_memory=0
  b=ENVIFIDToRaster(r_fid)
  meta=b.METADATA

  maskClass = envi_get_data(fid=fid,dims=dims,pos=0)
  x=HISTOGRAM(maskClass,min=min(maskClass))
  x=where(x gt 0);变成了序数,空白类255变成序数
  x=x+min(maskClass)
  s_fid=INDGEN(N_ELEMENTS(x))
  fileSV=[]
  ;del_fid=INDGEN(N_ELEMENTS(x))

  for i=0,N_ELEMENTS(x)-1 do begin
    tempFile = e.GetTemporaryFilename();
    maskData=envi_get_data(fid=c_fid,dims=dims,pos=0)
    ;将单个大类单抠出来
    maskData=(maskClass eq x[i])*maskData; + (maskClass ne x[i])*(SN+1)
    a=ENVIFIDToRaster(fid)
    newRaster = ENVIRaster(maskData, URI=tempFile,SPATIALREF=a.SPATIALREF,METADATA=meta)
    newRaster.Save
    Mfid = ENVIRasterToFID(newRaster)
    maskData=[]
    fileSV=[fileSV,tempFile]
    ;print,tempFile
    s_fid[i]=Mfid
  endfor
  ok=deleteFID(r_fid)

  saveFile='C:\temp\1.txt'
  openW,lun,saveFile,/GET_LUN
  PRINTF,lun,fileSV
  FREE_LUN,lun

  return, s_fid
end

;+
; :Description:generate class follow by the direction of core
;
;-
function generClass,allCluster,s_fid;s_fid是截取的各个大类的空间位置
  compile_opt idl2
  ;获取available，大类附属的小类
  setbackCore,core,std,extent,available
  e_obj=ENVI()
  new_fid=[];存放修改了的分类图像fid
  a=list()
  number=1
  for i=0,(size(available))[2]-1 do begin
    av=where(available[*,i] gt 0,count)
    a.Add,indgen(count-1)+number
    number+=(count-1)
  endfor
  ;因为最后一组会干扰计算，所以全部换成0
  a[-1] *= 0
  primary=[]
  foreach e, a do primary=[primary,e[0]]
  envi_file_query,s_fid[0],nl=nl,ns=ns,nb=nb,dims=dims
  for i=0,N_ELEMENTS(s_fid)-1 do begin
    Data=envi_get_data(fid=s_fid[i],dims=dims,pos=0)
    print,'ppp model',size(data,/DIMENSIONS)
    ;a[i]只对最主要的也就是第一个进行处理
    ;allCluster[i]，首先查看有无最主要的，其次对相关的部分聚类
    j=contain(a[i,0],allCluster[i])
    print,'j',j
    if j ne -1 then BEGIN;如果有首要核心，先处理首要核心
      foreach part,allCluster[i,j] do begin
        if part eq a[i,0] then continue
        Data=(Data eq part)*byte((a[i,0])[0])+(Data ne part)*Data
        print,'jJj model',size(data,/DIMENSIONS),part
      endforeach
    endif
    print,'now the process of k'
    for k=0,N_ELEMENTS(allCluster[i])-1 do begin
      if k eq j then CONTINUE
      ;key=???????????????
      key=0
      foreach e,primary do begin
        key=where(allCluster[i,k] eq e)
        if key ne -1 then break
      endforeach
      if key eq -1 then key = 0
      key=byte((allCluster[i,k,key])[0])

      foreach part,allCluster[i,k] do begin
        print,'kK-- model+',size(data,/DIMENSIONS)
        Data=(Data eq part)*key+(Data ne part)*Data
        print,'kKk model',size(data,/DIMENSIONS),part
      endforeach
    endfor
    ;存储data
    tempFile = e_obj.GetTemporaryFilename();
    newRaster = ENVIRaster(Data, URI=tempFile);,SPATIALREF=a.SPATIALREF,METADATA=meta)
    newRaster.Save
    afid = ENVIRasterToFID(newRaster)
    new_fid=[new_fid,afid]
    Data=[]
  endfor
  Data = bytarr(ns,nl)
  for i=0,N_ELEMENTS(new_fid)-1 do begin
    Data +=envi_get_data(fid=new_fid[i],dims=dims,pos=0)
  endfor
  tempFile = e_obj.GetTemporaryFilename();
  newRaster = ENVIRaster(Data, URI=tempFile);,SPATIALREF=a.SPATIALREF,METADATA=meta)
  newRaster.Save
  print,'save to new class:',tempFile
  return, tempFile
end

;+
; :Description:wheather or not a in list b
;
;-
function contain,a,b;a in b or not
  j=-1
  for i=0,N_ELEMENTS(b)-1 do begin
    if where(b[i] eq a) ne -1 then j=i
  endfor
  return, j
end

;+
; :Description:
;
;-
function qw,list1,saveFile
  compile_opt idl2
  ;saveFile='C:\temp\1.txt'
  openW,lun,saveFile,/GET_LUN
  ;for i=0,N_ELEMENTS(list1)-1 do begin
  PRINTF,lun,list1
  ;endfor
  FREE_LUN,lun

  return, 1
end



function rt,file;----File='C:\temp\1.txt'
  compile_opt idl2
  File='C:\temp\1.txt'
  ;nLines=file_lines(File)
  fs=fileread(file)
  list1=[]
  e=ENVI()
  foreach l,fs do begin
    raster = e.OpenRaster(l);
    fid = ENVIRasterToFID(raster)
    list1=[list1,fid]
  endforeach

  print,'o'
  return, list1
end