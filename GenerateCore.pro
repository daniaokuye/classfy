;+
; :Function: 使用isodata方法注入核心
;-报文生成函数
;
PRO GenerateCore,maskFile,logFile,abandon
;  RESOLVE_ROUTINE, ['nonWetland','staticClass','modeH',$
  ;    'isodata','modeALL','deleteFID'],$
  ;    /IS_FUNCTION, /NO_RECOMPILE
  path = 'F:\Data\test\L028_Index.tif'
  ;maskFile='F:\Data\IsoCopy\tpCode.tif'
  ;logFile = 'D:\360Downloads\5index\log.txt'
  ;abandon=0
  envi_open_file, path, r_fid=zzfid
  envi_file_query,zzfid,nl=nl,ns=ns,nb=nb,dims=dims
  print,'ns',ns,'nl:',nl,'nb:',nb,dims
  envi_open_file, maskFile, r_fid=mfid
  envi_file_query,mfid,nl=mnl,ns=mns,dims=mdims

  ;data = envi_get_data(fid=zzfid,dims=dims,pos=0)
  maskClass = envi_get_data(fid=mfid,dims=mdims,pos=0)
  ;print,zzfid,mfid,'dims of both',dims,mdims


  openW,lun,logFile,/GET_LUN
  PRINTF,lun,'The image was selected to generate the core for classification:'+STRING(10b)+$
    path
  ;FREE_LUN,lun
  ;classNUMs=max(maskClass);??????????????????????????;这有个恶心的7
  ;;现分类总数
  classNUMs=HISTOGRAM(maskClass,min=min(maskClass))
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(maskClass);maskdata是byte格式，所以binsize为1，现在变成DN值
  if abandon ne !NULL then begin
    foreach ab, abandon do begin
      ab=where(classNUMs eq ab)
      if ab ne -1 then classNUMs[ab]=-1
    endforeach
    ab=where(classNUMs ne -1)
    classNUMs=classNUMs[ab]
  endif

  e = ENVI()
  nameString=[]
  for i=0,N_ELEMENTS(classNUMs)-1 do begin
    ;openU,lun,logFile,/GET_LUN
    PRINTF,lun,format='(/,i-,a)',classnums[i],'is the value in calculating now'
    ;如果i大于类的总数，开始计算影像中未被分类的ISOData特征；
    print,format='($,i)',classnums[i]
    maskData=maskClass eq classnums[i]
    ;if i gt 7 then begin
    ;  maskData=nonWetland(zzfid,maskData)
    ;endif

    tempFile = e.GetTemporaryFilename();掩膜文件
    tempFile2 = e.GetTemporaryFilename();掩膜结果
    newRaster = ENVIRaster(maskData, URI=tempFile)
    newRaster.Save
    Mfid = ENVIRasterToFID(newRaster)

    ;    pos=indgen(nb)
    ;    ENVI_DOIT, 'ENVI_MASK_APPLY_DOIT', DIMS=dims, FID=zzfid,$;/IN_MEMORY,
    ;      M_FID=fid, M_POS=0, OUT_BNAME='mask',$
    ;      OUT_NAME=tempFile2, POS=pos, R_FID=rFid , VALUE=-32768;value是背景值

    isoFid = isodata(zzFid,Mfid);非监督分类结果
    IsSaveLog = staticClass(isoFid,zzfid,lun)

    ; Delete the output raster file when we're done with it
    ;DELfid=[Mfid,isoFid];[fid,rFid]
    ;ok=deleteFID(DELfid)
    print,'ok-1'
  endfor

  FREE_LUN,lun
  ;raster = data*mask
  ;isodataFile = isodata(mask*temporary(envi_get_data(fid=zzfid,dims=dims,pos=0)))
  ;delvar,raster
  ;    isodataFile='D:\360Downloads\test\test_isodata.tif'
  ;  IsSaveLog = staticClass(isodataFile,zzfid)
  print,'ok-2'
end

;+
; :得到一个maskdata，不属于任何一类又不包含边界:
;BackGloc是背景值取值
;-
function nonWetland,fid,maskData
  compile_opt idl2
  envi_file_query,fid,dims=dims
  data = envi_get_data(fid=fid,dims=dims,pos=0)
  maskData *= (data ne data[0,0])
  data=[]
  return, maskData
end

;+
; :统计生成报表:
;
;-
function staticClass,IFid,zzfid,lun
  envi_file_query,IFid,NUM_CLASSES =nums,CLASS_NAMES =clsname
  envi_file_query,zzfid,nb=nb,dims=dims
  PTR=(indgen(nums))[where(clsname ne 'Unclassified')]
  ENVI_DOIT, 'CLASS_STATS_DOIT', $
    CLASS_FID=Ifid, FID=zzfid, CLASS_PTR=PTR, POS=indgen(nb) , $
    CLASS_DIMS=dims, COMP_FLAG=2 , HIST=hist ,$;HIST:返回的hist为destiny密度，其总数和max-min相关 BINSIZE = (MAX – MIN) / (NBINS – 1).
    DMAX=dmax, DMIN=dmin, MEAN=dmean, STDV=dstev
  ;求众数
  dmode = modeH(hist,dmax,dmin)
  dcount = total(hist,1)
  result = [[[dmin]],[[dmax]],[[dmean]],[[dstev]],[[dmode]],[[dcount]]]
  result =TRANSPOSE(result,[2,0,1]);列行维分别是parameter、bands、classes。
  sz=size(result,/DIMENSIONS)
  if N_ELEMENTS(sz) gt 2 then sizeP=sz[2]-1 else sizeP=0
  for i=0,sizeP do begin;classes
    print,format='($,a)','*:'
    if result[-1,0,i] eq 0 then continue;count不为0时才打印
    PRINTF,lun,format='(/,6a+25,i)','min','max','mean','steddev','mode','pixels',i+1
    PRINTF,lun,format='(6(g,:,","))',result[*,*,i]
  endfor
  return, 1
end

;+
; :'CLASS_STATS_DOIT'中直方图中nbinsize:
;BINSIZE = (MAX – MIN) / (NBINS – 1)
;bins是直方图中的有效个数，
;doit统计中将每一个class下的band放到256个值下，
;但是不是拉伸到256这样的宽度下
;-
function modeH,hist,dmax,dmin
  ;首先统计每一列下的有效区间
  index = indgen(size(hist,/dimensions))
  index = (index mod (size(hist))[1]) +1;得到的是按列的空间复制
  histStat = (hist gt 0)*index;有效值区间，并得到其自动排序.但是现在还没处理0值
  ;打算使用min、max得到区间
  OMasked = total(histstat,1)/(total(histstat gt 0 , 1));用这个值代替histstat中的0值
  shape = size(OMasked,/DIMENSIONS);二维数组变成3维
  OM=[]
  cols=intarr((size(hist,/dimensions))[0])+1
  if N_ELEMENTS(shape) gt 1 then sizeOM=shape[1]-1 else sizeOM=0
  for i=0, sizeOM do begin
    temp=TRANSPOSE(omasked[*,i])##cols
    OM=[[[OM]],[[temp]]]
  endfor

  histstat += (hist eq 0)*OM;asked
  min=min(histStat,DIMENSION=1);最值同时也是序列值
  max=max(histStat,DIMENSION=1)
  nbins=max-min

  density=max(hist,DIMENSION=1,location);b是频率最大值，hist的列行维分别为nbins、bands、classes
  BINSIZE=round((dmax-dmin)/nbins)
  loc = (location mod (size(hist))[1])-min+1;这是各个列上的第n个。min+(max-min)/nbins*nth是它的DN值。
  mode = BINSIZE*loc+dmin
  return, mode
end


;+
; :isodata分类:
;
;-
function isodata,fid,M_fid
  compile_opt idl2
  envi_file_query,fid,dims=dims,nb=nb
  ;isodata的参数
  NUM_CLASSES = 15
  MIN_CLASSES = 3
  ITERATIONS = 10
  CHANGE_THRESH = .05
  ISO_MIN_PIXELS = 1000
  ISO_SPLIT_STD = 0.010
  ISO_MERGE_DIST = 0.0050
  ;ISO_SPLIT_SMULT = 1

  ISO_MERGE_PAIRS = 2

  out_bname='IsoData'
  e = ENVI()
  tempFile = e.GetTemporaryFilename()

  ENVI_DOIT,'class_doit',fid=fid,pos=indgen(nb),dims=dims,$
    out_bname=out_bname,out_name=tempFile,method=4,$
    r_fid=r_fid,M_FID=M_fid,M_POS=0,$
    CHANGE_THRESH = CHANGE_THRESH ,$
    NUM_CLASSES = NUM_CLASSES ,$
    ITERATIONS = ITERATIONS ,$
    ISO_MERGE_DIST = ISO_MERGE_DIST ,$
    ISO_MERGE_PAIRS = ISO_MERGE_PAIRS ,$
    ISO_MIN_PIXELS = ISO_MIN_PIXELS ,$
    ;ISO_SPLIT_SMULT = ISO_SPLIT_SMULT ,$
    ISO_SPLIT_STD = ISO_SPLIT_STD ,$
    MIN_CLASSES = MIN_CLASSES ,$
    in_memory=0

  return, r_fid
end

;+
; :众数:
;
;-
function modeALL,a
  his = histogram(a)
  idx = where(his eq max(his))
  return,a[idx[0]]
end
;+
; :删除fid及相关的文件:
;
;-
function deleteFID,fid
  for i=0,N_ELEMENTS(fid)-1 do begin
    raster = ENVIFIDToRaster(fid[i])
    path=raster.uri
    raster.close
    r_Pos = STRPOS(path,'.dat')
    fileName= STRMID(path,0,r_Pos)+'*'
    image_files=file_search(fileName,count=numfiles)
    for k=0,numfiles-1 do begin
      FILE_DELETE, image_files[k]
    endfor
  endfor
  return, 1
end
