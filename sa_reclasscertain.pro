;对特定的亚类集合重新寻找分类核心。

pro sa_reclassCertain
  ;realloc是据相似度的亚类集合
  sa_matrixCore,xx1,SIDmatrix
  realloc=CLUSTER(SIDmatrix)

  classPath = 'D:\360Downloads\5index\test2_gradClass.tif'
  logFile = 'D:\test\log_reC.txt';保存generate生成的core
  maskFile='F:\Data\IsoCopy\tpCode.tif'
  reSave= 'D:\test\FileData_reC.txt';保存setback得到的一系列参数

  envi_open_file, classPath, r_fid=zzfid
  envi_file_query,zzfid,nl=nl,ns=ns,nb=nb,dims=dims
  Data=envi_get_data(fid=zzfid,dims=dims,pos=0)
  e_obj=ENVI()
  raster=ENVIFIDToRaster(zzfid)

  newdata = bytarr(ns,nl)
  key=0b;key不会太多，应该10以内，约1-3个
  foreach a,realloc do begin
    if(N_ELEMENTS(a) le 3) then CONTINUE;3是普通类别拥有的亚类数目
    key += 1
    ;获取mask
    foreach part,a do begin
      ;print,'type:',(size(newdata))[-2]
      newdata += (Data eq part)*byte(key)
    endforeach
  endforeach
  Data=[]
  envi_open_file, maskFile, r_fid=mfid
  maskClass = envi_get_data(fid=mfid,dims=dims,pos=0)
  maskClass = (maskClass eq 2 or maskClass eq 3);times一定是10
  times=byte(ceil(key/10.0)*10)
  newdata = maskClass*newdata + (1b-maskClass)*(newdata+times);大类内外多增加一类
  newdata[where(newdata eq times)]=0
  ;maskClass=[]
  ;保存mask
  tempFile = e_obj.GetTemporaryFilename();
  newRaster = ENVIRaster(newdata, URI=tempFile,SPATIALREF=raster.SPATIALREF)
  newRaster.Save
  print,tempFile

  GenerateCore,tempFile,logFile,[0];至此，以上数据最终保存在logfile中



end

;+
; :Description:range是上面用于重组的亚类组合
;
;-
function reObtainCore,logFile,targets
  reSave= 'D:\test\FileData_reC.txt';保存setback得到的一系列参数
  setbackCore,logFile,reSave,hardCore,std,extent,ava1;得到特定区域的核心
  setbackCore,non1,non2,Core,std,extent,ava2;得到全部的核心
  print,format='(5(f8.3))',core
  useful=[ava2[1], ava2[2]]-1 ;由类2,3得出。序号为1,2。减一的目的是回归成序号
  newcore=[ [Core[*,useful]] , [hardCore]];用来寻找新核与大类间对应关系
  totalNum=(size(core,/DIMENSIONS))[1]
  core=[[core],[hardCore]];newcore[*,a]
  sa_matrixCore,newcore,SIDmatrix
  x=relocate(SIDmatrix,indgen(N_ELEMENTS(useful)));得到两行数组；0.旧-新；1、旧-旧
  ;x=[  6 ,  8,7 ,9,10 ,  8]
  a=x[0];x0是原23大类各亚类所对应新核心在newcore中的序号
  b=x[1];x1是原core之间最相近的对应关系，可认为是新核转变的暗道
  d=byte([[a[b]],[a]]-N_ELEMENTS(useful)+1);ava不是序号，是第几个
  ;
  ;+这儿需要思考当ava不止两个时怎么办，即不止两个大类需要重组时的办法
  ;根据x填入core，根据class12填入ava2
  ;////range是和ava2有关联的
  ;Core[*,useful]=0
  ;
  ;core原位置被替代，需要将相同的core减掉，并且调整ava2


  class1=[]
  class2=[]
  foreach i,ava1[0] do begin
    loc=where(d eq i) mod N_ELEMENTS(useful);因为x是按useful排列的，总数也同于它
    w=d[loc,*]
    w=w[UNIQ(w[sort(w)])]
    m=0
    foreach key, w do begin
      if where(ava1[1] eq key) ne -1 then m+=1
    endforeach
    if m eq 0 then class1=[class1,i] else class2=[class2,i]
    ;class1=[1,2];class2=[3]
  endforeach
  if class1 eq !NULL or class2 eq !NULL then begin
    Result_MESSAGE = DIALOG_MESSAGE('大类缺失亚条目',TITLE='ERROR',/ERROR)
  endif
  print,ava2
  ;对第2、3和第7（外）大类添加、调整
  ava2[1]=class1+totalNum
  ava2[2]=class2+totalNum
  ava2[7]=[ava2[7],ava1[1]+totalNum]
  print,ava2
  ;将core中设为0的出去，同时大类包含亚类的列表更新
  allIndex=indgen((size(core,/DIMENSIONS))[1])
  allIndex[targets-1]=-1
  allIndex=allIndex[where(allIndex gt -1)]
  core=Core[*,allIndex]

  ;调整ava列表
  for i=0,N_ELEMENTS(ava2)-1 do begin
    for j=0,N_ELEMENTS(ava2[i])-1 do begin
      ava2[i,j]=where(allIndex eq ava2[i,j]-1)+1
    endfor
    ava2[i]=(ava2[i])[where(ava2[i] gt 0)]
  endfor
  ;print,format='(5(f8.3))',core
  print,ava2
  ;现在得到了新的核和分组
  subarea,Core

  return, 1
end


