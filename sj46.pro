;+
; :Author: DN
; 为了给定一个度量值，衡量分类前后的斑块从属度和价值度
; 仿照kappa系数，全参与的方式，得到大类和聚类的关系度量
; classF设定为聚类值，truth为大类图值
; locate是如何对聚类cluster分类
;-
pro sj46,koppa,clusterF,TruthF,path,koppaFile,newLog,locate;,explore
  COMPILE_OPT IDL2
  if clusterF eq !NULL then clusterF='D:\360Downloads\June\envitempfileWedJun211034212017663_1.dat'
  if TruthF eq !NULL then TruthF='F:\Data\IsoCopy\tpCode.tif'
  if koppaFile eq !NULL then koppaFile = 'D:\360Downloads\June\koppa.txt'
  ;if clusterF eq !NULL then clusterF='D:\360Downloads\June\envitempfileWedJun211034212017663_1.dat'
  if path eq !NULL then path = 'F:\Data\test\L028_Index.tif'
  if newLog eq !NULL then newLog = 'D:\360Downloads\June\log2.txt'
  ;name = envi_pickfile(title='选择图像')

  envi_open_file, clusterF, r_fid=cfid
  envi_open_file, TruthF, r_fid=tfid;truth
  IF (cfid EQ -1 or tfid EQ -1) THEN return


  envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
  maskClass = envi_get_data(fid=tfid,dims=dims,pos=0)
  cluster=envi_get_data(fid=cfid,dims=dims,pos=0);聚类分类
  ;to obtain the valueable data of class
  classNUMs=HISTOGRAM(maskClass,min=min(maskClass));这是关键，按照min-max排列为一维，无值即为0--直方图
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(maskClass);maskdata是byte格式，所以binsize为1，现在变成DN值

  ;cluster聚类类别
  clusterNUMs=HISTOGRAM(cluster,min=min(cluster))
  clusterNUMs=where(clusterNUMs gt 0)
  clusterNUMs=clusterNUMs+min(clusterNUMs)


  coeffient=[];从属度 same/older
  valuable=[];价值度  same/itself

  foreach clus,clusterNUMs do begin

    clusterData=cluster eq clus
    countClu = total(clusterData);count
    cTemp=[]
    vTemp=[]
    countCls=0
    foreach cls,classNUMs  do begin
      maskData=maskClass eq cls
      new=clusterData+maskData
      _ = where(new eq 2,countNew);count
      _ = where(maskData eq 1,countCls);count
      _=[]
      new=[]
      maskData=[]
      cTemp=[cTemp,1.0*countNew/countCls];
      vTemp=[vTemp,1.0*countNew/countClu]
      print,format='($,a)',' : '

    endforeach
    coeffient=[[coeffient],[cTemp]]
    valuable=[[valuable],[vTemp]]
    clusterData=[]
    print,' clusterNUMs:',clus
  endforeach

  ;Save the matrix of coeffient & valuable
  openW,lun,koppaFile,/GET_LUN;,/append
  f=strcompress('('+String(fix(N_ELEMENTS(classNUMs)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'degree of membership : same/older'
  PRINTF,lun,format=f,coeffient
  PRINTF,lun,'degree of valueable : same/itself'
  PRINTF,lun,format=f,valuable
  FREE_LUN,lun
  print,'--koppaFile has been saved--'

  explore=list();其维度等于clusterNUMs的，为每一个聚类的类注释需要分解成几个
  for p=0,n_elements(clusterNUMs)-1 do explore.add,[];list中加载一串列表
  ;查看需要被分解，以及被分解成几个部分的聚类
  for i=0,n_elements(classNUMs)-1 do begin
    cLine= where(coeffient[i,*] ge 0.1);列--得到的是行,此处为10%，像元数的比值
    foreach line, cLine do begin
      partial=where(valuable[*,line] gt 0.1,count);每聚类在大类下的比重,查看它需要被分解成几部分
      ;explore重要的包括本类，即i.然后包括它重点所在的地方，最大不过n_elements(clusterNUMs)-1
      if count gt 0 then explore[line]=[ explore[line] , [i,partial]]
    endforeach
    ;主要针对255中的特别类，需要保留到下一循环
    ;    vLine= where(valuable[i,*] ge 0.3,count);列--得到的是行,此处为30%以上为单一的某个类时的情况
    ;    for vI=0, count-1 do begin
    ;      ;如果此行此列占比达到2%，就会记录下来，往下传（城镇，道路，等非湿地）
    ;      if coeffient[i,vLine[vI]] gt 0.02 then explore[vLine[vI]]=[ explore[vLine[vI]] , i]
    ;    endfor
  endfor

  ;整理explore的每一个列表，去除重复的部分
  for p=0,n_elements(clusterNUMs)-1 do begin
    cN=explore[p]
    if cN ne [] then begin
      cN=cN[sort(cN)]
      cN=cN[uniq(cN)]
      explore[p]=cN
    endif
  endfor

  ;kappa number calculating同时也增加一个确定类别的list-->locate
  ;however the last column is useless for it belongs to 255,or nonetype
  coeffient[[0,n_elements(classNUMs)-1],*]=0.0;0也将变成255
  valuable[[0,n_elements(classNUMs)-1],*]=0.0;[0,n_elements(classNUMs)-1]
  koppa=0
  e_Index=0;explore的下标号
  locate=intarr(n_elements(clusterNUMs))+255;长度与explore相同，空余的类设为255，与输入maskFile相同
  foreach ex,explore do begin;line
    numOfE=0
    kp=0
    foreach divied,ex do begin;column
      numOfE++
      kp+=coeffient[divied,e_Index] * valuable[divied,e_Index]
      print,format='($,a)',kp
    endforeach

    if numOfE ne 0 then kp/=numOfE
    koppa+=kp
    print,'   Kp till Now:',koppa,' with the index:',e_Index

    ;calculate the which class should the cluster be belonged
    n_ex = n_elements(ex)
    if n_ex ge 1 then begin
      index=where( valuable[*,e_Index] eq max(valuable[*,e_Index]))
      if valuable[index,e_Index] gt 0.1 then locate[e_Index]=index
    endif
    e_Index++
  endforeach
  ;得到聚类的炸开
  envi_open_file, path, r_fid=fid
  ok=blastCluster(explore,cfid,fid,newLog)

  print,'blaster Cluster done!'
  ;要关掉打开的影像文件
  foreach i_fid,[cfid,fid,tfid] do begin
    catch,error_status
    ;如果发生不能转出fid的错误，就跳过
    if error_status ne 0 then begin
      catch,/cancel
      continue
    endif
    raster = ENVIFIDToRaster(i_fid)
    raster.close
  endforeach

  ;ENVI_BATCH_EXIT
end

;+
; :explore method:
; 将explore中的元素炸开
; clusterF,path,newLog分别是聚类影像，原始影像，报表记录
; clusterF --> cfid
; path --> fid
;-
function blastCluster,explore,cfid,fid,newLog
  compile_opt idl2
  ;  if clusterF eq !NULL then clusterF='D:\360Downloads\June\envitempfileWedJun211034212017663_1.dat'
  ;  if path eq !NULL then path = 'F:\Data\test\L028_Index.tif'
  ;  if newLog eq !NULL then newLog = 'D:\360Downloads\June\log2.txt'
  openW,Lun,newLog,/GET_LUN,/append
  ;
  ;  envi_open_file, path, r_fid=fid
  ;  envi_open_file, clusterF, r_fid=cfid
  envi_file_query,fid,nl=nl,ns=ns,dims=dims
  envi_file_query,cfid,dims=cdims
  ;得到大类分类图以及他的类别数目序列
  ;ClassD = envi_get_data(fid=fid,dims=dims,pos=0)
  clusterD=envi_get_data(fid=cfid,dims=cdims,pos=0)
  ;NUMs=HISTOGRAM(ClassD,min=min(ClassD))
  ;NUMs=where(NUMs gt 0);变成了序数,空白类255变成序数
  ;NUMs=NUMs+min(ClassD);maskdata是byte格式，所以binsize为1，现在变成DN值
  e_Index=-1;用来指示cluster
  e = ENVI()
  print,'explore:',explore
  foreach ex,explore do begin
    e_Index++
    cN=explore[e_Index]
    n_cN = n_elements(cN)
    if n_cN ge 1 then begin;只要不是空就该保留，因为成分关系重大
      cluD=clusterD eq e_Index
      tempFile = e.GetTemporaryFilename();掩膜文件
      newRaster = ENVIRaster(cluD, URI=tempFile)
      newRaster.Save
      Mfid = ENVIRasterToFID(newRaster)
      if n_cN eq 1 then n_cN=2;对于0的处理，图像的黑边,影像被分成两类，黑边，整个图像（但图像的均值没什么意义才对）

      isoFid = isodata(Fid,Mfid,n_cN);非监督分类结果
      PRINTF,lun,format='(/,i-,a)',e_Index,'is the value in calculating now'
      if isoFid ne -1 then IsSaveLog = staticClass(isoFid,fid,lun)
      print,' with the index:',e_Index,' and its parts to be break',n_cN
      ; Delete the output raster file when we're done with it
      DELfid=[Mfid,isoFid];[fid,rFid]
      ok=deleteFID(DELfid);先删掉掩膜临时文件
      ;foreach ifid,DELfid do begin
      ;  if ifid ne -1 then begin
      ;    raster = ENVIFIDToRaster(ifid)
      ;    path=raster.uri
      ;    print,path
      ;    raster.close
      ;  endif
      ;endforeach

      ;continue
    endif
  endforeach
  FREE_LUN,lun

  return, 1
end


