;+
;计算各个层级的代号所对应的4/4影像中的代号
;主要是建立‘历史影像’和‘大簇’间的关系
;
;-
pro calentroy,newFile,maskFile,path,baseFile
  COMPILE_OPT IDL2
  ;如果有supply.tif这个数据，它就是额外的掩膜文件
  supply=baseFile+'supply.tif'
  newFile=supplyToMTKL(supply,newFile,path)

  ;默认按照[类别图 ,newFile]的组合方式
  stackFile=stack(newFile,maskFile)
  corr=stat(stackFile)
  code=obtainCode(corr,n_elements(newFile))
  tempFile=classSet(code,newFile,maskFile,supply);对通过的部分在次使用slope
  classFile= miniDistance(tempFile)
  print,'class file is: ',classFile
  print,'ok'
end


;+
; :Description:
;为组合而生成的code
;corr是按照[类别图 ,newFile]的组合成的相关矩阵，是空间相关系数
;N是聚簇的个数
;-
function obtainCode,corr,N
  compile_opt idl2
  zoneLine=N
  allN=n_elements(corr[*,0])
  zoneCol=allN-zoneLine
  corr=corr[0:zoneCol-1,allN-zoneLine:-1]

  code=list()
  for i=0,zoneCol-1 do code.add,[]
  for i=0,zoneLine-1 do begin
    tempCorr=corr[*,i]
    ;当存在相关系数大于0.1的时候
    if total(tempCorr gt 0.1) gt 0 then begin
      maxCorr=where(tempCorr eq max(tempCorr))
      tpI=maxCorr[0]
      code[tpI]=[code[tpI],i]
    endif
  endfor

  print,'code set as: ',code
  return, code
end


;+
; :Description:
;对历史分类图分割，产生独立的分类图像
;-maskList
function dividedClass,maskFile
  compile_opt idl2
  envi_open_file,maskFile, r_fid=fid
  envi_file_query,fid,dims=dims
  objRaster = ENVIFIDToRaster(fid)
  proj =objRaster.SPATIALREF

  maskClass = envi_get_data(fid=fid,dims=dims,pos=0)
  ;to obtain the valueable data of class
  classNUMs=HISTOGRAM(maskClass,min=min(maskClass));这是关键，按照min-max排列为一维，无值即为0--直方图
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(maskClass);maskdata是byte格式，所以binsize为1，现在变成DN值
  ;需要的那些不到200的编号
  classNUMs=classNUMs[where(classNUMs lt 200)]
  maskList=[]
  foreach nums,classNUMs do begin
    name=strtrim(string(nums),2)
    ;保存阈值文件
    maskData=maskClass eq nums
    tempFile=saveTemple(maskData,proj);data,name
    maskList=[maskList,tempFile]
  endforeach
  return,maskList
end




;+
;code=list()
;for i=0,6 do code.add,[];6是maxkClass的classNum数
;code[6]=[code[6],7,5];后面的数值代表fileList中的排序号
;code[4]=[code[4],6,4]
;code[1]=[code[1],3]
;code[0]=[code[0],0]
; :Description:本打算这是最后一个函数，根据多方归类准则，生成类别影像
; code是根据‘空间相关’得到的数据，没有被分开的数据将被给定一个200以上的值
;-
function classSet,code,fileList,maskFile,supply
  compile_opt idl2
  num=n_elements(code)
  unRegisted=200b;未注册部分从200开始编号
  codeList=[0];注册编号,0值肯定首先包含进去
  tacked=[]
  maskClass=[];save the intergert image

  envi_open_file, supply, r_fid=sfid
  if sfid ne -1 then begin;如果没有对应的掩膜文件，就是没必要有掩膜文件
    envi_file_query,sfid,nb=nb,dims=dims
    slope = envi_get_data(fid=sfid,dims=dims,pos=0)
  endif


  for i=0b,num-1 do begin
    this=code[i]
    tacked=[tacked,this]
    if n_elements(this) ge 1 then begin;非空时
      partCode=1b;组别编号
      foreach part,this do begin
        thresData=fileList[part]
        envi_open_file, thresData, r_fid=tfid
        envi_file_query,tfid,dims=dims
        Data=envi_get_data(fid=tfid,dims=dims,pos=0)
        ;在这儿加一步，对于通过编组的，用高程处理一下
        if sfid ne -1 then Data *= slope
        if maskClass eq [] then maskClass=((Data eq 1)*(i*10+partCode)) else $
          maskClass+=((Data eq 1)*(i*10+partCode))
        Data=[]
        codeList=[codeList,i*10+partCode];eq运算提供一个1
        ;异常值的处理，即重叠部分
        ok=unCompled(maskClass,codeList,i*10+partCode);重叠部分已剔除，所以直接加就好
        envi_file_mng,id = tfid,/remove;移除
        partCode++;不能提前，因为uncompled也用到了该变量
      endforeach
      print,format='($,a)','g'
      ;' registed is over'
    endif

    print,i
  endfor

  print,'tacked',tacked,' coeL',codeList
  ;未注册部分的处理
  for i=0,n_elements(fileList)-1 do begin
    if total(i eq tacked) eq 1 then continue;注册部分跳过
    thresData=fileList[i]
    envi_open_file, thresData, r_fid=tfid
    envi_file_query,tfid,dims=dims
    Data=envi_get_data(fid=tfid,dims=dims,pos=0)
    if maskClass eq [] then maskClass=((Data eq 1)*unRegisted) else $
      maskClass+=((Data eq 1)*unRegisted)

    codeList=[codeList,unRegisted];eq运算提供一个1
    ;异常值的处理，即重叠部分
    ok=unCompled(maskClass,codeList,unRegisted)
    print,format='($,a)','u'
    unRegisted++;不能提前，因为uncompled也用到了该变量
  endfor
  print,' Unregisted is over'

  ;envi_open_file,fileList[0], r_fid=mfid
  envi_open_file,maskFile, r_fid=mfid
  objRaster = ENVIFIDToRaster(mfid)
  proj =objRaster.SPATIALREF

  tempFile=saveTemple(maskClass,proj);project
  return,tempFile
end

;+
; :Description:异常值处理
;多部分影像合成分类图时，重叠部分找出的异常值，通过该方法归并到‘上一个的部分上’
;base 肯定是一个byte类型的数值
;-
function unCompled,maskClass,codeList,base
  compile_opt idl2
  classNUMs=HISTOGRAM(maskClass,min=min(maskClass));这是关键，按照min-max排列为一维，无值即为0--直方图
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(maskClass);maskdata是byte格式，所以binsize为1，现在变成DN值
  print,classNUMs
  foreach cln,classNUMs do begin
    if total(cln eq codeList) ne 1 then begin
      data=(maskClass eq cln)*base
      ;maskClass*=(maskClass ne cln)
      maskClass-=data;减掉重叠部分新增加的部分
      data=[]
    endif
  endforeach
  maskClass=BYTE(maskClass);理论上不能超过255，但也要小心看前面的处理
  return, 1
end




;+
; :Description:
;统计，computer statisics。envi中的函数，计算相关corrlection
;-
function stat,file
  compile_opt idl2
  envi_open_file,file, r_fid=fid
  envi_file_query,fid,nb=nb,dims=dims
  ;pos=[0,1]
  envi_doit, 'envi_stats_doit', dims=dims, fid=fid, pos=indgen(nb), $
    comp_flag=4, COV=cov,stdv=stdv,EVAL=eval, EVEC=evec
  ;corr=cov[0,1]/(stdv[0]*stdv[1])
  corr=cov
  base=stdv#transpose(stdv)
  for i=0,n_elements(cov)-1 do corr[i]=cov[i]/base[i]
  ;print,corr
  envi_file_mng,id = fid,/remove
  return, corr
end

;+
; :Description:
; layer stacking叠加
;-
function stack,thresholdFile,maskFile
  compile_opt idl2
  ;project of return file
  envi_open_file,maskFile, r_fid=mfid
  ;objRaster = ENVIFIDToRaster(mfid)
  ;proj2 =objRaster.SPATIALREF
  proj=ENVI_GET_PROJECTION(fid=mfid)
  ;envi_file_query,mfid,DATA_TYPE =dt,dims=dims

  ;class numbers and its respective file
  maskList=dividedClass(maskFile)
  stackList=[maskList,thresholdFile]

  tempFile=stackingExc(stackList,proj)
  return,tempFile
end

;+
; :Description:
;stackList是要组合的文件名，现在都是单一波段的组合
;以后可以试一试pos可不可以为、或者应该是list类型（for layer_stacking_doit)
;-
function stackingExc,stackList,proj,OUT_BNAME
  compile_opt idl2
  e=envi()
  tempFile = e.GetTemporaryFilename()
  num=n_elements(stackList)
  if OUT_BNAME eq !NULL then $
    OUT_BNAME = 'B '+strtrim(string(indgen(num)),2)
  d=[];dims矩阵
  dt=1;data type 1:byte
  thFidList=[]
  for i=0,n_elements(stackList)-1 do begin
    envi_open_file,stackList[i], r_fid=nfid
    envi_file_query,nfid,dims=dim,DATA_TYPE =dt
    d=[[d],[dim]]
    thFidList=[thFidList,nfid]
  endfor

  ;stacking
  ENVI_DOIT, 'ENVI_LAYER_STACKING_DOIT',dims=d,FID=thFidList,$;/EXCLUSIVE,
    OUT_DT=dt,OUT_BNAME=OUT_BNAME,OUT_NAME=tempFile ,OUT_PROJ=proj,$;,INTERP =0,,/IN_MEMORY
    OUT_PS=[30,30],pos=transpose(intarr(num));,R_FID=fidOUT_PS=proj.PIXEL_SIZE

  print,'stacking file is :',tempFile
  ;移除这些临时文件
  foreach ifid,thFidList  do envi_file_mng,id = ifid,/remove;, /DELETE

  return, tempFile
end





;+
; :最小值分类:
;暂时不用。
;现在用来将分类结果图，转成envi的分类图
;-
function miniDistance,file
  compile_opt idl2
  ;minmum distance
  envi_open_file, file, r_fid=fid
  envi_file_query,fid,nb=nb,dims=dims
  maskClass = envi_get_data(fid=fid,dims=dims,pos=0)
  classNUMs=numsClass(maskClass)
  len=n_elements(classNUMs)
  maskClass=[]
  className=strtrim(string(classNUMs),2)
  className=['Unclassified',className]
  fileName= STRMID(file,0,STRPOS(file,'.'))
  outputfile= fileName+'_class.tif'

  envi_doit, 'class_doit', fid=fid, pos=indgen(nb),$
    dims=dims, r_fid=r_fid, $
    out_bname='min', method=1, out_name=outputfile, $;
    mean=transpose(classNUMs), class_names=className, $
    ;lookup= bytarr(3,SN+1), in_memory=0
    lookup= byte(randomu(1,[3,len+1])*255), in_memory=0;,$
  ;
  envi_file_mng,id = fid,/remove
  return, outputfile;r_fid
end


;+
; :Description:
;
;-
function bond
  compile_opt idl2
  file=[]
  path = 'D:\360Downloads\July12\data\L120028_I.tif';影像文件，比如多指数文件
  envi_open_file,path, r_fid=fid
  miniFid=miniDistance(file)
  IsSaveLog = staticClass(miniFid,fid,lun)
  return, 1
end


;+
; :Description:外链
;根据p从n中找到关系,newOut/poentenOut
;-
function fari,n,p
  compile_opt idl2
  Log = 'D:\360Downloads\June\log4_wq.txt'

  num=n_elements(n);--1

  ;计算相似度，形成矩阵(据类号关系)
  outMatrix=[]
  len=0
  for line=0,num-1 do begin;n
    len=n_elements(n[line])
    inMatrix=[]
    for clo=0,num-1 do begin;p
      ratio=0.0
      ;diffO=0;用来验证是否属于另一项的子集，如是，则应化为一类
      ;diffI=0;同样作用
      for items=1,len-1 do begin;第一个肯定不同，所以不用比
        d_o=n_elements(n[line,items])
        d_i=n_elements(p[clo,items])
        degree=max([d_o,1]);非对称阵d_i;
        sames=0.0
        foreach s_o,n[line,items] do begin
          ;如果是[]，无论如何都会被跳过
          if total(s_o eq p[clo,items]) eq 1 then sames++
        endforeach
        r=sames/degree
        ;d=sames/d_o;相同的个数和它本身的个数
        ;if degree eq 0 then r=1.0
        ;if (d_o eq sames) then diffO++;或为空或为从属关系
        ;if (d_i eq sames) then diffI++;或为空或为从属关系
        ratio += r;最大值等于n_elements(history)-1,比如是7
      endfor
      ;if max([diffO,diffI]) eq len-1 then ratio=len-1
      inMatrix=[inMatrix,ratio]
    endfor
    outMatrix=[[outMatrix],[inMatrix]]
  endfor

  openW,lun,Log,/GET_LUN;,/append
  f=strcompress('('+String(fix(N_ELEMENTS(n)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'history similiarty:'
  printF,lun,format=f,outMatrix;————保存相似矩阵
  FREE_LUN,lun

  ;求相似
  threshold=floor((len-1)/2.0)
  similar=list()
  ;full=[]
  for i=0,num-1 do begin
    ;
    this=where(outMatrix[i,*] ge threshold)
    similar.add,this
    ;full=[full,this]
    ;full=full[sort(full)]
    ;full=full[uniq(full)]
  endfor

  return, similar
end




;+
; :Description:计算面积相互之间所占比例，聚簇和历史分类图间关系
;
;-
function sR,thresholdFile,maskFile
  compile_opt idl2
  Log = 'D:\360Downloads\June\log4_spat.txt'
  maskFile='F:\Data\IsoCopy\tpCode.tif'
  dims=[]
  thFidList=[]
  for i=0,n_elements(thresholdFile)-1 do begin
    envi_open_file,thresholdFile[i], r_fid=mfid
    if dims eq [] then envi_file_query,mfid,dims=dims
    thFidList=[thFidList,mfid]
  endfor
  ;maskData
  envi_open_file,maskFile, r_fid=fid
  maskData = envi_get_data(fid=fid,dims=dims,pos=0)
  classNUMs=HISTOGRAM(maskData,min=min(maskData))
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(maskData);maskdata是byte格式，所以binsize为1，现在变成DN值

  spaticalMat=[]
  foreach i,classNUMs do begin
    dataI = maskData eq i
    ;sumI=total(dataI)
    this=[]
    for j=0,n_elements(thresholdFile)-1 do begin
      dataJ = envi_get_data(fid=thFidList[j],dims=dims,pos=0)
      sumI=total(dataJ)
      dataJ=dataI+dataJ
      sumOL=total(dataJ eq 2)
      ratio=1.0*sumOL/sumI
      this=[this,ratio]
      print,format='($,a)','^'
    endfor
    print,i
    spaticalMat=[[spaticalMat],[this]]
  endforeach
  openW,lun,Log,/GET_LUN;,/append
  f=strcompress('('+String(fix(N_ELEMENTS(thresholdFile)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'history similiarty:'
  printF,lun,format=f,spaticalMat;————保存相似矩阵
  FREE_LUN,lun
  ;移除thFidList
  foreach ifid, thFidList do envi_file_mng,id = ifid,/remove

  return,1; spaticalMat
end



;+
; :Description:
;多个指数的计算
;-f2='C:\temp\b119029_flassh.dat'
function Bmath,fname
  COMPILE_OPT idl2
  ;ENVI,/RESTORE_BASE_SAVE_FILES
  ;ENVI_BATCH_INIT
  ;r_Pos = STRPOS(fname,'_')
  ;prefix=STRMID(fname,0,r_Pos)
  e=envi()
  ;prefix=e.GETPREFERENCE('temporary_directory')
  ;postfix='.dat'
  ENVI_OPEN_FILE,fname,r_fid=fid
  ENVI_FILE_QUERY,fid,dims=dims,ns=ns,nl=nl,nb=nb,data_type=dt

  fid=intarr(nb)+fid
  pos=INDGEN(nb)
  express=['1.0*(b1-(b2+b3+b4))/(b1+b2+b3+b4)',$;NWI
    '1.0*(b1-b2)/(b1+b2)',$;NDWI
    '1.0*(b1-b2)/(b1+b2)',$;NDVI
    '1.0*(b1-b2)/(b1+b2)',$;MNDWI
    '1.0*(b1-b2-b3)/(b1+b2+b3)'];EWI

  out_bname=['NWI','NDWI','NDVI','MNDWI','EWI']
  posAll=list()
  posAll.add,[1,4,5,6];NWI
  posAll.add,[2,4];NDWI
  posAll.add,[4,3];NDVI
  posAll.add,[2,5];MNDWI
  posAll.add,[2,4,5];EWI


  outList=[]

  proj=0
  FOR i = 0L, n_elements(express)-1 DO BEGIN
    tempFile = e.GetTemporaryFilename();prefix + OUT_BNAME[i] + postfix
    tfid=fid[posAll[i]]
    tpos=pos[posAll[i]]
    ENVI_DOIT,'math_doit',dims=dims,pos=tpos,fid=tfid,$
      exp=express[i],r_fid=r_fid,OUT_BNAME=OUT_BNAME[i],$
      out_name=tempFile;
    outList=[outList,tempFile]
    proj=ENVI_GET_PROJECTION(fid=r_fid)
  ENDFOR

  tempFile=stackingExc(outList,proj,out_bname)
  ;留意这个删除命令能否运行
  rfidList=[];按outList顺序得到fid代号
  for times=0,n_elements(outList)-1 do begin
    envi_open_file,outList[times], r_fid=fid
    rfidList=[rfidList,fid];一列fid文件
  endfor
  foreach ifid,rfidList  do envi_file_mng,id = ifid,/delete
  return, tempFile
end

;+
; :class the image according locate:
;
;-
function ClassImg,clusterF,locate
  compile_opt idl2
  envi_open_file, clusterF, r_fid=cfid
  envi_file_query,cfid,nl=nl,ns=ns,dims=dims
  clusterD=envi_get_data(fid=cfid,dims=dims,pos=0)
  cluD=0;用来保存影像
  e = ENVI()
  tempFile = e.GetTemporaryFilename();
  for i=0,n_elements(locate)-1 do begin
    if (size(cluD))[0] eq 0 then cluD=(clusterD eq i)*locate[i]
    if (size(cluD))[0] ne 0 then cluD+=(clusterD eq i)*locate[i]
  endfor

  newRaster = ENVIRaster(cluD, URI=tempFile)
  newRaster.Save
  newRaster.close
  return, tempFile
end



;+
; :Description:有时候需要用到一些辅助数据，如DEM，甚至云或者城市等等。
; 那么这个有必要了。
; 方法就是将bri_test中得到的MTKL数据（大簇）掩膜处理。
; newFile 是bri_test得到的数据路径名称
;-
function supplyToMTKL,supply,newFile,path
  compile_opt idl2
  envi_open_file, supply, r_fid=fid
  if fid eq -1 then return,newFile;如果没有对应的掩膜文件，就是没必要有掩膜文件

  envi_file_query,fid,nb=nb,dims=dims
  maskClass = envi_get_data(fid=fid,dims=dims,pos=0)
  objRaster = ENVIFIDToRaster(fid)
  proj =objRaster.SPATIALREF
  Masked=[]
  validF=[]

  index=0;记录自然增长的下标
  trans=[];newfile可能有不能打开的，所以用它来记录下标
  for i=0,n_elements(newFile)-1 do begin
    envi_open_file,newFile[i] , r_fid=ifid
    if ifid eq -1 then continue
    trans=[trans,i]
    newData = envi_get_data(fid=ifid,dims=dims,pos=0)
    mark=newData + maskClass
    mark=1.0*total(mark eq 2)/total(newData)
    if mark gt 0.4 then begin
      saved = newData * maskClass
      newData -= saved
      toValid=saveTemple(newData,proj);用硬盘换取内存压力
      validFile=[toValid,string(index)]
      tempFile=saveTemple(saved,proj)
      saved=[]
      validF=[[validF],[validFile]]
      ;print,tempFile,i
    endif else tempFile=newFile[i]
    Masked=[Masked,tempFile]
    newData=[]
    envi_file_mng,id = ifid,/remove;,/delete
    index++
  endfor
  print,Masked,trans

  back=validSNOW(path,maskClass,validF)
  rubbish=[]
  rubbish=[rubbish,Masked[back]]
  rubbish=[rubbish,transpose(validF[0,*])]
  newBack=trans[back];唯一一个需要转换的下标接口，和newfile相关联的地方
  if n_elements(back) gt 0 then Masked[back]=newFile[newBack]

  for i=0,n_elements(rubbish)-1 do begin
    envi_open_file, rubbish[i], r_fid=rfid
    envi_file_mng,id = rfid,/remove,/delete
  endfor
  return, Masked
end


;+
; :Description:
; 只用来辅助supplyTOMTKL函数，maskclass就是supply的影像数据
; 用来验证被切割掉的，即tovalid数据是否相关于积雪
;-
function validSNOW,path,maskClass,toValid
  compile_opt idl2
  envi_open_file, path, r_fid=fid
  envi_file_query,fid,nb=nb,dims=dims
  S3 = envi_get_data(fid=fid,dims=dims,pos=0);就是s3积雪指数
  S3=S3 gt 0.12;阈值化，这个数值可以小些，反正一会要被supply切割，
  S3 *= (1-maskClass)
  ;maskClass=[]
  result=[]
  for i=0,n_elements(toValid[0,*])-1 do begin
    envi_open_file,toValid[0,i] , r_fid=ifid
    if ifid eq -1 then continue
    newData = envi_get_data(fid=ifid,dims=dims,pos=0)
    mark=newData + S3
    mark=1.0*total(mark eq 2)/total(newData)
    if total(newData) eq 0 then mark=0.0
    print,mark,i
    if mark lt 0.4 then result=[result,fix(toValid[1,i])]
  endfor

  return, result
end








