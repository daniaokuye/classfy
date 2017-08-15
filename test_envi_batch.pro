;+
;所有参数都是临时用来测试用的
;-outterRelate是影像file长度-1的维度，是相邻两景聚类号对应关系的list列表
;matrixList也是列表，每一项包括两个相关的矩阵，是空间关系评价的矩阵
;allList是一个简单数据结构，因为同时只会用到上一幅和这一幅两幅影像
;注意！！！！！！spaticalDraw中有一个临时sav文件，目的是为了防止程序崩溃时，
;   继续上次运行点继续运行,savFile
;files是分类isodata图像
;
pro test_envi_Batch,file,maskFile,baseFile;,fidList,outterRelate,matrixList

  ;baseFile='D:\360Downloads\JULY17_120\log\'
  Log = baseFile+'log_testbatch.txt'
  savFile=baseFile+'temp.sav';Draw()函数的临时文件
  mtkl=baseFile+'mtkl.sav';保存mtkl file列表文件
  outter=baseFile+'outter.sav';保存相互关系outterRelate文件
  newF=baseFile+'newFile.sav'
  
  openW,lun,Log,/GET_LUN;,/append
  file=judgeClassNUm(file)
  ;开始for循环
  now=0
  ;fidList=[];按file顺序得到fid代号
  allLIst=list();数据结构，循环列表，表长为2，now & pre
  ;注意利用outterRelate的维度，它的维度是file-1。里面的维度是倒数n-1个file影像的聚类号数
  outterRelate=list();相关类（聚类号），list文件，包含n个list文件，n=n_elements(file)
  matrixList=list();包含多个内部的list，每个list包括coeffient、valueable;--A
  for i=0,1 do allList.add,[]
  for times=0,n_elements(file)-1 do begin
    clusterF=file[times]
    thisList=readData(clusterF)
    ;fidList=[fidList,fid];一列fid文件
    allList[now]=thisList
    ;allList.add,thisList;仅测试用
    if times lt 1 then begin
      now++
      continue
    endif

    resultMatrix=portionCal(now,allList);--B
    matrixList.add,resultMatrix;--C
    ;resultMatrix=matrixList[times-1];测试用，直接跳过matrix--对应行[A/B/C]注释掉
    innerRelate=relateC(now,allList,resultMatrix)
    outterRelate.add,innerRelate
    ok=saveMat(now,resultMatrix,allList,lun)
    print,'--test_envi_batch loop at -->',times

    now++
    now=now mod 2
  endfor
  FREE_LUN,lun


  ;print,'fidlist:'
  ;print,fidList,'outterRelate'
  foreach outR,outterRelate do print,outR,'innerRelate'
  mtklFile = spaticalDraw(file,outterRelate,savFile)
  print,'fileList:',mtklFile
  save,filename=mtkl,mtklFile
  save,filename=outter,outterRelate
  print,'start britest.pro'

  ;该方法得到聚类成“大簇”的文件，排列依旧按照mtklFile顺序,newFile为输出文件（簇）
  bri_test,outterRelate,file,mtklFile,maskFile,newFile,baseFile
  save,filename=newF,newFile
  
  ;得到分类图
  ;calentroy,newFile,maskFile
  print,'ok_test_envi_batch'

end



;+
; :Description:
; isodata聚类后的类别数目来判断是否需要这一影像参与mtkl
; 太少没有什么意义
;-
function judgeClassNUm,file
  compile_opt idl2
  numsAll=[];对应的类别数

  for i=0,n_elements(file)-1 do begin
    envi_open_file, file[i], r_fid=fid
    envi_file_query,fid,nb=nb,dims=dims
    maskClass = envi_get_data(fid=fid,dims=dims,pos=0)
    classNUMs=numsClass(maskClass)
    numsAll=[numsAll,n_elements(classNUMs)]
    envi_file_mng,id = fid,/remove
  endfor

  ;print,numsAll
  suitable=where(numsAll ge 15)
  return, file[suitable]
end








;将聚类类号，数据，文件名写入数据结构list中去
function readData,clusterF
  compile_opt idl2

  envi_open_file, clusterF, r_fid=cfid
  IF cfid EQ -1 THEN return,1
  envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
  data=envi_get_data(fid=cfid,dims=dims,pos=0)

  ;print,'now:',now,'n:',n_elements(allList)
  thisList=list()

  ;聚类类别，有DN值的跳跃histogram也可表示出来
  classNUMs=HISTOGRAM(data,min=min(data));这是关键，按照min-max排列为一维，无值即为0--直方图
  classNUMs=where(classNUMs gt 0);变成了序数,空白类255变成序数
  classNUMs=classNUMs+min(data);maskdata是byte格式，所以binsize为1，现在变成DN值

  thisList.add,data
  thisList.add,classNUMs
  thisList.add,clusterF
  data=[]
  envi_file_mng,id = cfid,/remove
  return, thisList
end



;+
; :calculation the portion of overlay over both classied image:
;输出相应的占比矩阵和自用率矩阵
;-
function portionCal,now,allList
  compile_opt idl2
  pre=(now + 1) mod 2;allList是循环列表,包含两项
  classNUMs = allList[now,1];列
  maskClass = allList[now,0]
  ;allList[now,0]=[];现在的不能置0，下次还要用
  clusterNUMs = allList[pre,1];行
  cluster = allList[pre,0]
  ;allList[pre,0]=[]

  print,'pre:',pre,' now:',now
  coeffient=[];从属度 same/older
  valueable=[]
  foreach clus,clusterNUMs do begin;行数

    clusterData=cluster eq clus
    countClu=total(clusterData)
    cTemp=[]
    vTemp=[]
    countCls=0
    foreach cls,classNUMs  do begin;列数
      maskData=maskClass eq cls
      countCls=total(maskData)
      new=clusterData+maskData
      countNew=total(new eq 2)
      new=[]
      maskData=[]
      cTemp=[cTemp,1.0*countNew/countClu];
      vTemp=[vTemp,1.0*countNew/countCls]
      print,format='($,a)',' : '

    endforeach
    coeffient=[[coeffient],[cTemp]]
    valueable=[[valueable],[vTemp]]
    clusterData=[]
    print,' clusterNUMs:',clus
  endforeach
  result=list()
  result.add,coeffient
  result.add,valueable
  return, result
end



;+
; :calculate the related cluster number:
;
;-
function relateC,now,allList,resultMatrix
  compile_opt idl2
  coeffient=resultMatrix[0]
  valueable=resultMatrix[1]
  pre=(now + 1) mod 2;allList是循环列表,包含两项
  classN = allList[now,1];列
  clusterN = allList[pre,1];行
  ;之所以用indgen，序号可能是不连续的。尤其是通过readData函数读入的。源自于isodata可能某项没有值
  classNUMs = indgen(n_elements(classN))
  clusterNUMs = indgen(n_elements(clusterN))
  ;obtain the relative cluster number
  innerRelate=list()
  ;行数clus,clusterNUMs,现在用列号了（优化）

  foreach cls,classNUMs do begin
    relate=[]
    ;符合条件的行号
    c=clusterNUMs[where(coeffient[cls,*] gt 0.1)];条件一
    portionC=coeffient[cls,c]*valueable[cls,c];
    v=clusterNUMs[where(valueable[cls,*] gt 0.1)];条件二
    portionV=coeffient[cls,v]*valueable[cls,v];
    if portionC ne [] or portionV ne [] then begin;条件三
      cI=c[where(portionC gt 0.1)];加入relate用，但先调到clusterN下
      vI=v[where(portionV gt 0.1)];同上
      if portionC ne [] then relate=[relate,clusterN[cI]]
      if portionV ne [] then relate=[relate,clusterN[vI]]
    endif else begin
      relate=[relate,clusterN[v]]
      if relate eq [] then relate=[relate,clusterN[c]]
    endelse
    relate=relate[sort(relate)]
    relate=relate[uniq(relate)]
    innerRelate.add,relate
  endforeach

  return, innerRelate
end


function saveMat,now,resultMatrix,allList,lun
  compile_opt idl2
  pre=(now + 1) mod 2;allList是循环列表,包含两项
  classNUMs = allList[now,1];列
  clusterNUMs = allList[pre,1];行
  ;Save the matrix of coeffient & valuable
  coeffient=resultMatrix[0]
  valueable=resultMatrix[1]
  f=strcompress('('+String(fix(N_ELEMENTS(classNUMs)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'col is --',allList[now,2],' nums is ',max(classNUMs)
  PRINTF,lun,'line is --',allList[pre,2],' nums is ',max(clusterNUMs)

  printF,lun,'coeffient'
  PRINTF,lun,format=f,coeffient
  printF,lun,'valueable'
  PRINTF,lun,format=f,valueable
  return,1
end

;
;关联空间位置。innerRelate的长度是clusterNUMs的列表list
;savFile='c:\temp\temp.sav'
;innerRelate
function spaticalDraw,file,outterRelate,savFile;,matrixList
  compile_opt idl2
  e = ENVI()
  fileList=[];保存相关聚类的文件
  experienceHistory=[];保存经历的文件，记录对应行号变化
  start=0;用来读取save的临时文件
  valid = 1
  while valid do begin
    ON_IOERROR, not_exist
    restore,savFile
    start=str.now
    fileList=str.fileList
    valid = 0
    not_exist:if valid ne 0 then begin
      str={now:0,fileList:0}
      save,filename=savFile,str
    endif
  endwhile
  if start eq 0 then fileList=[]

  for now=start,n_elements(file)-1 do begin
    envi_open_file,file[now], r_fid=cfid
    envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
    classData=envi_get_data(fid=cfid,dims=dims,pos=0)
    if now ne n_elements(file)-1 then $
      innerRelate=outterRelate[now] $
    else innerRelate=indgen(n_elements(outterRelate[now-1]))

    clusterData = [];没一类都会做存储
    clus=0;暂默认所有的聚类号是从0开始的序列
    NewFile=strarr(n_elements(innerRelate))

    foreach line,innerRelate do begin;行数
      maskData=[]
      foreach cls,line do begin;没有重合，所有最大值不超过1
        if maskData eq [] then maskData=(classData eq cls) else $
          maskData += (classData eq cls);所有的cls都同属于classData的互斥类
        print,format='($,a)','^'
      endforeach

      ;添加进来继承的数据
      if (n_elements(fileList) gt cls)  then begin;fileList have cotents
        if (fileList[cls] ne '') then begin;this one is valuable
          envi_open_file,fileList[cls], r_fid=tfid
          maskData += envi_get_data(fid=tfid,dims=dims,pos=0)
          ;raster = ENVIFIDToRaster(tfid)
          ;raster.close
        endif
      endif
      print,format='($,a)',fix(max(maskData))
      tempFile = e.GetTemporaryFilename()
      raster = ENVIRaster(maskData, URI=tempFile);clusterData
      raster.save
      raster.close
      maskData=[]
      NewFile[clus]=tempFile;重新给一个新文件

      print,' -->',clus
      clus++;随line对应行号变化
    endforeach
    delfile=fileList
    fileList=NewFile
    ;存储为临时文件
    str={now:now,fileList:fileList}
    save,filename=savFile,str
    print,'index of fid: ...',now+1
    foreach delF,delfile do begin
      envi_open_file,delF, r_fid=tfid
      ok=deleteFID(tfid)
      ;envi_file_mng,id = ifid,/remove;, /DELETE
    endforeach
    classData=[]
    envi_file_mng,id = cfid,/remove
  endfor
  ;print,fileList
  return, fileList
end


;+
; :the history of construction:
;倒序方式，分别列入,filenum是最后一组的一个聚类号，取值[0,max]
;对应于文件列表file。fileNum是最后一个文件的聚类号，而输出的“结构树”
;是剩下文件对应该聚类号的相应聚类号，倒序输出
;-
function hisExperience,outterRelate,fileNum
  compile_opt idl2
  ;traceFile=list();轨迹标记;外
  traceNum=list([]);中--先放一组数组，表示下标
  maxLimit=n_elements(outterRelate[-1])-1
  ;limit='['+strtrim(string(fix(0)),2)+','+strtrim(string(fix(maxLimit)),2)+']'
  if (fileNum lt 0) or (fileNum gt maxLimit) then return,0
  ;  begin
  ;    print,'input a suitable limit which is ' + limit
  ;    return,0
  ;  endif


  exist=[fileNum]
  for j=0,n_elements(outterRelate)-1 do begin;聚类影像的个数，8
    num=n_elements(outterRelate)-1-j
    this=list(num+1);“结构树”
    this.add,outterRelate[num,exist]

    exist=[]
    foreach next,this[1:-1] do $
      foreach nx,next do exist=[exist,nx]
    exist = exist[sort(exist)]
    exist = exist[uniq(exist)]
    traceNum[0]=[traceNum[0],num]
    traceNum.add,exist;this
  endfor
  ;traceFile.add,traceNum

  return, traceNum
end





