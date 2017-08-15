;+
; :Author: DN
;FileData -- 每隔5行（波段数）自定义一行：以区分大类（类别）如4000和小类[1-15]（亚类）
;             格式分别为：1（一个数字），数字，位置（0或6）
;category -- 类别，亚类，像元数，最小值，所在行。共5个数据
;available -- 每大类中数目像元超过最多者一半以上的亚类，注意其行数等于useful数据项个数
;core -- 影像各聚类中心组成波谱曲线;extent -- [min,max]
;alist -- 创建available的序列列表 0,1,2,3,4……
;alist -- 2107/6/13-alist的第一个数字变成大类了，如1000，1001
;partial -- 用来确定全部保留0还是保留部分1
;-
pro setbackCore,core,std,extent,alist,IsPartial,logFile,reSave
  compile_opt IDL2
  if logFile eq !NULL then logFile = 'D:\360Downloads\June25\Log_1.txt'
  if reSave eq !NULL then  reSave= 'D:\360Downloads\June\FileData2.txt'
  if IsPartial eq !NULL then IsPartial=1
  FileData=fileToMat(logFile,category,classes,numsOfSec)

  available=intarr(numsOfSec+1,N_ELEMENTS(classes))-1;创建类别相同的行数，列数为最大的亚类数+1，多一列为保存大类号
  ;print,'class:',classes
  for i=0,N_ELEMENTS(classes)-1 do begin
    ;if classes[i] mod 1000 eq 0 then CONTINUE;为了去除0所代表的背景值
    all = where(category[0,*] eq classes[i]);序数
    ;这是经验方法，如果某聚类- 像元数- 不足15k，应该抛弃
    if max(category[2,all]) lt 10000 then continue

    sort=sort(category[2,all]);面积/像元数 序数 in ascending order
    ;计算可保留序数
    keep=[classes[i],0];两维，第一个类别（大类，小类），第二个是像元数
    for j=N_ELEMENTS(sort)-1,0,-1 do begin
      line = all[sort[j]];所在行
      keep=[[keep],[category[1:2,line]]]
    endfor
    if IsPartial eq 1 then keep[1,*]/=keep[1,1];面积必须是最大的、有效的面积的一半以上。0.5,且读入时已为浮点数
    if IsPartial eq 0 then keep[1,1:-1] =1.0;即全部保留
    leaves=where(keep[1,*] gt 0.5);暂定0.5
    available[0:leaves[-1],i]=keep[0,0:leaves[-1]];多一列keep[0]保存着大类代号
  endfor
  useful=where(available[0,*] gt -1);意欲去掉无值的类
  available=available[*,useful]

  nb=category[-1,1]-category[-1,0];波段数目nb

  core=[];影像各聚类中心组成波谱曲线
  extent=[];
  std=[];
  for i=0,N_ELEMENTS(available[0,*])-1 do begin;(size(available,/DIMENSIONS))[1]
    gtZero = where(available[*,i] gt -1,count);j从1开始的
    for j=1,max(gtZero) do begin
      ;将所有有效亚类罗列起来，筛选结果放入indexLine
      indexLine = where(category[0,*] eq available[0,i] and category[1,*] eq available[j,i])
      line = category[-1,indexLine];fileData中对应行号
      core=[[core],[TRANSPOSE(FileData[2,line+1:line+nb-1])]];先用mean来代替core
      extent=[[extent],[TRANSPOSE(FileData[1,line+1:line+nb-1])]]; [[],[]]max
      extent=[[extent],[TRANSPOSE(FileData[0,line+1:line+nb-1])]];min
      std=[[std],[TRANSPOSE(FileData[3,line+1:line+nb-1])]]
    endfor
  endfor

  ;创建available的序列列表 0,1,2,3,4……
  alist=list()
  number=0
  for i=0,(size(available))[2]-1 do begin
    av=where(available[*,i] gt -1,count)
    ;if count then begin
    ;print ,format='($,a)',count
    alist.Add,[available[0,i],indgen(count-1)+number]
    number+=(count-1)
    ;endif
  endfor

  ;  ;统计表格的重新读入与编码
  openW,lun,reSave,/GET_LUN
  PRINTF,lun,'\n1.this is FileData one\n'
  PRINTF,lun,format='(6(g,:,","))',FileData


  ;所有亚类及其参数
  PRINTF,lun,'\n2.this is category one\n'
  PRINTF,lun,format='(5(g,:,","))',category


  PRINTF,lun,'\n3.this is available one\t'
  f=strcompress('('+String(numsOfSec+1)+'(i,:,","))',/REMOVE_ALL)
  PRINTF,lun,format=f,available;'(15(i,:,","))'


  PRINTF,lun,'\n4.this is return one\t','core(mean)'
  f=strcompress('('+String(fix(nb-1))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,format=f,core;'(7(g,:,","))'

  printf,lun,'\n5.extent(min-max)'
  PRINTF,lun,format=f,extent

  printf,lun,'\n6.std'
  PRINTF,lun,format=f,std

  FREE_LUN,lun
  print,"statistic's over"
end

;+
; :file turn to matirx:
; classes:大类数目、类别数目
; numOfSec：度，最大的亚类数
;-
function fileToMat,logfile,category,classes,numsOfSec
  compile_opt idl2
  nLines=file_lines(logFile) ;获取行
  ;  打开文件
  fs=fileread(logFile);
  FileData=make_array(6,nLines,/DOUBLE);将报表去除文字，压缩成数字格式
  index=0

  if (n_elements(fs) gt 0) then begin
    for n=0L,n_elements(fs)-1 DO BEGIN;n_elements(fs)-1 ,这儿可能需要的序列号很多
      LineStr=fs[n]
      rbool=StringToDoubleArray(LineStr,DoubleArray,Count,Pos);
      ;print,n,'ol',count
      ;if count gt 0 then print,'doub',DoubleArray
      ;if count eq 1 then print,'os:',pos
      if Count eq 1 then begin
        if pos eq 0 then coli=0 else coli=3
        FileData[coli:coli+2,index]=[1,DoubleArray,Pos]
        if coli eq 3 then index+=1
      endif else if Count gt 1 then begin
        FileData[*,index]=DoubleArray
        index+=1
      endif

    endfor
  endif
  ;print,'index',index-1
  FileData =FileData[*,0:index-1]
  category=[];用来存储各个亚类的面积，以便提取有效亚类;;
  MainCat=0
  classes=[];;
  numsOfSec=0;最大的亚类数;;
  number=0;临时变量
  for i=0,index-1 do begin
    if ARRAY_EQUAL(FileData[[0,2],i],[1,0]) then begin
      MainCat=FileData[1,i]
      classes=[[classes],MainCat];大类
      number=0
    endif
    if ARRAY_EQUAL(FileData[[3,5],i],[1,6]) then begin
      ;类别，亚类，像元数，最小值，所在行。共5个数据
      secCategory=[MainCat,FileData[4,i],FileData[5,i+1],FileData[0,i+1],i]
      category=[[category],[secCategory]]
      number+=1
      if number gt numsOfSec then numsOfSec=number
    endif
  endfor

  return, FileData
end

;+
; :pick certain column of file:
;返回的格式为：大类，小类序号，meanS
;-
function pickColumn,logFile,x
  compile_opt idl2
  contents=['min','max','mean','steddev','mode','pixels']
  index=where(contents eq x)
  FileData=fileToMat(logFile,category,classes,numsOfSec)
  nb=category[-1,1]-category[-1,0]-1;波段数目nb
  certain=[]
  class=0
  sz=size(FileData)
  for i=0,sz[2]-1 do begin
    if ARRAY_EQUAL(FileData[[0,2],i],[1,0]) then begin
      class=FileData[1,i];1此时是大类
    endif
    if ARRAY_EQUAL(FileData[[3,5],i],[1,6]) then begin
      ;4此时是小类
      certain=[[certain],[class,FileData[4,i],TRANSPOSE(FileData[index,i+1:i+nb])]]
    endif
  endfor
  return, certain
end




;StringToDoubleArray
;==============================================================
Function StringToDoubleArray,DblStr,DoubleArray,Count,pos
  ;用法IDL>Status=StringToDoubleArray(DblStr,DoubleArray,Count)
  ;DblStr要转换的字符串,字符串可以是'1.23 ,3.4;0.3  E2.2 afd 3er.7 '
  ;DoubleArray保存double类型的数组
  ;Count总共可以提取多少个double类型的数
  ;返回是否成功的标志
  Count=0
  ;先对/t处理一下
  Stab=STRSPLIT(DblStr, STRING(9b),/EXTRACT)
  str=STRJOIN(Stab, ' ')
  Si=STRSPLIT(Str,'[;,: ]',/EXTRACT,/REGEX)
  ;----------------------------------------------------------------
  nn = n_elements(Si)
  ;----------------------------------------------------------------
  if (nn GT 0) then begin
    DoubleValue=findgen(nn)
    LineSize=UINDGEN(nn)
    ;------------------------------------
    for n=0L,nn-1 DO BEGIN
      ;先去掉首尾空格
      Si[n] = STRTRIM(Si[n],2)
      if Si[n] eq '-NaN' or Si[n] eq 'NaN' then begin
        print,Si[n]
        Si[n]='0.0'
      endif
      valid=IsDoubleString(Si[n])
      if (valid GT 0) then begin
        DoubleValue[n]=fix(Si[n],type=5)
        LineSize[n]=1;
      endif else begin
        DoubleValue[n]=0
        LineSize[n]=0;
      endelse
    endfor
    ;------------------------------------
    count=total(LineSize);
    index=where(LineSize GT 0)
    if count eq 1 then pos = index else pos =-1
    ;------------------------------------
    if (count GT 0) then begin
      DoubleArray=DoubleValue[index]
    endif
    ;------------------------------------
    return,1
    ;------------------------------------
  endif else begin
    return,0
  end
  ;----------------------------------------------------------------
END

FUNCTION  fileread,file
  ;file = DIALOG_PICKFILE(FILTER='*.*');也可以省略
  ;2010-08-28
  ;wxp07@qq.com
  LineCount = FILE_LINES(file);
  if (LineCount gt 0) then begin
    StringArray = strarr(LineCount);
    OPENR, unit, file, /GET_LUN
    READF, unit, StringArray
    FREE_LUN, unit
    FileString=StringArray
  endif  else begin
    FileString=''
    LineCount=0
  endelse
  RETURN,FileString
END

;用法IDL>Status=IsDoubleString(dblstr)
;Status=0或1,成功为1,否则为0
;以下识别正确的话,肯定可以使用:
;IDL>DoubleValue=fix(dblstr,type=5)转换为double数字

function IsDoubleString,str
  ;--------------------------------------------------------------
  Status=1;假设可以转换

  ;---------------------------------------------
  ;查找并去掉末尾的非法字符
  pos = STREGEX(str, '([^0-9.eE+-]|[+-.][Ee]|[eE].)')
  if pos GT -1 then begin
    str=strmid(str,0,pos)
  endif
  ;---------------------------------------------
  ;转换为ASCII
  inputstr = byte(str)
  ;获取字符个数
  nn = n_elements(inputstr)
  ;---------------------------------------------
  ;判断第一个字符的合法性
  if nn GT 0 then begin
    FirstDoubleStr=byte('+-1234567890.')
    index=where(FirstDoubleStr eq inputstr[0],count)
    if count eq 0 then Status=0
  endif else begin
    ;如果字符的长度小于1,也不是合法字符
    Status=0
  endelse
  ;---------------------------------------------
  ;判断第二个字符的合法性
  if nn GT 1 then begin
    SecondDoubleStr=byte('1234567890.Ee')
    index=where(SecondDoubleStr eq inputstr[1],count)
    if count eq 0 then Status=0
  endif
  ;---------------------------------------------
  ;查找字符串中数字的总数,不能小于1
  num_total=0
  NumberStr=byte('1234567890')
  for n=0L,nn-1 DO BEGIN
    count=0;
    index=where(NumberStr eq inputstr[n],count)
    num_total=count+num_total
  endfor
  if num_total LT 1 then Status=0
  ;=================================
  ;后续处理
  ;if Status eq 1 then begin
  ;print,fix(dblstr,type=5)
  ;endif
  ;=================================
  ;---------------------------------------------
  return,Status
End