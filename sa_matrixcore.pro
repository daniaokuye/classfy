;+
; :Author: DN
;构建核心的彼此间相似度矩阵
;hardcore是核心，使用SAM算法是偷巧的方式，也可以用cos计算余弦相似，二者完全等价
;值越小，相似度越高，（未cos取值） 2017/6/13
;有价值的暂时为sa_matrixCore
;-
;
pro sa_matrixCore,hardCore,SIDmatrix,baseFile;realloc,
  compile_opt IDL2
  ;  RESOLVE_ROUTINE, ['relocate','cluster'],$
  ;    /IS_FUNCTION, /NO_RECOMPILE
  ;hardCore=hardCore[0:3,*]
  if hardCore eq !NULL then setbackCore,logFile,reSave,hardCore
  out_name = baseFile+'sa_testimg.tif';这俩是sid算法的输出项，临时之用
  ;ruleName = 'D:\360Downloads\June\testRule.tif'
  reSave= baseFile+'sa_matrix.txt'
  ;  envi_doit, 'class_doit', fid=fid, pos=0, dims=dims, $
  ;    out_bname='SID', out_name=out_name, method=3, $
  ;    mean=mean, r_fid=r_fid,rule_out_name=ruleName, $
  ;    lookup=lookup, class_names=class_names, $
  ;    in_memory=0, thresh=thresh
  ;OMasked = total(histstat,1)/(total(histstat gt 0 , 1));用这个值代替histstat中的0值
  ;  inout='D:\test\allInside.txt'
  ;  fs=fileread(inout);
  ;  hardcore=make_array(5,n_elements(fs),/DOUBLE)
  ;  for n=0L,n_elements(fs)-1 DO BEGIN
  ;    rbool=StringToDoubleArray(fs[n],dataArray,Count,Pos);
  ;    hardcore[*,n] = dataArray
  ;  endfor
  ;  fs=[]


  shape = size(hardCore,/DIMENSIONS);二维数组变成3维
  OM=[]
  cols=intarr(shape[1])+1
  if N_ELEMENTS(shape) gt 1 then sizeOM=shape[0]-1 else sizeOM=0
  for i=0, sizeOM do begin
    temp=TRANSPOSE(hardCore[i,*])##cols
    OM=[[[OM]],[[temp]]]
  endfor
  line=(size(hardCore,/DIMENSIONS))[1]

  ;  matrix=[]
  ;  for i=0,line-1 do begin
  ;    r = SAM(OM,hardCore[*,i])
  ;    matrix=[matrix,cos(r[0,*])]
  ;  endfor
  ;;SID
  e = ENVI()
  tempFile = e.GetTemporaryFilename()
  newRaster = ENVIRaster(OM, URI=tempFile)
  newRaster.Save
  fid = ENVIRasterToFID(newRaster)
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  ;thresh=replicate(0.5,num_classes)

  ;method is a classifiction method-SAM. spectral angle match
  envi_doit, 'class_doit', fid=fid, pos=indgen(nb),$
    dims=[-1,0,line-1,0,line-1], $
    out_bname='SID', method=3, r_fid=ofid, out_name=out_name,$;
    mean=hardCore, rule_fid=r_fid,class_names=indgen(line+1), $;
    ;lookup= bytarr(3,line+1),  $,rule_out_name=ruleName
    lookup= byte(randomu(1,[3,line+1])*255),rule_in_memory=1;, thresh=thresh
  ;注意：rule_out_name和rule_in_memory都设置的情况下，谁在前谁优先，互相冲突


  SIDmatrix=list()
  envi_file_query, r_fid, nl=nl,ns=ns, nb=nb
  for i=0,nb-1 do begin
    new = envi_get_data(fid = r_fid,dims =[-1,0,0,0,nl-1] ,pos = i)
    SIDmatrix.add,new
  endfor
  SIDmatrix=SIDmatrix.ToArray();多一个维度
  SIDmatrix=TRANSPOSE(SIDmatrix,[0,2,1])
  ;;SID

  openW,lun,reSave,/GET_LUN,/append
  PRINTF,lun,'\n1.this is FileData one\n'
  f=strcompress('('+String(line)+'(f5.2,:,","))',/REMOVE_ALL)
  ;print,'ffffffff:::',f
  ;PRINTF,lun,format=f,acos(matrix)
  PRINTF,lun,'SIDmatrix'
  PRINTF,lun,format=f,SIDmatrix;'(13(f5.2,:,","))'
  ;PRINTF,lun,'want a girl friend'
  ;PRINTF,lun,format=f,SIDmatrix-acos(matrix)

  ;  realloc=CLUSTER(SIDmatrix);relocate
  ;  foreach b,indgen(realloc.count()) do begin
  ;    PRINTF,lun,realloc[b]
  ;    ;print,realloc[b]
  ;  endforeach

  ;
  newRaster.close
  FILE_DELETE, tempFile

  FREE_LUN,lun

  ;移除r_fid
  foreach iFid,[ofid,r_fid] do begin
    raster = ENVIFIDToRaster(iFid)
    raster.close
  endforeach
  print,'matrixCore is over'

end


;+
; :function:用来寻找新旧核间的对应，range规定了新旧核
;-;matrix为acos，角度的矩阵--返回的为一列表list
function relocate,matrix,range
  range=indgen(6)
  sz=size(matrix,/DIMENSIONS)
  left=indgen(sz[1])
  left[range]=-1
  left = left[where(left ne -1)]

  ;寻找新旧核间的对应关系
  location=list()
  loc=[]
  foreach i, range do begin
    minI=byte(where(matrix[i,*] eq min(matrix[i,left])))
    loc=[loc, minI]
  endforeach
  location.add,loc
  loc=[]
  ;新旧核对应引导到大类的对应
  foreach i, range do begin
    matrix[i,i]=1
    minI=byte(where(matrix[i,*] eq min(matrix[i,range])))
    loc=[loc, minI]
  endforeach
  location.add,loc
  return, location
end

;+
; :function:找到所有核的相似，无相似则不在其列
;--返回的为一列表list
;realloc是重新分配类别的核心代号；
;-
function cluster,matrix
  compile_opt idl2
  thresh=matrix le 0.2
  sz=size(matrix,/DIMENSIONS)
  L=LIST()
  done=LIST()
  for i=0,sz[0]-1 do begin
    if done.Where(i) ne !NULL then continue
    x=where(thresh[i,*] eq 1,count);x是他自己包含的小弟序号
    temp=LIST()
    foreach a,x do temp.add,a
    if count gt 1 then begin
      foreach a,x do begin;对小弟的处理
        compare = 0
        if done.Where(a)  ne !NULL then begin
          ;temp基本上等于x，但排除了已经被挖走的
          ;挖走的顺序等于排列的先后顺序，这样好吗？/|/
          temp.remove,temp.where(a)
          CONTINUE
        endif
        if a eq i then continue
        y=where(thresh[a,*] eq 1,count);小弟a包含的朋友y序号
        foreach b,temp do begin;看是舍弃小弟的朋友，还是舍弃小弟/|/解决
          if where(y eq b) eq -1 then compare -=1 else compare +=1
        endforeach
        if compare lt 2 then temp.remove,temp.where(a)
      endforeach
    endif
    if temp.ToArray() ne !NULL then L.add,temp.ToArray()+1;此处有加1，matrix从0开始计数
    foreach a,temp do done.add,a
  endfor
  ;PRINT,L
  newL=bytarr(N_ELEMENTS(L))
  foreach a,L do newL[L.where(a)]=N_ELEMENTS(a)
  ;print,newL
  L.Remove,where(newL lt 2)

  return, L
end

