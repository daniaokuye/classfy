;+
; :Author: DN
; IsPartial取值[0,1]，0值即在setbackCore中为全保留，否则1为参考面积保留
; 在test_batch中，意思是同一类的不可取相似。
; classFile是最后的分类图像
;-
pro test_batch,logFile,reSave,IsPartial,classFile
  ;  COMPILE_OPT IDL2
  ;  ENVI,/RESTORE_BASE_SAVE_FILES
  ;  ENVI_BATCH_INIT
  ;  path_c='D:\modis2014\M_HKM_mulband.tif'
  ;  ENVI_OPEN_FILE,path_c,r_fid = fid               ;envi函数1
  ;  mapinfo=ENVI_GET_MAP_INFO(fid=fid)              ;envi函数2
  ;  state_1 = ENVI_GET_PROJECTION(fid = fid)        ;envi函数声明1
  ;  state_2 = ENVI_PROJ_CREATE(/geo)                ;envi函数声明2
  ;  aa=prj(fid)
  ;  print,'mapinfo'
  ;  ENVI_BATCH_EXIT
  if IsPartial eq !NULL then IsPartial=1
  ;setbackCore,core,std,extent,alist,IsPartial
  setbackCore,core,std,extent,alist,IsPartial,logFile,reSave
  lenAlist=(size(alist))[1]

  ;-----------------classNum-----------------------
  ;classnum得到所有的大类类号1001，1002系列中单一的代号，使用遍历的方式
  classNum=[]
  foreach L,alist do begin
    classNum=[classNum,L[0] mod 1000]
  endforeach
  classNum=classNum[sort(classNum)]
  classNum=classNum[uniq(classNum)]
  ;print,'classNum',classNum

  ;----------------sameClass------------------------
  startAlist=0
  sameClass=list(); 将4区间同一类下的所在alist中的行号记下来
  foreach similar, classNum do begin
    tempList=[]
    ;similar = alist[startAlist,0] mod 1000
    for i=startAlist, lenAlist-1 do begin
      similarTemp = alist[i,0] mod 1000
      if similar eq similarTemp then begin
        templist=[templist,i]
        ;endAlist=i
      endif
    endfor
    startAlist++
    sameClass.add,tempList
  endforeach
  ;print,'sameClass',sameClass
  ;----------------teamCore------------------------
  ;sameClass=sameClass[0:1]
  teamCore=list();将samelist所对应alist的行号--的聚类号记下
  foreach iList,sameClass do begin
    team=[]
    foreach j,iList do begin
      team=[team,alist[j,1:-1]]
    endforeach
    teamCore.add,team
  endforeach
  ;print,'teamCore',teamCore
  ;--------------outterList--------------------------
  ;将聚类的mean值按teamCore中的行号，从core中取出
  matrix=0
  outterList=list();仍是代号，不过是内外双层包裹着代表core的代号
  ;
  ;alert:--注意，这儿打算将按大类聚的方式，改成全部点一块聚的方式。
  ;         即所有的都放到一块，而不是按teamCore的集合方式聚集
  ;+++++++++++++++++++++++++++++++{--
  ;foreach t,teamCore do begin
  t=[]
  foreach temp,teamCore do t=[t,temp]

  sa_matrixCore,core[*,t],matrix;matrix是相关矩阵，越小越相近
  ;print,size(t)
  ;print,size(matrix);确认一下它是不是矩阵
  s=(size(matrix))[1]
  ;sample=s;行列数相同

  ;用来保存matrix中相似的代号，以便生成核心时用到
  coreList=list()
  full=[]
  for i=0,s-1 do begin
    if total(i eq full) ge 1 then continue
    coreN=where(matrix[*,i] lt 0.1,count);cos or sam,0 is the most similar one

    ;去掉coreN中已经存在full中的系数
    tempCN=[]
    foreach c,coreN do begin
      if total(c eq full) ge 1 then continue
      tempCN=[tempCN,c]
    endforeach
    coreN=tempCN

    if count ne 0 then begin
      full=[full,coreN]
      ;对coreN排序，最好能按相识度差异大小抽选??
      coreList.add,coreN
      full=full[sort(full)];full为列号集合，满序为indgen（s），现在排序
      full=full[uniq(full)];取唯一值
    endif
    if (size(full))[1] eq s then break
  endfor
  ;print,'coreList',coreList,'over'
  ;---
  ;获得相似的核心，先暂时使用随机抽取吧
  ;使用双层list，方便确认大类
  innerList=list()
  foreach coreN,coreList do begin
    innerList.add,t[coreN]
  endforeach
  outterList.add,innerList
  ;endforeach
  ;;+++++++++++++++++++++++++++++++--}
  ;print,'outterList',outterList

  ;对是否是IsPartial处理
  if IsPartial eq 0 then begin
    OL_Index=0;用来表示outterList当前所在项
    foreach OL,outterList do begin;outterList

      IL_Index=0;用来表示tempList当前的行
      tempList=list();innerList的替代者,其第一行是非同大类需要合，需要裂的依次列在后面
      foreach IL,OL do begin;OL-->innerlist;IL-->里面的列表
        tempList.add,[]
        if n_elements(IL) eq 1 then begin
          tempList[IL_Index]=IL;如果只有一项就跳过
          IL_Index++
          continue
        endif
        ;如不是，需要检查里面有无alist中同一行的，同属于一大类的
        divide=list()
        foreach AL,alist do divide.add,[];initiate divide
        for i=0,n_elements(IL)-1 do begin;targets of IL
          al_Index=0;al对应系数
          foreach AL,alist do begin;scan AL
            ;AL[0]是大类号
            if total(IL[i] eq AL[1:-1]) then begin
              ;这样，属于同一alist中类的将以IL系数的形式，列在divide列表中（list）
              divide[al_Index]=[divide[al_Index],IL[i]]
              break;找到了就确定下一个
            endif
            al_Index++
          endforeach;scan AL
        endfor;targets of IL
        ;至此，此IL中所有的项都做了相应归并，需要对divide做处理了
        ;newIL=list()
        jump=0;tempList需要增长的长度
        foreach DL,divide do begin;DL
          if DL eq [] then continue;为空则跳过
          ;如果只有一项直接添加当前列表中
          if n_elements(DL) eq 1 then tempList[IL_Index]=[tempList[IL_Index],DL]
          if n_elements(DL) gt 1 then begin
            ;如果有多项，第一项加入当前列表，其余项顺次添加
            tempList[IL_Index]=[tempList[IL_Index],DL[0]]
            for i=1,n_elements(DL)-1 do begin
              tempList.add,[DL[i]]
              jump++
            endfor
          endif
        endforeach;DL
        IL_Index+=(jump+1)
      endforeach;OL-->innerlist;IL-->里面的列表
      outterList[OL_Index]=tempList
      OL_Index++
    endforeach;outterList
  endif;IsPartial


  ;----------------------------------------
  ;每次都取得一样行数的hardcore，即现时分类类别数
  ;每次根据相似程度，对有相似类别的聚类做随机更换，共迭代totaltime次
  totalTime=1;5
  ruleFileList=[];个数为totalTime
  for times=1,totalTime do begin;6times
    ;对每一大类(class)中相同聚类(cluster)取随机一个，非相同聚类汇成此大类的核，此核用于分类
    ;然后做一个循环
    hardcore=[];核
    clsNums=[];每大类下的聚类的数目
    foreach cls,outterList do begin;每一大类
      i=0
      foreach clus,cls do begin;大类下相同聚类
        nums = N_Elements(clus);randomu创建[0,nums)的浮点数
        sub = clus[FIX(nums * RANDOMU(var, 1))];得到[0,nums-1]的整数
        hardcore=[[hardcore],[sub]];core[*,sub]
        i++
      endforeach
      clsNums=[clsNums,i]
    endforeach

    print,hardcore,'size',size(hardcore)
    ;此时hardcore的维度是一样的，即类别数一样，用它作为分类中心，调用sam方法，得到一个rule影像
    ;rule影像维度和hardcore一样
    sa_readDBF,core[*,hardCore],class_name,ruleName;调用Sam分类方法
    ruleFileList=[ruleFileList,ruleName]
    print,'over'
  endfor
  ;----------------------------------------
  ;将每一类的影像取最小值
  if N_Elements(ruleFileList) gt 1 then begin
    classFile=minify(ruleFileList);类的boost
  endif else begin
    classFile=class_name
  endelse


end


;+
; :function:选择每一类中的最小值
;  这儿需要两个 循环，一个是遍历所有的fileList，另外一个是遍历所有类别
;  同时需要注意内存，所以每次遍历都要清空内存，并且将文件存为临时文件
;  所以每一步直接先取最值，然后保存，然后再比较，也要至少打开两幅影像
;-;
function minify,fileList
  nums=N_Elements(fileList)

  e = ENVI()


  classList=[];保存tempfile
  nb=0
  while 1 do begin
    band=0
    ruleImg=[];保存影像
    for file=0,nums-1 do begin
      path=fileList[file]
      envi_open_file, path, r_fid=fid;要保证所有fid都可打开，因为连着外面的nb
      envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
      ruleImg =[ [ruleImg],[envi_get_data(fid=fid,dims=dims,pos=band)]]
      if (size(ruleImg))[0] eq 3 then ruleImg=min(ruleImg, DIMENSION=3);对应size，第4个即下标3为z维度
      ;关闭这个fid
      raster = ENVIFIDToRaster(fid)
      raster.close
    endfor
    tempFile = e.GetTemporaryFilename();这个应该该保存
    Raster = ENVIRaster(ruleImg, URI=tempFile)
    Raster.Save
    Raster.close
    band++
    ;保存tempfile
    classList=[classList,tempFile]
    if band eq nb then break;band++在前
  endwhile

  print,classList

  ruleImg=[];保存影像
  ;对classList处理，得到影像图。
  for band=0,nb-1 do begin
    path=classList[band]
    envi_open_file, path, r_fid=fid;要保证所有fid都可打开，因为连着外面的nb
    envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
    ruleImg =[ [ruleImg],[envi_get_data(fid=fid,dims=dims,pos=band)]]
    if (size(ruleImg))[0] eq 3 then begin

      ruleImg=min(ruleImg, DIMENSION=3);
    endif
  endfor


  return,1
end



