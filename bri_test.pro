;+
; :首先，确定一个阈值，先定为4，即一半
; 然后，据此确定聚类的主要来源。newSim为输出文件
;-fidList是iisodata后的分类文件，比如有8个
;蒙特卡洛File是test_evni_batch.pro得到的大簇文件
;-
pro bri_test,outterRelate,file,mtklFile,maskFile,newFile,baseFile

  newSav = baseFile+'newout.sav';newout
  poeSav = baseFile+'poe.sav';poentenOut
  thSav = baseFile+'th.sav'; thresholdFile
  matSav = baseFile+'mat.sav';spaticalMat

  fidlist = fileToFid(file)
  newOut=list()
  poentenOut=list()
  thresholdFile=[];阈值化影像保存
  for fileNum=0,n_elements(mtklFile)-1 do begin;n_elements(mtklFile)-1
    construct = updateHis(fidList,mtklFile,outterRelate,fileNum,thresholdFile)
    if size(construct, /TYPE) ne 11 then continue
    newOut.add,construct[0]
    poentenOut.add,construct[1]
  endfor

  save,filename=newSav,newOut
  save,filename=poeSav,poentenOut
  save,filename=thSav,thresholdFile;first time to save the varianation

  s = closer(newOut,baseFile);similar这个展示出那些需要聚合成更大的簇(一)
  c = conflictFinder(newOut,s);conflictList
  thresholdFile=reConstruct(newOut,fidList,thresholdFile)
  spaticalMat=spaticalRelativ(thresholdFile,baseFile);这个展示出那些需要聚合成更大的簇(二)

  save,filename=thSav,thresholdFile;second time to save the varianation
  save,filename=matSav,spaticalMat

  newSim=relationship(s,spaticalMat)
  ;dealConflict（）调用函数的顺序：ob/sub,howdivide,depart,departExecute
  outFile=dealConflict(c,newSim,fidList,thresholdFile,newOut);spaticalMat,
  ;移除fidList
  foreach ifid, fidList do envi_file_mng,id = ifid,/remove
  print,'removeConflictFile:',outFile

  ;相似
  spaticalMat=spaticalRelativ(outFile,baseFile)
  newSim=relationship(s,spaticalMat)

  ;maskFile='D:\360Downloads\July12\wetSHP\TP_119029.tif'
  envi_open_file,maskFile, r_fid=mfid
  objRaster = ENVIFIDToRaster(mfid)
  proj =objRaster.SPATIALREF

  newFile=crowedImage(newSim,outFile,proj)
  newFile=newFile[sort(newFile)]
  newFile=newFile[uniq(newFile)]
  print,'ok--bri_test',newFile
end



;+
; :Description:
; file to fid
;-
function fileToFid,file
  compile_opt idl2
  fidList=[];按file顺序得到fid代号
  for times=0,n_elements(file)-1 do begin
    envi_open_file,file[times], r_fid=fid
    fidList=[fidList,fid];一列fid文件
  endfor

  return, fidList
end


;+
; :Description:
;利用去除冲突后的影像和相关关系建立“簇”，输出为文件名
;文件中没有有效影像的将剔除
;-
function crowedImage,newSim,outFile,proj
  compile_opt idl2
  newFile=outFile
  for i=0,n_elements(newSim)-1 do begin
    if newFile[i] ne outFile[i] then continue;已经处理过了
    gather=newSim[i]
    Data=[]
    gatherFid=[]
    foreach items,gather do begin
      envi_open_file,outFile[i], r_fid=fid
      envi_file_query,fid,dims=dims
      if Data eq [] then Data=envi_get_data(fid=fid,dims=dims,pos=0) else $
        Data += envi_get_data(fid=fid,dims=dims,pos=0)
      gatherFid=[gatherFid,fid]
    endforeach
    Data=Data gt 0
    if total(Data) lt 1000 then begin
      newFile[[gather]]=''
    endif else begin
      name='crowed_'+strtrim(string(i),2)
      tempFile=saveTemple(Data,proj,name);这个需要保持成为带投影的文件
      newFile[[gather]]=tempFile
    endelse
    if total(newFile[[gather]] ne outFile[[gather]])gt 0 then begin
      foreach ffid,gatherFid do ENVI_FILE_MNG,ID=ffid,/DELETE
    endif
    print,format='($,a)',' cI '
  endfor
  return, newFile
end



;+
; :Description:根据similar和spaticalMat，即光谱来源相关和空间位置相关，
; 共同创建相似关系。;
;-
function relationship,similar,spaticalMat
  compile_opt idl2
  newSim=list()
  for i=0,n_elements(similar)-1 do begin
    line=similar[i]
    if n_elements(line) eq 1 then newSim.add,[line];如果只有一个，跳过不管
    if n_elements(line) gt 1 then begin
      ;检测是否需要剔除,过小0.5
      this=[]
      foreach items,line do begin
        if spaticalMat[items,i] lt 0.5 or $
          spaticalMat[i,items] lt 0.5 then continue
        this=[this,items]
      endforeach
      ;检测是否有需要添加的，比较大，0.6
      newOnes=where(spaticalMat[*,i] gt 0.6)
      foreach items,newOnes do begin
        if total(items eq this) eq 0 then begin;不在this中存在的项
          ;all in this is satisfied,for it does not meet the demand of origin of isodata classified image
          sum=0
          foreach friens,this do begin;相关项中也-都-满足0.6以上的条件
            if spaticalMat[items,friens] gt 0.6 then sum++
          endforeach
          if sum eq n_elements(this) then this=[this,items]
        endif
        ;完成两部操作
      endforeach
      newSim.add,[this]
    endif
  endfor
  return, newSim
end




;+
; :Description:根据优选代号重新构造
;
;-
function reConstruct,newOut,fidList,thresholdFile
  compile_opt idl2
  fileArr=[];保存返回的文件列表
  n=n_elements(fidList)-1
  ;ns and nl for null output
  ;envi_file_query,fidList[0],ns=ns,nl=nl
  
  for i=0,n_elements(newOut)-1 do begin
    inOut=newOut[i]
    fileIndex=0
    outPut=0;作为指针，保存文件
    foreach items,inOut do begin
      ;clusterF=file[n-fileIndex]
      ;按照update优化了的结构，再重新计算ＭＴＫＬ数据
      ;envi_open_file, clusterF, r_fid=cfid
      cfid=fidList[n-fileIndex]
      IF cfid EQ -1 THEN return,1
      envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
      Data=envi_get_data(fid=cfid,dims=dims,pos=0)
      foreach parts,items do begin
        if n_elements(outPut) eq 1 then outPut = Data eq parts else $
          outPut += Data eq parts
      endforeach
      Data=[]
      print,format='($,a)','s'
      fileIndex++
    endforeach
    print,i
    ;计算合适的阈值
    thresData=thresholdFile[i]
    envi_open_file, thresData, r_fid=tfid
    Data=envi_get_data(fid=tfid,dims=dims,pos=0)
    threshold=1;测试合适阈值
    feature=[];各个阈值下的特征
    maxTh=ceil(n_elements(fidList)/2.0)
    if max(outPut) gt 1 then begin
      for threshold=1,max(outPut)-1 do begin
        tempData=outPut ge threshold
        overlay=(tempData+Data) eq 2
        ;以下两个指数同步动作,lapNow为上升趋势，lapPre为下降趋势--单调,目标是取得二者都接近1的情况
        lapPre=1.0*total(overlay)/total(Data)
        lapNow=1.0*total(overlay)/total(tempData)
        tempData=[]
        overlay=2-(lapPre+lapNow)
        feature=[feature,overlay]
        print,lapPre,lapNow,threshold
      endfor
      Data=[]
      envi_file_mng,id = tfid,/remove;移除
      threshold=(where(feature eq min(feature)))[0]+1;feature从1开始
      outPut=outPut ge threshold
    endif
    if max(outPut) eq 0 then  outPut=Data*0
    Data=[]
    ;保存阈值文件
    print,feature,threshold
    tempFile=saveTemple(outPut);data
    fileArr=[fileArr,tempFile]
  endfor
  ok=delTemp(thresholdFile)
  return,fileArr
end


;+
; :Description:judge it's object or subject
;1是对target拆分出来，0则不是跳过
;a/ this 和target不可能有较大重叠
;b/ this不可分，target不可分都可能出现
;c/ 所以，按照即便不可分的情况出现，按照先后顺序处理也是公平的：
;       结合a/b重叠必然是源于similiar出现同源。即冲突情况下，必有一个可分
;       拆分不过是为了让相交的部分保持独立
;即便确定拆分，也没必要调整newout，即判定相关、互斥的信息源之一。也没必要做循环：
;因为如果继续某项发生冲突，最终的效果是在面积上来看的。保持一致的是对target处理，因为target终究
;会变成主体，可以考虑对dealConflict中的target解冻
;-
function objectOrSub,simliarList,this,target;spaticalMat
  compile_opt idl2

  ;0:不对target拆/;1:object客体，需要对target拆分
  result=0
  former = n_elements(simliarList[this])
  latter = n_elements(simliarList[target])
  if former gt latter then  result=1 else result=0

  ;对于不可分情况，target和i显然也是没有交集的（从来历上看），有必要切除
  IsmultiPart=1
  foreach ones,simliarList[this] do begin
    len=n_elements(ones)
    if len gt 1 then begin
      IsmultiPart=1
      break
    endif
    IsmultiPart=0
  endforeach
  ;如果本身不可分，而目标判定为不去分的情况下，避免this未来被拆分，对目标拆分
  if IsmultiPart eq 0 and result eq 0 then result=1
  return,result
end



;+
; :Description:
;根据dealConflict得到的拆分依据，对影像拆分。
;howDivid是集合了所有需要拆分的项，集中考虑进行拆分，而不是发现一项拆分一项
;-
function dePart,fidList,thresholdFile,howDivid
  compile_opt idl2
  outFile=thresholdFile;先照抄，这是要输出的文件
  homeSub=[];needed所对应的主体文件号
  needed=[];用来记录需要被拆分的文件号(thresholdFile)
  delFile=[];要删除文件下标
  for i=0,n_elements(howDivid)-1 do begin
    homeSub=[homeSub,howDivid[i,0,0]];和needed一一对应的主体
    needed=[needed,howDivid[i,0,1]];得到的是每一项中的target，参考compare（list，dealConflict-fun）
  endfor
  allItem=needed[sort(needed)]
  allItem=allItem[uniq(allItem)]
  ;处理所有独特的项allItem
  foreach needone,allItem do begin
    index=where(needone eq needed)
    ratio=[]
    ;得到最大的拆分比率
    foreach items,index do ratio=[ratio,howDivid[items,1,0]]
    ratio=(where( ratio eq max(ratio)))[0];此时ratio变成index的系数
    items=index[ratio]
    nums=howDivid[items,0,2]
    n=n_elements(fidList)-1
    divid=howDivid[items,3];分别为主体代表，客体代表--也是列表list一枚[2,3]
    ;thresholdFile[needone]
    tempF=departExecute(thresholdFile[needone],fidList[n-nums],divid)
    ;替换相应文件
    outFile[needone]=tempF
    delFile=[delFile,needone]
    print,needone
  endforeach
  ok=delTemp(thresholdFile[delFile])
  return, outFile
end

;+
; :Description:执行具体的分割动作
;clusterF:阈值化的簇影像文件；classF：聚类图像
;返回处理后的文件,object
;-
function departExecute,clusterF,classFid,divid
  compile_opt idl2
  envi_open_file, clusterF, r_fid=cfid
  IF cfid EQ -1 THEN return,1
  envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
  thresData=envi_get_data(fid=cfid,dims=dims,pos=0)
  ;envi_open_file,classF, r_fid=mfid
  mfid=classFid
  maskData=envi_get_data(fid=mfid,dims=dims,pos=0)
  data=maskData*0b
  foreach parts,divid do begin
    data+=maskData eq parts
  endforeach
  maskData=[]
  thresData=thresData*data
  tempF=saveTemple(thresData)
  thresData=[]
  data=[]
  return,tempF
end


;+
; :Description:根据冲突发现结果和相似结果，处理冲突
;conflictList
;simliarList
;file   分类图文件路径
;thresholdFile    阈值化结果文件名
;-spaticalMat,
function dealConflict,conflictList,simliarList,fidList,thresholdFile,newOut
  compile_opt idl2
  howDivid=list();怎么拆分，拆分的依据
  objects=[];列入客体的队伍
  for i=0,n_elements(conflictList)-1 do begin
    if conflictList[i] eq [] then continue;空值跳过
    print,'conflictList:',i
    foreach items,conflictList[i] do begin
      if total(items eq simliarList[i])eq 1 then continue;同类跳过
      ;if items lt i then continue;已经做过的跳过;-----这个要重写，和needDepart相关
      ;判断是否需要对它拆分;加了一个客体主体判读程序
      needDepart=objectOrSub(simliarList,i,items);(spaticalMat,i,items)
      ;如果i该项已经列入了客体,则不会再对现在的items进行拆分
      if total(i eq objects) eq 1 then needDepart=0
      print,format='($,a)',BYTE(needDepart)
      ;      full=[conflictList[items],i]
      ;      foreach part,full do begin
      ;        ;如果该item面对两类以上，就拆分,拆分时，只拆到现在重叠的部分.
      ;        if total(part eq simliarList[i])eq 0 then begin
      ;          needDepart=1
      ;          break
      ;        endif
      ;      endforeach
      ;如需要拆分，则对items做如下处理
      if needDepart eq 1 then begin
        objects=[objects,items];列入客体的队伍
        overlay=[];用于验证此拆分项合适否
        data=0;用来保存data数据
        foreach clusterF,thresholdFile[[i,items]] do begin
          ;一共有两项，最后一项是target
          envi_open_file, clusterF, r_fid=cfid
          IF cfid EQ -1 THEN return,1
          envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
          data=envi_get_data(fid=cfid,dims=dims,pos=0)
          if overlay eq [] then overlay=data else overlay += data;存在重叠
          ;raster = ENVIFIDToRaster(cfid)
          ;raster.close
          envi_file_mng,id = cfid,/remove
        endforeach
        overlay=overlay eq 2
        data=data-overlay
        if min(data) lt 0 then print,'min of data is little than 0'
        print,items
        ;返回结果结构参考函数中this（list）
        ;compare可能等于【】，目标不可分时,应跳过，等目标做主体，
        ;compare是包含所有经历可拆分项的list,后优选分离度最大一项
        compare = howDividePart(newOut[items],fidList,overlay,data)
        overlay=[]
        data=[]
        ;如果compare为空则跳过该项
        if n_elements(compare) gt 1 then begin
          ;要被分解的代号，和可以参考的分类图,以及当前的文件号i（此项没太大作用）
          compare[0]=[i,items,compare[0]]
          howDivid.add,compare
        endif
      endif else begin
        objects=[objects,i];列入客体的队伍
      endelse
    endforeach
  endfor
  if n_elements(howDivid) gt 0 then $
    outFile=dePart(fidList,thresholdFile,howDivid) else $
    outFile=thresholdFile
  print,'ik'
  return,outFile
end

;+
; :Description:divide the given part alongs its experice
;newOut[items]对所有经历可拆分项进行分析，找出可拆分比例，以便后续优选最佳拆分依据
;-
function howDividePart,now,fidList,overlay,TargetData
  compile_opt idl2

  ;now=newOut[items]
  n=n_elements(now)-1;file的个数下标
  compare=list()
  partialData=0
  for nums=0,n do begin;拆分依据newOut
    if n_elements(now[nums]) gt 1 then begin
      LineNow=now[nums]
      maxIndex=2.2;这和indexTO的因子相关
      while 1 do begin

        ;classFile=fidList[n-nums];newOut是倒序的
        ;envi_open_file, classFile, r_fid=cfid
        cfid=fidList[n-nums]
        envi_file_query,cfid,nb=nb,ns=ns,nl=nl,dims=dims
        data=envi_get_data(fid=cfid,dims=dims,pos=0)
        partialData=[]
        cost=data*0b;用来叠加计算
        for several=0,n_elements(LineNow)-1 do begin
          cost += data eq LineNow[several];他们linenow下的是互斥的
          ratio=(cost+overlay) eq 2;交集部分所能勾选出去的比例;lapO
          ratio=1.0*total(ratio)/total(overlay);面积对重叠区的可分性
          ;lapO=1.0*total(lapO)/total(cost);overlay<=cost有效性
          lapT=(cost+TargetData) eq 2;非交集的成本
          lapT=1.0*total(lapT)/total(TargetData);额外的扰动
          if several eq 0 then preT=0 else preT=partialData[1,several-1];上一个lapT
          if several eq 0 then preO=0 else preO=partialData[0,several-1];上一个ratio
          increaseT=lapT-preT;增幅，越小越好
          increaseO=ratio-preO;增幅，越da越好——2是加权因子，看重ratio所代表的意义
          indexTO=1.2*increaseO+(1-increaseT);两个数都是【0-1】之间,所以最大不过1+因子*1
          partialData =[[partialData] , [ratio,lapT,indexTO]];,lapO
        endfor
        ;raster = ENVIFIDToRaster(cfid)
        ;raster.close
        data=[]
        cost=[]
        ;sub=[0,indgen(n_elements(partialData[1,*])-1)]
        ;sub=partialData[1,*]-partialData[1,sub];得到的是每一个的增幅（扰动）
        a=partialData[2,*];indexTO
        LineNow=LineNow[REVERSE(SORT(a))];根据lap进行从大到小的倒序排列SORT(a)
        ;因为cost的计算是互斥的，所以ratio最大就为1
        if max(a) ge maxIndex then break;目的是为了有效选择合适的用于区分的类【LineNow】
        maxIndex = max(a);有些不止两次，需要再反转一次，但必须至少两次
      endwhile

      sub=indgen(n_elements(partialData[1,*]))+1;从第一个到最后一个,并重复最后一个
      sub=partialData[sub]-partialData[2,*];计算开始lap减少的系数,开始不划算的地方
      index=(where(sub le 0))[0]
      this=list()
      this.add,nums;用到的文件号,这个用的时候要用n-nums
      this.add,partialData[2,index];现在所能达到的最大的分离度
      this.add,LineNow[0:index];用到的类号
      if index eq n_elements(LineNow)-1 then divid=[] else divid=LineNow[index+1:-1]
      this.add,divid;所有的类号
      if index lt n_elements(LineNow)-1 then compare.add,this;不具可分性时，抛弃
    endif
  endfor
  ;即便本体被处理过，现在肯定还有重叠，因为重叠不具有传递性（根据newout经历来看）
  ;但如果目标也不可分，应该跳过compare
  ;选择分离度最大的一组
  divid=[]
  for i=0,n_elements(compare)-1 do divid=[divid,compare[i,1]];选择compare中分离度最大的
  if divid ne [] then begin;如果目标也不可分,compare直接被输出
    dx=where(divid eq max(divid))
    compare = compare[dx[0]] ;结构等价于this
  endif
  return,compare
end




;+
; :Description:根据更新类，找到有冲突的对象
;
;-
function conflictFinder,newOut,simliarList;,poentenOut
  compile_opt idl2
  conflictList=list()
  for i=0,n_elements(newOut)-1 do begin
    this=[]
    now=newOut[i]
    for j=0,n_elements(newOut)-1 do begin
      if j eq i then continue
      rightNow=newOut[j]
      Isbreak=0
      for k=0, n_elements(now)-1 do begin
        foreach new,now[k] do begin
          if total(new eq rightNow[k]) eq 1 then begin
            this=[this,j];[[this],[j,k]]
            Isbreak=1
            break
          endif
        endforeach
        if Isbreak eq 1 then break
      endfor
    endfor
    newThis=[]
    foreach part,this do begin
      if total(part eq simliarList[i]) eq 1 then continue
      newThis=[newThis,part]
    endforeach
    conflictList.add,newThis
  endfor

  return, conflictList
end


;+
; :根据最大化的聚类，优选直接贡献类别:首先阈值选择，然后计算history各元素的占比
; fidList是file的fid列表
; maskData是满足阈值的mtkl类
; clus是当前的聚类号
; outterRelate是file聚类图间的相互关系，相邻两项
;-
function updateHis,fidList,mtklFile,outterRelate,fileNum,thresholdFile
  compile_opt idl2
  ;验证fileNum范围
  maxLimit=n_elements(outterRelate[-1])-1
  limit='['+strtrim(string(fix(0)),2)+','+strtrim(string(fix(maxLimit)),2)+']'
  if (fileNum lt 0) or (fileNum gt maxLimit) then begin
    ;print,'input a suitable limit:' + limit
    return,0
  endif
  ;获取maskData
  threshold=ceil(n_elements(fidList)/2.0)
  envi_open_file,mtklFile[fileNum], r_fid=mfid
  envi_file_query,mfid,dims=dims
  maskData = envi_get_data(fid=mfid,dims=dims,pos=0)
  ;print,max(maskData)
  maskData=maskData ge threshold

  ;此时的clus对应了聚类号,相应文件为file中的最后一个
  traceNum=list()
  poentialNUm=list()
  hisIndex=0b
  history=hisExperience(outterRelate,fileNum)
  ;为file中最后一个文件留位置
  hisArr=[-1,history[0]];
  history[0]=fileNum
  ;print,'history',history
  foreach f,hisArr do begin
    print,format='($,a)','^'
    ;history[0]对应file的相应文件下标
    ffid=fidList[f]
    if dims eq [] then envi_file_query,ffid,dims=dims
    classData = envi_get_data(fid=ffid,dims=dims,pos=0)
    this=[]
    thisNot=[];构建连襟类，存在相关关系
    ;print,history[hisIndex]
    foreach clusNow,history[hisIndex] do begin
      clusNowData=(classData eq clusNow);1
      sumNow=total(clusNowData)
      clusNowData=clusNowData + maskData
      ;t=total(maskData)
      ratio=1.0*total(clusNowData eq 2)/sumNow;eq2/eq1
      ;print,t,sumNow
      if ratio gt 0.5 then this=[this,clusNow] else thisNot=[thisNot,clusNow]
      clusNowData=[]
    endforeach
    traceNum.add,this;对应位置的ratio
    poentialNUm.add,thisNot
    classData=[]
    ;因为是fidLIst，在本程序中多次应用，所以不要移除
    ;raster = ENVIFIDToRaster(ffid)
    ;raster.close
    ;
    hisIndex++
  endforeach
  result=list();返回值
  result.add,traceNum
  result.add,poentialNUm
  ;保存阈值文件
  tempFile=saveTemple(maskData);data
  maskData=[]
  ;移除用过了的thresholdFile
  envi_file_mng,id = mfid,/remove
  thresholdFile=[thresholdFile,tempFile]
  print,fileNum
  ;print,'traceNum',traceNum
  return,result
end


;+
; :保存临时的文件的函数，在测试函数时可临时关闭:
;name必须是string类型
;-
function saveTemple,data,proj,name
  compile_opt idl2
  e=envi()
  ;保存阈值图像
  if name eq !NULL then tempFile = e.GetTemporaryFilename()
  if name ne !NULL then begin
    tempDir = e.GETPREFERENCE('temporary_directory')
    tempFile= tempDir+name+'.dat'
  endif
  if proj eq !NULL then raster = ENVIRaster(data, URI=tempFile);clusterData
  if proj ne !NULL then raster = ENVIRaster(data, URI=tempFile,SPATIALREF=proj)
  raster.save
  raster.close

  return, tempFile
end

;+
; :Description:
;删除文件列表中的文件
;-
function delTemp,delfile
  compile_opt idl2
  foreach delF,delfile do begin
    envi_open_file,delF, r_fid=tfid
    ok=deleteFID(tfid)
  endforeach
  return, 1
end


;+
; :计算history的相似度:
; 计算筛选（updateHis）后的各个history的相识度
;
;num为最后一文件的最大类号num=n_elements(outterRelate[-1])-1
;-
function closer,traceFile,baseFile;,num,outterRelate
  compile_opt idl2
  Log = baseFile+'closerLog.txt'
  ;  traceFile=list();轨迹标记;外
  ;  num=n_elements(outterRelate[-1])
  ;  for i=0,num do begin
  ;    history=hisExperience(outterRelate,i);history 第一行是下标数组
  ;    if size(history, /TYPE) ne 11 then continue
  ;    traceFile.add,history
  ;  endfor
  num=n_elements(traceFile);--1

  ;计算相似度，形成矩阵(据类号关系)
  outMatrix=[]
  for outter=0,n_elements(traceFile)-1 do begin
    len=n_elements(traceFile[outter])
    inMatrix=[]
    for inner=0,n_elements(traceFile)-1 do begin
      ratio=0.0
      diffO=0;用来验证是否属于另一项的子集，如是，则应化为一类
      diffI=0;同样作用
      for items=1,len-1 do begin;第一个肯定不同，所以不用比
        d_o=n_elements(traceFile[outter,items])
        d_i=n_elements(traceFile[inner,items])
        degree=max([d_o,d_i]);非对称阵d_i;
        sames=0.0
        foreach s_o,traceFile[outter,items] do begin
          if total(s_o eq traceFile[inner,items]) eq 1 then sames++
          ;if (sames eq 0) and (d_o eq 0 or d_i eq 0) then diff++;如果这句从未被执行，说明属于另一项的子集
        endforeach
        r=sames/degree
        ;d=sames/d_o;相同的个数和它本身的个数
        if degree eq 0 then r=1.0
        if (d_o eq sames) then diffO++;或为空或为从属关系
        if (d_i eq sames) then diffI++;或为空或为从属关系
        ratio += r;最大值等于n_elements(history)-1,比如是7
      endfor
      if max([diffO,diffI]) eq len-1 then ratio=len-1
      inMatrix=[inMatrix,ratio]
    endfor
    outMatrix=[[outMatrix],[inMatrix]]
  endfor

  openW,lun,Log,/GET_LUN;,/append
  f=strcompress('('+String(fix(N_ELEMENTS(traceFile)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'history similiarty:'
  printF,lun,format=f,outMatrix;————保存相似矩阵
  FREE_LUN,lun

  ;求相似
  threshold=(n_elements(traceFile[0])/2.0)
  similar=list()
  ;full=[]
  for i=0,num-1 do begin
    ;if total(i eq full) eq 1 then continue
    this=where(outMatrix[i,*] ge threshold)
    similar.add,this
    ;full=[full,this]
    ;full=full[sort(full)]
    ;full=full[uniq(full)]
  endfor
  return,similar
end




;+
; :Description:计算面积相互之间所占比例，overlay/itself
;按行排列，两两相交的部分占“行头文件”的比例
;-
function spaticalRelativ,thresholdFile,baseFile
  compile_opt idl2
  Log = baseFile+'spatRltLog.txt'
  dims=[]
  thFidList=[]
  for i=0,n_elements(thresholdFile)-1 do begin
    envi_open_file,thresholdFile[i], r_fid=mfid
    if dims eq [] then envi_file_query,mfid,dims=dims
    thFidList=[thFidList,mfid]
  endfor

  spaticalMat=[]
  for i=0,n_elements(thresholdFile)-1 do begin
    dataI = envi_get_data(fid=thFidList[i],dims=dims,pos=0)
    sumI=total(dataI);行头文件
    this=[]
    for j=0,n_elements(thresholdFile)-1 do begin
      dataJ = envi_get_data(fid=thFidList[j],dims=dims,pos=0)
      dataJ=dataI+dataJ
      sumOL=total(dataJ eq 2)
      ratio=1.0*sumOL/sumI
      this=[this,ratio]
      print,format='($,a)','^'
    endfor
    print,i
    spaticalMat=[[spaticalMat],[this]]
  endfor
  openW,lun,Log,/GET_LUN;,/append
  f=strcompress('('+String(fix(N_ELEMENTS(thresholdFile)))+'(g,:,","))',/REMOVE_ALL)
  PRINTF,lun,'history similiarty:'
  printF,lun,format=f,spaticalMat;————保存相似矩阵
  FREE_LUN,lun
  ;移除thFidList
  foreach ifid, thFidList do envi_file_mng,id = ifid,/remove

  return, spaticalMat
end




