
pro classfy,data,hardCore,std,extent,result
  COMPILE_OPT IDL2
  ENVI,/RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT
  ;  hardCore = [[540.184094,700.8897,1008.021662,1267.864769,2066.874847,2882.352625,2622.021813],$;1
  ;    [692.481469,889.307856,1274.78381,1625.96863,2429.703448,3406.261622,3204.793726],$;2
  ;    [603.9,747.5,982.6,997.9,681.5,291,252.9],$;3
  ;    [771,961.6,1313,1562,2272.8,2449.7,2122.7],$;4
  ;    [458.52328,587.790516,841.532569,1004.340154,1994.709161,2584.728177,2123.945396],$;5
  ;    [880.096477,1102.876267,1508.84566,1876.746566,2759.855543,3572.120086,3266.54905],$;6
  ;    [252.508152,323.416209,445.724699,480.259562,548.397351,337.231641,258.254413],$;7
  ;    [675.949519,854.376136,1191.874814,1464.986125,2422.237747,3103.907365,2697.449468]$;8
  ;    ,[329.063894,432.875514,612.841611,610.108904,316.509771,185.174325,168.352566]$;9新增一枚
  ;    ,[263.408109,355.993702,516.028895,590.844718,887.342687,489.316117,379.100871]$;10
  ;    ]
  ;
  ;  extent = [[1098,1353,1702,1906,3477,3581,3483],[318,392,675,757,1355,1752,1604],$
  ;    [1582,1803,1929,2250,3435,4052,4907],[481,616,887,930,1991,2179,2027],$
  ;    [1404,1726,2319,2308,2260,1472,1295],[-49,45,76,87,0,-12,-14],$
  ;    [6764,7297,8882,10212,13622,13698,13715],[98,168,194,259,629,150,109],$
  ;    [1347,1522,1365,1716,3207,3611,4225],[138,249,56,290,850,1175,1333]$
  ;    ,[1963,2220,2683,3055,4455,5246,9151],[289,375,279,507,1082,1862,1697],$
  ;    [621,808,1224,1080,1408,1157,980],[-49,12,77,71,-26,-63,-38],$
  ;    [1963,2322,2431,2928,4439,4227,4789],[203,428,63,510,1289,1911,1632]$
  ;    ,[646,739,996,1124,853,463,469],[67,132,193,180,100,30,55]$
  ;    ,[683,844,1136,3822,4249,5706,5976],[-49,47,146,132,34,-23,18]$
  ;    ];max,min
  ;
  ;  std = [[62.7,76.8,105.8,140.5,269.6,233.4,283.9],$
  ;    [61.7,73.9,105.6,125.9,174.1,174.8,216.2],$
  ;    [236.7,272,334,345.9,356,157.2,122.4],$
  ;    [280.7,327.9,402,475.6,562.9,1123.1,1080.7],$
  ;    [75.658857,91.282232,117.207287,158.095632,276.728366,242.120154,296.923795]$
  ;    ,[154.436417,160.486879,181.47766,196.377617,225.118551,221.230583,285.044966],$
  ;    [84.569736,102.237815,146.190622,149.310578,225.150392,155.316177,103.977493],$
  ;    [125.33415,136.800838,161.004507,190.23313,275.854115,185.656553,260.972008]$
  ;    ,[74.743241,88.878772,123.053076,135.872646,89.974296,42.591881,38.735869 ]$
  ;    ,[66.868158,79.070988,101.848703,131.290584,354.658748,262.198874,177.004472]$
  ;    ]
  ;    ;;;
  ;setbackCoreavailable
  ;;;
  ;统计表格的重新读入与编码
  reSave= 'D:\360Downloads\hardcorenew.txt'
  openW,lun,reSave,/GET_LUN
  ;  inputfile = 'F:\Data\test\nlzz_\testField\overlap.tif'
  ;  outputfile = 'D:\360Downloads\test\oooo222.tif'
  ;  logfile = 'D:\360Downloads\test\log.txt'
  ;
  ;  envi_open_file, inputfile, r_fid = fid
  ;  IF (fid EQ -1 ) THEN return
  ;  envi_file_query,fid,nb = nb,ns = ns,$
  ;    data_type = dt,nl = nl,dims = dims,file_type  =  file_type
  ;  map_info = envi_get_map_info(fid = fid)
  ;  data=[]
  ;  for i=0,nb-1 do begin
  ;    new = envi_get_data(fid = fid,dims = dims,pos = i)
  ;    data = [[[data]],[[new]]]
  ;  endfor
  ;;;
  szData=size(data)
  ns=szData[1]
  nl=szData[2]
  nb=szData[3]
  ;;;;

  S= size(hardCore,/DIMENSIONS)
  colum=s[0] & line=s[1]
  ;设置一些重要参数
  ClassNum = 1 ;定义最初的类别数
  Ulimit = 0.5;U的限值
  gap_limit = 100;填空时的阈值

  ;创建几个数组
  RuleCenter = fltarr(nb,ClassNum)
  U=fltarr(ns,nl,ClassNum)

  print,'初始化'
  for i=0,line-1 do begin
    index = (where(total(total(U,1),1) eq 0))[0];没有就是-1
    if index eq -1 then begin
      Unew=fltarr(ns,nl)
      U=[[[U]],[[Unew]]]
      Cnew=fltarr(nb)
      RuleCenter=[[RuleCenter],[Cnew]]
    endif
    U[*,*,i]=setFloate(data,hardCore[*,i],extent[*,2*i:2*i+1],std[*,i])
    rulecenter[*,i]=hardCore[*,i]
  endfor
  print,'计算'
  ;设置循环用参数
  Max_Ite=3;循环次数
  times=0;循环次数
  errAll=0.1;限差
  err=1;差值
  DISTANCE=25;针对第一波段空隙的最大间距
  while(times lt Max_Ite and err gt errALL)do begin
    print,'填空'
    gap=intarr(ns,nl)+1
    index = (where(total(total(U,1),1) eq 0))[0];没有就是-1
    line_u = (size(U,/dimensions))[2]-1
    if index gt 0 then G_max = index-1 else G_max = line_u
    for i=0,G_max do begin;U的维度数
      zeroMask = where(U[*,*,i] lt Ulimit)
      gap[zeroMask]=0
      print,i,'gap',total(gap)
    endfor
    ;直方图统计整数;以便找到空间上逻辑不符点在灰度上的特征
    if mean(data[*,*,0]) lt 10 then TT=1000 else TT=1
    x=gap*data[*,*,0]*TT
    h=histogram(x);返回的h为destiny密度，其总数为max-min的范围
    ;plot,h
    ;x=s_file(inputfile,gap)
    location=where(h eq 0,count)
    ;
    gapList=[]
    for i=0,count-2 do begin
      temp=[0,0]
      j=i
      while(location[i+1] - location[i] eq 1)do begin
        i+=1
        ;if i gt 1695 then print,'i:',i
        if i eq count-1 then break
      endwhile
      if (j ne i) then begin
        ;得到灰度值,i是count的计数。这个数加上起步值便是灰度值
        temp[0]=location[j]/TT + min(gap*data[*,*,0])
        temp[1]=location[i]/TT + min(gap*data[*,*,0])
        ;判断他们的灰度值
        gapList=[[gapList],[temp]]
      endif
    endfor
    ;连接相近的空
    i=0
    if size(gapList,/N_DIMENSIONS) gt 1 then begin;如果gaplist维度大于0然后再往下
      while( (size(gapList,/dimensions))[1]-3 gt i ) do begin
        gapDis=gapList[0,i+1]-gapList[1,i]
        if gapDis lt DISTANCE then begin
          gapList[1,i] = gapList[1,i+1]
          gapList=[[gapList[*,0:i]],[gapList[*,i+2:-1]]]
        endif else begin
          i+=1
        endelse
      endwhile

      print,'对切片统计'
      clipCore = fltarr(nb);每个波段都有对应的中心
      clipStd = fltarr(nb)
      clipExtent = fltarr(nb,2)
      for i=0,(size(gapList,/dimensions))[1] - 1 do begin
        ;既然是统计值，就没必要保留矩阵格式，这样会带入不必要的0
        minzone=where(data[*,*,0] gt gapList[0,i] and data[*,*,0] lt gapList[1,i]);得到的是满足条件的序数
        print,i,'=',size(minzone,/N_Element)
        if (size(minzone,/N_Element)) lt gap_limit then continue
        for j=0,nb-1 do begin
          clip=data[minzone+j*ns*nl]
          clipCore[j]=mean(clip)
          clipExtent[j,*]=[[max(clip)],[min(clip)]]
          clipStd[j]=stddev(clip)
        endfor
        print,'now is ith:',i
        if min(clipStd) le 0 then continue
        index = (where(total(total(U,1),1) eq 0))[0];没有就是-1
        if index lt 0 then begin
          Unew=fltarr(ns,nl)
          U=[[[U]],[[Unew]]]
          Cnew=fltarr(nb)
          RuleCenter=[[RuleCenter],[Cnew]]
        endif
        ;第二次写入
        index = (where(total(total(U,1),1) eq 0))[0]
        U[*,*,index]=setFloate(data,clipCore,clipExtent,clipStd)
        rulecenter[*,index]=clipCore
      endfor

      print,'合并center'
      ;合并center:center有一系列点组成，查看新增的每个系列与自设center有无交叉
      ;如若无交叉、处于同样区间，则合并
      index = (where(total(total(U,1),1) eq 0))[0];这是第几列，但相对于共几列
      if index eq -1 then index = (size(U,/dimensions))[2]-1
      if index gt line-1 then begin;至少有一个填充中心
        character=intarr(2,index-line+1)
        hardC=hardCore[*,sort(hardCore[0,*])];安装hardcore的第一列排序
        for i=line,index do begin
          between=rulecenter[*,i]#TRANSPOSE(intarr(line)+1) gt hardC
          between=total(between,1)
          zeroLoc=(where(between eq 0))[0]
          if zeroLoc+1 then begin
            ;between存在0，即rulecter大于所有自定义中心时
            if(between[zeroloc-1] eq colum)then begin
              character[*,i-line]=[i,zeroloc];这种不产生交叉，是可能存在同类的，需要合并
            endif else begin
              character[*,i-line]=[i,-1];这种是要保留的
            endelse
          endif else begin
            ;此时所有的都不为0,zeroloc=-1
            character[*,i-line]=[i,line]
          endelse
        endfor
      endif
      ;插入新增的中心位置
      centerIndex=line
      for i=0,index-line do begin
        if character[1,i] eq -2 then continue;-2这是后来定义的
        if character[1,i] eq -1 then begin
          rulecenter[*,centerIndex]=RuleCenter[*,character[0,i]]
        endif else begin
          sametype=where(character[1,i] eq character[1,i],count)
          if(size(sametype,/N_ELEMENTS) eq 1) then begin
            rulecenter[*,centerIndex]=RuleCenter[*,character[0,sametype]]
            U[*,*,centerIndex]=U[*,*,character[0,sametype]]
          endif else begin
            rulecenter[*,centerIndex]=total(RuleCenter[*,character[0,sametype]],2)/count
            U[*,*,centerIndex]=total(U[*,*,character[0,sametype]],3)/count
          endelse
          character[1,sametype]=-2
        endelse
        centerIndex+=1
      endfor
      RuleCenter=RuleCenter[*,0:centerIndex-1];只保留前x行，多余的将舍弃
      U=U[*,*,0:centerIndex-1];U会被重新计算
      ;89行-210行是用来自动填充空隙
    endif else begin
      times=Max_Ite
    endelse

    oldU=U
    print,'adjust U'
    m=2.5;weight
    index = (where(total(total(U,1),1) eq 0))[0];这是第几列，但相对于共几列
    if index eq -1 then index = (size(U,/dimensions))[2]-1

    U *= 0.0
    for k=0,index do begin;调整U
      DisK = calEuclideanDis(data,RuleCenter[*,k])
      for i=0,index do begin
        DisI = calEuclideanDis(data,RuleCenter[*,i])
        mask = where(DisI eq 0)
        DisI[mask]=1
        U[*,*,k] += (DisK*1.0/DisI)^(1.0/(m-1))
      endfor
      U[*,*,k] = 1.0 / U[*,*,k]
    endfor
    ;计算err限差
    err=total(sqrt((oldU-U)^2))
    print,'err:',times,'th ',err
    ;调整新增的center
    print,'adjust center'
    for i=line, index do begin
      print,'need adjusted',RuleCenter[*,i]
      xx=total(U[*,*,i]^m)
      yy=[]
      for j=0,nb-1 do begin
        yy =[yy, total( U[*,*,i]^m * data[*,*,j] )]
      endfor
      RuleCenter[*,i]=yy/xx

    endfor
    print,'rule center:',RuleCenter[0,*]
    ;    ；
    PRINTF,lun,'this is :',times,'times'
    PRINTF,lun,format='(7(g,:,","))',RuleCenter
    ;
    times++
  endwhile
  FREE_LUN,lun

  print,'据u来对影像分类'
  UMax=max(U,dimension=3)
  index = (size(RuleCenter))[2]
  print,'index of ruleC',index
  result=0
  for i=0,index-1 do begin
    Mask = U[*,*,i] eq UMax
    if(total(Mask) lt 10) then continue;如果数目太少，则跳过
    locGT=where(Mask gt 0,count);如果属于多个类别，规定给后边的类别
    if count gt 0 then result[locGT]=0
    result += Mask*(i+1)
  endfor

  ;sv=sv_IMG(map_info,outputfile,result)

end
















