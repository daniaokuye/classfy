;所有的代码块连在一起
;.compile -v 'C:\Users\DN\Documents\IDL\classfy\all.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\classfy\sj50.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\classfy\GenerateCore.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\classfy\setbackCore.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\classfy\test_batch.pro'
;         .compile -v 'C:\Users\DN\Documents\IDL\classfy\sa_readdbf.pro'
;         .compile -v 'C:\Users\DN\Documents\IDL\classfy\sa_matrixcore.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\Default\bri_test.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\Default\test_envi_batch.pro'
;.compile -v 'C:\Users\DN\Documents\IDL\Default\calentroy.pro'

;
;all the procedure in processing
pro all
  baseFile='D:\update\process\L131037_6\'

  path = baseFile+'131037_6_I2.tif';影像文件，比如多指数文件
  maskFile=baseFile+'TP_131037.tif';历史分类图，用来记录往年的分类情况
  log = baseFile+'log\Log_';1.txt
  reSave= baseFile+'log\FileCore_';1.txt

  son_koppa=0;koppa系数
  times=1;次数
  clusterF=0;存储的分类影像
  cluList=[];保存所有的分类结果
  logFile = log + strtrim(string(times),2)+'.txt'
  locate=0;用来说明怎么分影像
  IsPartial=1;第一次全扫描，只保留部分
  cut='*******************************'
  koppa=0
  for blocks=1,4 do begin
    ;生成报表文件，存入logfile中。共4组，以取相似
    GenerateCore,blocks,path,maskFile,logFile
  endfor

  while 1 do begin
    if times gt 8 then break;次数不要太多了
    fileData=reSave+ strtrim(string(times),2)+'.txt'
    print,fileData
    ;setbackCore,core,std,extent,alist,IsPartial,logFile,reSave
    ;所有的滤过后的核心，对应std，范围[min,max]，alist是0开始的代号，是否可以保留部分以优选，读入报表文件，存入的报表文件
    ;IsPartial设为0将全部保留所有核心，这个参数到第二遍开始设为0，现在使用默认1。只用在了test_batch和setbackCore
    test_batch,path,logFile,fileData,IsPartial,clusterF,baseFile;将会得到一个分类图classFile，调用了sa_readDBF
    cluList=[cluList,clusterF]
    print,'index:',times,'koppa',koppa,cut,' &file:',clusterF
    times++
    koppaFile = log +'Koppa' + strtrim(string(times),2)+'.txt'
    logFile = log + strtrim(string(times),2)+'.txt'
    sj50,koppa,clusterF,maskFile,path,koppaFile,logFile,locate;,explore
    IsPartial=0
    print,'index:',times,'koppa',koppa,cut,'locate',locate
    print
    if abs(son_koppa-koppa) lt  0.01 then begin
      son_koppa=koppa
      break
    endif
  endwhile

  print,cluList
  print,cut,'start test_envi_Batch'
  test_envi_Batch,cluList,maskFile,baseFile


end


