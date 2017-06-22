;所有的代码块连在一起
;all the procedure in processing
pro all
  path = 'F:\Data\test\L028_Index.tif'
  maskFile='F:\Data\IsoCopy\tpCode.tif'
  log = 'D:\360Downloads\June\log';1.txt
  reSave= 'D:\360Downloads\June\FileData';1.txt
  ;koppaFile = 'D:\360Downloads\June\koppa.txt'
  newLog = 'D:\360Downloads\June\log2.txt'

  son_koppa=0;koppa系数
  times=1;次数
  clusterF=0;存储的分类影像
  logFile = log + strtrim(string(times),2)+'.txt'
  for blocks=1,4 do begin
    ;生成报表文件，存入logfile中。共4组，以取相似
    GenerateCore,blocks,maskFile,logFile
  endfor

  while 1 do begin
    fileData=reSave+ strtrim(string(times),2)+'.txt'
    ;setbackCore,core,std,extent,alist,IsPartial,logFile,reSave
    ;所有的滤过后的核心，对应std，范围[min,max]，alist是0开始的代号，是否可以保留部分以优选，读入报表文件，存入的报表文件
    ;IsPartial设为0将全部保留所有核心，这个参数到第二遍开始设为0，现在使用默认1。只用在了test_batch和setbackCore
    test_batch,logFile,fileData,IsPartial,clusterF;将会得到一个分类图classFile，调用了sa_readDBF
    times++
    koppaFile = log +'Koppa' + strtrim(string(times),2)+'.txt'
    logFile = log + strtrim(string(times),2)+'.txt'
    sj46,koppa,clusterF,maskFile,path,koppaFile,logFile;,explore
    if abs(son_koppa-koppa) lt  0.1 then begin
      son_koppa=koppa
      break
    endif

  endwhile
end