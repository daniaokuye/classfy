;
;+-1.首先，得到所有关于水的mask。将这个mask做膨胀
;距离先设为900m，30多个像元
;2.然后，对草、沙对象化，转变为roi
;如果对象内有一半（0.5）是和膨胀后的水无接触，
;就认定这些个对象不属于草、沙中的湿地成分
;——————这种假设是基于混合比例的假定——————
;3.将这种鉴定后的roi转回成栅格，根据相应的信息做相应调整，
;具体就是草、沙转成非湿地的草、沙
;-+
;

pro sa_reassign
  compile_opt idl2
  ;  classfile = 'D:\360Downloads\test2_grad.img';分类图
  ENVI, /restore_base_save_files ;加载核心save文件
  ENVI_batch_init
  classfile = 'F:\Data\IsoCopy\subClass.tif'
  outfile = 'D:\360Downloads\ad_class.tif'
  ;  outfile = 'D:\360Downloads\t_mor.tif'
  ;  img='F:\Data\IsoCopy\subImg.tif'
  ;  roiI = 'D:\360Downloads\t_roi2Img.tif'
  ;  envi_open_file, img, r_fid=img_fid
  envi_open_file, classfile, r_fid=fid  ;返回的fid 1

  ;  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;  evfFileS=['C:\temp\envitempfileWedNov232358232016992_1.evf']
  ;  ;roi_fid = evfToRoi(evfFileS,fid)
  ;  roi_fid = evfToRoi2(evfFileS,fid)
  ;  Filename = 'F:\Data\IsoCopy\subMaskR.roi'
  ;  ENVI_SAVE_ROIS, Filename, roi_fid
  ;  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims
  ;  envi_file_query,img_fid,nl=g_nl,ns=g_ns,nb=g_nb,dims=g_dims
  map_info = envi_get_map_info(fid = fid)
  data=envi_get_data(fid = fid,dims = dims,pos = 0)
  ;  setbackCore,available
  ;  sa_matrixCore,realloc
  ;marked = (data eq 1) + (data eq 2)*2+ (data eq 8)*8
  kernelSZ=25;kernel size
  augment=1.9;放大倍数，先故意用小树，改的时候容易些

  water=[1,2,8,12,9,10,11,13]
  identy=[3,4,5];,6,7]
  similar=14
  ;得到水的成分，然后膨胀用
  marked=intarr(ns,nl)
  for i=0, N_ELEMENTS(water)-1 do begin
    marked += (data eq water[i])
  endfor
  e = ENVI()
  tempFile = e.GetTemporaryFilename()
  MorFile = e.GetTemporaryFilename()
  ;  aufTempFile = e.GetTemporaryFilename()
  objRaster = ENVIFIDToRaster(fid)
  waterR = ENVIRaster(marked, URI=tempFile,SPATIALREF=objRaster.SPATIALREF)
  waterR.Save
  marked=[]
  print,'1',tempFile
  M_fid = ENVIRasterToFID(waterR)  ;返回的fid 2
  ENVI_Doit,'Morph_Doit', $
    FID = M_fid, $
    DIMS = dims, $
    POS = indgen(nb), $
    METHOD = 1, $
    GRAY = 1, $
    KERNEL = Fltarr(kernelSZ,kernelSZ) + 1, $
    VALUE = Fltarr(kernelSZ,kernelSZ), $
    OUT_NAME = MorFile, $
    R_FID = r_fid  ;返回的fid 3
  print,'2',MorFile
  ok=deleteFID(M_fid)
  ;  ENVI_Doit,'Morph_Doit', $
  ;    FID = M_fid, $
  ;    DIMS = dims, $
  ;    POS = indgen(nb), $
  ;    METHOD = 1, $
  ;    GRAY = 1, $
  ;    KERNEL = Fltarr(kernelSZ*augment,kernelSZ*augment) + 1, $;会自动采用floor（）的方式
  ;    VALUE = Fltarr(kernelSZ*augment,kernelSZ*augment), $
  ;    OUT_NAME = aufTempFile, $
  ;    R_FID = augR_fid  ;返回的fid 4
  ;  ;:the selected water adjacted area:
  ;  maskAug=envi_get_data(fid = augR_fid,dims = dims,pos = 0);属于水体周边的范围做成一个mask
  ;2.另外部分，非水体转成对象，验证距离
  ;因为这是知道分类结果的每一个值。而且不担心nondata值的影响
  nonWater=intarr(ns,nl)

  evfFileS=[];###
  for i=min(data),max(data) do begin
    if where(identy eq i) ne -1 then begin;where成立时：
      nonWater = (data eq i);*i+
      if(total(nonwater) eq 0) then continue

      ;nonWater = nonWater;*maskAug
      tempFile = e.GetTemporaryFilename();***************
      nonWaterR = ENVIRaster(nonWater, URI=tempFile,SPATIALREF=objRaster.SPATIALREF)
      nonWaterR.Save
      nonWater=[]
      print,'3',tempFile
      nW_fid = ENVIRasterToFID(nonWaterR)  ;返回的fid 4

      evfFile=e.GetTemporaryFilename('evf');*******************
      print,'4',evfFile
      ENVI_DOIT,'rtv_doit',$
        fid=nW_fid, pos=0, dims=dims, $
        IN_MEMORY = 0, $
        values=1,  $;l_name=l_name,
        out_names=evfFile
      evfFileS=[evfFileS,evfFile]
      print,'wow'
      ok=deleteFID(nW_fid)
    endif
  endfor

  ;roi_fid = evfToRoi(evfFileS,fid)
  roi_fid = evfToRoi2(evfFileS,fid)
  ;roi转存出，以查看你分布
  Filename = 'F:\Data\IsoCopy\subMaskR.roi'
  ENVI_SAVE_ROIS, Filename, roi_fid

  ;3.使用ENVI中roi的统计功能，对roi内膨胀体的成分分析
  num_classes = n_elements(roi_fid)
  class_values=intarr(num_classes)
  for j=0, num_classes-1 do begin
    ; get the statistics for each selected class
    roi_dims=[envi_get_roi_dims_ptr(roi_fid[j]),0,ns-1,0,nl-1]
    ;ENVI_GET_ROI_INFORMATION,roi_dims[0] , NPTS=npts;roi包含的像元数目
    ;if npts gt 0 then begin
    envi_doit, 'envi_stats_doit', fid=r_fid, pos=0, $
      dims=roi_dims, comp_flag=1, mean=c_mean
    ;endif else c_mean=0
    class_values[j]=round(c_mean)
  endfor


  print,'statics is over, then is  roi to img'
  ;4.roi转存回栅格数据，发生关系的地方在values这。
  tempFile = e.GetTemporaryFilename('tif')
  ENVI_DOIT, 'ENVI_ROI_TO_IMAGE_DOIT', CLASS_VALUES=class_values, $
    FID=fid, IN_MEMORY=0, OUT_NAME=tempFile , ROI_IDS=roi_fid, R_FID=Img_fid;[variable]
  print,'5',tempFile
  ;5.删除有关的roi，以及其他不用的raster

  print,'what is up'
  ;对身份确认完以后，开始进行转移
  transID=(envi_get_data(fid = Img_fid,dims = dims,pos = 0));框入的有机会变身
  result=bytarr(ns,nl)
  for i=min(data),max(data) do begin
    if where(identy eq i) ne -1 then begin;where成立时：
      ;变身的部分一部分(0)变成similar，一部分(1)变成他本身
      result +=(data eq i)*transID*i + (data eq i)*(transID eq 0)*similar;
    endif else begin;where不成立时：
      result +=(data eq i)*i
    endelse
  endfor
  ok = sv_IMG(map_info,outfile,result)
  print,'ok'
  envi_delete_rois, roi_fid,/all
end

;+
; :evf record to roi's:
;
;-
function record2roi,evf_FID,fid,i
  compile_opt idl2
  envi_file_query,fid,nl=nl,ns=ns
  record = ENVI_EVF_READ_RECORD(evf_FID, i,PARTS_PTR =ppt)
  if(N_ELEMENTS(ppt) gt 2) then begin
    ;现在隐藏了有hole的，注意。
    record=record[*,ppt[0]:ppt[1]-1]
  endif
  roi_ID = ENVI_CREATE_ROI(color=i+2,ns = ns,nl = nl)
  ;转换为文件坐标
  ENVI_CONVERT_FILE_COORDINATES,fid,xmap,ymap,record[0,*],record[1,*]
  ENVI_DEFINE_ROI, roi_id, /polygon, xpts=REFORM(xMap), ypts=REFORM(yMap)

  result=[i,roi_ID]
  return, result
end
;+
; :evf record to roi's:
;
;-
function record2roi2,evf,fid,sz
  compile_opt idl2
  ENVI, /restore_base_save_files ;加载核心save文件
  ENVI_batch_init
  ;return,LONARR(15000)
  ;print,'nowwwwwwwwww'
  evf_FID = ENVI_EVF_OPEN(evf)
  envi_evf_info, evf_FID, num_recs=num_recs
  roi_ids = LONARR(num_recs)
  ;return,num_recs
  ;envi_file_query,fid,nl=nl,ns=ns,dims=dims
  nl=sz[0]
  ns=sz[1]
  reSave= 'D:\360Downloads\evfFid.txt'
  openW,lun,reSave,/GET_LUN
  PRINTF,lun,'\n1.this is FileData one\n'
  FOR i=0,num_recs-1 DO BEGIN
    record = ENVI_EVF_READ_RECORD(evf_FID, i,PARTS_PTR =ppt)
    roi_ID = ENVI_CREATE_ROI(color=i+2,ns = ns,nl = nl)
    ;转换为文件坐标
    ENVI_CONVERT_FILE_COORDINATES,fid,xmap,ymap,record[0,*],record[1,*]
    if( N_ELEMENTS(xmap) lt 10 or N_ELEMENTS(ppt) gt 2 )then begin
      ENVI_DEFINE_ROI, roi_id, /POINT, xpts=REFORM(xMap), ypts=REFORM(yMap)
    endif else begin
      ENVI_DEFINE_ROI, roi_id, /polygon, xpts=REFORM(xMap), ypts=REFORM(yMap)
    endelse
    roi_ids[i] = roi_id
    if (i mod 100)eq 0 then print,format='($,a)',i
    ENVI_GET_ROI_INFORMATION, roi_id, NPTS=npts
    PRINTF,lun,format='(6(g,:,","))',[roi_id,N_ELEMENTS(xmap),npts]
  endfor
  envi_evf_close, evf_FID
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  FREE_LUN,lun
  ;ENVI_BATCH_EXIT
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  return,roi_ids
end
;+
; :evfToRoi:
;
;-
function evfToRoi,evfFileS,fid
  compile_opt idl2
  roiIDS=[]
  for j=0,N_ELEMENTS(evfFileS)-1 do begin

    evf_FID = ENVI_EVF_OPEN(evfFileS[j]) ;返回的矢量fid 1

    envi_evf_info, evf_FID, num_recs=num_recs
    roi_ids = LONARR(num_recs)
    ;读取各个记录的点数
    ;duo xian cheng
    threads=400;the number of the processes used
    k=0;k=times mod threads;由row得到k
    pos=bindgen(threads);统计进程闲置的个数
    p=objarr(threads);object of thread
    signal=intarr(threads);status of thread
    reserve=[]
    FOR i=0,num_recs-1 DO BEGIN
      ;;;;;;

      p[pos[k]] = OBJ_NEW('IDL_IDLBridge')
      p[pos[k]]->setvar,'ToRoi',record2roi(evf_FID,fid,i)
      p[pos[k]]->Execute,"ret=ToRoi",/nowait
      k+=1

      if( k ge N_ELEMENTS(pos) ) then  begin
        ;reserve=[reserve,pos]
        pos=-1
        while(ARRAY_EQUAL(pos,-1)) do begin
          ;for n=0,N_ELEMENTS(reserve)-1 do signal[reserve[n]]=p[reserve[n]]->Status()
          for n=0,threads-1 do signal[n]=p[n]->Status()
          pos=where(signal eq 0)
        endwhile

        ;reserve=bindgen(threads)
        ;reserve[pos]=-1
        ;reserve=where(reserve gt -1)

        for n=0,N_ELEMENTS(pos)-1 do begin
          thread_id = pos[n]
          x=p[thread_id]->Getvar('ret')
          roi_ids[x[0]] = x[1]
        endfor
        k=0
      endif
      ;如果到限制了，就清空。
      if i eq num_recs-1 then begin
        while(N_ELEMENTS(pos) ne threads) do begin
          ;for n=0,N_ELEMENTS(reserve)-1 do signal[reserve[n]]=p[reserve[n]]->Status()
          for n=0,threads-1 do signal[n]=p[n]->Status()
          pos=where(signal eq 0)
        endwhile
        for n=0,N_ELEMENTS(pos)-1 do begin
          thread_id = pos[n]
          x=p[thread_id]->Getvar('ret')
          roi_ids[x[0]] = x[1]
        endfor
      endif
      if (i mod 100)eq 0 then print,format='($,a)',i
      ;;;;;;;;;;

    endfor
    envi_evf_close, evf_FID
    roiIDS=[roiIDS,roi_ids]
  endfor

  return, roiIDS
end


function evfToRoi2,evfFileS,fid
  compile_opt idl2
  roiIDS=[]
  ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  threads=min([4,N_ELEMENTS(evfFileS)]);the number of the processes used
  k=0;k=times mod threads;由row得到k
  pos=bindgen(threads);统计进程闲置的个数
  p=objarr(threads);object of thread
  signal=intarr(threads);status of thread
  fun_pos = 'C:\Users\DN\Documents\IDL\Default\sa_reassign.pro'
  envi_file_query,fid,nl=nl,ns=ns
  sz=[nl,ns]
  ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  for j=0,N_ELEMENTS(evfFileS)-1 do begin

    evf= evfFileS[j] ;返回的矢量fid 1
    roi_ids=record2roi2(evf,fid,sz)
    roiIDS=[roiIDS,roi_ids]
    continue

    ;################%%%%%%%%%%%%%%%%
    p[pos[k]] = OBJ_NEW('IDL_IDLBridge')
    p[pos[k]]->setvar,'evf',evfFileS[j]
    p[pos[k]]->setvar,'fid',fid
    p[pos[k]]->setvar,'sz',sz
    p[pos[k]]->Execute,".compile "+"'"+fun_pos+"'"
    p[pos[k]]->Execute,"ret = record2roi2(evf,fid,sz)",/nowait
    ;p[pos[k]]->setvar,'ToRoi',record2roi2(evf_FID,fid)
    ;p[pos[k]]->Execute,"ret=ToRoi",/nowait
    k+=1
    if( k ge N_ELEMENTS(pos) or j eq N_ELEMENTS(evfFileS)-1) then  begin
      if k ge N_ELEMENTS(pos) then begin
        pos=-1
        while(ARRAY_EQUAL(pos,-1)) do begin;存在这样的点就执行
          for n=0,threads-1 do signal[n]=p[n]->Status()
          pos=where(signal eq 0)
        endwhile
      endif
      if j eq N_ELEMENTS(evfFileS)-1 then begin
        while (N_ELEMENTS(pos) lt threads) do begin;满足所有的点都存在
          for n=0,threads-1 do signal[n]=p[n]->Status()
          pos=where(signal eq 0)
        endwhile
      endif
      for n=0,N_ELEMENTS(pos)-1 do begin
        thread_id = pos[n]
        roi_ids=p[thread_id]->Getvar('ret')
        roiIDS=[roiIDS,roi_ids]
      endfor
      k=0
    endif
    ;    ;如果到限制了，就清空。
    ;    if  then begin
    ;
    ;      for n=0,N_ELEMENTS(pos)-1 do begin
    ;        thread_id = pos[n]
    ;        roi_ids=p[thread_id]->Getvar('ret')
    ;        roiIDS=[roiIDS,roi_ids]
    ;      endfor
    ;    endif
    ;#######################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  endfor

  return, roiIDS
end