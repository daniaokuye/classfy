;指数转成整数，然后stack

pro index2INTER
  path = 'F:\Data\test\L028_Index.tif'
  envi_open_file, path, r_fid=fid
  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims,BNAMES=bn

  ;data = envi_get_data(fid = fid,dims = dims,pos = 0)
  envi_doit, 'envi_stats_doit', fid=fid, pos=INDGEN(nb), $
    dims=dims, comp_flag=2, dmin=dmin, dmax=dmax, $
    mean=mean, stdv=stdv, hist=hist

  ;  min=min(histStat,DIMENSION=1);最值同时也是序列值
  ;  max=max(histStat,DIMENSION=1)
  ;  nbins=max-min
  minOne=[]
  maxOne=[]
  for i=0,nb-1 do begin
    span=where(hist[*,i] gt 0)
    minOne=[minOne,span[0]]
    maxOne=[maxOne,span[-1]]
  endfor
  print,minOne,'o',maxOne
  nbins=maxOne-minOne+1
  BINSIZE=(dmax-dmin)/nbins
  ;得到有有效数目的点位区间
  minThers=[]
  maxThers=[]
  for i=0,nb-1 do begin
    span=where(hist[*,i] ge 6)
    minThers=[minThers,span[0]]
    maxThers=[maxThers,span[-1]]
  endfor


  ;范围确定
  minThers = (minThers - minOne)* BINSIZE + dmin
  maxThers = (maxThers - minOne)* BINSIZE + dmin
  print,minThers,'o',maxThers
  e = ENVI()
  S_fid=intarr(nb);保存seg的各个返回fid
  for i=0,nb-1 do begin
    tempFile = e.GetTemporaryFilename()
    ENVI_Doit, 'Stretch_Doit',POS = i, $
      FID = fid,DIMS = dims,METHOD = 1,r_fid=temp, $
      I_MIN = minThers[i], I_MAX = maxThers[i],$
      RANGE_BY = 1,OUT_MIN = 0,OUT_MAX = 255, $
      OUT_DT = 1,OUT_BNAME=bn[i], OUT_NAME = tempFile
    S_fid[i]=temp
  endfor
  tempFile = e.GetTemporaryFilename()
  proj = ENVI_GET_PROJECTION(FID=fid, PIXEL_SIZE=ps)
  ENVI_DOIT, 'ENVI_LAYER_STACKING_DOIT',dims=dims,$
    fid=S_fid,OUT_DT=2,OUT_PROJ=proj,$;2-integer;1-byte
    OUT_PS=ps, POS=indgen(nb),$
    R_FID=v,OUT_NAME = tempFile
  ok=deleteFID(S_fid);第一步的文件删除，大类数目
  print,'oa'
end



