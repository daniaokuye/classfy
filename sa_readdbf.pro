pro sa_readDBF,hardCore,class_name,ruleName;r_fid,data
  compile_opt IDL2
  ;  outputfile = 'D:\360Downloads\test2___.tif'
  ;  envi_file_query,fid,nl=nl,ns=ns,nb=nb,dims=dims

  e = ENVI()
  ruleName = e.GetTemporaryFilename()
  class_name = e.GetTemporaryFilename()
  ;  newRaster = ENVIRaster(data, URI=tempFile)
  ;  newRaster.Save
  ;  fid = ENVIRasterToFID(newRaster)
  data = 'F:\Data\test\L028_Index.tif'
  ruleName = 'D:\360Downloads\June\testRule.tif'
  envi_open_file, data, r_fid=fid

  envi_file_query, fid, dims=dims, nb=nb
  SN= (size(hardCore,/DIMENSIONS))[1];行
  ;minmum distance
  ;  envi_doit, 'class_doit', fid=fid, pos=indgen(nb),$
  ;    dims=dims, r_fid=r_fid, $
  ;    out_bname='min', method=1, out_name=out_name, $;
  ;    mean=hardCore, class_names=indgen(SN+1), $
  ;    ;lookup= bytarr(3,SN+1), in_memory=0
  ;    lookup= byte(randomu(1,[3,SN+1])*255), in_memory=0;,$
  ;    ;


  ;Rfid=r_fid
  ;ok=deleteFID(fid)
  ;
  ; SAM or maybe the method can be replaced by maxmimum likelyhood
  envi_doit, 'class_doit', fid=fid, pos=indgen(nb),$
    dims=dims,r_fid=r_fid, $
    out_bname='SAM', method=3,  $;,
    mean=hardCore,class_names=indgen(SN+1),rule_fid=rule_fid, $;
    lookup= byte(randomu(1,[3,SN+1])*255),$
    rule_out_name=ruleName,out_name=class_name;,rule_in_memory=1, thresh=thresh
  ;移除r_fid
  foreach iFid,[r_fid,rule_fid] do begin
    raster = ENVIFIDToRaster(iFid)
    raster.close
  endforeach

  ;  envi_file_query,r_fid,dims=dims,nb=nb
  ;  ;dim=[0,0,col,0,line]
  ;  ;dim=[-1,0,end_s-star_s,0,end_l-star_l];/////////////
  ;  ;if ARRAY_EQUAL(dims,dim) then alert=DIALOG_MESSAGE('dims ne dim',TITLE='ERROR',/ERROR)
  ;  result = envi_get_data(fid = r_fid,dims =dims,pos = 0)
  ;  DELfid=[DELfid,r_Fid]


end