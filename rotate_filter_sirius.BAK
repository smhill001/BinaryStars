pro rotate_filter_sirius,casenumber,smoothing=smoothing

if keyword_set(smoothing) then smoothing=smoothing else smoothing=1

case casenumber of
 1: begin
    x=read_bmp("C:\Astronomy\Images\2007\12\23\Rigel_Extreme_Rotated.bmp")
    sz=size(x,/dimensions)
    xc=393
    yc=173
    end
 2: begin
    x=read_bmp("C:\Astronomy\Images\2008\02\19\Rigel_5500mm_60sec_Setttings1__0000thru0004_Stack300NoDarkDrizzleRGBShiftWaveletsGamma_LVProc_640x480_Rotated.bmp")
    sz=size(x,/dimensions)
    xc=319
    yc=262.5
    end
 3: begin
    x=read_bmp("C:\Astronomy\Images\2009\03\14\Rigel_R3_Extreme_Rotated.bmp")
    sz=size(x,/dimensions)
    xc=385
    yc=276
    end
 4: begin
    x=read_bmp("C:\Astronomy\Images\2009\03\14\Sirius_R4.bmp")
    sz=size(x,/dimensions)
    xc=324.2
    yc=247.2
    end
 5: begin
    x=read_bmp("C:\Astronomy\Images\2009\03\14\Sirius_R5_Extreme.bmp")
    sz=size(x,/dimensions)
    xc=295
    yc=265
    end
 6: begin
    x=read_bmp("C:\Astronomy\Images\2009\03\14\Sirius_R4R5_Stacked_Derotated.bmp")
    sz=size(x,/dimensions)
    xc=301.2
    yc=206.7
    end
 7: begin
    x=read_bmp("C:\Astronomy\Images\2009\07\07\Bu648Refocus3_60sec_Settings4__Rot4_Stack300_WaveletsStar_Drizzle0thru2_WaveletsGamma.bmp")
    sz=size(x,/dimensions)
    xc=551
    yc=515
    end
 8: begin
    x=read_bmp("C:\Astronomy\Images\2009\07\07\Bu648Refocus3_60sec_Settings4__Rot5_Stack300_WaveletsStar_Drizzle0thru2_WaveletsGamma.bmp")
    sz=size(x,/dimensions)
    xc=730
    yc=435.0
    end
 9: begin
    x=read_bmp("C:\Astronomy\Images\2009\07\07\BU648_0thru2_Rotated5.bmp")
    sz=size(x,/dimensions)
    xc=553
    yc=520.0
    end
 10: begin
    x=read_bmp("C:\Astronomy\Images\2010\01\15\SIRIUS_RED_10MS_STACK42_LOG.BMP")
    sz=size(x,/dimensions)
    xc=200.5
    yc=133.3
    end
endcase


if casenumber ne 10 then begin
  xsz=sz[1]
  ysz=sz[2]

  y=fltarr(xsz,ysz)
  y=float(x(0,*,*))+float(x(1,*,*))+float(x(2,*,*))
  y=reform(y,xsz,ysz)
endif else begin
  xsz=sz[0]
  ysz=sz[1]
  ysum=fltarr(xsz,ysz)
  y=x
endelse


for i=0,359 do begin
  yprime=rot(y,i,1,xc,yc,/pivot)
  ysum=ysum+yprime
endfor

yavg=ysum/360.

yratio=smooth(y,smoothing)/yavg
ydiff=smooth(y,smoothing)-yavg

window,0,xsize=640,ysize=480
tvscl,(yratio>0.3)<3
;tvscl,(yratio[xc-200:xc+200,yc-160:yc+160]>0.3)<3
window,1,xsize=640,ysize=480
tvscl,ydiff
;tvscl,ydiff[xc-200:xc+200,yc-160:yc+160]


end